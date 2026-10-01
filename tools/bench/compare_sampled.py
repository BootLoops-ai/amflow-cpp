#!/usr/bin/env python3
"""Compare sampled C++ benchmark output against a Mathematica reference JSON."""

from __future__ import annotations

import argparse
import json
import re
import sys
from dataclasses import dataclass
from decimal import Decimal, InvalidOperation, getcontext
from pathlib import Path
from typing import Any, Iterable


getcontext().prec = 120


@dataclass(frozen=True)
class IntegralKey:
    family: str
    indices: tuple[int, ...]

    def __str__(self) -> str:
        joined = ", ".join(str(x) for x in self.indices)
        return f"j[{self.family}, {joined}]"


@dataclass(frozen=True)
class ComplexValue:
    re: Decimal
    im: Decimal


@dataclass(frozen=True)
class ReferenceEntry:
    key: IntegralKey
    value: ComplexValue
    label: str
    # When set, this entry is one Laurent coefficient (eps^order) of an
    # `solve_integrals` result rather than a sampled value.  cpp_entries
    # then keys on (key, order) to look up the corresponding C++ coef.
    order: int | None = None


@dataclass(frozen=True)
class ComponentDiff:
    abs_error: Decimal
    rel_error: Decimal
    passed: bool


def parse_decimal(raw: Any) -> Decimal:
    text = str(raw).strip()
    if text.startswith("[") and text.endswith("]"):
        text = text[1:-1].strip()
    if "+/-" in text:
        text = text.split("+/-", 1)[0].strip()
    text = text.replace("`", "")
    try:
        return Decimal(text)
    except InvalidOperation as exc:
        raise ValueError(f"cannot parse decimal value {raw!r}") from exc


def parse_complex(obj: dict[str, Any]) -> ComplexValue:
    return ComplexValue(re=parse_decimal(obj.get("re", "0")),
                        im=parse_decimal(obj.get("im", "0")))


def parse_cpp_complex(obj: dict[str, Any]) -> ComplexValue:
    return ComplexValue(re=parse_decimal(obj.get("re", "0")),
                        im=parse_decimal(obj.get("im", "0")))


def parse_integral(text: str) -> IntegralKey:
    match = re.fullmatch(r"\s*j\[\s*([^,\]\s]+)\s*,\s*(.*?)\s*\]\s*", text)
    if not match:
        raise ValueError(f"cannot parse integral key {text!r}")
    family = match.group(1)
    rest = match.group(2).strip()
    indices = tuple(int(part.strip()) for part in rest.split(",") if part.strip())
    return IntegralKey(family=family, indices=indices)


def reference_entries(reference: dict[str, Any], include_auxiliary: bool) -> list[ReferenceEntry]:
    entries: list[ReferenceEntry] = []

    def add_entry(obj: dict[str, Any], value_key: str, label: str) -> None:
        # Laurent schema: entry has `coefficients` (list of {order, re, im}).
        # Emit one ReferenceEntry per Laurent term so cpp_entries can match
        # by (integral, order).
        if "coefficients" in obj:
            key = parse_integral(obj["integral"])
            for term in obj["coefficients"]:
                entries.append(ReferenceEntry(
                    key=key,
                    value=parse_complex(term),
                    label=label,
                    order=int(term["order"]),
                ))
            return
        # Sampled schema: entry has a single complex value at `value_key`.
        entries.append(ReferenceEntry(
            key=parse_integral(obj["integral"]),
            value=parse_complex(obj[value_key]),
            label=label,
        ))

    if "requested_target" in reference:
        add_entry(reference["requested_target"], "value", "requested_target")
    for obj in reference.get("requested_targets", []):
        add_entry(obj, "value", "requested_targets")
    if "top_integral" in reference:
        add_entry(reference["top_integral"], "value", "top_integral")
    for obj in reference.get("preferred_masters", []):
        add_entry(obj, "value", "preferred_masters")
    for obj in reference.get("masters_only", []):
        add_entry(obj, "mma", "masters_only")

    if include_auxiliary:
        if "scalar_master_from_same_run" in reference:
            add_entry(reference["scalar_master_from_same_run"], "value",
                      "scalar_master_from_same_run")
        for obj in reference.get("sampled_values_from_same_mma_run", []):
            add_entry(obj, "value", "sampled_values_from_same_mma_run")

    seen: set[tuple[IntegralKey, int | None]] = set()
    unique: list[ReferenceEntry] = []
    for entry in entries:
        ident = (entry.key, entry.order)
        if ident not in seen:
            unique.append(entry)
            seen.add(ident)
    return unique


def cpp_entries(
        cpp_output: dict[str, Any]
) -> dict[tuple[IntegralKey, int | None], ComplexValue]:
    """Index C++ output by (integral, order|None).

    Sampled outputs use order=None.  Laurent outputs use the per-term order
    from the C++ `coefficients` array.
    """
    out: dict[tuple[IntegralKey, int | None], ComplexValue] = {}
    for obj in cpp_output.get("result", []):
        integral = obj["integral"]
        key = IntegralKey(
            family=integral["family"],
            indices=tuple(int(x) for x in integral["indices"]),
        )
        if "coefficients" in obj:
            for term in obj["coefficients"]:
                order = int(term["order"])
                out[(key, order)] = parse_cpp_complex(term["value"])
            continue
        samples = obj.get("samples", [])
        if not samples:
            continue
        out[(key, None)] = parse_cpp_complex(samples[0]["value"])
    return out


def component_diff(expected: Decimal, actual: Decimal, rel_tol: Decimal,
                   abs_tol: Decimal) -> ComponentDiff:
    abs_error = abs(actual - expected)
    if expected == 0:
        rel_error = Decimal(0) if abs_error == 0 else Decimal("Infinity")
    else:
        rel_error = abs_error / abs(expected)
    return ComponentDiff(
        abs_error=abs_error,
        rel_error=rel_error,
        passed=(abs_error <= abs_tol or rel_error <= rel_tol),
    )


def decimal_for_report(value: Decimal) -> str:
    return f"{value:.3E}"


def compare(reference: dict[str, Any], cpp_output: dict[str, Any],
            rel_tol: Decimal, abs_tol: Decimal,
            include_auxiliary: bool) -> tuple[int, list[dict[str, Any]]]:
    refs = reference_entries(reference, include_auxiliary)
    actuals = cpp_entries(cpp_output)
    rows: list[dict[str, Any]] = []
    failures = 0

    for ref in refs:
        actual = actuals.get((ref.key, ref.order))
        if actual is None:
            rows.append({
                "integral": str(ref.key),
                "label": ref.label,
                "order": ref.order,
                "status": "missing",
            })
            failures += 1
            continue

        re_diff = component_diff(ref.value.re, actual.re, rel_tol, abs_tol)
        im_diff = component_diff(ref.value.im, actual.im, rel_tol, abs_tol)
        passed = re_diff.passed and im_diff.passed
        if not passed:
            failures += 1
        rows.append({
            "integral": str(ref.key),
            "label": ref.label,
            "order": ref.order,
            "status": "ok" if passed else "mismatch",
            "re_abs_error": str(re_diff.abs_error),
            "re_rel_error": str(re_diff.rel_error),
            "im_abs_error": str(im_diff.abs_error),
            "im_rel_error": str(im_diff.rel_error),
        })

    return failures, rows


def load_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def write_json_report(path: Path, rows: Iterable[dict[str, Any]],
                      failures: int) -> None:
    report = {
        "status": "ok" if failures == 0 else "failed",
        "failures": failures,
        "comparisons": list(rows),
    }
    with path.open("w", encoding="utf-8") as handle:
        json.dump(report, handle, indent=2, sort_keys=True)
        handle.write("\n")


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(
        description="Compare sampled C++ benchmark output against a Mathematica reference JSON.")
    parser.add_argument("reference", type=Path,
                        help="Mathematica sampled reference JSON")
    parser.add_argument("cpp_output", type=Path,
                        help="C++ amflow_cli sampled output JSON")
    parser.add_argument("--rel-tol", default="1e-25",
                        help="relative tolerance per real/imag component")
    parser.add_argument("--abs-tol", default="1e-40",
                        help="absolute tolerance per real/imag component")
    parser.add_argument("--include-auxiliary", action="store_true",
                        help="also require auxiliary values stored from the same MMA run")
    parser.add_argument("--json-out", type=Path,
                        help="optional machine-readable comparison report")
    args = parser.parse_args(argv)

    rel_tol = parse_decimal(args.rel_tol)
    abs_tol = parse_decimal(args.abs_tol)
    failures, rows = compare(
        reference=load_json(args.reference),
        cpp_output=load_json(args.cpp_output),
        rel_tol=rel_tol,
        abs_tol=abs_tol,
        include_auxiliary=args.include_auxiliary,
    )

    for row in rows:
        order_suffix = (
            f" eps^{row['order']:+d}" if row.get("order") is not None else ""
        )
        if row["status"] == "missing":
            print(f"MISSING {row['integral']}{order_suffix} ({row['label']})")
            continue
        print(
            f"{row['status'].upper():8} {row['integral']}{order_suffix} ({row['label']}): "
            f"re abs={decimal_for_report(Decimal(row['re_abs_error']))} "
            f"rel={decimal_for_report(Decimal(row['re_rel_error']))}; "
            f"im abs={decimal_for_report(Decimal(row['im_abs_error']))} "
            f"rel={decimal_for_report(Decimal(row['im_rel_error']))}")

    if args.json_out:
        write_json_report(args.json_out, rows, failures)

    print(f"Compared {len(rows)} integral(s); failures={failures}")
    return 0 if failures == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
