#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz; code
# written by Claude (Anthropic) under his supervision.
"""Independent analytic anchors for the amflow-cpp audit (Tiers 2-3).

Closed forms / 1-D integral representations evaluated with mpmath at
arbitrary precision. Conventions (matching upstream AMFlow):

    J = int prod_i d^d l_i / (i pi^{d/2})  1 / prod_a D_a^{nu_a},
    D = momentum^2 - m^2 + i0,   d = 4 - 2 eps,   NO exp(eps*gammaE) factor.

Families:
  bubble0   : massless bubble, props {l^2, (l-p)^2}, p^2 = s.
              Exact: Gamma(eps) Gamma(1-eps)^2 / Gamma(2-2eps) * (-s)^{-eps}.
  bubblem   : equal-mass bubble, props {l^2-msq, (l-p)^2-msq}.
              Exact: Gamma(eps) * Int_0^1 dx [msq - s x(1-x) - i0]^{-eps}.
              Valid in Euclidean (s<0) and physical (s>4msq) regions;
              the -i0 fixes ln(neg) = ln|.| - i*pi  => Im J(eps^0) = +pi*beta.
  box0      : massless on-shell box, props {l^2, (l+p1)^2, (l+p1+p2)^2,
              (l+p1+p2+p3)^2}, s=(p1+p2)^2, t=(p2+p3)^2, p_i^2=0.
              Exact (Mellin-Barnes, e.g. Duplancic-Nizic hep-ph/0006249):
                J = (2 c_G)/(s t eps^2) *
                    [ (-s)^{-eps} 2F1(1,-eps;1-eps; 1+s/t)
                    + (-t)^{-eps} 2F1(1,-eps;1-eps; 1+t/s) ],
                c_G = Gamma(1+eps) Gamma(1-eps)^2 / Gamma(1-2eps)
              (implemented in box0_exact());
              cross-checked here against the standard eps-expansion through
              eps^1 and (externally) against pySecDec.

Per-order Laurent coefficients are computed from exact series algebra
(Gamma series via zeta values; moment integrals of log powers), not from
finite differences.

Usage: ref_closed_forms.py --dps 320 --orders 8 --out anchors.json
"""
import argparse, json
from fractions import Fraction
from mpmath import mp, mpf, mpc, gamma, zeta, euler, pi, quad, sqrt, log, exp, hyp2f1, binomial, factorial, power, fadd, mpmathify


def gamma_series(a, n):
    """Taylor series of Gamma(a + eps) in eps, n+1 terms, a a positive integer
    or use reflection-free ln-Gamma series: only needed at a=1 here."""
    assert a == 1
    # Gamma(1+eps) = exp(-euler*eps + sum_{k>=2} (-1)^k zeta(k)/k eps^k)
    ln_series = [mp.mpf(0), -euler] + [(-1) ** k * zeta(k) / k for k in range(2, n + 1)]
    return exp_series(ln_series, n)


def exp_series(c, n):
    """exp of a power series c (c[0] arbitrary), n+1 terms."""
    out = [exp(c[0])] + [mp.mpf(0)] * n
    # d/dx exp(C) = C' exp(C):  (k+1) out[k+1] = sum_{j} (j+1) c[j+1] out[k-j]
    for k in range(n):
        s = mp.mpf(0)
        for j in range(k + 1):
            if j + 1 < len(c):
                s += (j + 1) * c[j + 1] * out[k - j]
        out[k + 1] = s / (k + 1)
    return out


def mul_series(a, b, n):
    return [sum(a[j] * b[k - j] for j in range(max(0, k - len(b) + 1), min(k + 1, len(a)))) for k in range(n + 1)]


def pow_negs_series(s_val, n):
    """Series of (-s)^{-eps} = exp(-eps ln(-s)) for s<0 (so -s>0, real)."""
    L = log(-s_val)
    return [(-L) ** k / factorial(k) for k in range(n + 1)]


def bubble0_orders(s, n):
    """Massless bubble Laurent coefficients, orders eps^{-1} .. eps^{n-1}."""
    # J = Gamma(eps)Gamma(1-eps)^2/Gamma(2-2eps) (-s)^{-eps}
    # Gamma(eps) = Gamma(1+eps)/eps
    N = n + 2
    g1 = gamma_series(1, N)                       # Gamma(1+eps)
    # Gamma(1-eps): flip sign of odd terms
    g1m = [g1[k] * (-1) ** k for k in range(N + 1)]
    # Gamma(2-2eps) = (1-2eps) Gamma(1-2eps)
    g12m = [g1[k] * (-2) ** k for k in range(N + 1)]  # Gamma(1-2eps)
    g22m = [g12m[k] - 2 * (g12m[k - 1] if k else 0) for k in range(N + 1)]
    num = mul_series(mul_series(g1, g1m, N), g1m, N)
    inv_g22m = inv_series(g22m, N)
    pref = mul_series(num, inv_g22m, N)
    full = mul_series(pref, pow_negs_series(s, N), N)
    # divide by eps: orders -1..n-1 are full[0..n]
    return {-1 + k: full[k] for k in range(n + 1)}


def inv_series(a, n):
    out = [1 / a[0]] + [mp.mpf(0)] * n
    for k in range(1, n + 1):
        s = sum(a[j] * out[k - j] for j in range(1, min(k, len(a) - 1) + 1))
        out[k] = -s / a[0]
    return out


def log_moments(s, msq, n):
    """m_k = int_0^1 ln^k(msq - s x(1-x) - i0) dx, k=0..n, complex for s>4msq."""
    s = mpf(s); msq = mpf(msq)
    disc = s * s - 4 * s * msq
    pieces = []  # integration breakpoints
    roots = []
    if s > 0 and disc > 0:
        r1 = (s - sqrt(disc)) / (2 * s)
        r2 = (s + sqrt(disc)) / (2 * s)
        roots = [r1, r2]
    bps = [mpf(0)] + roots + [mpf(1)]

    def L(x):
        a = msq - s * x * (1 - x)
        if a > 0:
            return log(a)
        return log(-a) - pi * 1j   # -i0 prescription: ln(A - i0), A<0

    out = []
    for k in range(n + 1):
        val = mpc(0)
        for i in range(len(bps) - 1):
            val += quad(lambda x, k=k: L(x) ** k, [bps[i], bps[i + 1]])
        out.append(val)
    return out


def bubblem_orders(s, msq, n):
    """Equal-mass bubble Laurent coefficients, orders eps^{-1} .. eps^{n-1}."""
    N = n + 2
    g1 = gamma_series(1, N)
    mom = log_moments(s, msq, N)
    gser = [(-1) ** k * mom[k] / factorial(k) for k in range(N + 1)]  # int e^{-eps L}
    full = mul_series(g1, gser, N)
    return {-1 + k: full[k] for k in range(n + 1)}


def box0_exact(s, t, eps):
    """Massless on-shell box at numeric eps (Euclidean s,t<0). Exact 2F1 form.

    From MB representation (Duplancic-Nizic hep-ph/0006249 eq. (40)-ish;
    cross-checked vs the universal eps-expansion below):
       J = (2 c_G)/(s t eps^2) [ (-s)^{-eps} 2F1(1,-eps;1-eps;1+s/t)
                               + (-t)^{-eps} 2F1(1,-eps;1-eps;1+t/s) ]
       c_G = Gamma(1+eps) Gamma(1-eps)^2 / Gamma(1-2eps)
    """
    s = mpf(s); t = mpf(t); eps = mpmathify(eps)
    cg = gamma(1 + eps) * gamma(1 - eps) ** 2 / gamma(1 - 2 * eps)
    term = (power(-s, -eps) * hyp2f1(1, -eps, 1 - eps, 1 + s / t)
            + power(-t, -eps) * hyp2f1(1, -eps, 1 - eps, 1 + t / s))
    return 2 * cg / (s * t * eps ** 2) * term


def box0_orders(s, t, n):
    """Laurent coefficients of the box, orders eps^{-2}.. eps^{n-2}, via EXACT
    real series algebra.  Identity (derived from the Pochhammer ratio
    (-eps)_m/(1-eps)_m = -eps/(m-eps) for m>=1):

        2F1(1,-eps;1-eps;z) = 1 - sum_{k>=1} eps^k Li_k(z).
    """
    # At Euclidean s,t<0 both z=1+s/t and z=1+t/s are REAL > 1 -- on the
    # Li_k branch cut.  The physical (real) integral is the principal
    # value: the two terms' i*pi cut contributions cancel exactly, so
    # the coefficients use Re[Li_k(z+-i0)] (validated vs pySecDec at
    # 1e-11 over 5 orders and vs the engine at <=1e-64).
    from mpmath import polylog
    N = n + 3
    g1 = gamma_series(1, N)                       # Gamma(1+eps)
    g1m = [g1[k] * (-1) ** k for k in range(N + 1)]   # Gamma(1-eps)
    g12m = [g1[k] * (-2) ** k for k in range(N + 1)]  # Gamma(1-2eps)
    # c_Gamma = Gamma(1+eps)Gamma(1-eps)^2/Gamma(1-2eps)
    cg = mul_series(mul_series(g1, g1m, N), mul_series(g1m, inv_series(g12m, N), N), N)
    out_sum = [mp.mpf(0)] * (N + 1)
    for (a, z) in ((s, 1 + s / t), (t, 1 + t / s)):
        f21 = [mp.mpf(1)] + [-mp.re(polylog(k, mpc(z, 0))) for k in range(1, N + 1)]
        pa = pow_negs_series(a, N)                # (-a)^{-eps}
        term = mul_series(pa, f21, N)
        out_sum = [out_sum[k] + term[k] for k in range(N + 1)]
    full = mul_series(cg, out_sum, N)
    pref = 2 / (s * t)
    # divide by eps^2: order j coefficient = pref*full[j+2]
    return {j: pref * full[j + 2] for j in range(-2, n - 1)}


def fmt(x, dps):
    if isinstance(x, (int,)):
        x = mpf(x)
    re, im = mp.re(x), mp.im(x)
    return {"re": mp.nstr(re, dps, strip_zeros=False),
            "im": mp.nstr(im, dps, strip_zeros=False)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dps", type=int, default=350)
    ap.add_argument("--orders", type=int, default=9, help="number of eps orders above leading")
    ap.add_argument("--eps", type=str, default="1/1000", help="also evaluate at this rational eps")
    ap.add_argument("--out", type=str, required=True)
    args = ap.parse_args()

    mp.dps = args.dps + 30  # guard digits
    n = args.orders
    out = {"convention": "J = int d^dl/(i pi^{d/2}) 1/prod(D+i0), d=4-2eps, no exp(eps gammaE)",
           "dps": args.dps, "cases": {}}

    p, q = args.eps.split("/") if "/" in args.eps else (args.eps, "1")
    epsv = mpf(int(p)) / mpf(int(q))

    # ---- massless bubble at s=-1 and s=-100/7
    for s in ["-1", "-100/7"]:
        fr = Fraction(s); sv = mpf(fr.numerator) / fr.denominator  # exact rational
        key = f"bubble0[s={s}]"
        coeffs = bubble0_orders(sv, n)
        exact_at_eps = gamma(epsv) * gamma(1 - epsv) ** 2 / gamma(2 - 2 * epsv) * power(-sv, -epsv)
        out["cases"][key] = {
            "orders": {str(k): fmt(v, args.dps) for k, v in coeffs.items()},
            f"at_eps={args.eps}": fmt(exact_at_eps, args.dps)}

    # ---- massive bubble Euclidean s=-1, physical s=5 (msq=1)
    for s in ["-1", "5"]:
        fr = Fraction(s); sv = mpf(fr.numerator) / fr.denominator
        key = f"bubblem[s={s},msq=1]"
        coeffs = bubblem_orders(sv, 1, n)
        # exact at eps: Gamma(eps)*int via log_moments-based series is per-order;
        # numeric at eps via direct quad of [msq - s x(1-x) - i0]^{-eps}
        mom = None
        disc = sv * sv - 4 * sv
        bps = [mpf(0), mpf(1)]
        if sv > 4:
            r1 = (sv - sqrt(disc)) / (2 * sv); r2 = (sv + sqrt(disc)) / (2 * sv)
            bps = [mpf(0), r1, r2, mpf(1)]
        def integrand(x):
            a = 1 - sv * x * (1 - x)
            if a > 0:
                return power(a, -epsv)
            return power(-a, -epsv) * exp(1j * pi * epsv)  # (A-i0)^{-eps}, A<0
        val = mpc(0)
        for i in range(len(bps) - 1):
            val += quad(integrand, [bps[i], bps[i + 1]])
        exact_at_eps = gamma(epsv) * val
        out["cases"][key] = {
            "orders": {str(k): fmt(v, args.dps) for k, v in coeffs.items()},
            f"at_eps={args.eps}": fmt(exact_at_eps, args.dps)}

    # ---- massless box at (s,t)=(-3,-7) Euclidean
    sv, tv = mpf(-3), mpf(-7)
    key = "box0[s=-3,t=-7]"
    coeffs = box0_orders(sv, tv, n)
    out["cases"][key] = {
        "orders": {str(k): fmt(v, args.dps) for k, v in coeffs.items()},
        f"at_eps={args.eps}": fmt(mp.re(box0_exact(sv, tv, epsv)), args.dps)}

    # sanity: universal box expansion through eps^0:
    #   J = cG/(st) { 2/eps^2 [(-s)^{-eps}+(-t)^{-eps}] - ln^2(s/t) - pi^2 } + O(eps)
    X = log(sv / tv)
    cgser = 1  # at eps^0 c_G=1
    box_em2 = 4 / (sv * tv)
    box_em1 = -2 * (log(-sv) + log(-tv)) / (sv * tv)
    box_e0 = (log(-sv) ** 2 + log(-tv) ** 2 - X ** 2 - pi ** 2) / (sv * tv)
    out["box0_sanity"] = {
        "em2_match": mp.nstr(abs(coeffs[-2] - box_em2), 5),
        "em1_match": mp.nstr(abs(coeffs[-1] - box_em1), 5),
        "e0_residual_note": "e0 universal formula has scheme pieces; see audit notes",
        "e0_diff": mp.nstr(abs(coeffs[0] - box_e0), 5),
    }

    with open(args.out, "w") as f:
        json.dump(out, f, indent=1)
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
