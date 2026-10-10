// SPDX-License-Identifier: MIT
// ode::sparse — implementation.
//

#include "amflow/ode/sparse.hpp"

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <ostream>
#include <stdexcept>
#include <utility>

#include <flint/acb.h>

namespace amflow::ode {

using numeric::AcbValue;
using numeric::acb_is_chop_zero;
using numeric::chop_pre;

namespace {

bool sparse_chop_enabled() {
    return std::getenv("AMFLOW_NO_SPARSE_CHOP") == nullptr;
}

}  // namespace

int sparse_chop_digits() {
    const char* env = std::getenv("AMFLOW_SPARSE_CHOP_DIGITS");
    if (env && *env) {
        char* end = nullptr;
        long v = std::strtol(env, &end, 10);
        if (end != env && *end == '\0') {
            if (v < 0) return 0;
            return static_cast<int>(v);
        }
    }
    // Default: chop_pre, but bumped to (working_pre - 40) for the
    // sparse layer.  In big sparse Gauss-eliminations (12-master
    // 240-order boundary-order systems are typical for multi-mass 3L
    // families) acb rounding noise compounds across divisions by
    // ill-conditioned pivots and can land a residual at ~10^-(wp/2)
    // to ~10^-(wp-30).  A 10^-20 cut keeps that noise as spurious
    // matrix entries, which re-route Gauss elimination's pivot choice
    // and inflate the `unsolved` set at the wrong columns.  Tying the
    // floor to `working_pre - 40` gives a 40-decimal-digit noise
    // budget while keeping ~80% of working precision available for
    // legitimate values.  Override via AMFLOW_SPARSE_CHOP_DIGITS.
    int chop_pre_val = chop_pre();
    int wp = numeric::working_pre();
    int floor_val = (wp > 60) ? (wp - 40) : chop_pre_val;
    return std::max(chop_pre_val, floor_val);
}

// ===========================================================================
//  Sparsify
// ===========================================================================

SparseRow sparsify_row(const std::vector<AcbValue>& dense) {
    SparseRow out;
    out.reserve(dense.size());
    for (long c = 0; c < static_cast<long>(dense.size()); ++c) {
        if (!dense[c].is_zero()) {
            out.emplace_back(c, dense[c].clone());
        }
    }
    return out;
}

SparseSystem sparsify_matrix(const std::vector<std::vector<AcbValue>>& dense) {
    SparseSystem out;
    out.reserve(dense.size());
    for (const auto& row : dense) out.push_back(sparsify_row(row));
    return out;
}

// ===========================================================================
//  Translate / Scalar / Plus
// ===========================================================================

void sparse_translate(SparseRow& row, long n) {
    for (auto& e : row) e.col += n;
}
void sparse_translate(SparseSystem& sys, long n) {
    for (auto& r : sys) sparse_translate(r, n);
}

void scalar_sparse(SparseRow& row, acb_srcptr s, long prec) {
    for (auto& e : row) acb_mul(e.val.raw(), e.val.raw(), s, prec);
}

SparseRow plus_sparse(const SparseRow& a,
                      const SparseRow& b,
                      long             prec) {
    SparseRow out;
    out.reserve(a.size() + b.size());

    std::size_t i = 0, j = 0;
    while (i < a.size() && j < b.size()) {
        if (a[i].col < b[j].col) {
            out.emplace_back(a[i].col, a[i].val.clone());
            ++i;
        } else if (a[i].col > b[j].col) {
            out.emplace_back(b[j].col, b[j].val.clone());
            ++j;
        } else {
            // Same column: add.
            AcbValue sum;
            acb_add(sum.raw(), a[i].val.raw(), b[j].val.raw(), prec);
            if (!sum.is_zero()) out.emplace_back(a[i].col, std::move(sum));
            ++i; ++j;
        }
    }
    for (; i < a.size(); ++i) out.emplace_back(a[i].col, a[i].val.clone());
    for (; j < b.size(); ++j) out.emplace_back(b[j].col, b[j].val.clone());
    return out;
}

// ===========================================================================
//  SplitSystem
// ===========================================================================

SplitResult split_system(const SparseSystem& sys, long m) {
    SplitResult r;
    r.dep.reserve(sys.size());
    r.indep.reserve(sys.size());
    for (std::size_t i = 0; i < sys.size(); ++i) {
        bool found = false;
        for (const auto& e : sys[i]) {
            if (e.col == m) { found = true; break; }
            if (e.col > m)  break;       // sorted ascending
        }
        if (found) r.dep.push_back(i);
        else       r.indep.push_back(i);
    }
    return r;
}

// ===========================================================================
//  ConstructMatrix  (.m line 605-622)
// ===========================================================================
//
//   table[[i+1, j+1]] := (totalorder + 1 - j) * dn[j - i] - an[j - i - 1]
//
//  where i in [0, total_order], j in [i, min(total_order+1, max(i+expd, i+1+expa))].
//
//  dn[k] = Coefficient[dx, eta, k] * Identity     (scalar diagonal block)
//  an[k] = Coefficient[ax, eta, k]                (full N x N block)
//
//  table[[i+1]] is then horizontally flattened, sparsified, and translated by
//  N * Range[0, totalorder] so that row block i lives in absolute columns
//  starting at N * i.

namespace {

inline acb_srcptr coef_or_zero(const std::vector<AcbValue>& exp_, long k) {
    if (k < 0 || k >= static_cast<long>(exp_.size())) return nullptr;
    return exp_[k].raw();
}

}  // namespace

ConstructedMatrix construct_matrix(const std::vector<AcbValue>& dxexp,
                                   const std::vector<std::vector<std::vector<AcbValue>>>& axexp,
                                   long total_order,
                                   long prec) {
    if (axexp.empty()) {
        throw std::invalid_argument("construct_matrix: ax must have at least one row");
    }
    long N = static_cast<long>(axexp.size());
    for (const auto& row : axexp) {
        if (static_cast<long>(row.size()) != N) {
            throw std::invalid_argument("construct_matrix: ax must be square");
        }
    }

    long expd = static_cast<long>(dxexp.size()) - 1;
    long expa = -1;
    for (const auto& row : axexp) {
        for (const auto& vec : row) {
            long deg = static_cast<long>(vec.size()) - 1;
            if (deg > expa) expa = deg;
        }
    }
    if (expa < 0) expa = 0;
    if (expd < 0) {
        throw std::invalid_argument("construct_matrix: dx is the zero polynomial");
    }

    long total_columns = (total_order + 2) * N;

    ConstructedMatrix cm;
    cm.block_size    = N;
    cm.total_order   = total_order;
    cm.total_columns = total_columns;
    cm.block_systems.resize(total_order + 1);

    // For every i (row block), build N rows.  MMA first drops zero N x N
    // blocks from the per-row-block table, then ArrayFlattens the survivors.
    // Local column offsets advance only when a whole block is kept.
    for (long i = 0; i <= total_order; ++i) {
        SparseSystem block_rows(N);

        long j_max = std::min(total_order + 1, std::max(i + expd, i + 1 + expa));
        long block_col = 0;
        for (long j = i; j <= j_max; ++j) {
            long k_dn = j - i;
            long k_an = j - i - 1;

            acb_srcptr dn_c = coef_or_zero(dxexp, k_dn);

            // Test: is this entire N×N sub-block zero?  If so skip and don't
            // advance block_col (matches MMA's drop-zero-block step).
            bool block_nonzero = false;
            for (long r = 0; r < N && !block_nonzero; ++r) {
                for (long c = 0; c < N; ++c) {
                    bool has_d = (dn_c != nullptr) && (r == c);
                    bool has_a = false;
                    if (k_an >= 0
                        && k_an < static_cast<long>(axexp[r][c].size())
                        && !axexp[r][c][k_an].is_zero()) {
                        has_a = true;
                    }
                    if (has_d || has_a) {
                        block_nonzero = true;
                        break;
                    }
                }
            }
            if (!block_nonzero) continue;

            for (long r = 0; r < N; ++r) {
                for (long c = 0; c < N; ++c) {
                    long col_local = block_col * N + c;

                    bool has_d = (dn_c != nullptr) && (r == c);
                    bool has_a = false;
                    acb_srcptr an_c = nullptr;
                    if (k_an >= 0
                        && k_an < static_cast<long>(axexp[r][c].size())
                        && !axexp[r][c][k_an].is_zero()) {
                        has_a = true;
                        an_c = axexp[r][c][k_an].raw();
                    }

                    if (!has_d && !has_a) continue;

                    AcbValue value;
                    if (has_d) {
                        // value = (totalorder + 1 - j) * dn_c
                        acb_mul_si(value.raw(), dn_c,
                                   static_cast<long>(total_order + 1 - j), prec);
                    }
                    if (has_a) {
                        // value -= an_c
                        acb_sub(value.raw(), value.raw(), an_c, prec);
                    }
                    // MMA's ConstructMatrix (DESolver.m:620-636) operates on
                    // exact rationals, so a `(k)*dn - an` cancellation gives
                    // an exact 0 and `Sparsify` discards the entry.  In acb
                    // arithmetic, the same cancellation yields a ball with
                    // midpoint 0 but nonzero radius (from set_fmpq rounding
                    // error), so `acb_is_zero` returns false and we keep a
                    // spurious nonzero entry.  Use `acb_contains_zero` so that
                    // any ball straddling zero is dropped — this matches
                    // MMA's exact-rational behavior at our working precision
                    // (legitimate nonzero values have |mid| >> rad).
                    if (acb_contains_zero(value.raw())) continue;

                    block_rows[r].emplace_back(col_local, std::move(value));
                }
            }
            ++block_col;
        }

        // Translate row block i by N * i so it lives in absolute columns
        // starting at N * i.  Mirrors MMA's
        //     MapThread[SparseTranslate, {sp, N * Range[0, totalorder]}]
        sparse_translate(block_rows, N * i);

        cm.block_systems[i] = std::move(block_rows);
    }

    return cm;
}

// ===========================================================================
//  ForwardSparseGaussian
// ===========================================================================

GaussianRowSet forward_sparse_gaussian(SparseSystem            sys,
                                       const std::vector<long>& vset,
                                       long                     prec) {
    GaussianRowSet result;
    result.reduced.reserve(vset.size());

    // MMA's exact-coefficient elimination has cancelled entries disappear
    // exactly; our acb path can leave ball-sized remnants (~1e-122) on
    // columns that should be zero.  Without chopping, downstream pivot
    // detection records the wrong leading column.  Keep the working rows
    // chopped after every step.
    if (sparse_chop_enabled()) chop_sparse(sys);

    for (long n : vset) {
        SplitResult sp = split_system(sys, n);
        if (sp.dep.empty()) continue;

        std::vector<SparseRow> dep_rows;
        dep_rows.reserve(sp.dep.size());
        for (auto idx : sp.dep) dep_rows.push_back(std::move(sys[idx]));

        SparseSystem indep_rows;
        indep_rows.reserve(sp.indep.size());
        for (auto idx : sp.indep) indep_rows.push_back(std::move(sys[idx]));

        SparseRow pivot = std::move(dep_rows[0]);

        AcbValue leading;
        for (const auto& e : pivot) {
            if (e.col == n) { leading = e.val.clone(); break; }
            if (e.col > n)  break;
        }
        if (leading.is_zero()) {
            // pivot did contain n (split confirmed it); leading=0 means
            // rounding hit zero.  Treat as non-pivotable.
            sys = std::move(indep_rows);
            for (std::size_t k = 1; k < dep_rows.size(); ++k) {
                sys.push_back(std::move(dep_rows[k]));
            }
            continue;
        }

        // ------------------------------------------------------------------
        // calcx00 audit: a pivot BALL that contains zero
        // cannot be certified nonzero; acb_inv on it returns an infinite
        // (indeterminate) ball and poisons every later row -> the NaN that
        // resurfaces downstream as the calcx00 "could not be rationalized"
        // abort.  The chop layer removes |mid| < 10^-(wp-40) entries, but an
        // entry with rad >= |mid| ABOVE that floor slips through — and since
        // the chop floor scales with working_pre, the race is precision-
        // INVARIANT (the observed goal-14/goal-25 identical failure).
        //   AMF_DBG_X00            : report such pivots ([X00-piv]).
        //   AMFLOW_X00_PIVOT_GUARD : treat the uncertifiable entry as a
        //     cancelled zero — drop it from the row and RETAIN the row for
        //     later columns (the exact-arithmetic analogue: a cancelled
        //     entry is exactly absent).  Default OFF: byte-identical.
        // ------------------------------------------------------------------
        if (acb_contains_zero(leading.raw())) {
            if (std::getenv("AMF_DBG_X00")) {
                char* ms = arb_get_str(acb_realref(leading.raw()), 8, 0);
                char* is = arb_get_str(acb_imagref(leading.raw()), 8, 0);
                std::fprintf(stderr,
                        "[X00-piv] pivot ball contains zero at col=%ld: "
                        "re=%s im=%s (guard %s)\n",
                        n, ms, is,
                        std::getenv("AMFLOW_X00_PIVOT_GUARD") ? "ON" : "off");
                flint_free(ms);
                flint_free(is);
            }
            const char* g = std::getenv("AMFLOW_X00_PIVOT_GUARD");
            if (g && *g == '1') {
                // Remove the uncertifiable head entry; keep the rest of the
                // row in play for later pivots.
                SparseRow trimmed;
                trimmed.reserve(pivot.size());
                for (auto& e : pivot) {
                    if (e.col != n) trimmed.push_back(std::move(e));
                }
                sys = std::move(indep_rows);
                if (!trimmed.empty()) sys.push_back(std::move(trimmed));
                for (std::size_t k = 1; k < dep_rows.size(); ++k) {
                    sys.push_back(std::move(dep_rows[k]));
                }
                continue;
            }
        }

        AcbValue inv_leading;
        acb_inv(inv_leading.raw(), leading.raw(), prec);
        scalar_sparse(pivot, inv_leading.raw(), prec);
        if (sparse_chop_enabled()) chop_sparse(pivot);

        std::vector<SparseRow> reduced_rest;
        reduced_rest.reserve(dep_rows.size() - 1);
        for (std::size_t k = 1; k < dep_rows.size(); ++k) {
            AcbValue factor;
            for (const auto& e : dep_rows[k]) {
                if (e.col == n) {
                    factor = e.val.clone();
                    acb_neg(factor.raw(), factor.raw());
                    break;
                }
                if (e.col > n) break;
            }
            // dep_rows[k] += factor * pivot
            SparseRow scaled = pivot;
            scalar_sparse(scaled, factor.raw(), prec);
            SparseRow combined = plus_sparse(dep_rows[k], scaled, prec);
            if (sparse_chop_enabled()) chop_sparse(combined);
            reduced_rest.push_back(std::move(combined));
        }

        sys = std::move(indep_rows);
        for (auto& r : reduced_rest) sys.push_back(std::move(r));
        // PERF: an unconditional `chop_sparse(sys)` here is a no-op —
        // `indep_rows` are not touched this iteration (already chopped
        // on entry / by a prior iteration) and every row in
        // `reduced_rest` was chopped at construction above.  On the
        // the 3-loop crossed-box 44-master CalcTaylor block this re-scan dominated
        // forward_sparse_gaussian after the acb_is_chop_zero fix
        // (chop_sparse(sys&) = 39% of samples).  Dropped.

        result.reduced.push_back(std::move(pivot));
    }

    // MMA also re-sorts each rset row by column.
    for (auto& r : result.reduced) {
        std::sort(r.begin(), r.end(),
                  [](const SparseEntry& a, const SparseEntry& b) { return a.col < b.col; });
    }

    result.residue = std::move(sys);
    return result;
}

// ===========================================================================
//  SparseGaussian  (.m line 660-677)
// ===========================================================================

SparseGaussianResult sparse_gaussian(const ConstructedMatrix& cm,
                                     const std::vector<AcbValue>& nh,
                                     long prec) {
    if (cm.block_systems.empty()) {
        throw std::invalid_argument("sparse_gaussian: empty constructed matrix");
    }
    long N = cm.block_size;

    // Append nh values as an extra RHS column at index `total_columns`,
    // visiting block_index then in-block row.
    std::vector<SparseSystem> systems = cm.block_systems;   // deep copy

    long k = 0;
    for (std::size_t i = 0; i < systems.size(); ++i) {
        for (long r = 0; r < N; ++r) {
            if (k >= static_cast<long>(nh.size())) break;
            if (!nh[k].is_zero()) {
                systems[i][r].emplace_back(cm.total_columns, nh[k].clone());
            }
            ++k;
        }
    }

    SparseGaussianResult result;
    SparseSystem residue;          // moves between blocks

    for (std::size_t i = 0; i < systems.size(); ++i) {
        std::vector<long> var;
        if (i + 1 < systems.size()) {
            // Variables of block i in 0-based: N*i .. N*(i+1)-1.
            var.reserve(N);
            for (long t = 0; t < N; ++t) var.push_back(N * static_cast<long>(i) + t);
        } else {
            // Last block: include both block i and the trailing block i+1.
            var.reserve(2 * N);
            for (long t = 0; t < N; ++t) var.push_back(N * static_cast<long>(i) + t);
            for (long t = 0; t < N; ++t) var.push_back(N * (static_cast<long>(i) + 1) + t);
        }

        SparseSystem combined = std::move(systems[i]);
        for (auto& r : residue) combined.push_back(std::move(r));

        auto step = forward_sparse_gaussian(std::move(combined), var, prec);
        for (auto& r : step.reduced) result.reduce.push_back(std::move(r));
        residue = std::move(step.residue);
    }

    result.residue = std::move(residue);
    return result;
}

// ===========================================================================
//  Chop / pretty-print
// ===========================================================================

void chop_sparse(SparseRow& row) {
    int chop = sparse_chop_digits();
    // calcx00 audit, AMFLOW_X00_PIVOT_GUARD=1: also chop
    // entries whose ball CONTAINS ZERO (cancellation remnants [~0 +/- r]
    // with r above the chop floor).  acb_is_chop_zero tests the ball's
    // upper bound, so a wide zero-containing remnant survives it, later
    // gets picked as a pivot, and acb_inv turns it into an infinite ball
    // (the calcx00 "could not be rationalized" NaN source).  This is the
    // ball-arithmetic translation of the MMA original's Chop[] semantics
    // in DESolver's machine-precision elimination: a value with no
    // certifiable sign is treated as the cancelled zero it is.  PROBE-PATH
    // knob: enable only for fp-probe runs (downstream consumes midpoints;
    // the eps-fit validates the result).  Default OFF: byte-identical.
    static const bool contains_zero_chop = [] {
        const char* g = std::getenv("AMFLOW_X00_PIVOT_GUARD");
        return g && *g == '1';
    }();
    SparseRow filtered;
    filtered.reserve(row.size());
    for (auto& e : row) {
        if (contains_zero_chop && acb_contains_zero(e.val.raw())) {
            continue;   // no certifiable sign -> cancelled zero
        }
        if (!acb_is_chop_zero(e.val.raw(), chop)) {
            filtered.push_back(std::move(e));
        }
    }
    row = std::move(filtered);
}

void chop_sparse(SparseSystem& sys) {
    for (auto& r : sys) chop_sparse(r);
}

std::ostream& operator<<(std::ostream& os, const SparseRow& row) {
    os << "{";
    for (std::size_t i = 0; i < row.size(); ++i) {
        if (i != 0) os << ", ";
        os << row[i].col << " -> " << row[i].val.to_string(8);
    }
    os << "}";
    return os;
}

}  // namespace amflow::ode
