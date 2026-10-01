// SPDX-License-Identifier: MIT
// ode::regular — implementation.
//

#include "amflow/ode/regular.hpp"

#include <algorithm>
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <flint/acb.h>
#include <flint/arb.h>
#include <flint/fmpz.h>
#include <flint/mag.h>

#include "amflow/numeric/log.hpp"
#include "amflow/ode/asy.hpp"   // map_rational_expansion

namespace amflow::ode {

using numeric::AcbValue;
using numeric::log_line;
using numeric::silent_mode;
using numeric::working_prec_bits;
using numeric::x_order;

// ===========================================================================
//  Internal helpers
// ===========================================================================

namespace {

bool debug_regular_enabled() {
    return numeric::log::trace_enabled("AMFLOW_DEBUG_REGULAR");
}

AcbValue midpoint_only(const AcbValue& src) {
    AcbValue out;
    arb_set_arf(acb_realref(out.raw()), arb_midref(acb_realref(src.raw())));
    arb_set_arf(acb_imagref(out.raw()), arb_midref(acb_imagref(src.raw())));
    return out;
}

void dump_regular_vector(const char* label, const std::vector<AcbValue>& values) {
    if (!debug_regular_enabled()) return;
    std::cerr << "[regular_stage] " << label
              << " size=" << values.size() << std::endl;
    for (std::size_t i = 0; i < values.size(); ++i) {
        std::cerr << "  [" << i << "] = " << values[i].to_string(20) << std::endl;
    }
}

void dump_regular_coeffs_head(const TaylorCoefficients& coeffs,
                              std::size_t head = 5) {
    if (!debug_regular_enabled()) return;
    std::cerr << "[regular_stage] coeff_head size=" << coeffs.size() << std::endl;
    for (std::size_t i = 0; i < coeffs.size(); ++i) {
        std::cerr << "  coeff[" << i << "]";
        std::size_t upto = std::min(head, coeffs[i].size());
        for (std::size_t k = 0; k < upto; ++k) {
            std::cerr << " c" << k << "=" << coeffs[i][k].to_string(12);
        }
        if (!coeffs[i].empty()) {
            std::size_t last = coeffs[i].size() - 1;
            std::cerr << " c" << last << "=" << coeffs[i][last].to_string(12);
        }
        std::cerr << std::endl;
    }
}

void dump_regular_coeff_summary(const TaylorCoefficients& coeffs, long prec) {
    if (!debug_regular_enabled()) return;
    std::cerr << "[regular_stage] coeff_summary size=" << coeffs.size() << std::endl;
    for (std::size_t i = 0; i < coeffs.size(); ++i) {
        double best_mid = -1.0;
        double best_rad = -1.0;
        std::size_t best_mid_idx = 0;
        std::size_t best_rad_idx = 0;
        arb_t mag;
        arb_init(mag);
        for (std::size_t k = 0; k < coeffs[i].size(); ++k) {
            acb_abs(mag, coeffs[i][k].raw(), prec);
            double mid = arf_get_d(arb_midref(mag), ARF_RND_NEAR);
            double rad = std::max(
                mag_get_d(arb_radref(acb_realref(coeffs[i][k].raw()))),
                mag_get_d(arb_radref(acb_imagref(coeffs[i][k].raw()))));
            if (mid > best_mid) { best_mid = mid; best_mid_idx = k; }
            if (rad > best_rad) { best_rad = rad; best_rad_idx = k; }
        }
        arb_clear(mag);
        std::cerr << "  coeff[" << i << "]"
                  << " max_mid@" << best_mid_idx << "="
                  << coeffs[i][best_mid_idx].to_string(12)
                  << " max_rad@" << best_rad_idx << "="
                  << coeffs[i][best_rad_idx].to_string(12)
                  << std::endl;
    }
}

std::vector<AcbValue> clone_vec(const std::vector<AcbValue>& v) {
    std::vector<AcbValue> out;
    out.reserve(v.size());
    for (const auto& x : v) out.push_back(x.clone());
    return out;
}

std::vector<std::vector<AcbValue>>
clone_mat(const std::vector<std::vector<AcbValue>>& m) {
    std::vector<std::vector<AcbValue>> out;
    out.reserve(m.size());
    for (const auto& row : m) out.push_back(clone_vec(row));
    return out;
}

std::vector<std::vector<std::vector<AcbValue>>>
clone_3d(const std::vector<std::vector<std::vector<AcbValue>>>& m) {
    std::vector<std::vector<std::vector<AcbValue>>> out;
    out.reserve(m.size());
    for (const auto& mat : m) out.push_back(clone_mat(mat));
    return out;
}

std::size_t max_len_1d(const std::vector<AcbValue>& v) { return v.size(); }

std::size_t max_len_2d(const std::vector<std::vector<std::vector<AcbValue>>>& m) {
    std::size_t mx = 0;
    for (const auto& row : m)
        for (const auto& vec : row)
            mx = std::max(mx, vec.size());
    return mx;
}

// Apply the Taylor-shift table to a single polynomial.  After expansion the
// resulting polynomial is truncated to (XOrder + 1) coefficients (matches
// .m line 807).
//
//   table[i][j] = Binomial(i, j) * x0^(i - j),  0 <= j <= i <= maxorder
//   new_poly[j] = sum_{i = j}^{Length[poly] - 1} table[i][j] * poly[i]
//                 for j in 0..min(Length[poly] - 1, XOrder)

std::vector<AcbValue>
shift_poly(const std::vector<AcbValue>& poly,
           const std::vector<std::vector<AcbValue>>& table,
           long xorder,
           long prec) {
    if (poly.empty()) return {};
    long L = static_cast<long>(poly.size());
    long jmax = std::min<long>(L - 1, xorder);

    std::vector<AcbValue> out(jmax + 1);
    AcbValue tmp;
    for (long j = 0; j <= jmax; ++j) {
        AcbValue acc;     // zero
        for (long i = j; i < L; ++i) {
            acb_mul(tmp.raw(), table[i][j].raw(), poly[i].raw(), prec);
            acb_add(acc.raw(), acc.raw(), tmp.raw(), prec);
        }
        out[j] = std::move(acc);
    }
    return out;
}

}  // namespace

// ===========================================================================
//  expand_nh_equations_num
// ===========================================================================

std::vector<BlockEquationNum>
expand_nh_equations_num(const std::vector<BlockEquationNum>& nheqn, long /*prec*/) {
    std::vector<BlockEquationNum> out;
    out.reserve(nheqn.size());
    for (const auto& eq : nheqn) {
        BlockEquationNum c;
        c.dxexp  = clone_vec(eq.dxexp);
        c.axexp  = clone_3d(eq.axexp);
        c.bxexpn = clone_3d(eq.bxexpn);
        c.bxexpd = clone_3d(eq.bxexpd);
        c.block  = eq.block;
        c.sub    = eq.sub;
        out.push_back(std::move(c));
    }
    return out;
}

std::vector<BlockEquationNum>
expand_nh_equations_num(const std::vector<BlockEquationNum>& nheqn,
                        acb_srcptr x0,
                        long prec) {
    if (acb_is_zero(x0)) return expand_nh_equations_num(nheqn, prec);

    // 1) Determine the maximum polynomial length across all dx, ax, bxn, bxd.
    std::size_t max_len = 0;
    for (const auto& eq : nheqn) {
        max_len = std::max(max_len, max_len_1d(eq.dxexp));
        max_len = std::max(max_len, max_len_2d(eq.axexp));
        max_len = std::max(max_len, max_len_2d(eq.bxexpn));
        max_len = std::max(max_len, max_len_2d(eq.bxexpd));
    }
    long maxorder = (max_len == 0) ? 0 : static_cast<long>(max_len) - 1;
    long xorder = static_cast<long>(x_order());

    // 2) table[i][j] = C(i, j) * x0^(i - j), 0 <= j <= i <= maxorder.
    std::vector<std::vector<AcbValue>> table(maxorder + 1);
    for (long i = 0; i <= maxorder; ++i) table[i].resize(i + 1);

    std::vector<AcbValue> x0_pow(maxorder + 1);
    x0_pow[0].set_one();
    for (long k = 1; k <= maxorder; ++k) {
        acb_mul(x0_pow[k].raw(), x0_pow[k - 1].raw(), x0, prec);
    }

    fmpz_t binom;
    fmpz_init(binom);
    for (long i = 0; i <= maxorder; ++i) {
        for (long j = 0; j <= i; ++j) {
            fmpz_bin_uiui(binom, static_cast<unsigned long>(i),
                                  static_cast<unsigned long>(j));
            acb_set_fmpz(table[i][j].raw(), binom);
            if (i != j) {
                acb_mul(table[i][j].raw(), table[i][j].raw(),
                        x0_pow[i - j].raw(), prec);
            }
        }
    }
    fmpz_clear(binom);

    // 3) Apply shift to every poly in every block equation.
    std::vector<BlockEquationNum> out;
    out.reserve(nheqn.size());

    auto shift_mat = [&](const std::vector<std::vector<std::vector<AcbValue>>>& src) {
        std::vector<std::vector<std::vector<AcbValue>>> dst(src.size());
        for (std::size_t i = 0; i < src.size(); ++i) {
            dst[i].resize(src[i].size());
            for (std::size_t j = 0; j < src[i].size(); ++j) {
                dst[i][j] = shift_poly(src[i][j], table, xorder, prec);
            }
        }
        return dst;
    };

    for (const auto& eq : nheqn) {
        BlockEquationNum c;
        c.dxexp  = shift_poly(eq.dxexp, table, xorder, prec);
        c.axexp  = shift_mat(eq.axexp);
        c.bxexpn = shift_mat(eq.bxexpn);
        c.bxexpd = shift_mat(eq.bxexpd);
        c.block  = eq.block;
        c.sub    = eq.sub;
        out.push_back(std::move(c));
    }
    return out;
}

// ===========================================================================
//  calcx1x2
// ===========================================================================
//
//  Reference (.m line 816-837):
//    Initialise f0[i, 0] = bc[i].
//    For each block k:
//      maxd = Length[dxexp] - 1
//      maxa = max length of axexp inner - 1   (or -1 if all empty)
//      mata[m] = PickMat[axexp, m]    (constant matrix per m)
//      nh = MapRationalExpansion[bxexpn, bxexpd, f0[sub, 0..XOrder]]
//      For n in 0..XOrder-1:
//          vec = 1/((n+1) * dxexp[0]) *
//                  ( nh[All, n+1]
//                  + sum_{m=0}^{min(n, maxa)} mata[m] . f0[block, n-m]
//                  + sum_{m=1}^{min(n, maxd)} -(n-m+1) * dxexp[m] * f0[block, n-m+1] )
//          f0[block[i], n+1] = vec[i]

TaylorCoefficients
calcx1x2(const std::vector<BlockEquationNum>& nheqn,
         const std::vector<AcbValue>& bc,
         acb_srcptr x0,
         long prec) {
    auto nheqnx0 = expand_nh_equations_num(nheqn, x0, prec);

    long xorder = static_cast<long>(x_order());
    std::size_t n_int = bc.size();

    TaylorCoefficients f0(n_int);
    for (auto& row : f0) row.resize(xorder + 1);
    for (std::size_t i = 0; i < n_int; ++i) {
        f0[i][0] = bc[i].clone();
    }

    AcbValue tmp, tmp2;
    for (const auto& eq : nheqnx0) {
        const auto& dxexp = eq.dxexp;
        const auto& axexp = eq.axexp;
        const auto& bxexpn = eq.bxexpn;
        const auto& bxexpd = eq.bxexpd;
        const auto& block = eq.block;
        const auto& sub   = eq.sub;
        const std::size_t Nblk = block.size();

        if (dxexp.empty() || dxexp[0].is_zero()) {
            throw std::runtime_error(
                "calcx1x2: dx[0] is zero at the requested point -- the matrix "
                "is not regular here.  Did you intend to call calc_zero instead?");
        }

        long maxd = static_cast<long>(dxexp.size()) - 1;
        long maxa = -1;
        for (const auto& row : axexp) {
            for (const auto& vec : row) {
                long L = static_cast<long>(vec.size()) - 1;
                if (L > maxa) maxa = L;
            }
        }

        // Precompute mata[m] for m in 0..maxa.
        std::vector<std::vector<std::vector<AcbValue>>> mata(maxa < 0 ? 0 : (maxa + 1));
        for (long m = 0; m <= maxa; ++m) {
            mata[m].resize(Nblk);
            for (std::size_t i = 0; i < Nblk; ++i) {
                mata[m][i].resize(Nblk);
                for (std::size_t j = 0; j < Nblk; ++j) {
                    if (m < static_cast<long>(axexp[i][j].size())
                        && !axexp[i][j][m].is_zero()) {
                        mata[m][i][j] = axexp[i][j][m].clone();
                    }
                }
            }
        }

        // nh = MapRationalExpansion(bxn, bxd, f0[sub, *])
        std::vector<std::vector<AcbValue>> sub_expansions(sub.size());
        for (std::size_t k = 0; k < sub.size(); ++k) {
            sub_expansions[k] = clone_vec(f0[sub[k]]);
        }
        std::vector<std::vector<AcbValue>> nh =
            map_rational_expansion(bxexpn, bxexpd, sub_expansions, prec);

        for (std::size_t i = 0; i < nh.size(); ++i) {
            if (nh[i].empty()) {
                nh[i].resize(xorder + 1);
            }
        }

        AcbValue inv_factor;

        for (long n = 0; n < xorder; ++n) {
            AcbValue factor;
            acb_mul_si(factor.raw(), dxexp[0].raw(), n + 1, prec);
            acb_inv(inv_factor.raw(), factor.raw(), prec);

            std::vector<AcbValue> vec(Nblk);
            for (std::size_t i = 0; i < Nblk; ++i) {
                AcbValue acc;     // zero

                // (1) nh[i][n]
                if (n < static_cast<long>(nh[i].size())) {
                    acb_add(acc.raw(), acc.raw(), nh[i][n].raw(), prec);
                }

                // (2) sum_{m=0}^{min(n, maxa)} (mata[m] . f0[block, n-m])[i]
                long m_top_a = std::min<long>(n, maxa);
                for (long m = 0; m <= m_top_a; ++m) {
                    AcbValue dot;     // zero
                    for (std::size_t j = 0; j < Nblk; ++j) {
                        if (mata[m][i][j].is_zero()) continue;
                        acb_mul(tmp.raw(), mata[m][i][j].raw(),
                                f0[block[j]][n - m].raw(), prec);
                        acb_add(dot.raw(), dot.raw(), tmp.raw(), prec);
                    }
                    acb_add(acc.raw(), acc.raw(), dot.raw(), prec);
                }

                // (3) sum_{m=1}^{min(n, maxd)} -(n - m + 1) * dxexp[m] * f0[block[i], n - m + 1]
                long m_top_d = std::min<long>(n, maxd);
                for (long m = 1; m <= m_top_d; ++m) {
                    if (dxexp[m].is_zero()) continue;
                    acb_mul_si(tmp.raw(), dxexp[m].raw(), n - m + 1, prec);
                    acb_mul(tmp2.raw(), tmp.raw(),
                            f0[block[i]][n - m + 1].raw(), prec);
                    acb_sub(acc.raw(), acc.raw(), tmp2.raw(), prec);
                }

                acb_mul(vec[i].raw(), inv_factor.raw(), acc.raw(), prec);
            }

            for (std::size_t i = 0; i < Nblk; ++i) {
                f0[block[i]][n + 1] = std::move(vec[i]);
            }
        }
    }

    return f0;
}

// ===========================================================================
//  evaluate_taylor
// ===========================================================================

std::vector<AcbValue>
evaluate_taylor(const TaylorCoefficients& coeffs, acb_srcptr dh, long prec) {
    std::vector<AcbValue> out(coeffs.size());
    AcbValue dh_mid;
    arb_set_arf(acb_realref(dh_mid.raw()), arb_midref(acb_realref(dh)));
    arb_set_arf(acb_imagref(dh_mid.raw()), arb_midref(acb_imagref(dh)));

    const bool tail_bound = xorder_tail_bound_enabled();

    for (std::size_t i = 0; i < coeffs.size(); ++i) {
        const auto& c = coeffs[i];
        AcbValue val;     // zero
        // Midpoint-only Horner: mirrors Mathematica's point arithmetic and
        // avoids catastrophic interval blow-up on cancellation-heavy regular
        // contours.
        // CAVEAT: midpoint-only zeroes the Arb radius, so the eta-series
        // truncation residue is silently dropped.  With
        // AMFLOW_XORDER_TAIL_BOUND=1 we re-inject a geometric tail estimate
        // *after* the Horner sum so the propagated ball's radius honestly
        // reflects the x_order cut.  See series_tail_mag() / asy.hpp.
        for (long k = static_cast<long>(c.size()) - 1; k >= 0; --k) {
            acb_mul(val.raw(), val.raw(), dh_mid.raw(), prec);
            val = midpoint_only(val);
            AcbValue coeff_mid = midpoint_only(c[k]);
            acb_add(val.raw(), val.raw(), coeff_mid.raw(), prec);
            val = midpoint_only(val);
        }
        if (tail_bound) {
            // c[0] is the previous step's bc value; its radius already
            // carries the accumulated tail from earlier inf_to_regular /
            // calc_run steps (calcx1x2 propagates radii through the
            // recursion).  Preserve it across the midpoint-only Horner by
            // re-injecting max(prev_rad, this_step_tail).
            mag_t tail, prev_rad;
            mag_init(tail);
            mag_init(prev_rad);
            if (!c.empty()) {
                mag_max(prev_rad, arb_radref(acb_realref(c[0].raw())),
                                  arb_radref(acb_imagref(c[0].raw())));
            }
            series_tail_mag(tail, c, dh);
            mag_max(tail, tail, prev_rad);
            acb_add_error_mag(val.raw(), tail);
            mag_clear(prev_rad);
            mag_clear(tail);
        }
        out[i] = std::move(val);
    }
    return out;
}

// ===========================================================================
//  calc_run
// ===========================================================================

std::vector<AcbValue>
calc_run(const std::vector<BlockEquationNum>& nheqn,
         const std::vector<AcbValue>& bc,
         const std::vector<AcbValue>& run,
         long prec) {
    if (run.empty()) return {};
    if (run.size() == 1) {
        std::vector<AcbValue> out;
        out.reserve(bc.size());
        for (const auto& v : bc) out.push_back(v.clone());
        return out;
    }

    std::vector<AcbValue> bcr;
    bcr.reserve(bc.size());
    for (const auto& v : bc) bcr.push_back(v.clone());

    AcbValue dh;
    for (std::size_t i = 0; i + 1 < run.size(); ++i) {
        if (!silent_mode()) {
            log_line(std::string("CalcRun: current regular point -> ")
                     + run[i].to_string(20));
        }
        if (debug_regular_enabled()) {
            std::cerr << "[regular_stage] step=" << i
                      << " from " << run[i].to_string(20)
                      << " to " << run[i + 1].to_string(20) << std::endl;
            dump_regular_vector("before_step", bcr);
        }

        TaylorCoefficients rule = calcx1x2(nheqn, bcr, run[i].raw(), prec);
        dump_regular_coeffs_head(rule);
        dump_regular_coeff_summary(rule, prec);
        acb_sub(dh.raw(), run[i + 1].raw(), run[i].raw(), prec);
        bcr = evaluate_taylor(rule, dh.raw(), prec);
        dump_regular_vector("after_step", bcr);
    }
    return bcr;
}

std::vector<AcbValue>
calc_run(const numeric::RationalMatrix& de,
         const std::vector<AcbValue>& bc,
         const std::vector<AcbValue>& run,
         long prec) {
    auto nheq  = nh_equations(de, EquationMode::Regular);
    auto nheqn = nh_equations_num(nheq, prec);
    return calc_run(nheqn, bc, run, prec);
}

// ===========================================================================
//  calcx1x2_multi / calc_run_multi  (fp-probe column batching)
// ===========================================================================
//
//  Shares the column-independent work of calcx1x2 across all columns:
//    - expand_nh_equations_num(nheqn, x0)   (the per-step re-expansion)
//    - mata[m] pick-matrix assembly per block
//  The per-column recursion body is the byte-for-byte operation sequence of
//  calcx1x2, so each returned TaylorCoefficients is bit-identical to a
//  single-column calcx1x2 call on the same bc.

std::vector<TaylorCoefficients>
calcx1x2_multi(const std::vector<BlockEquationNum>& nheqn,
               const std::vector<std::vector<AcbValue>>& bc_cols,
               acb_srcptr x0,
               long prec) {
    auto nheqnx0 = expand_nh_equations_num(nheqn, x0, prec);

    long xorder = static_cast<long>(x_order());
    const std::size_t ncols = bc_cols.size();
    if (ncols == 0) return {};
    std::size_t n_int = bc_cols[0].size();

    std::vector<TaylorCoefficients> f0s(ncols);
    for (std::size_t c = 0; c < ncols; ++c) {
        if (bc_cols[c].size() != n_int) {
            throw std::invalid_argument(
                "calcx1x2_multi: ragged bc column sizes");
        }
        f0s[c].resize(n_int);
        for (auto& row : f0s[c]) row.resize(xorder + 1);
        for (std::size_t i = 0; i < n_int; ++i) {
            f0s[c][i][0] = bc_cols[c][i].clone();
        }
    }

    AcbValue tmp, tmp2;
    for (const auto& eq : nheqnx0) {
        const auto& dxexp = eq.dxexp;
        const auto& axexp = eq.axexp;
        const auto& bxexpn = eq.bxexpn;
        const auto& bxexpd = eq.bxexpd;
        const auto& block = eq.block;
        const auto& sub   = eq.sub;
        const std::size_t Nblk = block.size();

        if (dxexp.empty() || dxexp[0].is_zero()) {
            throw std::runtime_error(
                "calcx1x2_multi: dx[0] is zero at the requested point -- the "
                "matrix is not regular here.");
        }

        long maxd = static_cast<long>(dxexp.size()) - 1;
        long maxa = -1;
        for (const auto& row : axexp) {
            for (const auto& vec : row) {
                long L = static_cast<long>(vec.size()) - 1;
                if (L > maxa) maxa = L;
            }
        }

        // Column-independent: precompute mata[m] ONCE for all columns.
        std::vector<std::vector<std::vector<AcbValue>>> mata(maxa < 0 ? 0 : (maxa + 1));
        for (long m = 0; m <= maxa; ++m) {
            mata[m].resize(Nblk);
            for (std::size_t i = 0; i < Nblk; ++i) {
                mata[m][i].resize(Nblk);
                for (std::size_t j = 0; j < Nblk; ++j) {
                    if (m < static_cast<long>(axexp[i][j].size())
                        && !axexp[i][j][m].is_zero()) {
                        mata[m][i][j] = axexp[i][j][m].clone();
                    }
                }
            }
        }

        for (std::size_t c = 0; c < ncols; ++c) {
            TaylorCoefficients& f0 = f0s[c];

            // nh = MapRationalExpansion(bxn, bxd, f0[sub, *])
            std::vector<std::vector<AcbValue>> sub_expansions(sub.size());
            for (std::size_t k = 0; k < sub.size(); ++k) {
                sub_expansions[k] = clone_vec(f0[sub[k]]);
            }
            std::vector<std::vector<AcbValue>> nh =
                map_rational_expansion(bxexpn, bxexpd, sub_expansions, prec);

            for (std::size_t i = 0; i < nh.size(); ++i) {
                if (nh[i].empty()) {
                    nh[i].resize(xorder + 1);
                }
            }

            AcbValue inv_factor;

            for (long n = 0; n < xorder; ++n) {
                AcbValue factor;
                acb_mul_si(factor.raw(), dxexp[0].raw(), n + 1, prec);
                acb_inv(inv_factor.raw(), factor.raw(), prec);

                std::vector<AcbValue> vec(Nblk);
                for (std::size_t i = 0; i < Nblk; ++i) {
                    AcbValue acc;     // zero

                    // (1) nh[i][n]
                    if (n < static_cast<long>(nh[i].size())) {
                        acb_add(acc.raw(), acc.raw(), nh[i][n].raw(), prec);
                    }

                    // (2) sum_{m=0}^{min(n, maxa)} (mata[m] . f0[block, n-m])[i]
                    long m_top_a = std::min<long>(n, maxa);
                    for (long m = 0; m <= m_top_a; ++m) {
                        AcbValue dot;     // zero
                        for (std::size_t j = 0; j < Nblk; ++j) {
                            if (mata[m][i][j].is_zero()) continue;
                            acb_mul(tmp.raw(), mata[m][i][j].raw(),
                                    f0[block[j]][n - m].raw(), prec);
                            acb_add(dot.raw(), dot.raw(), tmp.raw(), prec);
                        }
                        acb_add(acc.raw(), acc.raw(), dot.raw(), prec);
                    }

                    // (3) sum_{m=1}^{min(n, maxd)} -(n-m+1)*dxexp[m]*f0[block[i], n-m+1]
                    long m_top_d = std::min<long>(n, maxd);
                    for (long m = 1; m <= m_top_d; ++m) {
                        if (dxexp[m].is_zero()) continue;
                        acb_mul_si(tmp.raw(), dxexp[m].raw(), n - m + 1, prec);
                        acb_mul(tmp2.raw(), tmp.raw(),
                                f0[block[i]][n - m + 1].raw(), prec);
                        acb_sub(acc.raw(), acc.raw(), tmp2.raw(), prec);
                    }

                    acb_mul(vec[i].raw(), inv_factor.raw(), acc.raw(), prec);
                }

                for (std::size_t i = 0; i < Nblk; ++i) {
                    f0[block[i]][n + 1] = std::move(vec[i]);
                }
            }
        }
    }

    return f0s;
}

std::vector<std::vector<AcbValue>>
calc_run_multi(const std::vector<BlockEquationNum>& nheqn,
               const std::vector<std::vector<AcbValue>>& bc_cols,
               const std::vector<AcbValue>& run,
               long prec) {
    const std::size_t ncols = bc_cols.size();
    if (run.empty()) return std::vector<std::vector<AcbValue>>(ncols);

    std::vector<std::vector<AcbValue>> bcr(ncols);
    for (std::size_t c = 0; c < ncols; ++c) {
        bcr[c].reserve(bc_cols[c].size());
        for (const auto& v : bc_cols[c]) bcr[c].push_back(v.clone());
    }
    if (run.size() == 1) return bcr;

    AcbValue dh;
    for (std::size_t i = 0; i + 1 < run.size(); ++i) {
        if (!silent_mode()) {
            log_line(std::string("CalcRun: current regular point -> ")
                     + run[i].to_string(20)
                     + " [batched " + std::to_string(ncols) + " cols]");
        }
        auto rules = calcx1x2_multi(nheqn, bcr, run[i].raw(), prec);
        acb_sub(dh.raw(), run[i + 1].raw(), run[i].raw(), prec);
        for (std::size_t c = 0; c < ncols; ++c) {
            bcr[c] = evaluate_taylor(rules[c], dh.raw(), prec);
        }
        if (std::getenv("AMF_DBG_X00")) {
            // calcx00 audit: localize the FIRST contour
            // point after which a transported value goes non-finite.
            long bad = 0;
            for (std::size_t c = 0; c < ncols && bad < 4; ++c) {
                for (std::size_t m = 0; m < bcr[c].size(); ++m) {
                    if (!acb_is_finite(bcr[c][m].raw())) {
                        std::fprintf(stderr,
                                "[X00-run] NON-FINITE after step %zu -> %zu: "
                                "col=%zu master=%zu\n",
                                i, i + 1, c, m);
                        ++bad;
                        break;
                    }
                }
            }
        }
    }
    return bcr;
}

}  // namespace amflow::ode
