// SPDX-License-Identifier: MIT
// ode::asy — implementation.
//

#include "amflow/ode/asy.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <stdexcept>
#include <utility>
#include <vector>

#include <flint/acb.h>
#include <flint/arb.h>
#include <flint/fmpq.h>
#include <flint/fmpq_poly.h>
#include <flint/mag.h>

#include "amflow/ode/blocks.hpp"   // for amflow::ode::to_num

namespace amflow::ode {

using numeric::AcbValue;
using numeric::FmpqPoly;
using numeric::RationalFunction;
using numeric::RationalMatrix;
using numeric::acb_is_chop_zero;
using numeric::chop_pre;
using numeric::x_order;

// ===========================================================================
//  Internal helpers
// ===========================================================================

namespace {

// Allocate (rows x cols) matrix of zero AcbValue.
std::vector<std::vector<AcbValue>>
zeros_matrix(std::size_t rows, std::size_t cols) {
    std::vector<std::vector<AcbValue>> out(rows);
    for (auto& row : out) row.resize(cols);
    return out;
}

// Deep-copy a 2D matrix of AcbValue.
std::vector<std::vector<AcbValue>>
clone_matrix(const std::vector<std::vector<AcbValue>>& m) {
    std::vector<std::vector<AcbValue>> out(m.size());
    for (std::size_t i = 0; i < m.size(); ++i) {
        out[i].reserve(m[i].size());
        for (const auto& v : m[i]) out[i].push_back(v.clone());
    }
    return out;
}

}  // namespace

// ===========================================================================
//  ToPS
// ===========================================================================

PowerSeries to_power_series(const RationalFunction& r, long prec) {
    PowerSeries ps;
    if (r.is_zero()) return ps;   // empty coeffs marks ToPS[0]

    // exponent = Exponent[Numerator, eta, Min] - Exponent[Denominator, eta, Min]
    long vn = r.numerator().valuation();
    long vd = r.denominator().valuation();
    long offset = (vn < 0 ? 0 : vn) - (vd < 0 ? 0 : vd);

    // newpoly = r * eta^(-offset) ; this should land back as a polynomial.
    RationalFunction shifted = r;
    shifted.multiply_by_eta_power(-offset);

    if (!shifted.is_polynomial()) {
        throw std::invalid_argument(
            "to_power_series: input is not reducible to a polynomial after "
            "factoring out the leading eta-power (unexpected pole away from 0)");
    }

    ps.leading_offset = offset;
    ps.coeffs = to_num(shifted.numerator(), prec);
    return ps;
}

std::vector<std::vector<PowerSeries>>
to_power_series_matrix(const RationalMatrix& m, long prec) {
    std::vector<std::vector<PowerSeries>> out(m.rows());
    for (std::size_t i = 0; i < m.rows(); ++i) {
        out[i].reserve(m.cols());
        for (std::size_t j = 0; j < m.cols(); ++j) {
            out[i].push_back(to_power_series(m(i, j), prec));
        }
    }
    return out;
}

// ===========================================================================
//  Truncated power-series arithmetic
// ===========================================================================
//
//  Mathematica:
//    TimesPoly[f, poly] :=
//      Table[
//        If[r <= Length[poly],
//          f[[;; r]] . Reverse[poly[[;; r]]],
//          f[[r - Length[poly] + 1 ;; r]] . Reverse[poly]
//        ],
//        {r, Length[f]}]
//
//  For 0-based output index n, the result is
//      sum_{k = max(0, n - |poly| + 1)}^{n}  f[k] * poly[n - k]
//  -- truncated polynomial multiplication.

std::vector<AcbValue>
times_poly(const std::vector<AcbValue>& f,
           const std::vector<AcbValue>& poly,
           long prec) {
    long Nf = static_cast<long>(f.size());
    long Np = static_cast<long>(poly.size());
    std::vector<AcbValue> out(Nf);

    AcbValue prod;
    for (long n = 0; n < Nf; ++n) {
        AcbValue acc;     // zero
        long k_lo = std::max<long>(0, n - Np + 1);
        long k_hi = n;
        for (long k = k_lo; k <= k_hi; ++k) {
            acb_mul(prod.raw(), f[k].raw(), poly[n - k].raw(), prec);
            acb_add(acc.raw(), acc.raw(), prod.raw(), prec);
        }
        out[n] = std::move(acc);
    }
    return out;
}

// InversePoly[f, poly]: solve  poly * out = f  (truncated to length |f|).
//
//   coe[n] = ( f[n] - sum_{j} coe[n-j] * poly[j] ) / poly[0]
//
//   where j ranges over 1 .. min(n, |poly|-1).

std::vector<AcbValue>
inverse_poly(const std::vector<AcbValue>& f,
             const std::vector<AcbValue>& poly,
             long prec) {
    long Nf = static_cast<long>(f.size());
    long Np = static_cast<long>(poly.size());
    if (Np == 0 || poly[0].is_zero()) {
        throw std::domain_error("inverse_poly: leading polynomial coefficient is zero");
    }

    AcbValue inv0;
    acb_inv(inv0.raw(), poly[0].raw(), prec);

    std::vector<AcbValue> coe(Nf);
    AcbValue tmp;
    for (long n = 0; n < Nf; ++n) {
        AcbValue acc;     // zero
        long j_max = std::min<long>(n, Np - 1);
        for (long j = 1; j <= j_max; ++j) {
            acb_mul(tmp.raw(), coe[n - j].raw(), poly[j].raw(), prec);
            acb_add(acc.raw(), acc.raw(), tmp.raw(), prec);
        }
        AcbValue diff;
        acb_sub(diff.raw(), f[n].raw(), acc.raw(), prec);
        acb_mul(coe[n].raw(), diff.raw(), inv0.raw(), prec);
    }
    return coe;
}

std::vector<AcbValue>
rational_expansion(const std::vector<AcbValue>& num,
                   const std::vector<AcbValue>& den,
                   const std::vector<AcbValue>& expansion,
                   long prec) {
    auto times = times_poly(expansion, num, prec);
    return inverse_poly(times, den, prec);
}

std::vector<std::vector<AcbValue>>
map_rational_expansion(const std::vector<std::vector<std::vector<AcbValue>>>& nums,
                       const std::vector<std::vector<std::vector<AcbValue>>>& dens,
                       const std::vector<std::vector<AcbValue>>& expansions,
                       long prec) {
    if (nums.size() != dens.size()) {
        throw std::invalid_argument("map_rational_expansion: nums/dens size mismatch");
    }

    std::vector<std::vector<AcbValue>> out(nums.size());
    for (std::size_t i = 0; i < nums.size(); ++i) {
        const auto& num_row = nums[i];
        const auto& den_row = dens[i];
        if (num_row.size() != den_row.size() || num_row.size() != expansions.size()) {
            throw std::invalid_argument("map_rational_expansion: inner-row size mismatch");
        }

        // Determine the truncation length from the first nonempty expansion.
        std::size_t L = 0;
        for (const auto& e : expansions) {
            if (!e.empty()) { L = e.size(); break; }
        }
        if (L == 0) {
            // All expansions empty -> output is empty (matches Mathematica
            // returning 0 for Total[{}]).
            out[i].clear();
            continue;
        }

        std::vector<AcbValue> acc(L);    // zero
        for (std::size_t k = 0; k < num_row.size(); ++k) {
            if (den_row[k].empty()) continue;
            if (expansions[k].empty()) continue;
            auto term = rational_expansion(num_row[k], den_row[k], expansions[k], prec);
            if (term.size() != L) {
                throw std::invalid_argument(
                    "map_rational_expansion: inconsistent expansion lengths");
            }
            for (std::size_t n = 0; n < L; ++n) {
                acb_add(acc[n].raw(), acc[n].raw(), term[n].raw(), prec);
            }
        }
        out[i] = std::move(acc);
    }
    return out;
}

// ===========================================================================
//  Per-mu helpers
// ===========================================================================

void extend_expansion(std::vector<std::vector<AcbValue>>& exp, std::size_t n) {
    if (exp.size() >= n) return;
    long row_len = static_cast<long>(x_order()) + 1;
    while (exp.size() < n) {
        exp.emplace_back(row_len);    // zero-initialised AcbValues
    }
}

std::vector<std::vector<AcbValue>>
plus_expansion(const std::vector<std::vector<std::vector<AcbValue>>>& exps,
               long prec) {
    if (exps.empty()) return {};
    std::size_t maxrows = 0;
    for (const auto& e : exps) maxrows = std::max(maxrows, e.size());
    if (maxrows == 0) return {};

    long row_len = static_cast<long>(x_order()) + 1;
    auto out = zeros_matrix(maxrows, static_cast<std::size_t>(row_len));

    for (const auto& e : exps) {
        for (std::size_t k = 0; k < e.size(); ++k) {
            for (std::size_t n = 0; n < e[k].size() && static_cast<long>(n) < row_len; ++n) {
                acb_add(out[k][n].raw(), out[k][n].raw(), e[k][n].raw(), prec);
            }
        }
    }
    return out;
}

std::vector<std::vector<AcbValue>>
rescale_expansion(const std::vector<std::vector<AcbValue>>& exp,
                  std::size_t order) {
    long row_len = static_cast<long>(x_order()) + 1;
    std::vector<std::vector<AcbValue>> out(exp.size());
    for (std::size_t k = 0; k < exp.size(); ++k) {
        out[k].resize(static_cast<std::size_t>(row_len));     // zeros
        for (std::size_t n = 0; n + order < static_cast<std::size_t>(row_len)
                                && n < exp[k].size(); ++n) {
            out[k][n + order] = exp[k][n].clone();
        }
    }
    return out;
}

// ===========================================================================
//  Series-truncation tail bound  (bug U1)
// ===========================================================================

bool xorder_tail_bound_enabled() {
    static const bool enabled = []() -> bool {
        const char* v = std::getenv("AMFLOW_XORDER_TAIL_BOUND");
        return v != nullptr && *v != '\0' && *v != '0';
    }();
    return enabled;
}

void series_tail_mag(mag_t out,
                     const std::vector<AcbValue>& c,
                     acb_srcptr x) {
    mag_zero(out);
    const long N = static_cast<long>(c.size()) - 1;
    if (N < 2) return;

    mag_t mN, mNm1, mx, mxN, ratio, geo, last3;
    mag_init(mN);   mag_init(mNm1); mag_init(last3);
    mag_init(mx);   mag_init(mxN);  mag_init(ratio); mag_init(geo);

    // Work from MIDPOINTS only: the c[k] coefficients can carry an
    // accumulated tail radius from the previous step (propagated through
    // the calcx1x2 recursion), which would otherwise feed back and blow up.
    auto mid = [](const AcbValue& z) -> AcbValue {
        AcbValue m;
        arf_set(arb_midref(acb_realref(m.raw())), arb_midref(acb_realref(z.raw())));
        arf_set(arb_midref(acb_imagref(m.raw())), arb_midref(acb_imagref(z.raw())));
        return m;
    };
    AcbValue cN   = mid(c[N]);
    AcbValue cNm1 = mid(c[N - 1]);
    AcbValue cNm2 = mid(c[N - 2]);

    acb_get_mag(mN,   cN.raw());
    acb_get_mag(mNm1, cNm1.raw());
    acb_get_mag(mx, x);

    // Last-3-terms partial sum |c_{N-2} x^{N-2} + c_{N-1} x^{N-1} + c_N x^N|
    // computed as a *complex number* (Horner) so sign/phase cancellation is
    // respected.  The naive cmax*|x|^N bound over-estimates by many orders
    // of magnitude on the cancellation-heavy regular-contour Taylor series
    // (kite@w=5: ~e+165 on an O(e9) value); the partial-sum probe tracks the
    // observed ~32-digit bias.
    {
        const long prec = numeric::working_prec_bits();
        AcbValue p;        // = c_N
        acb_set(p.raw(), cN.raw());
        acb_mul(p.raw(), p.raw(), x, prec);
        acb_add(p.raw(), p.raw(), cNm1.raw(), prec);
        acb_mul(p.raw(), p.raw(), x, prec);
        acb_add(p.raw(), p.raw(), cNm2.raw(), prec);
        // p = c_{N-2} + c_{N-1} x + c_N x^2;  scale by x^{N-2}
        acb_get_mag(last3, p.raw());
        mag_pow_ui(mxN, mx, static_cast<unsigned long>(N - 2));
        mag_mul(last3, last3, mxN);
    }
    if (mag_is_zero(last3)) {
        mag_clear(mN); mag_clear(mNm1); mag_clear(last3);
        mag_clear(mx); mag_clear(mxN);  mag_clear(ratio); mag_clear(geo);
        return;
    }

    // r = |c_N|/|c_{N-1}| * |x|   (fallback r = |x| if c_{N-1} == 0)
    if (mag_is_zero(mNm1)) {
        mag_set(ratio, mx);
    } else {
        mag_div(ratio, mN, mNm1);
        mag_mul(ratio, ratio, mx);
    }
    // Clamp ratio at 0.9 if >= 1 (so 1/(1-r) stays finite; the geometric
    // factor is heuristic anyway -- the exact constant is irrelevant once
    // the honest radius is within a factor ~10 of the true bias).
    double r_d = mag_get_d(ratio);
    if (!(r_d < 1.0)) r_d = 0.9;       // also catches NaN/Inf
    mag_set_d(geo, 1.0 / (1.0 - r_d));

    // out = |last-3-terms partial| * 1/(1-r)
    mag_mul(out, last3, geo);

    mag_clear(mN); mag_clear(mNm1); mag_clear(last3);
    mag_clear(mx); mag_clear(mxN);  mag_clear(ratio); mag_clear(geo);
}

// ===========================================================================
//  Evaluation
// ===========================================================================
//
//  EvaluateExpansion: Horner along eta first, then a direct sum over k with
//  precomputed log^k(x0).

AcbValue evaluate_expansion(const std::vector<std::vector<AcbValue>>& exp,
                            acb_srcptr x0,
                            long prec) {
    if (exp.empty()) return AcbValue();

    std::size_t Krows = exp.size();
    std::size_t Ncols = exp[0].size();

    // accum[k] = sum_n exp[k][n] * x0^n   (Horner on x0)
    std::vector<AcbValue> accum(Krows);
    if (Ncols == 0) return AcbValue();
    const bool tail_bound = xorder_tail_bound_enabled();
    for (std::size_t k = 0; k < Krows; ++k) {
        AcbValue val;            // zero
        for (long n = static_cast<long>(Ncols) - 1; n >= 0; --n) {
            // val = val * x0 + exp[k][n]
            acb_mul(val.raw(), val.raw(), x0, prec);
            acb_add(val.raw(), val.raw(), exp[k][n].raw(), prec);
        }
        // [bug U1] eta-series truncation residue -> Arb radius.
        // Without this the InfToRegular matching value claims working_pre
        // accuracy even when x_order is too small for |1/x0|.
        if (tail_bound) {
            mag_t tail;
            mag_init(tail);
            series_tail_mag(tail, exp[k], x0);
            acb_add_error_mag(val.raw(), tail);
            mag_clear(tail);
        }
        accum[k] = std::move(val);
    }

    if (Krows == 1) return std::move(accum[0]);

    // Sum over k:  result = sum_k accum[k] * log^k(x0)
    AcbValue logx;
    acb_log(logx.raw(), x0, prec);

    AcbValue result;     // zero
    AcbValue logpow;
    logpow.set_one();
    AcbValue tmp;
    for (std::size_t k = 0; k < Krows; ++k) {
        if (k > 0) acb_mul(logpow.raw(), logpow.raw(), logx.raw(), prec);
        acb_mul(tmp.raw(), accum[k].raw(), logpow.raw(), prec);
        acb_add(result.raw(), result.raw(), tmp.raw(), prec);
    }
    return result;
}

AcbValue evaluate_asy_term(acb_srcptr mu,
                           const std::vector<std::vector<AcbValue>>& exp,
                           acb_srcptr x0,
                           long prec) {
    AcbValue inner = evaluate_expansion(exp, x0, prec);
    if (inner.is_zero()) return inner;

    // Result = x0^mu * inner
    AcbValue x0_mu;
    acb_pow(x0_mu.raw(), x0, mu, prec);

    AcbValue out;
    acb_mul(out.raw(), x0_mu.raw(), inner.raw(), prec);
    return out;
}

// ===========================================================================
//  AsyTerm helpers
// ===========================================================================

AsyTerm AsyTerm::clone() const {
    return AsyTerm(mu.clone(), clone_matrix(exp));
}

AsyExpansion clone_asy(const AsyExpansion& asy) {
    AsyExpansion out;
    out.reserve(asy.size());
    for (const auto& t : asy) out.push_back(t.clone());
    return out;
}

std::vector<AsyExpansion>
clone_asy_list(const std::vector<AsyExpansion>& list) {
    std::vector<AsyExpansion> out;
    out.reserve(list.size());
    for (const auto& a : list) out.push_back(clone_asy(a));
    return out;
}

AcbValue evaluate_asy_expansion(const AsyExpansion& asy,
                                acb_srcptr           x0,
                                long                  prec) {
    AcbValue total;     // zero
    for (const auto& t : asy) {
        AcbValue v = evaluate_asy_term(t.mu.raw(), t.exp, x0, prec);
        acb_add(total.raw(), total.raw(), v.raw(), prec);
    }
    return total;
}

// ===========================================================================
//  Rule-set algebra
// ===========================================================================
//
//  Two mu values "differ by integer" if (mu1 - mu2) is real and within
//  10^-chop_pre of an integer.  We use that to group asymptotic terms into
//  regions.

bool mu_differ_by_integer(acb_srcptr mu1, acb_srcptr mu2, long prec) {
    AcbValue diff;
    acb_sub(diff.raw(), mu1, mu2, prec);

    // Use the RESONANCE tolerance, not the coarse chop_pre: deciding that two
    // indicial exponents "differ by an integer" groups them into one Frobenius
    // region/block and is exactly the eps=0 resonance merge.  An eps-deformed
    // pair {mu, mu + a*eps} differs by a*eps, which is NOT an integer for
    // eps != 0; with the coarse chop_pre (~10^-chop_pre) any |a*eps| below that
    // threshold is rounded to the integer 0 and the two sectors collapse,
    // silently discarding the genuine 1/(a*eps) pole and capping the eps->0
    // extrapolation (the calcx00 eps=0 resonance ceiling).  The resonance
    // tolerance (~0.8*working_pre) sits at the arithmetic noise floor, so a
    // genuine integer separation (residual identically 0) is still caught while
    // an eps-scale separation stays distinct.  (See numeric::resonance_chop_pre.)
    int chop = numeric::resonance_chop_pre();
    if (chop <= 0) chop = 20;
    const double eps_double = std::pow(10.0, -chop);

    // Imag part of (mu1 - mu2) should be (numerically) zero.
    arf_t mid_im;
    arf_init(mid_im);
    arf_abs(mid_im, arb_midref(acb_imagref(diff.raw())));
    bool im_small = (arf_cmp_d(mid_im, eps_double) < 0);
    arf_clear(mid_im);
    if (!im_small) return false;

    // Real part: |re - round(re)| should be (numerically) zero.
    arb_t re;
    arb_init(re);
    arb_set(re, acb_realref(diff.raw()));

    fmpz_t nearest;
    fmpz_init(nearest);
    arf_get_fmpz(nearest, arb_midref(re), ARF_RND_NEAR);

    arb_t residual;
    arb_init(residual);
    arb_set_fmpz(residual, nearest);
    arb_sub(residual, re, residual, prec);

    arf_t mid_res;
    arf_init(mid_res);
    arf_abs(mid_res, arb_midref(residual));
    bool re_int = (arf_cmp_d(mid_res, eps_double) < 0);
    arf_clear(mid_res);

    arb_clear(residual);
    arb_clear(re);
    fmpz_clear(nearest);
    return re_int;
}

// PlusRuleS: combine all terms in a single region, returning a single AsyTerm
// whose mu equals the *minimum integer-shift representative* of the region.
//
// Mathematica behaviour (line 935):
//    min = Min[Keys[rules0] - Keys[rules0][[1]]] + Keys[rules0][[1]]
// i.e. compute minimum offset from the first key (in the real sense), then
// add the first key back.  This avoids comparing complex numbers wrongly.

AsyTerm plus_rule_set(const AsyExpansion& region, long prec) {
    if (region.empty()) {
        return AsyTerm(AcbValue(), {});
    }

    // base = mu[0]
    AcbValue base = region[0].mu.clone();

    // For each i, compute mu[i] - base.  All differences are *real integers*
    // by hypothesis; we store them as long ints.
    std::vector<long> offsets(region.size());
    for (std::size_t i = 0; i < region.size(); ++i) {
        AcbValue diff;
        acb_sub(diff.raw(), region[i].mu.raw(), base.raw(), prec);
        fmpz_t z;
        fmpz_init(z);
        arf_get_fmpz(z, arb_midref(acb_realref(diff.raw())), ARF_RND_NEAR);
        offsets[i] = fmpz_get_si(z);
        fmpz_clear(z);
    }

    long min_off = offsets[0];
    for (long o : offsets) if (o < min_off) min_off = o;

    // min_mu = base + min_off
    AcbValue min_mu;
    acb_set(min_mu.raw(), base.raw());
    AcbValue add;
    add.set_si(min_off);
    acb_add(min_mu.raw(), min_mu.raw(), add.raw(), prec);

    // For each term, rescale_expansion(term.exp, offsets[i] - min_off) and sum.
    std::vector<std::vector<std::vector<AcbValue>>> rescaled;
    rescaled.reserve(region.size());
    for (std::size_t i = 0; i < region.size(); ++i) {
        long shift = offsets[i] - min_off;
        rescaled.push_back(rescale_expansion(region[i].exp, static_cast<std::size_t>(shift)));
    }
    auto combined = plus_expansion(rescaled, prec);

    return AsyTerm(std::move(min_mu), std::move(combined));
}

AsyExpansion union_rule_set(const AsyExpansion& asy, long prec) {
    AsyExpansion out;

    // Group terms by integer-equivalence on mu.
    std::vector<bool> done(asy.size(), false);
    for (std::size_t i = 0; i < asy.size(); ++i) {
        if (done[i]) continue;
        AsyExpansion group;
        group.push_back(asy[i].clone());
        done[i] = true;
        for (std::size_t j = i + 1; j < asy.size(); ++j) {
            if (done[j]) continue;
            if (mu_differ_by_integer(asy[i].mu.raw(), asy[j].mu.raw(), prec)) {
                group.push_back(asy[j].clone());
                done[j] = true;
            }
        }
        out.push_back(plus_rule_set(group, prec));
    }
    return out;
}

// PSTimesRuleS:
//   For each (i in 1..|coeffs|, j in 1..|asy|), produce
//     (mu_j + leading_offset + (i - 1)) -> (exp_j scaled by coeffs[i-1])

AsyExpansion ps_times_rule_set(const PowerSeries& ps,
                               const AsyExpansion& asy,
                               long prec) {
    if (ps.is_zero()) return {};

    AsyExpansion out;
    out.reserve(ps.coeffs.size() * asy.size());

    for (std::size_t i = 0; i < ps.coeffs.size(); ++i) {
        const auto& c = ps.coeffs[i];
        if (c.is_zero()) continue;
        for (std::size_t j = 0; j < asy.size(); ++j) {
            // mu_new = mu_j + ps.leading_offset + i
            AcbValue mu_new = asy[j].mu.clone();
            AcbValue add;
            add.set_si(ps.leading_offset + static_cast<long>(i));
            acb_add(mu_new.raw(), mu_new.raw(), add.raw(), prec);

            // exp_new = c * exp_j   element-wise
            std::vector<std::vector<AcbValue>> exp_new(asy[j].exp.size());
            for (std::size_t k = 0; k < asy[j].exp.size(); ++k) {
                exp_new[k].resize(asy[j].exp[k].size());
                for (std::size_t n = 0; n < asy[j].exp[k].size(); ++n) {
                    acb_mul(exp_new[k][n].raw(), c.raw(),
                            asy[j].exp[k][n].raw(), prec);
                }
            }

            out.emplace_back(std::move(mu_new), std::move(exp_new));
        }
    }
    return out;
}

std::vector<AsyExpansion>
ps_map_rule_set(const std::vector<std::vector<PowerSeries>>& psmap,
                const std::vector<AsyExpansion>& asy_list,
                long prec) {
    std::size_t M = psmap.size();
    std::size_t K = asy_list.size();

    std::vector<AsyExpansion> out(M);
    for (std::size_t i = 0; i < M; ++i) {
        if (psmap[i].size() != K) {
            throw std::invalid_argument(
                "ps_map_rule_set: psmap row size does not match asy_list size");
        }
        AsyExpansion concat;
        for (std::size_t j = 0; j < K; ++j) {
            auto term = ps_times_rule_set(psmap[i][j], asy_list[j], prec);
            for (auto& t : term) concat.push_back(std::move(t));
        }
        out[i] = union_rule_set(concat, prec);
    }
    return out;
}

// ===========================================================================
//  TimesAsyExp / PlusAsyExp
// ===========================================================================
//
//  TimesAsyExp[rational, asyexp]:
//    de       = Denominator[Together[rational]]
//    exponent = Exponent[de, eta, Min]                  -- eta-order pulled from den
//    de       = ToNum[de / eta^exponent // Together]    -- de coefficients (constant term nonzero)
//    num      = ToNum[Numerator[Together[rational]]]
//    Each AsyTerm (mu, exp) -> (mu - exponent,
//                                RationalExpansion(num, de, exp[k]) for each k)

AsyExpansion times_asy_exp(const RationalFunction& rational,
                           const AsyExpansion& asy,
                           long prec) {
    if (rational.is_zero()) return {};

    const FmpqPoly& num_poly = rational.numerator();
    const FmpqPoly& den_poly = rational.denominator();
    long exponent = den_poly.valuation();    // >= 0

    // de_shifted = den / eta^exponent
    FmpqPoly de_shifted = den_poly.clone();
    if (exponent > 0) {
        FmpqPoly tmp;
        long L = de_shifted.length();
        fmpq_t c;
        fmpq_init(c);
        for (long k = exponent; k < L; ++k) {
            de_shifted.coeff(k, c);
            tmp.set_coeff_fmpq(k - exponent, c);
        }
        fmpq_clear(c);
        de_shifted = std::move(tmp);
    }
    auto de_num = to_num(de_shifted, prec);
    auto num_num = to_num(num_poly, prec);

    AsyExpansion out;
    out.reserve(asy.size());
    for (const auto& term : asy) {
        AcbValue new_mu = term.mu.clone();
        AcbValue sub;
        sub.set_si(exponent);
        acb_sub(new_mu.raw(), new_mu.raw(), sub.raw(), prec);

        std::vector<std::vector<AcbValue>> new_exp(term.exp.size());
        for (std::size_t k = 0; k < term.exp.size(); ++k) {
            if (term.exp[k].empty()) {
                new_exp[k].clear();
                continue;
            }
            new_exp[k] = rational_expansion(num_num, de_num, term.exp[k], prec);
        }
        out.emplace_back(std::move(new_mu), std::move(new_exp));
    }
    return out;
}

AsyExpansion plus_asy_exp(const std::vector<AsyExpansion>& list, long prec) {
    AsyExpansion concat;
    for (const auto& a : list) {
        for (const auto& t : a) concat.push_back(t.clone());
    }
    return union_rule_set(concat, prec);
}

// ===========================================================================
//  FindLogPower
// ===========================================================================

long find_log_power(const std::vector<std::vector<AcbValue>>& exp,
                    int chop_digits) {
    long k = static_cast<long>(exp.size());
    while (k >= 1) {
        bool all_zero = true;
        for (const auto& v : exp[k - 1]) {
            if (!acb_is_chop_zero(v.raw(), chop_digits)) { all_zero = false; break; }
        }
        if (all_zero) --k; else break;
    }
    return k - 1;
}

}  // namespace amflow::ode
