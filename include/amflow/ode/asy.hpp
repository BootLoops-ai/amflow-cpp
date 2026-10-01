// SPDX-License-Identifier: MIT
// ode::asy — asymptotic-expansion algebra (Layer 4).
//
//
// Mirrors MMA DESolver.m lines 161-184 and 920-979:
//
//   EvaluateExpansion / EvaluateAsymptoticExpansion
//   TimesPoly / InversePoly / RationalExpansion / MapRationalExpansion
//   ExtendExpansion / PlusExpansion / RescaleExpansion
//   PlusRuleS / UnionRuleS / PSTimesRuleS / PSMapRuleS
//   ToPS / TimesAsyExp / PlusAsyExp
//
// What is being represented
// -------------------------
//
//   In a neighbourhood of eta = 0 every integral has a Frobenius-style
//   expansion:
//
//        I(eta) = sum_mu  eta^mu * sum_k log^k(eta) * sum_n c_{mu,k,n} eta^n
//
//   * `mu`    : exponent, a complex number (epsilon already numeric).
//   * `k`     : log power, 0 <= k <= K_mu.
//   * `n`     : integer eta-offset, 0 <= n <= XOrder.
//
//   In Mathematica this is encoded as a list of rules
//      { mu1 -> matrix1, mu2 -> matrix2, ... }
//   where matrix is 2D: rows = log power (k), columns = eta order (n).
//
//   We mirror that with `AsyTerm { mu, exp[k][n] }` and
//   `AsyExpansion = vector<AsyTerm>`.
//
// Power-series convention
// -----------------------
//
//   "PowerSeries" in this layer is `(coefficients, leading_offset)`:
//      f(eta) = eta^leading_offset * sum_n coefficients[n] eta^n
//   This matches Mathematica's ToPS output `{coef_list, exp}` exactly.
//
// Truncation conventions
// ----------------------
//
//   * `times_poly`/`inverse_poly` return a list with the same length as their
//     truncated input `f`.
//   * Asymptotic expansions are truncated to (XOrder + 1) eta orders by
//     convention (`exp[k]` has length XOrder + 1).  Some helpers grow the
//     log-power dimension dynamically; the eta-order dimension is global.

#ifndef AMFLOW_ODE_ASY_HPP
#define AMFLOW_ODE_ASY_HPP

#include <cstddef>
#include <iosfwd>
#include <utility>
#include <vector>

#include "amflow/numeric/acb_wrap.hpp"
#include "amflow/numeric/options.hpp"
#include "amflow/numeric/rational.hpp"

namespace amflow::ode {

// ---------------------------------------------------------------------------
//  PowerSeries  (output of ToPS)
// ---------------------------------------------------------------------------
//
//  Represents     eta^leading_offset * (coeffs[0] + coeffs[1] eta + ...).
//
//  An empty PowerSeries (coeffs.empty()) is used as the canonical "zero"
//  output of ToPS[0], distinct from a series that happens to start with a
//  zero coefficient.

struct PowerSeries {
    std::vector<numeric::AcbValue> coeffs;        // ordinary polynomial coefficients
    long                           leading_offset = 0;

    bool is_zero() const { return coeffs.empty(); }
};

// Convert a rational function r(eta) into a PowerSeries by factoring out the
// common eta-power.  Throws std::invalid_argument if the result is not a
// polynomial after factoring (i.e. the function still has poles at eta = 0
// or elsewhere -- ToPS in DESolver.m assumes the input has only an eta=0
// singularity, with the matrix T from NormalizeMat being the typical input).
PowerSeries to_power_series(const numeric::RationalFunction& r,
                            long prec = numeric::working_prec_bits());

// Apply ToPS element-wise to a matrix.  The result is a 2D PowerSeries grid.
std::vector<std::vector<PowerSeries>>
to_power_series_matrix(const numeric::RationalMatrix& m,
                       long prec = numeric::working_prec_bits());

// ---------------------------------------------------------------------------
//  Truncated power-series arithmetic
// ---------------------------------------------------------------------------
//
//   times_poly(f, p) :
//      Multiply the truncated series `f` by the polynomial `p`.  The result
//      has the same length as `f`.  `p` is given as a coefficient list
//      (low-order first).
//
//   inverse_poly(f, p) :
//      Divide the truncated series `f` by the polynomial `p`.  Same length
//      as `f`.  Throws if p[0] is zero.
//
//   rational_expansion(num, den, expansion) :
//      Returns InversePoly(TimesPoly(expansion, num), den) -- i.e. the
//      truncated series of (num / den) * expansion.

std::vector<numeric::AcbValue>
times_poly(const std::vector<numeric::AcbValue>& f,
           const std::vector<numeric::AcbValue>& poly,
           long prec = numeric::working_prec_bits());

std::vector<numeric::AcbValue>
inverse_poly(const std::vector<numeric::AcbValue>& f,
             const std::vector<numeric::AcbValue>& poly,
             long prec = numeric::working_prec_bits());

std::vector<numeric::AcbValue>
rational_expansion(const std::vector<numeric::AcbValue>& num,
                   const std::vector<numeric::AcbValue>& den,
                   const std::vector<numeric::AcbValue>& expansion,
                   long prec = numeric::working_prec_bits());

// MapRationalExpansion :
//   For each i, compute  rational_expansion(nums[i][k], dens[i][k], expansions[k])
//   and sum over k.  Result is a vector of truncated series, one per i.
std::vector<std::vector<numeric::AcbValue>>
map_rational_expansion(const std::vector<std::vector<std::vector<numeric::AcbValue>>>& nums,
                       const std::vector<std::vector<std::vector<numeric::AcbValue>>>& dens,
                       const std::vector<std::vector<numeric::AcbValue>>& expansions,
                       long prec = numeric::working_prec_bits());

// ---------------------------------------------------------------------------
//  Per-mu expansion utilities  (Extend / Plus / Rescale)
// ---------------------------------------------------------------------------
//
//  An "expansion matrix" is `vector<vector<AcbValue>>` where rows are log
//  powers and columns are eta orders.  These helpers do not touch the mu
//  exponent.

// Pad with zero rows so that exp.size() == n.  Each new row has length
// (XOrder + 1).  No-op if exp.size() >= n.
void extend_expansion(std::vector<std::vector<numeric::AcbValue>>& exp,
                      std::size_t n);

// Sum a list of expansion matrices, padding rows to the maximum length.
std::vector<std::vector<numeric::AcbValue>>
plus_expansion(const std::vector<std::vector<std::vector<numeric::AcbValue>>>& exps,
               long prec = numeric::working_prec_bits());

// Shift each row to the right by `order` columns; truncate to (XOrder + 1).
//   exp -> exp * eta^order  (logically, in the eta-order grid)
std::vector<std::vector<numeric::AcbValue>>
rescale_expansion(const std::vector<std::vector<numeric::AcbValue>>& exp,
                  std::size_t order);

// ---------------------------------------------------------------------------
//  Evaluation
// ---------------------------------------------------------------------------
//
//  evaluate_expansion(exp, x0)
//     = sum_{k, n} exp[k][n] * log^k(x0) * x0^n
//
//  evaluate_asy_term(mu, exp, x0)
//     = x0^mu * evaluate_expansion(exp, x0)

numeric::AcbValue
evaluate_expansion(const std::vector<std::vector<numeric::AcbValue>>& exp,
                   acb_srcptr x0,
                   long prec = numeric::working_prec_bits());

numeric::AcbValue
evaluate_asy_term(acb_srcptr mu,
                  const std::vector<std::vector<numeric::AcbValue>>& exp,
                  acb_srcptr x0,
                  long prec = numeric::working_prec_bits());

// ---------------------------------------------------------------------------
//  AsyTerm / AsyExpansion
// ---------------------------------------------------------------------------

struct AsyTerm {
    numeric::AcbValue                            mu;     // exponent
    std::vector<std::vector<numeric::AcbValue>>  exp;    // exp[log_power][eta_order]

    AsyTerm() = default;
    AsyTerm(numeric::AcbValue m,
            std::vector<std::vector<numeric::AcbValue>> e)
        : mu(std::move(m)), exp(std::move(e)) {}

    AsyTerm(AsyTerm&&) noexcept = default;
    AsyTerm& operator=(AsyTerm&&) noexcept = default;

    // Disable accidental copy; use clone() if needed.
    AsyTerm(const AsyTerm&)            = delete;
    AsyTerm& operator=(const AsyTerm&) = delete;

    AsyTerm clone() const;

    numeric::AcbValue evaluate(acb_srcptr x0,
                               long prec = numeric::working_prec_bits()) const {
        return evaluate_asy_term(mu.raw(), exp, x0, prec);
    }
};

using AsyExpansion = std::vector<AsyTerm>;

// Apply EvaluateAsymptoticExpansion to each term and sum.
numeric::AcbValue
evaluate_asy_expansion(const AsyExpansion& asy,
                       acb_srcptr           x0,
                       long                  prec = numeric::working_prec_bits());

// ---------------------------------------------------------------------------
//  Series-truncation tail bound  (bug U1: x_order silent-truncation fix)
// ---------------------------------------------------------------------------
//
// The eta-flow Taylor sums  Sum_{n=0}^{N} c_n x^n  at every matching /
// propagation point are truncated at N = x_order (+ extra_x_order for
// CalcTaylor) but the dropped tail  Sum_{n>N} c_n x^n  was NOT added to the
// returned Arb radius -- the midpoint-only Horner in evaluate_taylor()
// actively *zeroes* it.  The returned ball therefore claimed ~10^-working_pre
// while the actual value carried a ~|c_N x^N|/(1-r) bias (kite@w=5,
// x_order=100: ~32 true digits vs ~341 claimed).
//
// `series_tail_mag` writes a geometric estimate of the tail into `out`:
//
//   tail ~ |c_{N-2} x^{N-2} + c_{N-1} x^{N-1} + c_N x^N| / (1 - r),
//   r    = |c_N|/|c_{N-1}| * |x|,  clamped to 0.9 if r >= 1.
//
// (The last-3-terms partial is summed as a complex number so sign/phase
// cancellation is respected; the naive max|c_k|*|x|^N bound over-estimated
// the cancellation-heavy regular-contour series by ~250 orders on the kite.)
//
// Gated behind env var AMFLOW_XORDER_TAIL_BOUND=1 (default OFF -- does not
// change behaviour of jobs already running on the old expectation).

bool xorder_tail_bound_enabled();

void series_tail_mag(mag_t out,
                     const std::vector<numeric::AcbValue>& c,
                     acb_srcptr x);

// Deep-copy helper for a list of AsyTerm.
AsyExpansion clone_asy(const AsyExpansion& asy);

// Make a list-of-lists clone.
std::vector<AsyExpansion> clone_asy_list(const std::vector<AsyExpansion>& list);

// ---------------------------------------------------------------------------
//  Rule-set algebra  (PlusRuleS / UnionRuleS / PSTimesRuleS / PSMapRuleS)
// ---------------------------------------------------------------------------

// Test: do two mu values differ by a real integer (within numerical chop)?
bool mu_differ_by_integer(acb_srcptr mu1, acb_srcptr mu2,
                          long prec = numeric::working_prec_bits());

// Single-region combine (PlusRuleS):
//   * Find the minimum mu (in the integer-offset sense -- equivalent to
//     "subtract the leftmost mu from each, take the smallest result").
//   * Shift every term's mu to that minimum and rescale_expansion accordingly.
//   * plus_expansion the rescaled matrices.
// Returns the combined term `(min_mu, combined_exp)`.
AsyTerm plus_rule_set(const AsyExpansion& region,
                      long prec = numeric::working_prec_bits());

// Multi-region union (UnionRuleS):
//   * Group terms by "mu differ by integer".
//   * Apply plus_rule_set within each group.
AsyExpansion union_rule_set(const AsyExpansion& asy,
                            long prec = numeric::working_prec_bits());

// Power-series times an asymptotic expansion (PSTimesRuleS):
//   ps = (coeffs, leading_offset).  Each input AsyTerm
//        (mu_j, exp_j) becomes terms (mu_j + leading_offset + i, exp_j * coeffs[i])
//   for i = 0 .. coeffs.size() - 1.  Result is concatenated, NOT yet unioned.
AsyExpansion ps_times_rule_set(const PowerSeries& ps,
                               const AsyExpansion& asy,
                               long prec = numeric::working_prec_bits());

// Matrix variant (PSMapRuleS):
//   psmap is an M x K grid of PowerSeries (the result of to_power_series_matrix
//   on the variable-change matrix T).  asy_list is a list of K AsyExpansions
//   (one per source integral).  For each output index i (0 <= i < M),
//   compute  union_rule_set( concat over k of ps_times_rule_set(psmap[i][k], asy_list[k]) ).
std::vector<AsyExpansion>
ps_map_rule_set(const std::vector<std::vector<PowerSeries>>& psmap,
                const std::vector<AsyExpansion>& asy_list,
                long prec = numeric::working_prec_bits());

// ---------------------------------------------------------------------------
//  RationalFunction times AsyExpansion  (TimesAsyExp / PlusAsyExp)
// ---------------------------------------------------------------------------

// TimesAsyExp[rational, asy]:
//   Treat `rational` as an explicit eta-power factor times a polynomial /
//   polynomial fraction.  Pull out the common eta-offset of the
//   denominator, evaluate (poly/poly) on each row of each AsyTerm via
//   rational_expansion, and shift mu by  (-denominator_eta_offset).
AsyExpansion times_asy_exp(const numeric::RationalFunction& rational,
                           const AsyExpansion& asy,
                           long prec = numeric::working_prec_bits());

// PlusAsyExp[list]:
//   Concatenate a list of AsyExpansions and apply union_rule_set.
AsyExpansion plus_asy_exp(const std::vector<AsyExpansion>& list,
                          long prec = numeric::working_prec_bits());

// ---------------------------------------------------------------------------
//  Diagnostics
// ---------------------------------------------------------------------------

// FindLogPower[exp] = largest k such that row exp[k] has any nonzero entry,
// or -1 if all rows are entirely zero.
long find_log_power(const std::vector<std::vector<numeric::AcbValue>>& exp,
                    int chop_digits = numeric::chop_pre());

}  // namespace amflow::ode

#endif  // AMFLOW_ODE_ASY_HPP
