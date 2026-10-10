// SPDX-License-Identifier: MIT
// numeric::acb_wrap — implementation.
//

#include "amflow/numeric/acb_wrap.hpp"

#include <cmath>
#include <ostream>
#include <string>

#include "amflow/numeric/log.hpp"

#include <flint/acb.h>
#include <flint/arb.h>
#include <flint/arf.h>
#include <flint/fmpq.h>
#include <flint/fmpz.h>

namespace amflow::numeric {

// ---------------------------------------------------------------------------
//  Simplest rational within tolerance.  See acb_wrap.hpp.
// ---------------------------------------------------------------------------

bool arb_to_fmpq_simplest(fmpq_t out, arb_srcptr x, int digits) {
    if (!arb_is_finite(x)) return false;
    if (digits < 0) digits = 0;

    // A midpoint below 2^-bits(digits) <= 10^-digits rationalizes to 0;
    // answer before materializing it, so a midpoint with a huge negative
    // exponent never allocates a denominator of that many bits.
    if (arf_cmpabs_2exp_si(arb_midref(x), -decimal_digits_to_bits(digits)) < 0) {
        fmpq_zero(out);
        return true;
    }

    fmpq_t mid, tol, lo, hi;
    fmpq_init(mid);
    fmpq_init(tol);
    fmpq_init(lo);
    fmpq_init(hi);

    // Both ingredients are exact: an arf midpoint is a dyadic rational and
    // 10^-digits is formed from an integer power.  The radius is ignored
    // (see the header comment).
    arf_get_fmpq(mid, arb_midref(x));
    {
        fmpz_t pow10;
        fmpz_init(pow10);
        fmpz_set_ui(pow10, 10);
        fmpz_pow_ui(pow10, pow10, static_cast<unsigned long>(digits));
        fmpq_one(tol);
        fmpq_div_fmpz(tol, tol, pow10);
        fmpz_clear(pow10);
    }

    fmpq_sub(lo, mid, tol);
    fmpq_add(hi, mid, tol);
    fmpq_simplest_between(out, lo, hi);

    fmpq_clear(mid);
    fmpq_clear(tol);
    fmpq_clear(lo);
    fmpq_clear(hi);
    return true;
}

bool acb_real_to_fmpq_simplest(fmpq_t out, acb_srcptr x, int digits) {
    arf_t thr;
    arf_init(thr);
    arf_one(thr);
    arf_mul_2exp_si(thr, thr, -decimal_digits_to_bits(digits));
    const bool im_zero = (arf_cmpabs(arb_midref(acb_imagref(x)), thr) < 0);
    arf_clear(thr);
    if (!im_zero) return false;
    return arb_to_fmpq_simplest(out, acb_realref(x), digits);
}

// ---------------------------------------------------------------------------
//  Chop predicate.  See comment in acb_wrap.hpp / numeric.md §8c.
// ---------------------------------------------------------------------------

bool acb_is_chop_zero(acb_srcptr x, int chop_digits) {
    if (chop_digits <= 0) return false;

    // PERF: a chop_digits>300 branch that recomputes
    // (1/10)^chop_digits via arb_pow_ui on EVERY call is quadratic waste.  Under
    // sparse_gaussian's chop_sparse this fired O(10^8) times per
    // CalcTaylor on a 44-master 3-loop block at working_pre=390 →
    // sparse_chop_digits=350, and perf attributed ~95% of CalcTaylor
    // wall time to arb_pow_fmpz_binexp under this function (mpn_sqr/mul
    // _basecase 67% of total samples).  A chop predicate only needs an
    // order-of-magnitude threshold; compare |mid| against 2^e directly.
    //
    //   |x| < 10^-d  ⇔  |x| < 2^{-d·log2(10)}  ≈  2^{-ceil(3.3219·d)}
    //
    // arf_cmpabs_2exp_si is O(1) (exponent compare, no limb arithmetic,
    // no allocation).  The ≤1-bit slack vs the exact 10^-d threshold is
    // immaterial for a chop predicate — the d≤300 path already used a
    // 53-bit-mantissa double approximation to 10^-d.
    // Inline d→bits without decimal_digits_to_bits()'s 53-bit floor so
    // small chop_digits (unused in practice, but be defensive) keep the
    // intended 10^-d magnitude.
    const slong neg_bits = -static_cast<slong>(
        std::ceil(static_cast<double>(chop_digits) * 3.321928094887362));
    arf_srcptr mid_re = arb_midref(acb_realref(x));
    arf_srcptr mid_im = arb_midref(acb_imagref(x));
    return (arf_cmpabs_2exp_si(mid_re, neg_bits) < 0)
        && (arf_cmpabs_2exp_si(mid_im, neg_bits) < 0);
}

bool AcbValue::is_chop_zero(int chop_digits) const {
    return acb_is_chop_zero(handle_, chop_digits);
}

// ---------------------------------------------------------------------------
//  Pretty printing
// ---------------------------------------------------------------------------

std::string AcbValue::to_string(long display_digits) const {
    // FLINT 3 dropped acb_get_str.  Build "(re + im*I)" by formatting the
    // two arb halves independently.
    //
    // AMFLOW_DEBUG_PRECISE=1 sets ARB_STR_MORE|ARB_STR_NO_RADIUS, which
    // prints midpoint digits beyond the ball radius (diagnostic only —
    // untrusted digits).
    int arb_flags = 0;
    if (::amflow::numeric::log::trace_enabled("AMFLOW_DEBUG_PRECISE")) arb_flags = 1 | 2;
    char* re = arb_get_str(acb_realref(handle_), display_digits, arb_flags);
    char* im = arb_get_str(acb_imagref(handle_), display_digits, arb_flags);
    std::string out = "(";
    out += (re != nullptr) ? re : "<null>";
    out += " + ";
    out += (im != nullptr) ? im : "<null>";
    out += "*I)";
    if (re != nullptr) flint_free(re);
    if (im != nullptr) flint_free(im);
    return out;
}

std::ostream& operator<<(std::ostream& os, const AcbValue& v) {
    return os << v.to_string(20);
}

}  // namespace amflow::numeric
