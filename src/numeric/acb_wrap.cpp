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

namespace amflow::numeric {

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
