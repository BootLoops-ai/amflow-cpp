// SPDX-License-Identifier: MIT
// numeric::Rational — FmpqPoly + RationalComplex + RationalFunction.
//
// amflow::numeric.

#include "amflow/numeric/rational.hpp"

#include <ostream>
#include <stdexcept>
#include <string>

#include <flint/acb.h>
#include <flint/acb_poly.h>
#include <flint/arb.h>
#include <flint/fmpq.h>
#include <flint/fmpq_poly.h>

namespace amflow::numeric {

// ===========================================================================
//  Free helpers
// ===========================================================================

long fmpq_poly_valuation(const fmpq_poly_t p) {
    if (fmpq_poly_is_zero(p)) return -1;
    long len = fmpq_poly_length(p);
    fmpq_t coeff;
    fmpq_init(coeff);
    long v = -1;
    for (long i = 0; i < len; ++i) {
        fmpq_poly_get_coeff_fmpq(coeff, p, i);
        if (!fmpq_is_zero(coeff)) {
            v = i;
            break;
        }
    }
    fmpq_clear(coeff);
    return v;
}

void fmpq_poly_make_denominator_monic_pair(fmpq_poly_t num, fmpq_poly_t den) {
    if (fmpq_poly_is_zero(den)) return;
    long d = fmpq_poly_degree(den);
    fmpq_t lead;
    fmpq_init(lead);
    fmpq_poly_get_coeff_fmpq(lead, den, d);
    if (fmpq_is_zero(lead)) {
        fmpq_clear(lead);
        return;
    }
    if (!fmpq_is_one(lead)) {
        fmpq_t inv;
        fmpq_init(inv);
        fmpq_inv(inv, lead);
        fmpq_poly_scalar_mul_fmpq(num, num, inv);
        fmpq_poly_scalar_mul_fmpq(den, den, inv);
        fmpq_clear(inv);
    }
    fmpq_clear(lead);
}

// ===========================================================================
//  FmpqPoly
// ===========================================================================

long FmpqPoly::valuation() const { return fmpq_poly_valuation(handle_); }

std::string FmpqPoly::to_string(const char* var) const {
    char* s = fmpq_poly_get_str_pretty(handle_, var);
    std::string out = (s != nullptr) ? s : "0";
    if (s != nullptr) flint_free(s);
    return out;
}

std::ostream& operator<<(std::ostream& os, const FmpqPoly& p) {
    return os << p.to_string();
}

// ===========================================================================
//  RationalComplex
// ===========================================================================

RationalComplex::RationalComplex() {
    fmpq_init(re_);
    fmpq_init(im_);
}
RationalComplex::~RationalComplex() {
    fmpq_clear(re_);
    fmpq_clear(im_);
}

RationalComplex::RationalComplex(const RationalComplex& other) {
    fmpq_init(re_);
    fmpq_init(im_);
    fmpq_set(re_, other.re_);
    fmpq_set(im_, other.im_);
}

RationalComplex::RationalComplex(RationalComplex&& other) noexcept {
    fmpq_init(re_);
    fmpq_init(im_);
    fmpq_swap(re_, other.re_);
    fmpq_swap(im_, other.im_);
}

RationalComplex& RationalComplex::operator=(const RationalComplex& other) {
    if (this != &other) {
        fmpq_set(re_, other.re_);
        fmpq_set(im_, other.im_);
    }
    return *this;
}

RationalComplex& RationalComplex::operator=(RationalComplex&& other) noexcept {
    if (this != &other) {
        fmpq_swap(re_, other.re_);
        fmpq_swap(im_, other.im_);
    }
    return *this;
}

RationalComplex RationalComplex::from_fractions(long re_num, long re_den,
                                                long im_num, long im_den) {
    if (re_den == 0 || im_den == 0) {
        throw std::invalid_argument("RationalComplex::from_fractions: zero denominator");
    }
    RationalComplex z;
    fmpq_set_si(z.re_, re_num, re_den);
    fmpq_set_si(z.im_, im_num, im_den);
    return z;
}

bool RationalComplex::is_real() const { return fmpq_is_zero(im_); }
bool RationalComplex::is_zero() const { return fmpq_is_zero(re_) && fmpq_is_zero(im_); }

AcbValue RationalComplex::to_acb(long prec) const {
    AcbValue out;
    arb_set_fmpq(acb_realref(out.raw()), re_, prec);
    arb_set_fmpq(acb_imagref(out.raw()), im_, prec);
    return out;
}

std::string RationalComplex::to_string() const {
    char* sr = fmpq_get_str(nullptr, 10, re_);
    char* si = fmpq_get_str(nullptr, 10, im_);
    std::string out;
    out += (sr != nullptr) ? sr : "0";
    out += " + ";
    out += (si != nullptr) ? si : "0";
    out += "*I";
    if (sr != nullptr) flint_free(sr);
    if (si != nullptr) flint_free(si);
    return out;
}

bool RationalComplex::operator==(const RationalComplex& other) const {
    return fmpq_equal(re_, other.re_) && fmpq_equal(im_, other.im_);
}

std::ostream& operator<<(std::ostream& os, const RationalComplex& z) {
    return os << z.to_string();
}

// ===========================================================================
//  RationalFunction
// ===========================================================================

RationalFunction::RationalFunction() {
    num_.set_zero();
    den_.set_one();
}

RationalFunction::RationalFunction(FmpqPoly num, FmpqPoly den)
    : num_(std::move(num)), den_(std::move(den)) {
    if (den_.is_zero()) {
        throw std::invalid_argument("RationalFunction: denominator must be nonzero");
    }
    reduce();
}

RationalFunction RationalFunction::from_polynomial(FmpqPoly num) {
    FmpqPoly den;
    den.set_one();
    return RationalFunction(std::move(num), std::move(den));
}

RationalFunction RationalFunction::from_fmpq(const fmpq_t q) {
    FmpqPoly num;
    num.set_fmpq(q);
    FmpqPoly den;
    den.set_one();
    return RationalFunction(std::move(num), std::move(den));
}

RationalFunction RationalFunction::from_si(long n) {
    FmpqPoly num;
    num.set_si(n);
    FmpqPoly den;
    den.set_one();
    return RationalFunction(std::move(num), std::move(den));
}

RationalFunction RationalFunction::from_si_si(long p, long q) {
    if (q == 0) throw std::invalid_argument("from_si_si: zero denominator");
    fmpq_t r;
    fmpq_init(r);
    fmpq_set_si(r, p, q);
    auto out = from_fmpq(r);
    fmpq_clear(r);
    return out;
}

RationalFunction RationalFunction::monomial(long k) {
    fmpq_t one;
    fmpq_init(one);
    fmpq_one(one);
    auto out = monomial(k, one);
    fmpq_clear(one);
    return out;
}

RationalFunction RationalFunction::monomial(long k, const fmpq_t c) {
    FmpqPoly num, den;
    if (k >= 0) {
        num.set_coeff_fmpq(k, c);
        den.set_one();
    } else {
        num.set_fmpq(c);
        den.set_coeff_si(-k, 1);
    }
    return RationalFunction(std::move(num), std::move(den));
}

RationalFunction::RationalFunction(const RationalFunction& other)
    : num_(other.num_.clone()), den_(other.den_.clone()) {}

RationalFunction& RationalFunction::operator=(const RationalFunction& other) {
    if (this != &other) {
        num_ = other.num_.clone();
        den_ = other.den_.clone();
    }
    return *this;
}

bool RationalFunction::is_polynomial() const {
    return den_.degree() == 0 && !den_.is_zero();
}

bool RationalFunction::is_constant() const {
    return num_.degree() <= 0 && den_.degree() == 0;
}

long RationalFunction::net_valuation() const {
    if (is_zero()) return 0;
    long vn = num_.valuation();
    long vd = den_.valuation();
    if (vn < 0) return 0;
    return vn - vd;
}

void RationalFunction::reduce() {
    // 1. Cancel common polynomial factor (fmpq_poly_gcd returns a monic gcd).
    fmpq_poly_t g;
    fmpq_poly_init(g);
    fmpq_poly_gcd(g, num_.raw(), den_.raw());

    if (!fmpq_poly_is_zero(g) && fmpq_poly_degree(g) > 0) {
        fmpq_poly_div(num_.raw(), num_.raw(), g);
        fmpq_poly_div(den_.raw(), den_.raw(), g);
    }
    fmpq_poly_clear(g);

    // 2. Canonicalise: 0/1 if num is zero; else make denominator monic.
    if (num_.is_zero()) {
        den_.set_one();
        return;
    }
    fmpq_poly_make_denominator_monic_pair(num_.raw(), den_.raw());
}

RationalFunction RationalFunction::operator-() const {
    RationalFunction out = *this;
    fmpq_poly_neg(out.num_.raw(), out.num_.raw());
    return out;
}

// (a/b) + (c/d) = (a*d + b*c) / (b*d)
RationalFunction& RationalFunction::operator+=(const RationalFunction& other) {
    fmpq_poly_t ad, bc, bd;
    fmpq_poly_init(ad);
    fmpq_poly_init(bc);
    fmpq_poly_init(bd);

    fmpq_poly_mul(ad, num_.raw(),       other.den_.raw());
    fmpq_poly_mul(bc, den_.raw(),       other.num_.raw());
    fmpq_poly_mul(bd, den_.raw(),       other.den_.raw());
    fmpq_poly_add(num_.raw(), ad, bc);
    fmpq_poly_set(den_.raw(), bd);

    fmpq_poly_clear(ad);
    fmpq_poly_clear(bc);
    fmpq_poly_clear(bd);
    reduce();
    return *this;
}

RationalFunction& RationalFunction::operator-=(const RationalFunction& other) {
    RationalFunction neg_other = -other;
    return *this += neg_other;
}

RationalFunction& RationalFunction::operator*=(const RationalFunction& other) {
    fmpq_poly_mul(num_.raw(), num_.raw(), other.num_.raw());
    fmpq_poly_mul(den_.raw(), den_.raw(), other.den_.raw());
    reduce();
    return *this;
}

RationalFunction& RationalFunction::operator/=(const RationalFunction& other) {
    if (other.is_zero()) throw std::domain_error("RationalFunction divide by zero");
    fmpq_poly_mul(num_.raw(), num_.raw(), other.den_.raw());
    fmpq_poly_mul(den_.raw(), den_.raw(), other.num_.raw());
    reduce();
    return *this;
}

RationalFunction& RationalFunction::multiply_by_polynomial(const FmpqPoly& p) {
    fmpq_poly_mul(num_.raw(), num_.raw(), p.raw());
    reduce();
    return *this;
}

RationalFunction& RationalFunction::multiply_by_eta_power(long k) {
    if (k == 0) return *this;
    fmpq_poly_t shift;
    fmpq_poly_init(shift);
    if (k > 0) {
        fmpq_poly_set_coeff_si(shift, k, 1);
        fmpq_poly_mul(num_.raw(), num_.raw(), shift);
    } else {
        fmpq_poly_set_coeff_si(shift, -k, 1);
        fmpq_poly_mul(den_.raw(), den_.raw(), shift);
    }
    fmpq_poly_clear(shift);
    reduce();
    return *this;
}

RationalFunction& RationalFunction::multiply_by_fmpq(const fmpq_t q) {
    fmpq_poly_scalar_mul_fmpq(num_.raw(), num_.raw(), q);
    reduce();
    return *this;
}

RationalFunction RationalFunction::reciprocal() const {
    if (is_zero()) throw std::domain_error("RationalFunction reciprocal of zero");
    return RationalFunction(den_.clone(), num_.clone());
}

RationalFunction RationalFunction::derivative() const {
    // d/dη (n/d) = (n' d - n d') / d^2
    FmpqPoly nprime, dprime, lhs, rhs, top, bottom;
    fmpq_poly_derivative(nprime.raw(), num_.raw());
    fmpq_poly_derivative(dprime.raw(), den_.raw());

    fmpq_poly_mul(lhs.raw(), nprime.raw(), den_.raw());
    fmpq_poly_mul(rhs.raw(), num_.raw(),   dprime.raw());
    fmpq_poly_sub(top.raw(), lhs.raw(),    rhs.raw());

    fmpq_poly_mul(bottom.raw(), den_.raw(), den_.raw());

    return RationalFunction(std::move(top), std::move(bottom));
}

AcbValue RationalFunction::evaluate(acb_srcptr x, long prec) const {
    // FLINT 3.x dropped fmpq_poly_evaluate_acb.  Convert to acb_poly first.
    acb_poly_t pn, pd;
    acb_poly_init(pn);
    acb_poly_init(pd);
    acb_poly_set_fmpq_poly(pn, num_.raw(), prec);
    acb_poly_set_fmpq_poly(pd, den_.raw(), prec);

    AcbValue numv, denv, out;
    acb_poly_evaluate(numv.raw(), pn, x, prec);
    acb_poly_evaluate(denv.raw(), pd, x, prec);

    acb_poly_clear(pn);
    acb_poly_clear(pd);

    if (acb_contains_zero(denv.raw())) {
        // After canonical reduction the denominator vanishes only at true
        // poles.  contains_zero is conservative: "ball touches zero" is
        // treated as indeterminate.
        throw std::domain_error("RationalFunction::evaluate: denominator vanishes at the requested point");
    }
    acb_div(out.raw(), numv.raw(), denv.raw(), prec);
    return out;
}

AcbValue RationalFunction::evaluate(const RationalComplex& x, long prec) const {
    AcbValue ax = x.to_acb(prec);
    return evaluate(ax.raw(), prec);
}

// (a(η)/b(η)) /. η -> 1/η
//   = a(1/η) / b(1/η)
//
// p(1/η) = sum_i a_i η^-i = η^-d * sum_i a_i η^(d-i) = η^-d * reverse(p)
// so
//   p(1/η) / q(1/η) = (η^-dp * rev(p)) / (η^-dq * rev(q))
//                   = η^(dq - dp) * rev(p) / rev(q)
RationalFunction RationalFunction::substitute_inverse_eta() const {
    if (is_zero()) return *this;

    long dn = num_.degree();
    long dd = den_.degree();
    if (dn < 0) dn = 0;
    if (dd < 0) dd = 0;

    FmpqPoly rev_num, rev_den;
    {
        fmpq_t c;
        fmpq_init(c);
        for (long i = 0; i <= dn; ++i) {
            fmpq_poly_get_coeff_fmpq(c, num_.raw(), i);
            if (!fmpq_is_zero(c)) rev_num.set_coeff_fmpq(dn - i, c);
        }
        for (long i = 0; i <= dd; ++i) {
            fmpq_poly_get_coeff_fmpq(c, den_.raw(), i);
            if (!fmpq_is_zero(c)) rev_den.set_coeff_fmpq(dd - i, c);
        }
        fmpq_clear(c);
    }

    RationalFunction out(std::move(rev_num), std::move(rev_den));
    long shift = dd - dn;
    if (shift != 0) out.multiply_by_eta_power(shift);
    return out;
}

std::string RationalFunction::to_string(const char* var) const {
    if (is_zero()) return "0";
    if (den_.is_one()) return num_.to_string(var);
    std::string out = "(";
    out += num_.to_string(var);
    out += ") / (";
    out += den_.to_string(var);
    out += ")";
    return out;
}

bool RationalFunction::operator==(const RationalFunction& other) const {
    return fmpq_poly_equal(num_.raw(), other.num_.raw())
        && fmpq_poly_equal(den_.raw(), other.den_.raw());
}

RationalFunction operator+(RationalFunction a, const RationalFunction& b) { a += b; return a; }
RationalFunction operator-(RationalFunction a, const RationalFunction& b) { a -= b; return a; }
RationalFunction operator*(RationalFunction a, const RationalFunction& b) { a *= b; return a; }
RationalFunction operator/(RationalFunction a, const RationalFunction& b) { a /= b; return a; }

std::ostream& operator<<(std::ostream& os, const RationalFunction& r) {
    return os << r.to_string();
}

}  // namespace amflow::numeric
