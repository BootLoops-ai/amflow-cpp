// SPDX-License-Identifier: MIT
// algebra::mpoly — implementation.
//
// amflow::algebra.

#include "amflow/algebra/mpoly.hpp"

#include <algorithm>
#include <ostream>
#include <stdexcept>
#include <utility>

#include <flint/acb.h>
#include <flint/fmpq.h>
#include <flint/fmpz.h>
#include <flint/fmpz_mpoly.h>
#include <flint/fmpz_mpoly_q.h>

namespace amflow::algebra {

// ===========================================================================
//  MpolyContext
// ===========================================================================

MpolyContext::MpolyContext(std::vector<std::string> var_names)
    : var_names_(std::move(var_names)),
      var_cstr_(),
      ctx_(nullptr) {
    if (var_names_.empty()) {
        throw std::invalid_argument("MpolyContext: at least one variable required");
    }
    // Deduplicate check.
    std::vector<std::string> sorted = var_names_;
    std::sort(sorted.begin(), sorted.end());
    for (std::size_t i = 1; i < sorted.size(); ++i) {
        if (sorted[i] == sorted[i - 1]) {
            throw std::invalid_argument(
                "MpolyContext: duplicate variable name '" + sorted[i] + "'");
        }
    }

    var_cstr_.reserve(var_names_.size());
    for (const auto& n : var_names_) var_cstr_.push_back(n.c_str());

    ctx_ = static_cast<fmpz_mpoly_ctx_struct*>(
        flint_malloc(sizeof(fmpz_mpoly_ctx_struct)));
    fmpz_mpoly_ctx_init(ctx_, static_cast<long>(var_names_.size()), ORD_LEX);
}

MpolyContext::~MpolyContext() {
    if (ctx_) {
        fmpz_mpoly_ctx_clear(ctx_);
        flint_free(ctx_);
        ctx_ = nullptr;
    }
}

long MpolyContext::var_index(const std::string& name) const {
    for (std::size_t i = 0; i < var_names_.size(); ++i) {
        if (var_names_[i] == name) return static_cast<long>(i);
    }
    return -1;
}

const std::string& MpolyContext::var_name(long idx) const {
    if (idx < 0 || idx >= static_cast<long>(var_names_.size())) {
        throw std::out_of_range("MpolyContext::var_name: index out of range");
    }
    return var_names_[idx];
}

// ===========================================================================
//  Internal helpers
// ===========================================================================

namespace {

fmpz_mpoly_struct* alloc_mpoly() {
    return static_cast<fmpz_mpoly_struct*>(flint_malloc(sizeof(fmpz_mpoly_struct)));
}
fmpz_mpoly_q_struct* alloc_mfrac() {
    return static_cast<fmpz_mpoly_q_struct*>(flint_malloc(sizeof(fmpz_mpoly_q_struct)));
}

}  // namespace

// ===========================================================================
//  Mpoly
// ===========================================================================

Mpoly::Mpoly() noexcept
    : ctx_(nullptr), p_(nullptr) {}

Mpoly::Mpoly(std::shared_ptr<MpolyContext> ctx)
    : ctx_(std::move(ctx)), p_(alloc_mpoly()) {
    fmpz_mpoly_init(p_, ctx_->raw());
}

Mpoly::Mpoly(Mpoly&& other) noexcept
    : ctx_(std::move(other.ctx_)), p_(other.p_) {
    other.p_ = nullptr;
}

Mpoly& Mpoly::operator=(Mpoly&& other) noexcept {
    if (this != &other) {
        if (p_ && ctx_) {
            fmpz_mpoly_clear(p_, ctx_->raw());
            flint_free(p_);
        }
        ctx_ = std::move(other.ctx_);
        p_   = other.p_;
        other.p_ = nullptr;
    }
    return *this;
}

Mpoly::~Mpoly() {
    if (p_ && ctx_) {
        fmpz_mpoly_clear(p_, ctx_->raw());
        flint_free(p_);
        p_ = nullptr;
    }
}

Mpoly Mpoly::clone() const {
    Mpoly out(ctx_);
    fmpz_mpoly_set(out.p_, p_, ctx_->raw());
    return out;
}

Mpoly Mpoly::zero(std::shared_ptr<MpolyContext> ctx) {
    return Mpoly(std::move(ctx));
}

Mpoly Mpoly::one(std::shared_ptr<MpolyContext> ctx) {
    Mpoly out(std::move(ctx));
    fmpz_mpoly_set_si(out.p_, 1, out.ctx_->raw());
    return out;
}

Mpoly Mpoly::constant(std::shared_ptr<MpolyContext> ctx, long c) {
    Mpoly out(std::move(ctx));
    fmpz_mpoly_set_si(out.p_, c, out.ctx_->raw());
    return out;
}

Mpoly Mpoly::constant(std::shared_ptr<MpolyContext> ctx, const fmpz_t c) {
    Mpoly out(std::move(ctx));
    fmpz_mpoly_set_fmpz(out.p_, c, out.ctx_->raw());
    return out;
}

Mpoly Mpoly::variable(std::shared_ptr<MpolyContext> ctx, long var_idx) {
    if (var_idx < 0 || var_idx >= ctx->n_vars()) {
        throw std::out_of_range("Mpoly::variable: index out of range");
    }
    Mpoly out(std::move(ctx));
    fmpz_mpoly_gen(out.p_, var_idx, out.ctx_->raw());
    return out;
}

Mpoly Mpoly::variable(std::shared_ptr<MpolyContext> ctx, const std::string& name) {
    long idx = ctx->var_index(name);
    if (idx < 0) {
        throw std::invalid_argument("Mpoly::variable: unknown variable '" + name + "'");
    }
    return variable(std::move(ctx), idx);
}

Mpoly Mpoly::monomial(std::shared_ptr<MpolyContext> ctx,
                      const fmpz_t c,
                      const std::vector<unsigned long>& exponents) {
    if (static_cast<long>(exponents.size()) != ctx->n_vars()) {
        throw std::invalid_argument("Mpoly::monomial: exponents size mismatch");
    }
    Mpoly out(std::move(ctx));
    if (fmpz_is_zero(c)) return out;
    fmpz_mpoly_set_coeff_fmpz_ui(out.p_, c, exponents.data(), out.ctx_->raw());
    return out;
}

Mpoly Mpoly::from_string(std::shared_ptr<MpolyContext> ctx,
                         const std::string& expr,
                         fmpz_t out_scale) {
    // Use the rational-function parser to accept fractional coefficients.
    fmpz_mpoly_q_t q;
    fmpz_mpoly_q_init(q, ctx->raw());
    int ok = fmpz_mpoly_q_set_str_pretty(
        q, expr.c_str(),
        const_cast<const char**>(ctx->var_names_cstr()),
        ctx->raw());
    if (ok != 0) {
        fmpz_mpoly_q_clear(q, ctx->raw());
        throw std::invalid_argument(
            "Mpoly::from_string: parse error in '" + expr + "'");
    }

    const fmpz_mpoly_struct* num = fmpz_mpoly_q_numref(q);
    const fmpz_mpoly_struct* den = fmpz_mpoly_q_denref(q);

    if (!fmpz_mpoly_is_fmpz(den, ctx->raw())) {
        fmpz_mpoly_q_clear(q, ctx->raw());
        throw std::invalid_argument(
            "Mpoly::from_string: input has non-trivial denominator (use Mfrac::from_string)");
    }
    fmpz_t den_scalar;
    fmpz_init(den_scalar);
    if (fmpz_mpoly_length(den, ctx->raw()) > 0) {
        fmpz_mpoly_get_term_coeff_fmpz(den_scalar, den, 0, ctx->raw());
    } else {
        fmpz_one(den_scalar);
    }

    Mpoly out(ctx);
    fmpz_mpoly_set(out.p_, num, ctx->raw());

    if (out_scale != nullptr) {
        fmpz_set(out_scale, den_scalar);
    } else {
        if (!fmpz_is_one(den_scalar)) {
            fmpz_clear(den_scalar);
            fmpz_mpoly_q_clear(q, ctx->raw());
            throw std::invalid_argument(
                "Mpoly::from_string: input has fractional coefficient; "
                "pass out_scale != nullptr to receive the common denominator");
        }
    }

    fmpz_clear(den_scalar);
    fmpz_mpoly_q_clear(q, ctx->raw());
    return out;
}

bool Mpoly::is_zero() const noexcept {
    return fmpz_mpoly_is_zero(p_, ctx_->raw()) != 0;
}
bool Mpoly::is_one() const noexcept {
    return fmpz_mpoly_is_one(p_, ctx_->raw()) != 0;
}
bool Mpoly::is_constant() const noexcept {
    return fmpz_mpoly_is_fmpz(p_, ctx_->raw()) != 0;
}

long Mpoly::total_degree() const noexcept {
    if (is_zero()) return -1;
    return fmpz_mpoly_total_degree_si(p_, ctx_->raw());
}

long Mpoly::degree(long var_idx) const noexcept {
    if (is_zero()) return -1;
    return fmpz_mpoly_degree_si(p_, var_idx, ctx_->raw());
}

void Mpoly::check_same_ctx(const Mpoly& other) const {
    if (ctx_.get() != other.ctx_.get()) {
        throw std::invalid_argument("Mpoly: context mismatch");
    }
}

Mpoly Mpoly::operator-() const {
    Mpoly out(ctx_);
    fmpz_mpoly_neg(out.p_, p_, ctx_->raw());
    return out;
}

Mpoly& Mpoly::operator+=(const Mpoly& other) {
    check_same_ctx(other);
    fmpz_mpoly_add(p_, p_, other.p_, ctx_->raw());
    return *this;
}

Mpoly& Mpoly::operator-=(const Mpoly& other) {
    check_same_ctx(other);
    fmpz_mpoly_sub(p_, p_, other.p_, ctx_->raw());
    return *this;
}

Mpoly& Mpoly::operator*=(const Mpoly& other) {
    check_same_ctx(other);
    fmpz_mpoly_mul(p_, p_, other.p_, ctx_->raw());
    return *this;
}

Mpoly& Mpoly::multiply_by_si(long c) {
    fmpz_mpoly_scalar_mul_si(p_, p_, c, ctx_->raw());
    return *this;
}

Mpoly& Mpoly::multiply_by_fmpz(const fmpz_t c) {
    fmpz_mpoly_scalar_mul_fmpz(p_, p_, c, ctx_->raw());
    return *this;
}

Mpoly Mpoly::derivative(long var_idx) const {
    if (var_idx < 0 || var_idx >= ctx_->n_vars()) {
        throw std::out_of_range("Mpoly::derivative: var_idx out of range");
    }
    Mpoly out(ctx_);
    fmpz_mpoly_derivative(out.p_, p_, var_idx, ctx_->raw());
    return out;
}

Mpoly Mpoly::coeff_of(const std::vector<long>& vars,
                      const std::vector<unsigned long>& exps) const {
    if (vars.size() != exps.size()) {
        throw std::invalid_argument("Mpoly::coeff_of: vars and exps must have same length");
    }
    for (long v : vars) {
        if (v < 0 || v >= ctx_->n_vars()) {
            throw std::out_of_range("Mpoly::coeff_of: variable index out of range");
        }
    }
    {
        std::vector<long> sorted = vars;
        std::sort(sorted.begin(), sorted.end());
        for (std::size_t i = 1; i < sorted.size(); ++i) {
            if (sorted[i] == sorted[i - 1]) {
                throw std::invalid_argument(
                    "Mpoly::coeff_of: duplicate variable in vars");
            }
        }
    }
    Mpoly out(ctx_);
    fmpz_mpoly_get_coeff_vars_ui(out.p_, p_,
                                 vars.data(),
                                 exps.data(),
                                 (slong)vars.size(),
                                 ctx_->raw());
    return out;
}

Mpoly Mpoly::coeff_of(long var, unsigned long power) const {
    return coeff_of(std::vector<long>{var}, std::vector<unsigned long>{power});
}

Mpoly Mpoly::substitute(long var_idx, const fmpz_t value) const {
    if (var_idx < 0 || var_idx >= ctx_->n_vars()) {
        throw std::out_of_range("Mpoly::substitute: var_idx out of range");
    }
    Mpoly out(ctx_);
    int ok = fmpz_mpoly_evaluate_one_fmpz(out.p_, p_, var_idx, value, ctx_->raw());
    if (ok == 0) {
        throw std::runtime_error(
            "Mpoly::substitute: FLINT failed (coefficient blowup?)");
    }
    return out;
}

Mpoly Mpoly::substitute(long var_idx, long value) const {
    fmpz_t v;
    fmpz_init(v);
    fmpz_set_si(v, value);
    Mpoly out = substitute(var_idx, v);
    fmpz_clear(v);
    return out;
}

int Mpoly::evaluate(fmpz * const * values, fmpz_t out) const {
    return fmpz_mpoly_evaluate_all_fmpz(out, p_, values, ctx_->raw());
}

void Mpoly::evaluate_acb(acb_ptr out, acb_srcptr values, long prec) const {
    // fmpz_mpoly_evaluate_acb is only declared in public FLINT headers from
    // 3.4.0 on (before that the symbol exists in libflint but has no public
    // prototype, so the call does not compile against FLINT 3.0-3.3).
    // Go through fmpz_mpoly_q_evaluate_acb instead, which is public in all
    // FLINT 3.x: with a denominator of one it takes a fast path that runs
    // the exact same polynomial evaluator on the numerator, so the result
    // (value and error bounds at the given precision) is identical.
    fmpz_mpoly_q_t q;
    fmpz_mpoly_q_init(q, ctx_->raw());
    fmpz_mpoly_set(fmpz_mpoly_q_numref(q), p_, ctx_->raw());
    fmpz_mpoly_one(fmpz_mpoly_q_denref(q), ctx_->raw());
    fmpz_mpoly_q_evaluate_acb(out, q, values, prec, ctx_->raw());
    fmpz_mpoly_q_clear(q, ctx_->raw());
}

std::string Mpoly::to_string() const {
    char* s = fmpz_mpoly_get_str_pretty(
        p_,
        const_cast<const char**>(ctx_->var_names_cstr()),
        ctx_->raw());
    std::string out = s ? s : "0";
    if (s) flint_free(s);
    return out;
}

bool Mpoly::operator==(const Mpoly& other) const {
    if (ctx_.get() != other.ctx_.get()) return false;
    return fmpz_mpoly_equal(p_, other.p_, ctx_->raw()) != 0;
}

Mpoly operator+(const Mpoly& a, const Mpoly& b) { Mpoly t = a.clone(); t += b; return t; }
Mpoly operator-(const Mpoly& a, const Mpoly& b) { Mpoly t = a.clone(); t -= b; return t; }
Mpoly operator*(const Mpoly& a, const Mpoly& b) { Mpoly t = a.clone(); t *= b; return t; }

std::ostream& operator<<(std::ostream& os, const Mpoly& p) {
    return os << p.to_string();
}

// ===========================================================================
//  Mfrac
// ===========================================================================

Mfrac::Mfrac() noexcept
    : ctx_(nullptr), q_(nullptr) {}

Mfrac::Mfrac(std::shared_ptr<MpolyContext> ctx)
    : ctx_(std::move(ctx)), q_(alloc_mfrac()) {
    fmpz_mpoly_q_init(q_, ctx_->raw());
}

Mfrac::Mfrac(Mpoly num, Mpoly den)
    : ctx_(num.ctx()), q_(alloc_mfrac()) {
    if (num.ctx().get() != den.ctx().get()) {
        flint_free(q_);
        throw std::invalid_argument("Mfrac(Mpoly, Mpoly): context mismatch");
    }
    if (den.is_zero()) {
        flint_free(q_);
        throw std::invalid_argument("Mfrac(Mpoly, Mpoly): denominator is zero");
    }
    fmpz_mpoly_q_init(q_, ctx_->raw());
    fmpz_mpoly_set(fmpz_mpoly_q_numref(q_), num.raw(), ctx_->raw());
    fmpz_mpoly_set(fmpz_mpoly_q_denref(q_), den.raw(), ctx_->raw());
    fmpz_mpoly_q_canonicalise(q_, ctx_->raw());
}

Mfrac::Mfrac(Mfrac&& other) noexcept
    : ctx_(std::move(other.ctx_)), q_(other.q_) {
    other.q_ = nullptr;
}

Mfrac& Mfrac::operator=(Mfrac&& other) noexcept {
    if (this != &other) {
        if (q_ && ctx_) {
            fmpz_mpoly_q_clear(q_, ctx_->raw());
            flint_free(q_);
        }
        ctx_ = std::move(other.ctx_);
        q_   = other.q_;
        other.q_ = nullptr;
    }
    return *this;
}

Mfrac::~Mfrac() {
    if (q_ && ctx_) {
        fmpz_mpoly_q_clear(q_, ctx_->raw());
        flint_free(q_);
        q_ = nullptr;
    }
}

Mfrac Mfrac::clone() const {
    Mfrac out(ctx_);
    fmpz_mpoly_q_set(out.q_, q_, ctx_->raw());
    return out;
}

Mfrac Mfrac::zero(std::shared_ptr<MpolyContext> ctx) {
    return Mfrac(std::move(ctx));
}

Mfrac Mfrac::one(std::shared_ptr<MpolyContext> ctx) {
    Mfrac out(std::move(ctx));
    fmpz_mpoly_q_set_si(out.q_, 1, out.ctx_->raw());
    return out;
}

Mfrac Mfrac::from_mpoly(Mpoly p) {
    Mfrac out(p.ctx());
    fmpz_mpoly_set(fmpz_mpoly_q_numref(out.q_), p.raw(), out.ctx_->raw());
    return out;
}

Mfrac Mfrac::from_si(std::shared_ptr<MpolyContext> ctx, long c) {
    Mfrac out(std::move(ctx));
    fmpz_mpoly_q_set_si(out.q_, c, out.ctx_->raw());
    return out;
}

Mfrac Mfrac::from_fmpq(std::shared_ptr<MpolyContext> ctx, const fmpq_t q) {
    Mfrac out(std::move(ctx));
    fmpz_mpoly_q_set_fmpq(out.q_, q, out.ctx_->raw());
    return out;
}

Mpoly Mfrac::numerator() const {
    Mpoly out(ctx_);
    fmpz_mpoly_set(out.raw(), fmpz_mpoly_q_numref(q_), ctx_->raw());
    return out;
}

Mpoly Mfrac::denominator() const {
    Mpoly out(ctx_);
    fmpz_mpoly_set(out.raw(), fmpz_mpoly_q_denref(q_), ctx_->raw());
    return out;
}

bool Mfrac::is_zero() const noexcept {
    return fmpz_mpoly_is_zero(fmpz_mpoly_q_numref(q_), ctx_->raw()) != 0;
}

bool Mfrac::is_one() const noexcept {
    return fmpz_mpoly_equal(fmpz_mpoly_q_numref(q_), fmpz_mpoly_q_denref(q_),
                             ctx_->raw())
        && !is_zero();
}

bool Mfrac::is_polynomial() const noexcept {
    const fmpz_mpoly_struct* den = fmpz_mpoly_q_denref(q_);
    return fmpz_mpoly_is_fmpz(den, ctx_->raw()) != 0;
}

void Mfrac::check_same_ctx(const Mfrac& other) const {
    if (ctx_.get() != other.ctx_.get()) {
        throw std::invalid_argument("Mfrac: context mismatch");
    }
}

Mfrac Mfrac::operator-() const {
    Mfrac out(ctx_);
    fmpz_mpoly_q_neg(out.q_, q_, ctx_->raw());
    return out;
}

Mfrac& Mfrac::operator+=(const Mfrac& other) {
    check_same_ctx(other);
    fmpz_mpoly_q_add(q_, q_, other.q_, ctx_->raw());
    return *this;
}
Mfrac& Mfrac::operator-=(const Mfrac& other) {
    check_same_ctx(other);
    fmpz_mpoly_q_sub(q_, q_, other.q_, ctx_->raw());
    return *this;
}
Mfrac& Mfrac::operator*=(const Mfrac& other) {
    check_same_ctx(other);
    fmpz_mpoly_q_mul(q_, q_, other.q_, ctx_->raw());
    return *this;
}
Mfrac& Mfrac::operator/=(const Mfrac& other) {
    check_same_ctx(other);
    if (other.is_zero()) throw std::domain_error("Mfrac /= zero");
    fmpz_mpoly_q_div(q_, q_, other.q_, ctx_->raw());
    return *this;
}

Mfrac& Mfrac::multiply_by_mpoly(const Mpoly& p) {
    if (ctx_.get() != p.ctx().get()) {
        throw std::invalid_argument("Mfrac::multiply_by_mpoly: ctx mismatch");
    }
    fmpz_mpoly_mul(fmpz_mpoly_q_numref(q_),
                   fmpz_mpoly_q_numref(q_), p.raw(), ctx_->raw());
    fmpz_mpoly_q_canonicalise(q_, ctx_->raw());
    return *this;
}

Mfrac Mfrac::derivative(long var_idx) const {
    // Quotient rule: d(n/d) = (n' d - n d') / d^2.  Done manually.
    if (var_idx < 0 || var_idx >= ctx_->n_vars()) {
        throw std::out_of_range("Mfrac::derivative: var_idx out of range");
    }
    Mfrac out(ctx_);
    fmpz_mpoly_struct* num     = fmpz_mpoly_q_numref(q_);
    fmpz_mpoly_struct* den     = fmpz_mpoly_q_denref(q_);
    fmpz_mpoly_struct* out_num = fmpz_mpoly_q_numref(out.q_);
    fmpz_mpoly_struct* out_den = fmpz_mpoly_q_denref(out.q_);

    fmpz_mpoly_t nprime, dprime, lhs, rhs;
    fmpz_mpoly_init(nprime, ctx_->raw());
    fmpz_mpoly_init(dprime, ctx_->raw());
    fmpz_mpoly_init(lhs,    ctx_->raw());
    fmpz_mpoly_init(rhs,    ctx_->raw());

    fmpz_mpoly_derivative(nprime, num, var_idx, ctx_->raw());
    fmpz_mpoly_derivative(dprime, den, var_idx, ctx_->raw());
    fmpz_mpoly_mul(lhs, nprime, den, ctx_->raw());
    fmpz_mpoly_mul(rhs, num,    dprime, ctx_->raw());
    fmpz_mpoly_sub(out_num, lhs, rhs, ctx_->raw());
    fmpz_mpoly_mul(out_den, den, den, ctx_->raw());
    fmpz_mpoly_q_canonicalise(out.q_, ctx_->raw());

    fmpz_mpoly_clear(nprime, ctx_->raw());
    fmpz_mpoly_clear(dprime, ctx_->raw());
    fmpz_mpoly_clear(lhs,    ctx_->raw());
    fmpz_mpoly_clear(rhs,    ctx_->raw());
    return out;
}

Mfrac Mfrac::substitute(long var_idx, const fmpq_t value) const {
    if (var_idx < 0 || var_idx >= ctx_->n_vars()) {
        throw std::out_of_range("Mfrac::substitute: var_idx out of range");
    }

    fmpz_t qn, qd;
    fmpz_init(qn);
    fmpz_init(qd);
    fmpz_set(qn, fmpq_numref(value));
    fmpz_set(qd, fmpq_denref(value));

    fmpz_mpoly_struct* num = fmpz_mpoly_q_numref(q_);
    fmpz_mpoly_struct* den = fmpz_mpoly_q_denref(q_);
    long dn = fmpz_mpoly_degree_si(num, var_idx, ctx_->raw());
    long dd = fmpz_mpoly_degree_si(den, var_idx, ctx_->raw());
    if (dn < 0) dn = 0;
    if (dd < 0) dd = 0;
    long d_max = std::max(dn, dd);

    // For each monomial of `in`: if exponent of var_idx is k, multiply by
    // qd^(d_max - k).  After substituting var_idx -> qn, the result equals
    // (in at var_idx = qn/qd) * qd^d_max.  Identical numerator/denominator
    // multipliers cancel during canonicalisation, so the final fraction is
    // unchanged.
    auto scale_var = [&](fmpz_mpoly_t out, const fmpz_mpoly_struct* in, long deg_max) {
        fmpz_mpoly_zero(out, ctx_->raw());
        long len = fmpz_mpoly_length(in, ctx_->raw());
        for (long i = 0; i < len; ++i) {
            ulong* exp = static_cast<ulong*>(
                flint_malloc(sizeof(ulong) * ctx_->n_vars()));
            fmpz_mpoly_get_term_exp_ui(exp, in, i, ctx_->raw());
            long k = static_cast<long>(exp[var_idx]);

            fmpz_t coeff;
            fmpz_init(coeff);
            fmpz_mpoly_get_term_coeff_fmpz(coeff, in, i, ctx_->raw());

            fmpz_t scale;
            fmpz_init(scale);
            fmpz_pow_ui(scale, qd, static_cast<ulong>(deg_max - k));
            fmpz_mul(coeff, coeff, scale);
            fmpz_clear(scale);

            fmpz_mpoly_t term;
            fmpz_mpoly_init(term, ctx_->raw());
            fmpz_mpoly_set_coeff_fmpz_ui(term, coeff, exp, ctx_->raw());
            fmpz_mpoly_add(out, out, term, ctx_->raw());
            fmpz_mpoly_clear(term, ctx_->raw());

            flint_free(exp);
            fmpz_clear(coeff);
        }
    };

    fmpz_mpoly_t num_scaled, den_scaled, num_eval, den_eval;
    fmpz_mpoly_init(num_scaled, ctx_->raw());
    fmpz_mpoly_init(den_scaled, ctx_->raw());
    fmpz_mpoly_init(num_eval, ctx_->raw());
    fmpz_mpoly_init(den_eval, ctx_->raw());

    scale_var(num_scaled, num, d_max);
    scale_var(den_scaled, den, d_max);
    fmpz_mpoly_evaluate_one_fmpz(num_eval, num_scaled, var_idx, qn, ctx_->raw());
    fmpz_mpoly_evaluate_one_fmpz(den_eval, den_scaled, var_idx, qn, ctx_->raw());

    Mfrac out(ctx_);
    fmpz_mpoly_set(fmpz_mpoly_q_numref(out.q_), num_eval, ctx_->raw());
    fmpz_mpoly_set(fmpz_mpoly_q_denref(out.q_), den_eval, ctx_->raw());
    fmpz_mpoly_q_canonicalise(out.q_, ctx_->raw());

    fmpz_mpoly_clear(num_scaled, ctx_->raw());
    fmpz_mpoly_clear(den_scaled, ctx_->raw());
    fmpz_mpoly_clear(num_eval,   ctx_->raw());
    fmpz_mpoly_clear(den_eval,   ctx_->raw());
    fmpz_clear(qn);
    fmpz_clear(qd);
    return out;
}

Mfrac Mfrac::substitute(long var_idx, long value) const {
    fmpq_t v;
    fmpq_init(v);
    fmpq_set_si(v, value, 1);
    Mfrac out = substitute(var_idx, v);
    fmpq_clear(v);
    return out;
}

bool Mfrac::evaluate(fmpq * const * values, fmpq_t out) const {
    Mfrac tmp = clone();
    for (long i = 0; i < ctx_->n_vars(); ++i) {
        tmp = tmp.substitute(i, values[i]);
    }
    if (tmp.is_zero()) {
        fmpq_zero(out);
        return true;
    }
    if (!fmpz_mpoly_is_fmpz(fmpz_mpoly_q_numref(tmp.q_), ctx_->raw())
        || !fmpz_mpoly_is_fmpz(fmpz_mpoly_q_denref(tmp.q_), ctx_->raw())) {
        return false;
    }
    fmpz_t n, d;
    fmpz_init(n); fmpz_init(d);
    if (fmpz_mpoly_length(fmpz_mpoly_q_numref(tmp.q_), ctx_->raw()) > 0) {
        fmpz_mpoly_get_term_coeff_fmpz(n,
            fmpz_mpoly_q_numref(tmp.q_), 0, ctx_->raw());
    } else {
        fmpz_zero(n);
    }
    if (fmpz_mpoly_length(fmpz_mpoly_q_denref(tmp.q_), ctx_->raw()) > 0) {
        fmpz_mpoly_get_term_coeff_fmpz(d,
            fmpz_mpoly_q_denref(tmp.q_), 0, ctx_->raw());
    } else {
        fmpz_one(d);
    }
    if (fmpz_is_zero(d)) {
        fmpz_clear(n); fmpz_clear(d);
        return false;
    }
    fmpq_set_fmpz_frac(out, n, d);
    fmpz_clear(n); fmpz_clear(d);
    return true;
}

void Mfrac::evaluate_acb(acb_ptr out, acb_srcptr values, long prec) const {
    fmpz_mpoly_q_evaluate_acb(out, q_, values, prec, ctx_->raw());
}

std::string Mfrac::to_string() const {
    Mpoly n = numerator();
    Mpoly d = denominator();
    if (d.is_one()) return n.to_string();
    return "(" + n.to_string() + ")/(" + d.to_string() + ")";
}

bool Mfrac::operator==(const Mfrac& other) const {
    if (ctx_.get() != other.ctx_.get()) return false;
    return fmpz_mpoly_q_equal(q_, other.q_, ctx_->raw()) != 0;
}

Mfrac operator+(const Mfrac& a, const Mfrac& b) { Mfrac t = a.clone(); t += b; return t; }
Mfrac operator-(const Mfrac& a, const Mfrac& b) { Mfrac t = a.clone(); t -= b; return t; }
Mfrac operator*(const Mfrac& a, const Mfrac& b) { Mfrac t = a.clone(); t *= b; return t; }
Mfrac operator/(const Mfrac& a, const Mfrac& b) { Mfrac t = a.clone(); t /= b; return t; }

std::ostream& operator<<(std::ostream& os, const Mfrac& f) {
    return os << f.to_string();
}

}  // namespace amflow::algebra
