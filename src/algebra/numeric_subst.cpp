// SPDX-License-Identifier: MIT
// amflow::algebra::numeric_subst implementation.  See header for
// rationale.

#include "amflow/algebra/numeric_subst.hpp"

#include <stdexcept>
#include <vector>

#include <flint/fmpq.h>
#include <flint/fmpz.h>
#include <flint/fmpz_mpoly.h>

namespace amflow::algebra {

namespace {

class FmpqHolder {
public:
    FmpqHolder() { fmpq_init(q_); }
    FmpqHolder(const FmpqHolder& o) { fmpq_init(q_); fmpq_set(q_, o.q_); }
    FmpqHolder(FmpqHolder&& o) noexcept { fmpq_init(q_); fmpq_swap(q_, o.q_); }
    FmpqHolder& operator=(const FmpqHolder& o) {
        if (this != &o) fmpq_set(q_, o.q_);
        return *this;
    }
    FmpqHolder& operator=(FmpqHolder&& o) noexcept {
        if (this != &o) fmpq_swap(q_, o.q_);
        return *this;
    }
    ~FmpqHolder() { fmpq_clear(q_); }
    fmpq_t& raw() { return q_; }
    const fmpq_t& raw() const { return q_; }
private:
    fmpq_t q_;
};

}  // namespace

bool mpoly_has_var(const Mpoly& p, long v) {
    auto ctx = p.ctx();
    long len = fmpz_mpoly_length(p.raw(), ctx->raw());
    std::vector<unsigned long> exp((std::size_t)ctx->n_vars());
    for (long t = 0; t < len; ++t) {
        fmpz_mpoly_get_term_exp_ui(exp.data(), p.raw(), t, ctx->raw());
        if (exp[(std::size_t)v] > 0) return true;
    }
    return false;
}

bool parse_rational_string(const std::string& s, fmpq_t v) {
    // First try flint's native parser, which accepts arbitrary-precision
    // "p/q" and "p" (no decimals).  The previous std::stol path overflowed
    // on >18-digit numerators (e.g. the FD-stencil values in a 3-loop
    // maxcut diffeq driver, h = 10⁻²⁰).
    if (fmpq_set_str(v, s.c_str(), 10) == 0) {
        fmpq_canonicalise(v);
        return true;
    }
    std::size_t slash = s.find('/');
    try {
        if (slash != std::string::npos) {
            long p = std::stol(s.substr(0, slash));
            long q = std::stol(s.substr(slash + 1));
            fmpq_set_si(v, p, q);
            return true;
        }
        std::size_t pos;
        long p = std::stol(s, &pos);
        if (pos == s.size()) {
            fmpq_set_si(v, p, 1);
            return true;
        }
        std::size_t dot = s.find('.');
        if (dot == std::string::npos) return false;
        std::string sint = s.substr(0, dot) + s.substr(dot + 1);
        long denom = 1;
        for (std::size_t i = dot + 1; i < s.size(); ++i) denom *= 10;
        long num = std::stol(sint);
        fmpq_set_si(v, num, denom);
        return true;
    } catch (...) { return false; }
}

Mfrac substitute_fc_vars(
        const Mfrac& src,
        const std::map<std::string, std::string>& numeric_values,
        const std::set<std::string>& keep_names) {
    auto ctx = src.ctx();
    Mfrac cur = src.clone();
    for (long v = 0; v < ctx->n_vars(); ++v) {
        const std::string& name = ctx->var_name(v);
        if (keep_names.count(name)) continue;
        if (!mpoly_has_var(cur.numerator(), v)
                && !mpoly_has_var(cur.denominator(), v)) {
            continue;
        }
        auto it = numeric_values.find(name);
        if (it == numeric_values.end()) {
            throw std::runtime_error(
                "substitute_fc_vars: no numeric value for variable '"
                + name + "' (and it is not in the keep set).  "
                "Caller bug: numeric_values must cover every variable "
                "that is not in keep_names.");
        }
        FmpqHolder h;
        if (!parse_rational_string(it->second, h.raw())) {
            throw std::runtime_error(
                "substitute_fc_vars: cannot parse numeric value '"
                + it->second + "' for variable '" + name + "'");
        }
        cur = cur.substitute(v, h.raw());
    }
    return cur;
}

Mfrac mfrac_to_ctx_lenient(const Mfrac& src,
                            const std::shared_ptr<MpolyContext>& dst_ctx) {
    if (src.ctx().get() == dst_ctx.get()) return src.clone();

    std::map<std::string, long> name_to_dst;
    for (long i = 0; i < dst_ctx->n_vars(); ++i) {
        name_to_dst[dst_ctx->var_name(i)] = i;
    }
    std::vector<long> src_to_dst((std::size_t)src.ctx()->n_vars(), -1);
    for (long i = 0; i < src.ctx()->n_vars(); ++i) {
        auto it = name_to_dst.find(src.ctx()->var_name(i));
        if (it != name_to_dst.end()) {
            src_to_dst[(std::size_t)i] = it->second;
        }
    }

    auto walk = [&](const Mpoly& p) -> Mpoly {
        Mpoly out(dst_ctx);
        long len = fmpz_mpoly_length(p.raw(), src.ctx()->raw());
        std::vector<unsigned long> exp((std::size_t)src.ctx()->n_vars());
        std::vector<unsigned long> dexp((std::size_t)dst_ctx->n_vars(), 0);
        fmpz_t coeff;
        fmpz_init(coeff);
        for (long t = 0; t < len; ++t) {
            fmpz_mpoly_get_term_exp_ui(exp.data(), p.raw(), t,
                                        src.ctx()->raw());
            fmpz_mpoly_get_term_coeff_fmpz(coeff, p.raw(), t,
                                            src.ctx()->raw());
            std::fill(dexp.begin(), dexp.end(), 0);
            for (long i = 0; i < src.ctx()->n_vars(); ++i) {
                long pi = src_to_dst[(std::size_t)i];
                if (pi < 0) {
                    if (exp[(std::size_t)i] != 0) {
                        fmpz_clear(coeff);
                        throw std::runtime_error(
                            "mfrac_to_ctx_lenient: variable '"
                            + src.ctx()->var_name(i)
                            + "' has non-zero exponent but is not in "
                              "destination ctx (caller forgot to "
                              "substitute numeric values?)");
                    }
                    continue;
                }
                dexp[(std::size_t)pi] = exp[(std::size_t)i];
            }
            fmpz_mpoly_set_coeff_fmpz_ui(out.raw(), coeff,
                                          dexp.data(), dst_ctx->raw());
        }
        fmpz_clear(coeff);
        return out;
    };
    return Mfrac(walk(src.numerator()), walk(src.denominator()));
}

Mfrac substitute_and_narrow(
        const Mfrac& src,
        const std::map<std::string, std::string>& numeric_values,
        const std::shared_ptr<MpolyContext>& dst_ctx) {
    std::set<std::string> keep_names;
    for (long i = 0; i < dst_ctx->n_vars(); ++i) {
        keep_names.insert(dst_ctx->var_name(i));
    }
    auto subbed = substitute_fc_vars(src, numeric_values, keep_names);
    return mfrac_to_ctx_lenient(subbed, dst_ctx);
}

}  // namespace amflow::algebra
