// SPDX-License-Identifier: MIT
// qft::boundary — implementation.
//

#include "amflow/qft/boundary.hpp"

#include <algorithm>
#include <cstring>
#include <map>
#include <set>
#include <stdexcept>
#include <utility>

#include <flint/flint.h>
#include <flint/fmpq.h>
#include <flint/fmpz.h>
#include <flint/fmpz_mpoly.h>
#include <flint/fmpz_mpoly_factor.h>
#include <flint/fmpz_mpoly_q.h>

#include "amflow/qft/family_uf.hpp"

namespace amflow::qft {

using algebra::Mfrac;
using algebra::Mpoly;
using algebra::MpolyContext;

namespace {

bool mpoly_is_int(const Mpoly& p) {
    auto ctx = p.ctx();
    return fmpz_mpoly_is_fmpz(p.raw(), ctx->raw()) != 0;
}

void mpoly_get_int(fmpz_t out, const Mpoly& p) {
    auto ctx = p.ctx();
    fmpz_mpoly_get_fmpz(out, p.raw(), ctx->raw());
}

bool mfrac_is_rational(const Mfrac& f, fmpq_t out) {
    Mpoly num = f.numerator();
    Mpoly den = f.denominator();
    if (!mpoly_is_int(num) || !mpoly_is_int(den)) return false;
    fmpz_t n, d;
    fmpz_init(n);
    fmpz_init(d);
    mpoly_get_int(n, num);
    mpoly_get_int(d, den);
    fmpz_set(fmpq_numref(out), n);
    fmpz_set(fmpq_denref(out), d);
    fmpq_canonicalise(out);
    fmpz_clear(n);
    fmpz_clear(d);
    return true;
}

bool mpoly_mentions_var(const Mpoly& p, long var) {
    auto ctx = p.ctx();
    long len = fmpz_mpoly_length(p.raw(), ctx->raw());
    std::vector<unsigned long> exp((std::size_t)ctx->n_vars());
    for (long t = 0; t < len; ++t) {
        fmpz_mpoly_get_term_exp_ui(exp.data(), p.raw(), t, ctx->raw());
        if (exp[(std::size_t)var] > 0) return true;
    }
    return false;
}

long mpoly_max_deg_in(const Mpoly& p, long var) {
    auto ctx = p.ctx();
    long len = fmpz_mpoly_length(p.raw(), ctx->raw());
    std::vector<unsigned long> exp((std::size_t)ctx->n_vars());
    long m = 0;
    for (long t = 0; t < len; ++t) {
        fmpz_mpoly_get_term_exp_ui(exp.data(), p.raw(), t, ctx->raw());
        if ((long)exp[(std::size_t)var] > m) m = (long)exp[(std::size_t)var];
    }
    return m;
}

Mpoly mpoly_pow(const Mpoly& base, long exp) {
    if (exp < 0) {
        throw std::invalid_argument("mpoly_pow: negative exponent");
    }
    auto ctx = base.ctx();
    Mpoly out = Mpoly::one(ctx);
    for (long i = 0; i < exp; ++i) out *= base;
    return out;
}

Mfrac rational_mfrac(const std::shared_ptr<MpolyContext>& ctx,
                     long num, long den) {
    fmpq_t q;
    fmpq_init(q);
    fmpq_set_si(q, num, den);
    fmpq_canonicalise(q);
    Mfrac out = Mfrac::from_fmpq(ctx, q);
    fmpq_clear(q);
    return out;
}

std::vector<Mfrac> multiply_truncated_series(const std::vector<Mfrac>& lhs,
                                             const std::vector<Mfrac>& rhs) {
    if (lhs.empty()) return {};
    if (rhs.empty()) return {};
    auto ctx = lhs[0].ctx();
    std::vector<Mfrac> out;
    out.reserve(lhs.size());
    for (std::size_t i = 0; i < lhs.size(); ++i) {
        out.push_back(Mfrac::zero(ctx));
    }
    for (std::size_t i = 0; i < lhs.size(); ++i) {
        for (std::size_t j = 0; i + j < lhs.size() && j < rhs.size(); ++j) {
            out[i + j] += lhs[i] * rhs[j];
        }
    }
    return out;
}

std::vector<Mfrac> shifted_inverse_power_series(const Mpoly& shift,
                                                long exp,
                                                long max_order) {
    if (exp <= 0) {
        throw std::invalid_argument(
            "shifted_inverse_power_series: exponent must be positive");
    }
    auto ctx = shift.ctx();
    std::vector<Mfrac> coeffs;
    coeffs.reserve((std::size_t)(max_order + 1));

    Mpoly shift_pow = mpoly_pow(shift, exp);
    coeffs.emplace_back(Mpoly::one(ctx), std::move(shift_pow));

    if (max_order == 0) return coeffs;

    Mfrac inv_shift(Mpoly::one(ctx), shift.clone());
    for (long m = 0; m < max_order; ++m) {
        Mfrac next = coeffs[(std::size_t)m].clone();
        next *= rational_mfrac(ctx, -(exp + m), m + 1);
        next *= inv_shift;
        coeffs.push_back(std::move(next));
    }
    return coeffs;
}

int mfrac_compare_max(const Mfrac& a, const Mfrac& b) {
    Mfrac d = a.clone();
    d -= b;
    Mpoly dn = d.numerator();
    Mpoly dd = d.denominator();
    auto ctx = d.ctx();
    if (dn.is_zero()) return 0;

    long n_vars = ctx->n_vars();

    auto eval_q = [&](const Mpoly& p, fmpq_t out) {
        fmpq_t accum, term, eps_pow, base;
        fmpq_init(accum); fmpq_init(term);
        fmpq_init(eps_pow); fmpq_init(base);
        long eps_var = ctx->var_index("eps");
        if (eps_var < 0) eps_var = ctx->var_index("__amf_eps");

        std::vector<unsigned long> exp((std::size_t)n_vars);
        fmpz_t coeff;
        fmpz_init(coeff);
        long len = fmpz_mpoly_length(p.raw(), ctx->raw());
        for (long t = 0; t < len; ++t) {
            fmpz_mpoly_get_term_exp_ui(exp.data(), p.raw(), t, ctx->raw());
            fmpz_mpoly_get_term_coeff_fmpz(coeff, p.raw(), t, ctx->raw());
            bool drop = false;
            for (long v = 0; v < n_vars; ++v) {
                if (v == eps_var) continue;
                if (exp[(std::size_t)v] != 0) { drop = true; break; }
            }
            if (drop) continue;
            fmpq_set_fmpz(term, coeff);
            if (eps_var >= 0 && exp[(std::size_t)eps_var] > 0) {
                fmpq_set_si(eps_pow, 1, 1);
                fmpq_set_si(base,    1, 100);
                for (unsigned long k = 0;
                        k < exp[(std::size_t)eps_var]; ++k) {
                    fmpq_mul(eps_pow, eps_pow, base);
                }
                fmpq_mul(term, term, eps_pow);
            }
            fmpq_add(accum, accum, term);
        }
        fmpz_clear(coeff);
        fmpq_set(out, accum);
        fmpq_clear(accum); fmpq_clear(term);
        fmpq_clear(eps_pow); fmpq_clear(base);
    };

    fmpq_t qn, qd;
    fmpq_init(qn); fmpq_init(qd);
    eval_q(dn, qn);
    eval_q(dd, qd);
    if (fmpq_is_zero(qd)) {
        std::vector<unsigned long> z((std::size_t)n_vars, 0);
        fmpz_t cn, cd;
        fmpz_init(cn); fmpz_init(cd);
        fmpz_mpoly_get_coeff_fmpz_ui(cn, dn.raw(), z.data(), ctx->raw());
        fmpz_mpoly_get_coeff_fmpz_ui(cd, dd.raw(), z.data(), ctx->raw());
        int sgn = fmpz_sgn(cn) * fmpz_sgn(cd);
        fmpz_clear(cn); fmpz_clear(cd);
        fmpq_clear(qn); fmpq_clear(qd);
        return sgn;
    }
    fmpq_div(qn, qn, qd);
    int sgn = fmpq_sgn(qn);
    fmpq_clear(qn); fmpq_clear(qd);
    return sgn;
}

}  // namespace

std::vector<std::vector<Mfrac>>
boundary_pattern(const std::vector<std::vector<Mfrac>>& powers) {
    if (powers.empty()) return {};

    std::size_t n_regions = powers.size();
    std::size_t n_ints    = powers.front().size();
    for (const auto& v : powers) {
        if (v.size() != n_ints) {
            throw std::invalid_argument(
                "boundary_pattern: powers rows have inconsistent length");
        }
    }

    if (n_ints == 0) return {};

    struct Group {
        std::size_t std_region;
        std::vector<std::size_t> regions;
    };
    std::vector<Group> groups;
    for (std::size_t r = 0; r < n_regions; ++r) {
        bool placed = false;
        for (auto& g : groups) {
            Mfrac diff = powers[r][0].clone();
            diff -= powers[g.std_region][0];
            fmpq_t q;
            fmpq_init(q);
            bool ok = mfrac_is_rational(diff, q);
            bool is_int = ok && fmpz_is_one(fmpq_denref(q));
            fmpq_clear(q);
            if (is_int) {
                g.regions.push_back(r);
                placed = true;
                break;
            }
        }
        if (!placed) {
            Group g;
            g.std_region = r;
            g.regions.push_back(r);
            groups.push_back(std::move(g));
        }
    }

    std::vector<std::vector<Mfrac>> result;
    result.reserve(groups.size());
    for (const auto& g : groups) {
        std::vector<Mfrac> per_int;
        per_int.reserve(n_ints);
        const Mfrac& std_value = powers[g.std_region][0];

        for (std::size_t i = 0; i < n_ints; ++i) {
            Mfrac max_off = powers[g.regions.front()][i].clone();
            max_off -= std_value;
            for (std::size_t r_idx = 1; r_idx < g.regions.size(); ++r_idx) {
                std::size_t r = g.regions[r_idx];
                Mfrac off = powers[r][i].clone();
                off -= std_value;
                if (mfrac_compare_max(off, max_off) > 0) {
                    max_off = std::move(off);
                }
            }
            Mfrac p = std_value.clone();
            p += max_off;
            per_int.push_back(std::move(p));
        }
        result.push_back(std::move(per_int));
    }
    return result;
}

std::vector<Mfrac>
apart_one_var(const Mpoly& poly, long var) {
    auto ctx = poly.ctx();
    if (poly.is_zero()) {
        throw std::invalid_argument("apart_one_var: poly is zero");
    }

    fmpz_mpoly_factor_t fac;
    fmpz_mpoly_factor_init(fac, ctx->raw());
    int ok = fmpz_mpoly_factor(fac, poly.raw(), ctx->raw());
    if (!ok) {
        fmpz_mpoly_factor_clear(fac, ctx->raw());
        throw std::runtime_error("apart_one_var: factorization failed");
    }

    fmpz_t constant;
    fmpz_init(constant);
    fmpz_mpoly_factor_get_constant_fmpz(constant, fac, ctx->raw());

    long n_fac = (long)fmpz_mpoly_factor_length(fac, ctx->raw());
    struct VarFactor { Mpoly alpha; Mpoly leading; long exp; };
    std::vector<VarFactor> var_factors;
    Mpoly const_in_var(ctx);
    {
        Mpoly tmp(ctx);
        fmpz_mpoly_set_fmpz(tmp.raw(), constant, ctx->raw());
        const_in_var = std::move(tmp);
    }

    for (long i = 0; i < n_fac; ++i) {
        Mpoly base(ctx);
        fmpz_mpoly_factor_get_base(base.raw(), fac, i, ctx->raw());
        long e = fmpz_mpoly_factor_get_exp_si(fac, i, ctx->raw());

        if (!mpoly_mentions_var(base, var)) {
            Mpoly p = base.clone();
            for (long k = 1; k < e; ++k) p *= base;
            const_in_var *= p;
            continue;
        }

        long deg = mpoly_max_deg_in(base, var);
        if (deg != 1) {
            fmpz_clear(constant);
            fmpz_mpoly_factor_clear(fac, ctx->raw());
            throw std::runtime_error(
                "apart_one_var: factor not linear in var (degree="
                + std::to_string(deg) + ")");
        }
        Mpoly leading = base.coeff_of(var, 1);
        Mpoly rest = base.coeff_of(var, 0);
        if (!mpoly_is_int(leading)) {
            fmpz_clear(constant);
            fmpz_mpoly_factor_clear(fac, ctx->raw());
            throw std::runtime_error(
                "apart_one_var: factor's leading coefficient depends on other variables");
        }
        fmpz_t lead_int;
        fmpz_init(lead_int);
        mpoly_get_int(lead_int, leading);
        if (!fmpz_is_one(lead_int) && !fmpz_equal_si(lead_int, -1)) {
            fmpz_clear(lead_int);
            fmpz_clear(constant);
            fmpz_mpoly_factor_clear(fac, ctx->raw());
            throw std::runtime_error(
                "apart_one_var: factor's leading coefficient is not ±1");
        }
        Mpoly alpha(ctx);
        if (fmpz_is_one(lead_int)) {
            alpha = -rest;
        } else {
            alpha = rest.clone();
        }
        VarFactor vf;
        vf.alpha   = std::move(alpha);
        vf.leading = leading.clone();
        vf.exp     = e;
        var_factors.push_back(std::move(vf));
        fmpz_clear(lead_int);
    }

    fmpz_clear(constant);
    fmpz_mpoly_factor_clear(fac, ctx->raw());

    Mpoly prefactor = const_in_var.clone();
    for (auto& vf : var_factors) {
        prefactor *= mpoly_pow(vf.leading, vf.exp);
    }

    std::vector<Mfrac> result;
    if (var_factors.empty()) {
        Mpoly one_p = Mpoly::one(ctx);
        Mfrac r(std::move(one_p), prefactor.clone());
        result.push_back(std::move(r));
        return result;
    }

    Mpoly x_minus_alpha_base = Mpoly::variable(ctx, var);
    for (std::size_t i = 0; i < var_factors.size(); ++i) {
        long max_order = var_factors[i].exp - 1;
        std::vector<Mfrac> local_series;
        local_series.reserve((std::size_t)(max_order + 1));
        for (long m = 0; m <= max_order; ++m) {
            local_series.push_back(Mfrac::zero(ctx));
        }
        local_series[0] = Mfrac::one(ctx);

        for (std::size_t k = 0; k < var_factors.size(); ++k) {
            if (k == i) continue;

            Mpoly diff = var_factors[i].alpha.clone();
            diff -= var_factors[k].alpha;
            if (diff.is_zero()) {
                throw std::runtime_error(
                    "apart_one_var: repeated root normalization failed");
            }
            auto factor_series =
                shifted_inverse_power_series(diff, var_factors[k].exp, max_order);
            local_series = multiply_truncated_series(local_series, factor_series);
        }

        Mpoly minus_alpha = -var_factors[i].alpha;
        Mpoly x_minus_alpha = x_minus_alpha_base.clone() + minus_alpha;
        for (long j = var_factors[i].exp; j >= 1; --j) {
            long m = var_factors[i].exp - j;
            Mfrac coeff = local_series[(std::size_t)m].clone();
            if (coeff.is_zero()) continue;
            Mpoly denom = prefactor.clone();
            denom *= mpoly_pow(x_minus_alpha, j);
            coeff /= Mfrac::from_mpoly(std::move(denom));
            result.push_back(std::move(coeff));
        }
    }

    return result;
}

std::vector<std::vector<Mfrac>>
apart_rationals(const std::vector<Mfrac>& rationals,
                const DListContext& dctx) {
    std::vector<std::vector<Mfrac>> out;
    out.reserve(rationals.size());

    for (const auto& r : rationals) {
        std::vector<Mfrac> pieces;
        pieces.push_back(r.clone());
        for (long k = 0; k < dctx.n_d; ++k) {
            long var = dctx.first_d_var + k;
            std::vector<Mfrac> next;
            for (auto& p : pieces) {
                Mpoly num = p.numerator();
                Mpoly den = p.denominator();
                auto pfd = apart_one_var(den, var);
                for (auto& fr : pfd) {
                    Mfrac contribution = std::move(fr);
                    contribution *= Mfrac::from_mpoly(num.clone());
                    next.push_back(std::move(contribution));
                }
            }
            pieces = std::move(next);
        }
        out.push_back(std::move(pieces));
    }
    return out;
}

namespace {

struct LaurentSeries {
    std::shared_ptr<MpolyContext> ctx;
    long max_order;
    std::vector<Mfrac> coeff;

    LaurentSeries() : max_order(-1) {}
    LaurentSeries(std::shared_ptr<MpolyContext> c, long N)
        : ctx(std::move(c)), max_order(N), coeff() {
        coeff.reserve((std::size_t)(N + 1));
        for (long i = 0; i <= N; ++i) coeff.push_back(Mfrac::zero(ctx));
    }

    LaurentSeries clone() const {
        LaurentSeries s(ctx, max_order);
        for (long i = 0; i <= max_order; ++i) s.coeff[(std::size_t)i] = coeff[(std::size_t)i].clone();
        return s;
    }

    void scale_by(const Mfrac& f) {
        for (auto& c : coeff) c *= f;
    }
};

LaurentSeries laurent_mul(const LaurentSeries& a, const LaurentSeries& b) {
    if (a.ctx.get() != b.ctx.get())
        throw std::invalid_argument("laurent_mul: ctx mismatch");
    long N = std::min(a.max_order, b.max_order);
    LaurentSeries out(a.ctx, N);
    for (long k = 0; k <= N; ++k) {
        Mfrac acc = Mfrac::zero(a.ctx);
        for (long j = 0; j <= k; ++j) {
            Mfrac term = a.coeff[(std::size_t)j].clone();
            term *= b.coeff[(std::size_t)(k - j)];
            acc += term;
        }
        out.coeff[(std::size_t)k] = std::move(acc);
    }
    return out;
}

LaurentSeries laurent_invert(const LaurentSeries& a) {
    LaurentSeries out(a.ctx, a.max_order);
    if (a.coeff[0].is_zero())
        throw std::invalid_argument("laurent_invert: leading term is zero");
    Mfrac inv_a0 = Mfrac::one(a.ctx);
    inv_a0 /= a.coeff[0];
    out.coeff[0] = inv_a0.clone();
    for (long k = 1; k <= a.max_order; ++k) {
        Mfrac acc = Mfrac::zero(a.ctx);
        for (long j = 1; j <= k; ++j) {
            Mfrac term = a.coeff[(std::size_t)j].clone();
            term *= out.coeff[(std::size_t)(k - j)];
            acc += term;
        }
        Mfrac neg_inv_a0 = -inv_a0;
        acc *= neg_inv_a0;
        out.coeff[(std::size_t)k] = std::move(acc);
    }
    return out;
}

LaurentSeries laurent_pow(const LaurentSeries& a, long p) {
    if (p == 0) {
        LaurentSeries out(a.ctx, a.max_order);
        out.coeff[0] = Mfrac::one(a.ctx);
        return out;
    }
    if (p < 0) {
        return laurent_pow(laurent_invert(a), -p);
    }
    LaurentSeries out = a.clone();
    for (long k = 1; k < p; ++k) out = laurent_mul(out, a);
    return out;
}

Mfrac coeff_half_eta(const Mfrac& f, long half_eta_var, long net_he) {
    Mpoly num = f.numerator();
    Mpoly den = f.denominator();
    auto ctx = num.ctx();
    long n_vars = ctx->n_vars();

    long len_d = fmpz_mpoly_length(den.raw(), ctx->raw());
    long den_he_max = 0;
    long den_he_min = -1;
    std::vector<unsigned long> dexp((std::size_t)n_vars);
    for (long t = 0; t < len_d; ++t) {
        fmpz_mpoly_get_term_exp_ui(dexp.data(), den.raw(), t, ctx->raw());
        long d = (long)dexp[(std::size_t)half_eta_var];
        if (d > den_he_max) den_he_max = d;
        if (den_he_min < 0 || d < den_he_min) den_he_min = d;
    }
    if (den_he_min < 0) den_he_min = 0;
    if (den_he_max != den_he_min) {
        throw std::runtime_error(
            "coeff_half_eta: denominator is not a monomial in half_eta");
    }
    long target_num_he = net_he + den_he_min;
    if (target_num_he < 0) {
        return Mfrac::zero(ctx);
    }

    long len_n = fmpz_mpoly_length(num.raw(), ctx->raw());
    std::vector<unsigned long> nexp((std::size_t)n_vars);
    Mpoly out_num(ctx);
    fmpz_t coeff;
    fmpz_init(coeff);
    for (long t = 0; t < len_n; ++t) {
        fmpz_mpoly_get_term_exp_ui(nexp.data(), num.raw(), t, ctx->raw());
        if ((long)nexp[(std::size_t)half_eta_var] != target_num_he) continue;
        fmpz_mpoly_get_term_coeff_fmpz(coeff, num.raw(), t, ctx->raw());
        std::vector<unsigned long> e2 = nexp;
        e2[(std::size_t)half_eta_var] = 0;
        fmpz_mpoly_set_coeff_fmpz_ui(out_num.raw(), coeff,
                                      e2.data(), ctx->raw());
    }
    fmpz_clear(coeff);

    Mpoly out_den(ctx);
    fmpz_t dcoeff;
    fmpz_init(dcoeff);
    for (long t = 0; t < len_d; ++t) {
        fmpz_mpoly_get_term_exp_ui(dexp.data(), den.raw(), t, ctx->raw());
        fmpz_mpoly_get_term_coeff_fmpz(dcoeff, den.raw(), t, ctx->raw());
        std::vector<unsigned long> e2 = dexp;
        e2[(std::size_t)half_eta_var] = 0;
        fmpz_mpoly_set_coeff_fmpz_ui(out_den.raw(), dcoeff,
                                      e2.data(), ctx->raw());
    }
    fmpz_clear(dcoeff);

    return Mfrac(std::move(out_num), std::move(out_den));
}

Mfrac lift_mfrac_prefix(const Mfrac& src,
                          const std::shared_ptr<MpolyContext>& dst) {
    auto src_ctx = src.ctx();
    if (src_ctx.get() == dst.get()) return src.clone();

    std::vector<Mpoly> gens;
    gens.reserve((std::size_t)src_ctx->n_vars());
    for (long i = 0; i < src_ctx->n_vars(); ++i) {
        gens.push_back(Mpoly::variable(dst, i));
    }
    auto compose = [&](const Mpoly& p) -> Mpoly {
        std::vector<fmpz_mpoly_struct*> ptrs(gens.size());
        for (std::size_t i = 0; i < gens.size(); ++i) {
            ptrs[i] = const_cast<fmpz_mpoly_struct*>(gens[i].raw());
        }
        Mpoly out(dst);
        int ok = fmpz_mpoly_compose_fmpz_mpoly(
            out.raw(), p.raw(), ptrs.data(),
            src_ctx->raw(), dst->raw());
        if (ok == 0) throw std::runtime_error("lift_mfrac_prefix: compose failed");
        return out;
    };
    return Mfrac(compose(src.numerator()), compose(src.denominator()));
}

}  // namespace

BoundaryIntegrandsResult
boundary_integrands(const FamilyConfig&             fc,
                    const RegionContext&            rctx,
                    const std::vector<JIntegral>&   integrals,
                    const std::vector<long>&        border,
                    const Region&                   region) {
    if (integrals.empty()) {
        BoundaryIntegrandsResult r;
        r.dctx = nullptr;
        r.first_d_var = 0;
        return r;
    }
    if (border.size() != integrals.size()) {
        throw std::invalid_argument(
            "boundary_integrands: border length != integrals length");
    }

    auto rule = region_rule(rctx, region.transform, region.scale);
    long fc_eta_var = fc.ctx->var_index("eta");

    std::size_t N = fc.propagators_after_conservation.size();
    std::vector<Mfrac> fullde;
    fullde.reserve(N);
    for (std::size_t i = 0; i < N; ++i) {
        Mfrac raw = apply_region_rule(
            fc, rctx, fc.propagators_after_conservation[i], rule);
        if (fc_eta_var >= 0) {
            auto sub_eta = [&](const Mpoly& p) -> Mpoly {
                Mpoly out(rctx.ctx);
                long len = fmpz_mpoly_length(p.raw(), rctx.ctx->raw());
                std::vector<unsigned long> exp((std::size_t)rctx.ctx->n_vars());
                fmpz_t coeff;
                fmpz_init(coeff);
                for (long t = 0; t < len; ++t) {
                    fmpz_mpoly_get_term_exp_ui(exp.data(), p.raw(), t,
                                                rctx.ctx->raw());
                    fmpz_mpoly_get_term_coeff_fmpz(coeff, p.raw(), t,
                                                    rctx.ctx->raw());
                    long ek = (long)exp[(std::size_t)fc_eta_var];
                    std::vector<unsigned long> e2 = exp;
                    e2[(std::size_t)fc_eta_var] = 0;
                    e2[(std::size_t)rctx.half_eta_var] += (unsigned long)(2 * ek);
                    fmpz_mpoly_set_coeff_fmpz_ui(out.raw(), coeff,
                                                  e2.data(), rctx.ctx->raw());
                }
                fmpz_clear(coeff);
                return out;
            };
            Mpoly num = sub_eta(raw.numerator());
            Mpoly den = sub_eta(raw.denominator());
            raw = Mfrac(std::move(num), std::move(den));
        }
        Mfrac reduced = apply_uf_replacement(fc, rctx.ctx, raw);
        fullde.push_back(std::move(reduced));
    }

    std::vector<long> factor_he(N, 0);
    std::vector<Mfrac> stripped;
    stripped.reserve(N);
    for (std::size_t i = 0; i < N; ++i) {
        Mpoly num = fullde[i].numerator();
        long len = fmpz_mpoly_length(num.raw(), rctx.ctx->raw());
        // Linear (worldline) propagator fix: a propagator linear in the
        // loop momentum (e.g. worldline 1/(2u·k)) scales as he^1, not he^2,
        // under l -> he*l.  The previous `has_he ? 2 : 0` hardcoded the
        // quadratic case and made the leading Laurent coefficient zero,
        // tripping `laurent_invert`.  Track the ACTUAL max half_eta degree.
        long max_he = 0;
        std::vector<unsigned long> exp((std::size_t)rctx.ctx->n_vars());
        for (long t = 0; t < len; ++t) {
            fmpz_mpoly_get_term_exp_ui(exp.data(), num.raw(), t,
                                        rctx.ctx->raw());
            long e = (long)exp[(std::size_t)rctx.half_eta_var];
            if (e > max_he) max_he = e;
        }
        long fac_he = max_he;
        factor_he[i] = fac_he;

        if (fac_he == 0) {
            stripped.push_back(fullde[i].clone());
        } else {
            Mpoly he = Mpoly::variable(rctx.ctx, rctx.half_eta_var);
            Mpoly he_pow = he.clone();
            for (long k = 1; k < fac_he; ++k) he_pow = he_pow * he;
            Mpoly new_den = fullde[i].denominator() * he_pow;
            stripped.push_back(Mfrac(fullde[i].numerator(),
                                       std::move(new_den)));
        }
    }

    std::vector<Mfrac> expde;
    expde.reserve(N);
    for (std::size_t i = 0; i < N; ++i) {
        expde.push_back(coeff_half_eta(stripped[i],
                                          rctx.half_eta_var,
                                          0));
    }

    std::vector<Mpoly> expde_fc;
    expde_fc.reserve(N);
    for (std::size_t i = 0; i < N; ++i) {
        Mpoly num = expde[i].numerator();
        Mpoly den = expde[i].denominator();
        if (mpoly_mentions_var(den, rctx.eta_var) ||
            mpoly_mentions_var(den, rctx.half_eta_var)) {
            throw std::runtime_error(
                "boundary_integrands: expde denominator depends on eta");
        }
        auto project = [&](const Mpoly& p) -> Mpoly {
            Mpoly out(fc.ctx);
            long len = fmpz_mpoly_length(p.raw(), rctx.ctx->raw());
            std::vector<unsigned long> rexp((std::size_t)rctx.ctx->n_vars());
            std::vector<unsigned long> fexp((std::size_t)fc.ctx->n_vars(), 0);
            fmpz_t coeff;
            fmpz_init(coeff);
            for (long t = 0; t < len; ++t) {
                fmpz_mpoly_get_term_exp_ui(rexp.data(), p.raw(), t,
                                            rctx.ctx->raw());
                fmpz_mpoly_get_term_coeff_fmpz(coeff, p.raw(), t,
                                                rctx.ctx->raw());
                for (long j = 0; j < fc.ctx->n_vars(); ++j) {
                    fexp[(std::size_t)j] = rexp[(std::size_t)j];
                }
                fmpz_mpoly_set_coeff_fmpz_ui(out.raw(), coeff,
                                              fexp.data(), fc.ctx->raw());
            }
            fmpz_clear(coeff);
            return out;
        };
        Mpoly num_fc = project(num);
        Mpoly den_fc = project(den);
        if (!mpoly_is_int(den_fc)) {
            throw std::runtime_error(
                "boundary_integrands: expde denominator is not constant after projection");
        }
        expde_fc.push_back(std::move(num_fc));
    }

    std::vector<std::size_t> top_positions;
    for (std::size_t i = 0; i < N; ++i) {
        bool any_pos = false;
        for (const auto& integ : integrals) {
            if (integ.indices()[i] > 0) { any_pos = true; break; }
        }
        if (any_pos) top_positions.push_back(i);
    }

    std::vector<Mpoly> top_expde;
    top_expde.reserve(top_positions.size());
    for (std::size_t i : top_positions) {
        top_expde.push_back(expde_fc[i].clone());
    }

    auto completede_raw = to_complete_explicit(fc, top_expde);
    auto dctx = make_dlist_context(fc, completede_raw.size());
    auto sptoD = sp_list_to_dlist_symbol(fc, dctx, completede_raw);

    std::vector<std::string> merged_names;
    merged_names.reserve((std::size_t)dctx.ctx->n_vars() + 2);
    for (long i = 0; i < dctx.ctx->n_vars(); ++i) {
        merged_names.push_back(dctx.ctx->var_name(i));
    }
    merged_names.emplace_back("__amf_eta");
    merged_names.emplace_back("__amf_half_eta");
    auto mctx = std::make_shared<MpolyContext>(std::move(merged_names));
    long mctx_eta_var      = dctx.ctx->n_vars();
    long mctx_half_eta_var = dctx.ctx->n_vars() + 1;

    std::vector<Mfrac> sptoD_in_mctx;
    sptoD_in_mctx.reserve(sptoD.size());
    for (const auto& s : sptoD) sptoD_in_mctx.push_back(lift_mfrac_prefix(s, mctx));

    auto rctx_to_mctx = [&](const Mfrac& src) -> Mfrac {
        std::vector<Mpoly> gens;
        gens.reserve((std::size_t)rctx.ctx->n_vars());
        long ng = rctx.ctx->n_vars();
        for (long i = 0; i < ng; ++i) {
            if (i < fc.ctx->n_vars()) {
                gens.push_back(Mpoly::variable(mctx, i));
            } else if (i == rctx.eta_var) {
                gens.push_back(Mpoly::variable(mctx, mctx_eta_var));
            } else if (i == rctx.half_eta_var) {
                gens.push_back(Mpoly::variable(mctx, mctx_half_eta_var));
            } else {
                throw std::runtime_error("rctx_to_mctx: unexpected rctx var");
            }
        }
        auto compose = [&](const Mpoly& p) -> Mpoly {
            std::vector<fmpz_mpoly_struct*> ptrs(gens.size());
            for (std::size_t i = 0; i < gens.size(); ++i) {
                ptrs[i] = const_cast<fmpz_mpoly_struct*>(gens[i].raw());
            }
            Mpoly out(mctx);
            int ok = fmpz_mpoly_compose_fmpz_mpoly(
                out.raw(), p.raw(), ptrs.data(),
                rctx.ctx->raw(), mctx->raw());
            if (ok == 0) throw std::runtime_error("rctx_to_mctx: compose failed");
            return out;
        };
        return Mfrac(compose(src.numerator()), compose(src.denominator()));
    };

    auto substitute_splist_in_mfrac = [&](const Mfrac& f) -> Mfrac {
        auto walk = [&](const Mpoly& p) -> Mfrac {
            Mfrac acc = Mfrac::zero(mctx);
            long len = fmpz_mpoly_length(p.raw(), mctx->raw());
            std::vector<unsigned long> exp((std::size_t)mctx->n_vars());
            fmpz_t coeff;
            fmpz_init(coeff);
            long n_loop = (long)fc.n_loops();
            long n_red  = (long)fc.n_red_legs();
            for (long t = 0; t < len; ++t) {
                fmpz_mpoly_get_term_exp_ui(exp.data(), p.raw(), t,
                                            mctx->raw());
                fmpz_mpoly_get_term_coeff_fmpz(coeff, p.raw(), t,
                                                mctx->raw());

                std::vector<long> loop_slots;
                std::vector<long> leg_slots;
                for (long i = 0; i < n_loop; ++i) {
                    long e = (long)exp[(std::size_t)i];
                    for (long k = 0; k < e; ++k) loop_slots.push_back(i);
                }
                for (long i = 0; i < n_red; ++i) {
                    long e = (long)exp[(std::size_t)(n_loop + i)];
                    for (long k = 0; k < e; ++k) leg_slots.push_back(i);
                }

                long n_loop_total = (long)loop_slots.size();
                long n_leg_total  = (long)leg_slots.size();
                if (n_leg_total > n_loop_total) {
                    throw std::runtime_error(
                        "substitute_splist: malformed monomial (more legs than loops)");
                }
                long n_loop_loop_pairs = (n_loop_total - n_leg_total) / 2;
                long n_loop_leg_pairs  = n_leg_total;
                if ((n_loop_total - n_leg_total) % 2 != 0) {
                    throw std::runtime_error(
                        "substitute_splist: monomial cannot be paired");
                }

                Mfrac mfac = Mfrac::one(mctx);

                for (long pair = 0; pair < n_loop_loop_pairs; ++pair) {
                    long a = loop_slots[(std::size_t)(2 * pair)];
                    long b = loop_slots[(std::size_t)(2 * pair + 1)];
                    long ai = std::min(a, b), bi = std::max(a, b);
                    long idx = -1;
                    for (long k = 0; k < (long)fc.sp_list.size(); ++k) {
                        std::vector<unsigned long> sp_exp((std::size_t)fc.ctx->n_vars(), 0);
                        fmpz_mpoly_get_term_exp_ui(
                            sp_exp.data(), fc.sp_list[(std::size_t)k].raw(),
                            0, fc.ctx->raw());
                        bool match = (long)sp_exp[(std::size_t)ai] >= 1
                                  && (long)sp_exp[(std::size_t)bi] >= 1;
                        if (ai == bi) match = (long)sp_exp[(std::size_t)ai] >= 2;
                        if (match) { idx = k; break; }
                    }
                    if (idx < 0) {
                        throw std::runtime_error(
                            "substitute_splist: cannot find loop-loop sp_list entry");
                    }
                    mfac *= sptoD_in_mctx[(std::size_t)idx];
                }
                long loop_offset = 2 * n_loop_loop_pairs;
                for (long lp = 0; lp < n_loop_leg_pairs; ++lp) {
                    long a = loop_slots[(std::size_t)(loop_offset + lp)];
                    long b = leg_slots[(std::size_t)lp];
                    long idx = -1;
                    for (long k = 0; k < (long)fc.sp_list.size(); ++k) {
                        std::vector<unsigned long> sp_exp((std::size_t)fc.ctx->n_vars(), 0);
                        fmpz_mpoly_get_term_exp_ui(
                            sp_exp.data(), fc.sp_list[(std::size_t)k].raw(),
                            0, fc.ctx->raw());
                        bool match = (long)sp_exp[(std::size_t)a] >= 1
                                  && (long)sp_exp[(std::size_t)(n_loop + b)] >= 1;
                        if (match) { idx = k; break; }
                    }
                    if (idx < 0) {
                        throw std::runtime_error(
                            "substitute_splist: cannot find loop-leg sp_list entry");
                    }
                    mfac *= sptoD_in_mctx[(std::size_t)idx];
                }

                std::vector<unsigned long> rest_exp((std::size_t)mctx->n_vars(), 0);
                for (long i = n_loop + n_red; i < mctx->n_vars(); ++i) {
                    rest_exp[(std::size_t)i] = exp[(std::size_t)i];
                }
                Mpoly rest_mono = Mpoly::zero(mctx);
                fmpz_mpoly_set_coeff_fmpz_ui(rest_mono.raw(), coeff,
                                              rest_exp.data(), mctx->raw());
                mfac *= Mfrac::from_mpoly(std::move(rest_mono));

                acc += mfac;
            }
            fmpz_clear(coeff);
            return acc;
        };

        Mfrac new_num = walk(f.numerator());
        Mfrac new_den = walk(f.denominator());
        if (new_den.is_zero()) {
            throw std::runtime_error(
                "substitute_splist: denominator vanished after substitution");
        }
        return Mfrac(new_num.numerator() * new_den.denominator(),
                      new_den.numerator() * new_num.denominator());
    };

    std::vector<Mfrac> stripped_in_dctx;
    stripped_in_dctx.reserve(N);
    for (std::size_t i = 0; i < N; ++i) {
        Mfrac in_mctx = rctx_to_mctx(stripped[i]);
        Mfrac in_dctx_mctx = substitute_splist_in_mfrac(in_mctx);
        stripped_in_dctx.push_back(std::move(in_dctx_mctx));
    }

    long max_border = 0;
    for (long b : border) if (b > max_border) max_border = b;

    long max_order = 2 * max_border;

    std::vector<LaurentSeries> base;
    base.reserve(N);
    for (std::size_t i = 0; i < N; ++i) {
        LaurentSeries s(mctx, max_order);
        // Region engine: fill ALL base Laurent
        // coefficients down to he^-max_order instead of the previous
        // hardcoded he^0..he^-2 window.  For the sqrt(eta) regions this is
        // bit-identical (quadratic propagators strip he^2, so coefficients
        // beyond he^-2 are structurally zero and coeff_half_eta returns
        // the same zero the constructor already placed); for the l ~ eta
        // regions (AMFLOW_EIKONAL_REGIONS=1) quadratics strip he^4 and
        // mixed regions produce he^-1/he^-3 terms, all captured here.
        for (long j = 0; j <= max_order; ++j) {
            s.coeff[(std::size_t)j] = coeff_half_eta(stripped_in_dctx[i],
                                         mctx_half_eta_var, -j);
        }
        base.push_back(std::move(s));
    }

    BoundaryIntegrandsResult result;
    result.completede = squared_denominators(fc, completede_raw);
    result.dctx       = mctx;
    result.first_d_var = dctx.first_d_var;

    result.integrands.reserve(integrals.size());
    for (std::size_t k = 0; k < integrals.size(); ++k) {
        const auto& integ = integrals[k];
        if ((long)integ.indices().size() != (long)N) {
            throw std::invalid_argument(
                "boundary_integrands: integral indices length mismatch");
        }
        LaurentSeries acc(mctx, max_order);
        acc.coeff[0] = Mfrac::one(mctx);
        for (std::size_t i = 0; i < N; ++i) {
            long p = integ.indices()[i];
            if (p == 0) continue;
            LaurentSeries factor = laurent_pow(base[i], -p);
            acc = laurent_mul(acc, factor);
        }
        std::vector<Mfrac> terms;
        long b = border[k];
        if (b >= 0) {
            terms.reserve((std::size_t)(b + 1));
            for (long j = 0; j <= b; ++j) {
                long order_x = 2 * j;
                if (order_x > max_order) break;
                terms.push_back(acc.coeff[(std::size_t)order_x].clone());
            }
        }
        result.integrands.push_back(std::move(terms));
    }

    return result;
}

std::vector<LaportaTerm>
laporta_integrals(const Mfrac& term, const DListContext& dctx) {
    Mpoly num = term.numerator();
    Mpoly den = term.denominator();
    auto ctx = num.ctx();

    long len_den = fmpz_mpoly_length(den.raw(), ctx->raw());
    if (len_den != 1) {
        throw std::runtime_error(
            "laporta_integrals: denominator is not a monomial of DListSymbol");
    }

    std::vector<unsigned long> den_exp((std::size_t)ctx->n_vars());
    fmpz_mpoly_get_term_exp_ui(den_exp.data(), den.raw(),
                                0, ctx->raw());
    fmpz_t den_coeff;
    fmpz_init(den_coeff);
    fmpz_mpoly_get_term_coeff_fmpz(den_coeff, den.raw(),
                                    0, ctx->raw());

    long first_d = dctx.first_d_var;
    long n_d     = dctx.n_d;

    std::vector<long> den_d_exp((std::size_t)n_d, 0);
    for (long k = 0; k < n_d; ++k) {
        den_d_exp[(std::size_t)k] = (long)den_exp[(std::size_t)(first_d + k)];
    }

    long len_num = fmpz_mpoly_length(num.raw(), ctx->raw());
    std::vector<LaportaTerm> out;
    std::vector<unsigned long> num_exp((std::size_t)ctx->n_vars());
    fmpz_t num_coeff;
    fmpz_init(num_coeff);
    for (long t = 0; t < len_num; ++t) {
        fmpz_mpoly_get_term_exp_ui(num_exp.data(), num.raw(), t,
                                    ctx->raw());
        fmpz_mpoly_get_term_coeff_fmpz(num_coeff, num.raw(), t,
                                        ctx->raw());

        std::vector<long> indices((std::size_t)n_d);
        for (long k = 0; k < n_d; ++k) {
            indices[(std::size_t)k] = den_d_exp[(std::size_t)k]
                - (long)num_exp[(std::size_t)(first_d + k)];
        }
        std::vector<unsigned long> rest_exp = num_exp;
        for (long k = 0; k < n_d; ++k) rest_exp[(std::size_t)(first_d + k)] = 0;
        Mpoly rest_mono = Mpoly::zero(ctx);
        fmpz_mpoly_set_coeff_fmpz_ui(rest_mono.raw(), num_coeff,
                                      rest_exp.data(), ctx->raw());

        Mfrac coef_mfrac = Mfrac::from_mpoly(std::move(rest_mono));
        Mpoly den_const(ctx);
        fmpz_mpoly_set_fmpz(den_const.raw(), den_coeff, ctx->raw());
        coef_mfrac /= Mfrac::from_mpoly(std::move(den_const));

        LaportaTerm lt;
        lt.indices = std::move(indices);
        lt.coeff   = std::move(coef_mfrac);
        out.push_back(std::move(lt));
    }
    fmpz_clear(den_coeff);
    fmpz_clear(num_coeff);
    return out;
}

std::vector<BoundaryFamily>
boundary_integrals(const FamilyConfig& fc,
                   const RegionContext& /*rctx*/,
                   const BoundaryIntegrandsResult& bi) {
    if (bi.integrands.empty()) return {};

    DListContext dctx;
    dctx.ctx          = bi.dctx;
    dctx.first_d_var  = bi.first_d_var;
    dctx.n_d          = (long)bi.completede.size();

    std::vector<std::size_t> partition;
    partition.reserve(bi.integrands.size());
    std::vector<Mfrac> flat;
    for (const auto& v : bi.integrands) {
        partition.push_back(v.size());
        for (const auto& f : v) flat.push_back(f.clone());
    }

    auto pfd_per_flat = apart_rationals(flat, dctx);

    auto analyze_family = [&](const Mfrac& term) -> std::vector<Mpoly> {
        Mpoly den = term.denominator();
        auto ctx = den.ctx();
        std::vector<Mpoly> sig;
        sig.reserve((std::size_t)dctx.n_d);
        for (long k = 0; k < dctx.n_d; ++k) sig.emplace_back(ctx);

        fmpz_mpoly_factor_t fac;
        fmpz_mpoly_factor_init(fac, ctx->raw());
        int ok = fmpz_mpoly_factor(fac, den.raw(), ctx->raw());
        if (!ok) {
            fmpz_mpoly_factor_clear(fac, ctx->raw());
            throw std::runtime_error("boundary_integrals: factor failed");
        }
        long n_fac = (long)fmpz_mpoly_factor_length(fac, ctx->raw());
        for (long i = 0; i < n_fac; ++i) {
            Mpoly base(ctx);
            fmpz_mpoly_factor_get_base(base.raw(), fac, i, ctx->raw());

            long which_d = -1;
            for (long k = 0; k < dctx.n_d; ++k) {
                long var = dctx.first_d_var + k;
                long len_b = fmpz_mpoly_length(base.raw(), ctx->raw());
                std::vector<unsigned long> exp((std::size_t)ctx->n_vars());
                bool depends = false;
                for (long t = 0; t < len_b && !depends; ++t) {
                    fmpz_mpoly_get_term_exp_ui(exp.data(), base.raw(),
                                                t, ctx->raw());
                    if (exp[(std::size_t)var] > 0) depends = true;
                }
                if (depends) {
                    if (which_d != -1) {
                        fmpz_mpoly_factor_clear(fac, ctx->raw());
                        throw std::runtime_error(
                            "boundary_integrals: factor depends on multiple D variables");
                    }
                    which_d = k;
                }
            }
            if (which_d < 0) continue;

            long var = dctx.first_d_var + which_d;
            Mpoly leading = base.coeff_of(var, 1);
            Mpoly rest    = base.coeff_of(var, 0);
            if (!fmpz_mpoly_is_fmpz(leading.raw(), ctx->raw())) {
                fmpz_mpoly_factor_clear(fac, ctx->raw());
                throw std::runtime_error(
                    "boundary_integrals: D-factor leading coefficient not constant");
            }
            fmpz_t lead_int;
            fmpz_init(lead_int);
            fmpz_mpoly_get_fmpz(lead_int, leading.raw(), ctx->raw());
            Mpoly alpha(ctx);
            if (fmpz_is_one(lead_int)) {
                alpha = -rest;
            } else if (fmpz_equal_si(lead_int, -1)) {
                alpha = rest.clone();
            } else {
                fmpz_clear(lead_int);
                fmpz_mpoly_factor_clear(fac, ctx->raw());
                throw std::runtime_error(
                    "boundary_integrals: D-factor leading coefficient not ±1");
            }
            fmpz_clear(lead_int);

            if (!sig[(std::size_t)which_d].is_zero()
                    && !(sig[(std::size_t)which_d] - alpha).is_zero()) {
                fmpz_mpoly_factor_clear(fac, ctx->raw());
                throw std::runtime_error(
                    "boundary_integrals: inconsistent D-factor offsets");
            }
            sig[(std::size_t)which_d] = std::move(alpha);
        }
        fmpz_mpoly_factor_clear(fac, ctx->raw());
        return sig;
    };

    auto shift_d_vars = [&](const Mfrac& src,
                              const std::vector<Mpoly>& alpha) -> Mfrac {
        std::vector<Mpoly> gens;
        gens.reserve((std::size_t)dctx.ctx->n_vars());
        for (long i = 0; i < dctx.ctx->n_vars(); ++i) {
            Mpoly g = Mpoly::variable(dctx.ctx, i);
            for (long k = 0; k < dctx.n_d; ++k) {
                if (i == dctx.first_d_var + k
                        && !alpha[(std::size_t)k].is_zero()) {
                    g += alpha[(std::size_t)k];
                }
            }
            gens.push_back(std::move(g));
        }
        std::vector<fmpz_mpoly_struct*> ptrs(gens.size());
        for (std::size_t i = 0; i < gens.size(); ++i) {
            ptrs[i] = const_cast<fmpz_mpoly_struct*>(gens[i].raw());
        }
        Mpoly num_out(dctx.ctx);
        Mpoly den_out(dctx.ctx);
        int rc1 = fmpz_mpoly_compose_fmpz_mpoly(
            num_out.raw(), src.numerator().raw(), ptrs.data(),
            dctx.ctx->raw(), dctx.ctx->raw());
        int rc2 = fmpz_mpoly_compose_fmpz_mpoly(
            den_out.raw(), src.denominator().raw(), ptrs.data(),
            dctx.ctx->raw(), dctx.ctx->raw());
        if (!rc1 || !rc2) {
            throw std::runtime_error("shift_d_vars: compose failed");
        }
        return Mfrac(std::move(num_out), std::move(den_out));
    };

    auto sig_key = [](const std::vector<Mpoly>& sig) -> std::string {
        std::string k;
        for (const auto& a : sig) { k += a.to_string(); k += "|"; }
        return k;
    };

    auto alpha_to_fc = [&](const Mpoly& a) -> Mpoly {
        Mpoly out(fc.ctx);
        long len = fmpz_mpoly_length(a.raw(), dctx.ctx->raw());
        std::vector<unsigned long> dexp((std::size_t)dctx.ctx->n_vars());
        std::vector<unsigned long> fexp((std::size_t)fc.ctx->n_vars(), 0);
        fmpz_t coeff;
        fmpz_init(coeff);
        for (long t = 0; t < len; ++t) {
            fmpz_mpoly_get_term_exp_ui(dexp.data(), a.raw(), t,
                                        dctx.ctx->raw());
            fmpz_mpoly_get_term_coeff_fmpz(coeff, a.raw(), t,
                                            dctx.ctx->raw());
            for (long k = 0; k < dctx.n_d; ++k) {
                if (dexp[(std::size_t)(dctx.first_d_var + k)] > 0) {
                    fmpz_clear(coeff);
                    throw std::runtime_error(
                        "alpha_to_fc: alpha unexpectedly depends on D");
                }
            }
            for (long j = 0; j < fc.ctx->n_vars(); ++j) {
                fexp[(std::size_t)j] = dexp[(std::size_t)j];
            }
            fmpz_mpoly_set_coeff_fmpz_ui(out.raw(), coeff,
                                          fexp.data(), fc.ctx->raw());
        }
        fmpz_clear(coeff);
        return out;
    };

    struct SigSummedEntry {
        std::vector<Mpoly> sig;
        std::string        sig_k;
        Mfrac              summed;
    };
    std::vector<std::vector<SigSummedEntry>> flat_grouped(pfd_per_flat.size());
    for (std::size_t f = 0; f < pfd_per_flat.size(); ++f) {
        std::map<std::string, std::size_t> key_to_idx;
        for (std::size_t j = 0; j < pfd_per_flat[f].size(); ++j) {
            auto sig = analyze_family(pfd_per_flat[f][j]);
            std::string key = sig_key(sig);
            auto it = key_to_idx.find(key);
            if (it == key_to_idx.end()) {
                SigSummedEntry e{std::move(sig), key,
                                  pfd_per_flat[f][j].clone()};
                key_to_idx[key] = flat_grouped[f].size();
                flat_grouped[f].push_back(std::move(e));
            } else {
                flat_grouped[f][it->second].summed += pfd_per_flat[f][j];
            }
        }
    }

    std::map<std::string, std::size_t> family_index;
    std::vector<BoundaryFamily> out;

    auto ensure_family = [&](const std::vector<Mpoly>& sig,
                              const std::string&       key) -> BoundaryFamily& {
        auto it = family_index.find(key);
        if (it != family_index.end()) return out[it->second];
        BoundaryFamily fam;
        fam.prop.reserve(bi.completede.size());
        for (std::size_t k = 0; k < bi.completede.size(); ++k) {
            Mpoly p = bi.completede[k].clone();
            if (k < sig.size() && !sig[k].is_zero()) {
                Mpoly a_fc = alpha_to_fc(sig[k]);
                p -= a_fc;
            }
            fam.prop.push_back(std::move(p));
        }
        fam.terms.resize(bi.integrands.size());
        for (std::size_t i = 0; i < bi.integrands.size(); ++i) {
            fam.terms[i].resize(bi.integrands[i].size());
        }
        family_index[key] = out.size();
        out.push_back(std::move(fam));
        return out.back();
    };

    std::size_t flat_idx = 0;
    for (std::size_t input_i = 0; input_i < bi.integrands.size(); ++input_i) {
        for (std::size_t order_j = 0;
                 order_j < bi.integrands[input_i].size();
                 ++order_j) {
            for (auto& entry : flat_grouped[flat_idx]) {
                ensure_family(entry.sig, entry.sig_k);
                Mfrac shifted = shift_d_vars(entry.summed, entry.sig);
                auto lt_vec = laporta_integrals(shifted, dctx);
                std::size_t fam_idx = family_index[entry.sig_k];
                out[fam_idx].terms[input_i][order_j] = std::move(lt_vec);
            }
            ++flat_idx;
        }
    }

    return out;
}

}  // namespace amflow::qft
