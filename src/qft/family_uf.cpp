// SPDX-License-Identifier: MIT
// qft::family_uf — implementation.
//

#include "amflow/qft/family_uf.hpp"

#include <stdexcept>
#include <utility>

#include <flint/fmpz.h>
#include <flint/fmpz_mpoly.h>
#include <flint/fmpz_mpoly_q.h>

namespace amflow::qft {

using algebra::Mfrac;
using algebra::Mpoly;
using algebra::MpolyContext;
using algebra::MpolyMatrix;

namespace {

std::shared_ptr<MpolyContext> build_uf_ctx(const FamilyConfig& fc,
                                            std::size_t n_x) {
    std::vector<std::string> names;
    names.reserve(static_cast<std::size_t>(fc.ctx->n_vars()) + n_x);
    for (long i = 0; i < fc.ctx->n_vars(); ++i) {
        names.push_back(fc.ctx->var_name(i));
    }
    for (std::size_t i = 0; i < n_x; ++i) {
        names.push_back("__feyn_x" + std::to_string(i + 1));
    }
    return std::make_shared<MpolyContext>(std::move(names));
}

Mpoly compose_to(const Mpoly& src,
                 const std::shared_ptr<MpolyContext>& dst_ctx,
                 const std::vector<Mpoly>& gen_map) {
    auto src_ctx = src.ctx();
    if (static_cast<long>(gen_map.size()) != src_ctx->n_vars()) {
        throw std::invalid_argument("compose_to: gen_map size mismatch");
    }
    std::vector<fmpz_mpoly_struct*> ptrs(gen_map.size());
    for (std::size_t i = 0; i < gen_map.size(); ++i) {
        if (gen_map[i].ctx().get() != dst_ctx.get()) {
            throw std::invalid_argument("compose_to: gen ctx mismatch");
        }
        ptrs[i] = const_cast<fmpz_mpoly_struct*>(gen_map[i].raw());
    }
    Mpoly out(dst_ctx);
    int ok = fmpz_mpoly_compose_fmpz_mpoly(
        out.raw(), src.raw(), ptrs.data(),
        src_ctx->raw(), dst_ctx->raw());
    if (ok == 0) {
        throw std::runtime_error("compose_to: FLINT compose failed");
    }
    return out;
}

std::vector<Mpoly>
make_lift_gen_map(const std::shared_ptr<MpolyContext>& uf_ctx,
                  long n_src_vars) {
    std::vector<Mpoly> gens;
    gens.reserve(static_cast<std::size_t>(n_src_vars));
    for (long i = 0; i < n_src_vars; ++i) {
        gens.push_back(Mpoly::variable(uf_ctx, i));
    }
    return gens;
}

Mpoly bilinear_coef(const Mpoly& poly, long var0, long var1) {
    if (var0 == var1) {
        return poly.coeff_of(var0, 2);
    }
    return poly.coeff_of(std::vector<long>{var0, var1},
                         std::vector<unsigned long>{1, 1});
}

}  // namespace

Mfrac apply_uf_replacement(const FamilyConfig& fc,
                            const std::shared_ptr<MpolyContext>& uf_ctx,
                            const Mfrac& p) {
    auto local_leg_pair_key = [](const std::string& a, const std::string& b)
        -> std::string {
        if (a <= b) return a + "*" + b;
        return b + "*" + a;
    };
    auto local_lift_mfrac = [&](const Mfrac& src,
                                 const std::vector<Mpoly>& gens) -> Mfrac {
        std::vector<fmpz_mpoly_struct*> ptrs(gens.size());
        for (std::size_t i = 0; i < gens.size(); ++i) {
            ptrs[i] = const_cast<fmpz_mpoly_struct*>(gens[i].raw());
        }
        Mpoly num_in(uf_ctx);
        Mpoly den_in(uf_ctx);
        int ok1 = fmpz_mpoly_compose_fmpz_mpoly(
            num_in.raw(), src.numerator().raw(), ptrs.data(),
            src.ctx()->raw(), uf_ctx->raw());
        int ok2 = fmpz_mpoly_compose_fmpz_mpoly(
            den_in.raw(), src.denominator().raw(), ptrs.data(),
            src.ctx()->raw(), uf_ctx->raw());
        if (!ok1 || !ok2) {
            throw std::runtime_error("apply_uf_replacement: lift failed");
        }
        return Mfrac(std::move(num_in), std::move(den_in));
    };
    auto local_make_lift_gens = [&](long n_src) -> std::vector<Mpoly> {
        std::vector<Mpoly> gs;
        gs.reserve((std::size_t)n_src);
        for (long i = 0; i < n_src; ++i) {
            gs.push_back(Mpoly::variable(uf_ctx, i));
        }
        return gs;
    };
    long n_loop  = static_cast<long>(fc.n_loops());
    long n_red   = static_cast<long>(fc.n_red_legs());
    if (fc.reduced_replacement.empty()) return p.clone();

    std::vector<Mpoly> lift_gens = local_make_lift_gens(fc.ctx->n_vars());
    std::map<std::string, Mfrac> lifted_rules;
    for (const auto& [k, v] : fc.reduced_replacement) {
        lifted_rules.emplace(k, local_lift_mfrac(v, lift_gens));
    }

    auto walk_mpoly = [&](const Mpoly& src) -> Mfrac {
        Mfrac acc = Mfrac::zero(uf_ctx);
        long len = fmpz_mpoly_length(src.raw(), uf_ctx->raw());
        std::vector<unsigned long> exp((std::size_t)uf_ctx->n_vars());
        fmpz_t coeff;
        fmpz_init(coeff);
        for (long t = 0; t < len; ++t) {
            fmpz_mpoly_get_term_coeff_fmpz(coeff, src.raw(), t, uf_ctx->raw());
            fmpz_mpoly_get_term_exp_ui(exp.data(), src.raw(), t, uf_ctx->raw());

            std::vector<long> leg_powers;
            for (long i = 0; i < n_red; ++i) {
                long e = static_cast<long>(exp[(std::size_t)(n_loop + i)]);
                for (long k = 0; k < e; ++k) leg_powers.push_back(i);
            }

            Mfrac leg_replacement = Mfrac::one(uf_ctx);
            std::size_t paired = 0;
            for (; paired + 1 < leg_powers.size(); paired += 2) {
                std::string key = local_leg_pair_key(
                    fc.reduced_legs[(std::size_t)leg_powers[paired]],
                    fc.reduced_legs[(std::size_t)leg_powers[paired + 1]]);
                auto it = lifted_rules.find(key);
                if (it == lifted_rules.end()) {
                    Mpoly mono = Mpoly::variable(
                        uf_ctx, n_loop + leg_powers[paired]);
                    mono *= Mpoly::variable(
                        uf_ctx, n_loop + leg_powers[paired + 1]);
                    leg_replacement *= Mfrac::from_mpoly(std::move(mono));
                } else {
                    leg_replacement *= it->second;
                }
            }
            for (; paired < leg_powers.size(); ++paired) {
                Mpoly mono = Mpoly::variable(
                    uf_ctx, n_loop + leg_powers[paired]);
                leg_replacement *= Mfrac::from_mpoly(std::move(mono));
            }

            std::vector<unsigned long> rest_exp = exp;
            for (long i = 0; i < n_red; ++i) {
                rest_exp[(std::size_t)(n_loop + i)] = 0;
            }
            Mpoly rest_mono = Mpoly::zero(uf_ctx);
            fmpz_mpoly_set_coeff_fmpz_ui(rest_mono.raw(), coeff,
                                          rest_exp.data(), uf_ctx->raw());

            Mfrac contrib = Mfrac::from_mpoly(std::move(rest_mono));
            contrib *= leg_replacement;
            acc += contrib;
        }
        fmpz_clear(coeff);
        return acc;
    };

    Mfrac new_num = walk_mpoly(p.numerator());
    Mfrac new_den = walk_mpoly(p.denominator());
    if (new_den.is_zero()) {
        throw std::runtime_error(
            "apply_uf_replacement: denominator vanished after substitution");
    }
    return Mfrac(new_num.numerator() * new_den.denominator(),
                  new_den.numerator() * new_num.denominator());
}

namespace {

// "Mass term of D":  D == momentum^2 + mass; mass = coe3 - coe2^2 / (4*coe1).
Mfrac mass_term_uf(const FamilyConfig& fc,
                    const std::shared_ptr<MpolyContext>& uf_ctx,
                    const std::vector<Mpoly>& lift_gens,
                    const Mpoly& den) {
    long n_loop = (long)fc.n_loops();

    long chosen_loop = -1;
    for (long j = 0; j < n_loop; ++j) {
        Mpoly c2 = den.coeff_of(j, 2);
        if (!c2.is_zero()) {
            chosen_loop = j;
            break;
        }
    }
    if (chosen_loop < 0) {
        Mpoly lifted = compose_to(den, uf_ctx, lift_gens);
        return Mfrac::from_mpoly(std::move(lifted));
    }

    Mpoly coe1 = den.coeff_of(chosen_loop, 2);
    Mpoly coe2 = den.coeff_of(chosen_loop, 1);
    Mpoly coe3 = den.coeff_of(chosen_loop, 0);

    Mpoly coe1_uf = compose_to(coe1, uf_ctx, lift_gens);
    Mpoly coe2_uf = compose_to(coe2, uf_ctx, lift_gens);
    Mpoly coe3_uf = compose_to(coe3, uf_ctx, lift_gens);

    Mpoly four_coe1(uf_ctx);
    {
        fmpz_t four;
        fmpz_init_set_ui(four, 4);
        fmpz_mpoly_scalar_mul_fmpz(four_coe1.raw(), coe1_uf.raw(),
                                    four, uf_ctx->raw());
        fmpz_clear(four);
    }
    Mpoly num = four_coe1 * coe3_uf - coe2_uf * coe2_uf;
    Mfrac mass_raw(std::move(num), std::move(four_coe1));

    return apply_uf_replacement(fc, uf_ctx, mass_raw);
}

}  // namespace

ABCResult evaluate_abc(const FamilyConfig& fc,
                       const std::vector<Mpoly>& denominators) {
    const std::size_t n = denominators.size();
    if (n == 0) {
        throw std::invalid_argument("evaluate_abc: empty denominator list");
    }
    for (const auto& d : denominators) {
        if (d.ctx().get() != fc.ctx.get()) {
            throw std::invalid_argument("evaluate_abc: denominator ctx mismatch");
        }
    }

    const long l = (long)fc.n_loops();
    const long e = (long)fc.n_red_legs();

    auto uf_ctx = build_uf_ctx(fc, n);
    const long first_x = fc.ctx->n_vars();

    auto lift_gens = make_lift_gen_map(uf_ctx, fc.ctx->n_vars());

    MpolyMatrix Ax = MpolyMatrix::zeros(uf_ctx,
                                          (std::size_t)l, (std::size_t)l);

    std::vector<Mpoly> Bx;
    Bx.reserve((std::size_t)l);
    for (long j = 0; j < l; ++j) Bx.emplace_back(uf_ctx);

    Mpoly Cx(uf_ctx);

    for (std::size_t i = 0; i < n; ++i) {
        const Mpoly& de = denominators[i];

        Mpoly x_i = Mpoly::variable(uf_ctx, first_x + (long)i);

        for (long j = 0; j < l; ++j) {
            for (long k = 0; k < l; ++k) {
                Mpoly coef = bilinear_coef(de, j, k);
                if (coef.is_zero()) continue;

                Mpoly lifted = compose_to(coef, uf_ctx, lift_gens);
                Mpoly scaled(uf_ctx);
                if (j == k) {
                    fmpz_t two;
                    fmpz_init_set_ui(two, 2);
                    fmpz_mpoly_scalar_mul_fmpz(scaled.raw(), lifted.raw(),
                                                two, uf_ctx->raw());
                    fmpz_clear(two);
                } else {
                    scaled = std::move(lifted);
                }
                Mpoly term = x_i * scaled;
                Ax(j, k) += term;
            }
        }

        for (long j = 0; j < l; ++j) {
            Mpoly bpre_j(uf_ctx);
            for (long k = 0; k < e; ++k) {
                long leg_var = l + k;
                Mpoly coef = bilinear_coef(de, j, leg_var);
                if (coef.is_zero()) continue;

                Mpoly lifted_coef = compose_to(coef, uf_ctx, lift_gens);
                Mpoly leg_in_uf   = Mpoly::variable(uf_ctx, leg_var);
                bpre_j += lifted_coef * leg_in_uf;
            }
            if (!bpre_j.is_zero()) {
                Bx[(std::size_t)j] += x_i * bpre_j;
            }
        }

        Mpoly Ci = de.clone();
        fmpz_t zero;
        fmpz_init(zero);
        for (long j = 0; j < l; ++j) {
            Ci = Ci.substitute(j, zero);
        }
        fmpz_clear(zero);
        Mpoly Ci_uf = compose_to(Ci, uf_ctx, lift_gens);
        Cx += x_i * Ci_uf;
    }

    ABCResult out;
    out.uf_ctx     = uf_ctx;
    out.first_x_var = first_x;
    out.Ax         = std::move(Ax);
    out.Bx         = std::move(Bx);
    out.Cx         = std::move(Cx);
    return out;
}

UFResult evaluate_uf(const FamilyConfig& fc,
                     const std::vector<Mpoly>& denominators) {
    ABCResult abc = evaluate_abc(fc, denominators);
    const long l = (long)fc.n_loops();
    const std::size_t n = denominators.size();

    UFResult out;
    out.uf_ctx     = abc.uf_ctx;
    out.first_x_var = abc.first_x_var;

    Mpoly det_A_int = abc.Ax.det();

    {
        Mpoly two_pow_l(abc.uf_ctx);
        {
            fmpz_t v;
            fmpz_init_set_ui(v, 1);
            fmpz_mul_2exp(v, v, (ulong)l);
            fmpz_mpoly_set_fmpz(two_pow_l.raw(), v, abc.uf_ctx->raw());
            fmpz_clear(v);
        }
        out.u = Mfrac(det_A_int.clone(), std::move(two_pow_l));
    }

    {
        auto lift_gens = make_lift_gen_map(abc.uf_ctx, fc.ctx->n_vars());
        Mfrac f0(Mfrac::zero(abc.uf_ctx));
        for (std::size_t i = 0; i < n; ++i) {
            Mfrac mt = mass_term_uf(fc, abc.uf_ctx, lift_gens, denominators[i]);
            if (mt.is_zero()) continue;
            Mpoly x_i = Mpoly::variable(abc.uf_ctx,
                                          abc.first_x_var + (long)i);
            mt *= Mfrac::from_mpoly(std::move(x_i));
            f0 += mt;
        }
        out.f0 = -f0;
    }

    if (det_A_int.is_zero()) {
        out.degenerate = true;
        out.f = Mfrac::zero(abc.uf_ctx);
        out.u = Mfrac::zero(abc.uf_ctx);
        return out;
    }
    out.degenerate = false;

    MpolyMatrix adj_A = abc.Ax.adjugate();

    Mpoly BAB(abc.uf_ctx);
    for (long j = 0; j < l; ++j) {
        for (long k = 0; k < l; ++k) {
            BAB += abc.Bx[(std::size_t)j] * adj_A(j, k) *
                   abc.Bx[(std::size_t)k];
        }
    }

    Mpoly detC = det_A_int * abc.Cx;

    Mpoly numerator(abc.uf_ctx);
    {
        fmpz_t two;
        fmpz_init_set_ui(two, 2);
        Mpoly two_detC(abc.uf_ctx);
        fmpz_mpoly_scalar_mul_fmpz(two_detC.raw(), detC.raw(),
                                    two, abc.uf_ctx->raw());
        fmpz_clear(two);
        numerator = BAB - two_detC;
    }
    Mpoly two_pow_lp1(abc.uf_ctx);
    {
        fmpz_t v;
        fmpz_init_set_ui(v, 1);
        fmpz_mul_2exp(v, v, (ulong)(l + 1));
        fmpz_mpoly_set_fmpz(two_pow_lp1.raw(), v, abc.uf_ctx->raw());
        fmpz_clear(v);
    }
    Mfrac f_raw(std::move(numerator), std::move(two_pow_lp1));

    Mfrac f_after = apply_uf_replacement(fc, abc.uf_ctx, f_raw);

    Mfrac uf0 = out.u.clone();
    uf0 *= out.f0;
    out.f = std::move(f_after);
    out.f -= uf0;

    return out;
}

}  // namespace amflow::qft
