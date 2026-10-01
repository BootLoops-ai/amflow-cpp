// SPDX-License-Identifier: MIT
// qft::region — region utilities (basic).
//

#include "amflow/qft/region.hpp"

#include <set>
#include <stdexcept>

#include <flint/fmpz.h>
#include <flint/fmpz_mpoly.h>
#include <flint/fmpz_mpoly_q.h>

namespace amflow::qft {

using algebra::Mfrac;
using algebra::Mpoly;
using algebra::MpolyContext;
using algebra::MpolyMatrix;

namespace {
Mfrac do_apply_region_rule(const RegionContext& rctx,
                           const Mfrac& expr,
                           const std::vector<Mfrac>& rule,
                           long n_loop);
}  // namespace

RegionContext make_region_context(const FamilyConfig& fc) {
    std::vector<std::string> names;
    names.reserve(static_cast<std::size_t>(fc.ctx->n_vars()) + 2);
    for (long i = 0; i < fc.ctx->n_vars(); ++i) {
        names.push_back(fc.ctx->var_name(i));
    }
    names.emplace_back("__amf_eta");
    names.emplace_back("__amf_half_eta");
    auto rc_ctx = std::make_shared<MpolyContext>(std::move(names));

    RegionContext rctx;
    rctx.ctx          = rc_ctx;
    rctx.eta_var      = fc.ctx->n_vars();
    rctx.half_eta_var = fc.ctx->n_vars() + 1;

    rctx.lift_gens.reserve(static_cast<std::size_t>(fc.ctx->n_vars()));
    for (long i = 0; i < fc.ctx->n_vars(); ++i) {
        rctx.lift_gens.push_back(Mpoly::variable(rc_ctx, i));
    }
    return rctx;
}

Mpoly lift_to_region(const Mpoly& src, const RegionContext& rctx) {
    auto src_ctx = src.ctx();
    if (src_ctx->n_vars() != static_cast<long>(rctx.lift_gens.size())) {
        throw std::invalid_argument(
            "lift_to_region: src ctx variable count mismatch");
    }
    std::vector<fmpz_mpoly_struct*> ptrs(rctx.lift_gens.size());
    for (std::size_t i = 0; i < rctx.lift_gens.size(); ++i) {
        ptrs[i] = const_cast<fmpz_mpoly_struct*>(rctx.lift_gens[i].raw());
    }
    Mpoly out(rctx.ctx);
    int ok = fmpz_mpoly_compose_fmpz_mpoly(
        out.raw(), src.raw(), ptrs.data(),
        src_ctx->raw(), rctx.ctx->raw());
    if (ok == 0) {
        throw std::runtime_error("lift_to_region: FLINT compose failed");
    }
    return out;
}

std::vector<Mpoly>
branch_momenta(const FamilyConfig& fc,
               const std::vector<Mpoly>& denominators) {
    std::vector<Mpoly> out;
    out.reserve(denominators.size());

    long n_loop = static_cast<long>(fc.n_loops());
    long n_red  = static_cast<long>(fc.n_red_legs());

    for (const auto& d : denominators) {
        if (d.ctx().get() != fc.ctx.get()) {
            throw std::invalid_argument(
                "branch_momenta: denominator ctx mismatch");
        }

        long chosen_loop = -1;
        for (long j = 0; j < n_loop; ++j) {
            Mpoly c2 = d.coeff_of(j, 2);
            if (!c2.is_zero()) {
                chosen_loop = j;
                break;
            }
        }
        if (chosen_loop < 0) {
            out.emplace_back(fc.ctx);
            continue;
        }

        Mpoly coe1 = d.coeff_of(chosen_loop, 2);
        Mpoly coe2 = d.coeff_of(chosen_loop, 1);

        if (!fmpz_mpoly_is_fmpz(coe1.raw(), fc.ctx->raw())) {
            throw std::runtime_error(
                "branch_momenta: leading loop-squared coefficient must be 1");
        }
        {
            fmpz_t coe1_z;
            fmpz_init(coe1_z);
            fmpz_mpoly_get_fmpz(coe1_z, coe1.raw(), fc.ctx->raw());
            const int eq_one = fmpz_is_one(coe1_z);
            fmpz_clear(coe1_z);
            if (!eq_one) {
                throw std::runtime_error(
                    "branch_momenta: leading loop-squared coefficient must be 1");
            }
        }

        Mpoly mom2 = Mpoly::variable(fc.ctx, chosen_loop);
        {
            Mpoly coe2_half(fc.ctx);
            long len = fmpz_mpoly_length(coe2.raw(), fc.ctx->raw());
            std::vector<unsigned long> exp((std::size_t)fc.ctx->n_vars());
            fmpz_t coeff, q, r, two;
            fmpz_init(coeff);
            fmpz_init(q); fmpz_init(r);
            fmpz_init_set_ui(two, 2);
            bool indivisible = false;
            for (long t = 0; t < len; ++t) {
                fmpz_mpoly_get_term_exp_ui(exp.data(), coe2.raw(), t,
                                            fc.ctx->raw());
                fmpz_mpoly_get_term_coeff_fmpz(coeff, coe2.raw(), t,
                                                fc.ctx->raw());
                fmpz_tdiv_qr(q, r, coeff, two);
                if (!fmpz_is_zero(r)) { indivisible = true; break; }
                fmpz_mpoly_set_coeff_fmpz_ui(coe2_half.raw(), q,
                                              exp.data(), fc.ctx->raw());
            }
            fmpz_clear(coeff);
            fmpz_clear(q); fmpz_clear(r);
            fmpz_clear(two);
            if (indivisible) {
                throw std::runtime_error(
                    "branch_momenta: coe2 coefficient is not even -- "
                    "non-physical propagator?");
            }
            mom2 += coe2_half;
        }

        {
            fmpz_t zero;
            fmpz_init(zero);
            for (long k = 0; k < n_red; ++k) {
                Mpoly tmp(fc.ctx);
                fmpz_mpoly_evaluate_one_fmpz(tmp.raw(), mom2.raw(),
                                              n_loop + k, zero,
                                              fc.ctx->raw());
                mom2 = std::move(tmp);
            }
            fmpz_clear(zero);
        }

        out.push_back(std::move(mom2));
    }
    return out;
}

LoopTransform
branch_to_loop(const FamilyConfig& fc,
               const RegionContext& rctx,
               const std::vector<Mpoly>& candidates) {
    LoopTransform out;
    out.ok = false;

    long L = static_cast<long>(fc.n_loops());
    if (static_cast<long>(candidates.size()) != L) return out;

    MpolyMatrix A(rctx.ctx, (std::size_t)L, (std::size_t)L);
    for (long i = 0; i < L; ++i) {
        for (long j = 0; j < L; ++j) {
            Mpoly co = candidates[(std::size_t)i].coeff_of(j, 1);
            A(i, j) = lift_to_region(co, rctx);
        }
    }

    Mpoly det_A = A.det();
    if (det_A.is_zero()) return out;

    // Defensive Jacobian-is-unit check (mirrors AMFlow.m:731-732, where
    // upstream multiplies BoundaryIntegrands by `Abs[Det[...]]^(4-2eps)`).
    // The C++ port omits the explicit Jacobian factor, which is correct
    // *only* when |det_A| == 1.  By construction (`branch_momenta`
    // enforces each branch propagator's leading loop coefficient = 1),
    // the matrix `A` is a permutation matrix and `det_A` is the constant
    // ±1.  If a future change relaxes that precondition, we want to
    // catch it here rather than silently produce wrong boundary
    // integrands.
    if (!fmpz_mpoly_is_fmpz(det_A.raw(), rctx.ctx->raw())) {
        throw std::runtime_error(
            "branch_to_loop: loop-redefinition determinant is not a "
            "constant; the missing |Det|^(4-2eps) Jacobian factor "
            "(AMFlow.m:731-732) would produce incorrect boundary "
            "integrands.  This indicates a violation of the unit-leading-"
            "loop-coefficient precondition in branch_momenta.");
    }
    {
        fmpz_t det_z;
        fmpz_init(det_z);
        fmpz_mpoly_get_fmpz(det_z, det_A.raw(), rctx.ctx->raw());
        fmpz_t abs_det;
        fmpz_init(abs_det);
        fmpz_abs(abs_det, det_z);
        const int abs_one = fmpz_is_one(abs_det);
        fmpz_clear(det_z);
        fmpz_clear(abs_det);
        if (!abs_one) {
            throw std::runtime_error(
                "branch_to_loop: |det_A| != 1; the missing "
                "|Det|^(4-2eps) Jacobian factor (AMFlow.m:731-732) "
                "would produce incorrect boundary integrands.");
        }
    }

    MpolyMatrix adj = A.adjugate();

    out.map.reserve(static_cast<std::size_t>(L));
    for (long i = 0; i < L; ++i) {
        Mpoly num(rctx.ctx);
        for (long j = 0; j < L; ++j) {
            Mpoly term = adj(i, j) * Mpoly::variable(rctx.ctx, j);
            num += term;
        }
        out.map.push_back(Mfrac(std::move(num), det_A.clone()));
    }
    out.ok = true;
    return out;
}

std::vector<Mfrac>
region_rule(const RegionContext& rctx,
            const LoopTransform& tr,
            const std::vector<int>& scale) {
    if (!tr.ok) {
        throw std::invalid_argument("region_rule: invalid LoopTransform");
    }
    if (scale.size() != tr.map.size()) {
        throw std::invalid_argument("region_rule: scale size mismatch");
    }
    std::vector<Mfrac> out;
    out.reserve(tr.map.size());

    Mpoly half_eta = Mpoly::variable(rctx.ctx, rctx.half_eta_var);
    std::vector<Mfrac> scaled_loops;
    scaled_loops.reserve(tr.map.size());
    for (std::size_t k = 0; k < tr.map.size(); ++k) {
        // Eikonal region engine: scale
        // entries are half_eta degrees, sigma_k = scale_k/2 in eta units:
        //   0 -> l ~ 1,  1 -> l ~ he = sqrt(eta),  2 -> l ~ he^2 = eta.
        // Bit 2 is only ever generated by find_all_region when
        // AMFLOW_EIKONAL_REGIONS is enabled.
        if (scale[k] != 0 && scale[k] != 1 && scale[k] != 2) {
            throw std::invalid_argument(
                "region_rule: scale entries must be 0, 1 or 2");
        }
        Mfrac loop = Mfrac::from_mpoly(Mpoly::variable(rctx.ctx, (long)k));
        for (int s = 0; s < scale[k]; ++s) {
            loop *= Mfrac::from_mpoly(half_eta.clone());
        }
        scaled_loops.push_back(std::move(loop));
    }

    for (std::size_t k = 0; k < tr.map.size(); ++k) {
        out.push_back(
            do_apply_region_rule(rctx, tr.map[k], scaled_loops, (long)tr.map.size()));
    }
    return out;
}

namespace {

Mfrac do_apply_region_rule(const RegionContext& rctx,
                            const Mfrac& expr,
                            const std::vector<Mfrac>& rule,
                            long n_loop) {
    if ((long)rule.size() != n_loop) {
        throw std::invalid_argument(
            "apply_region_rule: rule size != n_loops");
    }

    long n_total = rctx.ctx->n_vars();
    std::vector<Mfrac> full_map;
    full_map.reserve(static_cast<std::size_t>(n_total));
    for (long i = 0; i < n_loop; ++i) {
        full_map.push_back(rule[(std::size_t)i].clone());
    }
    for (long i = n_loop; i < n_total; ++i) {
        full_map.push_back(Mfrac::from_mpoly(
            Mpoly::variable(rctx.ctx, i)));
    }

    auto sub_in_mpoly = [&](const Mpoly& src) -> Mfrac {
        Mfrac acc = Mfrac::zero(rctx.ctx);
        long len = fmpz_mpoly_length(src.raw(), rctx.ctx->raw());
        std::vector<unsigned long> exp((std::size_t)n_total);
        fmpz_t coeff;
        fmpz_init(coeff);
        for (long t = 0; t < len; ++t) {
            fmpz_mpoly_get_term_coeff_fmpz(coeff, src.raw(), t,
                                            rctx.ctx->raw());
            fmpz_mpoly_get_term_exp_ui(exp.data(), src.raw(), t,
                                        rctx.ctx->raw());
            Mfrac term = Mfrac::from_mpoly(
                Mpoly::constant(rctx.ctx, 1));
            for (long v = 0; v < n_total; ++v) {
                if (exp[(std::size_t)v] == 0) continue;
                Mfrac base = full_map[(std::size_t)v].clone();
                Mfrac power = Mfrac::from_mpoly(
                    Mpoly::constant(rctx.ctx, 1));
                for (unsigned long p = 0; p < exp[(std::size_t)v]; ++p) {
                    power *= base;
                }
                term *= power;
            }
            Mpoly c_mpoly = Mpoly::constant(rctx.ctx, coeff);
            term *= Mfrac::from_mpoly(std::move(c_mpoly));
            acc += term;
        }
        fmpz_clear(coeff);
        return acc;
    };

    Mfrac new_num = sub_in_mpoly(expr.numerator());
    Mfrac new_den = sub_in_mpoly(expr.denominator());
    if (new_den.is_zero()) {
        throw std::runtime_error(
            "apply_region_rule: denominator vanished after substitution");
    }
    return Mfrac(new_num.numerator() * new_den.denominator(),
                  new_den.numerator() * new_num.denominator());
}

}  // namespace

Mfrac apply_region_rule(const FamilyConfig& fc,
                         const RegionContext& rctx,
                         const Mfrac& expr,
                         const std::vector<Mfrac>& rule) {
    Mfrac lifted = (expr.ctx().get() == rctx.ctx.get())
        ? expr.clone()
        : Mfrac(lift_to_region(expr.numerator(), rctx),
                lift_to_region(expr.denominator(), rctx));
    return do_apply_region_rule(rctx, lifted, rule, (long)fc.n_loops());
}

Mfrac apply_region_rule(const FamilyConfig& fc,
                         const RegionContext& rctx,
                         const Mpoly& expr,
                         const std::vector<Mfrac>& rule) {
    Mpoly lifted = (expr.ctx().get() == rctx.ctx.get())
        ? expr.clone()
        : lift_to_region(expr, rctx);
    return do_apply_region_rule(rctx,
                                  Mfrac::from_mpoly(std::move(lifted)),
                                  rule, (long)fc.n_loops());
}

std::vector<int>
branch_scale(const FamilyConfig& fc,
             const std::vector<Mfrac>& transformed_branches,
             const std::vector<int>& scale) {
    long L = static_cast<long>(fc.n_loops());
    if ((long)scale.size() != L) {
        throw std::invalid_argument("branch_scale: scale size != n_loops");
    }
    std::vector<int> out;
    out.reserve(transformed_branches.size());

    // Eikonal region engine: the branch
    // class is the he-degree of the branch momentum under the region
    // scaling, i.e. max_k { scale[k] : loop k appears in the branch }.
    // For scale entries restricted to {0,1} this reproduces the previous
    // boolean "touches a large loop" semantics exactly; with the eikonal
    // engine enabled it distinguishes l ~ sqrt(eta) (1) from l ~ eta (2)
    // branches so that bscale grouping keeps the sigma families separate.
    for (const auto& br : transformed_branches) {
        Mpoly num = br.numerator();
        Mpoly den = br.denominator();
        auto ctx = num.ctx();
        long n_total = ctx->n_vars();
        std::vector<unsigned long> exp((std::size_t)n_total);
        int cls = 0;
        auto walk = [&](const Mpoly& p) {
            long len = fmpz_mpoly_length(p.raw(), ctx->raw());
            for (long t = 0; t < len; ++t) {
                fmpz_mpoly_get_term_exp_ui(exp.data(), p.raw(), t,
                                            ctx->raw());
                for (long k = 0; k < L; ++k) {
                    if (exp[(std::size_t)k] > 0
                            && scale[(std::size_t)k] > cls) {
                        cls = scale[(std::size_t)k];
                    }
                }
            }
        };
        walk(num);
        walk(den);
        out.push_back(cls);
    }
    return out;
}

}  // namespace amflow::qft
