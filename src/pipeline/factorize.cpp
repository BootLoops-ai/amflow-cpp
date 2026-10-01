// SPDX-License-Identifier: MIT
// amflow::factorize — implementation.
//

#include "amflow/pipeline/factorize.hpp"

#include <algorithm>
#include <map>
#include <set>
#include <stdexcept>
#include <utility>

#include <flint/flint.h>
#include <flint/fmpq.h>
#include <flint/fmpz.h>
#include <flint/fmpz_mpoly.h>

#include "amflow/qft/dense_q.hpp"
#include "amflow/qft/family_uf.hpp"
#include "amflow/qft/topology.hpp"

namespace amflow::pipeline {

using algebra::Mpoly;

namespace {

Mpoly apply_loop_redef(const qft::FamilyConfig& fc,
                        const Mpoly& expr,
                        const std::vector<std::pair<long, Mpoly>>& redef) {
    auto ctx = fc.ctx;
    std::vector<Mpoly> subs;
    subs.reserve((std::size_t)ctx->n_vars());
    for (long v = 0; v < ctx->n_vars(); ++v) {
        subs.push_back(Mpoly::variable(ctx, v));
    }
    for (const auto& [v, mp] : redef) {
        subs[(std::size_t)v] = mp.clone();
    }
    std::vector<fmpz_mpoly_struct*> ptrs(subs.size());
    for (std::size_t i = 0; i < subs.size(); ++i) {
        ptrs[i] = const_cast<fmpz_mpoly_struct*>(subs[i].raw());
    }
    Mpoly out(ctx);
    int rc = fmpz_mpoly_compose_fmpz_mpoly(out.raw(), expr.raw(),
                                              ptrs.data(),
                                              ctx->raw(), ctx->raw());
    if (rc == 0) {
        throw std::runtime_error(
            "factorize_family: compose for loop redef failed");
    }
    return out;
}

}  // namespace

std::vector<FactorizedComponent>
factorize_family(const qft::FamilyConfig& fc,
                 const std::vector<Mpoly>& denominators,
                 const std::vector<std::vector<long>>& patts) {
    auto info = qft::analyze_topology(fc, denominators);
    if (info.empty()) return {};

    long n_x = (long)denominators.size();
    long L   = (long)fc.n_loops();

    qft::UFResult uf = qft::evaluate_uf(fc, denominators);
    long first_x = uf.first_x_var;

    std::map<long, long> var_to_pos;
    for (long k = 0; k < n_x; ++k) {
        var_to_pos[first_x + k] = k;
    }

    std::vector<FactorizedComponent> out;
    out.reserve(info.size());

    for (auto& comp : info) {
        FactorizedComponent fcomp;

        auto u0_num = comp.u0.numerator();
        auto ctx_uf = comp.u0.ctx();
        long len = fmpz_mpoly_length(u0_num.raw(), ctx_uf->raw());
        if (len == 0) {
            throw std::runtime_error(
                "factorize_family: empty u0 component");
        }
        std::vector<std::set<long>> mono_vars;
        mono_vars.reserve((std::size_t)len);
        for (long t = 0; t < len; ++t) {
            std::vector<unsigned long> exp_t((std::size_t)ctx_uf->n_vars());
            fmpz_mpoly_get_term_exp_ui(exp_t.data(), u0_num.raw(),
                                        t, ctx_uf->raw());
            std::set<long> s;
            for (long v = first_x; v < first_x + n_x; ++v) {
                if (exp_t[(std::size_t)v] > 0) s.insert(v);
            }
            mono_vars.push_back(std::move(s));
        }

        long count_one  = 0;
        long count_zero = 0;
        long pos_one    = -1;
        for (std::size_t k = 0; k < comp.mass.size(); ++k) {
            if (comp.mass[k].is_one()) {
                ++count_one;
                if (pos_one < 0) pos_one = (long)k;
            } else if (comp.mass[k].is_zero()) {
                ++count_zero;
            }
        }
        if (comp.vacQ
                && count_one == 1
                && count_zero == (long)comp.mass.size() - 1
                && pos_one >= 0) {
            long target_var = comp.var[(std::size_t)pos_one];

            std::stable_sort(
                mono_vars.begin(), mono_vars.end(),
                [target_var](const std::set<long>& a,
                              const std::set<long>& b) {
                    bool a_has = a.count(target_var) > 0;
                    bool b_has = b.count(target_var) > 0;
                    return a_has && !b_has;
                });
        }

        std::vector<long> tobeloop;
        for (long v : mono_vars.front()) {
            auto it = var_to_pos.find(v);
            if (it != var_to_pos.end()) tobeloop.push_back(it->second);
        }
        std::sort(tobeloop.begin(), tobeloop.end());

        long n_tbl = (long)tobeloop.size();

        std::vector<fmpq*> mat_flat((std::size_t)(n_tbl * L), nullptr);
        for (long i = 0; i < n_tbl * L; ++i) {
            mat_flat[(std::size_t)i] = (fmpq*)flint_malloc(sizeof(fmpq));
            fmpq_init(mat_flat[(std::size_t)i]);
        }
        for (long ti = 0; ti < n_tbl; ++ti) {
            long k = tobeloop[(std::size_t)ti];
            const Mpoly& d = denominators[(std::size_t)k];
            long chosen = -1;
            for (long j = 0; j < L; ++j) {
                if (!d.coeff_of(j, 2).is_zero()) { chosen = j; break; }
            }
            if (chosen < 0) continue;
            // Mirror upstream `ToSquare` (`AMFlow.m:437-444`): the per-row
            // entry MMA produces is `Coefficient[mom, Loop[j]]` where
            // `mom = Sqrt[alpha]*l_chosen + beta_j/(2*Sqrt[alpha])` and
            // `alpha = Coefficient[d, l_chosen, 2]`,
            // `beta_j = Coefficient[d, l_chosen, 1]` (a polynomial in
            // other loops/legs).  Under the AMFlow propagator convention
            // alpha == 1 (no metric factor); we assert that and then take
            // `Coefficient[mom, l_j] = (coef of l_chosen*l_j in d) / 2`.
            std::vector<unsigned long> z(
                (std::size_t)fc.ctx->n_vars(), 0);
            {
                Mpoly chosen_sq = d.coeff_of(chosen, 2);
                fmpz_t alpha;
                fmpz_init(alpha);
                fmpz_mpoly_get_coeff_fmpz_ui(alpha, chosen_sq.raw(),
                                              z.data(),
                                              fc.ctx->raw());
                if (!fmpz_is_one(alpha)) {
                    fmpz_clear(alpha);
                    throw std::runtime_error(
                        "factorize_family: chosen loop's squared "
                        "coefficient is not 1; AMFlow propagator "
                        "convention violated (alpha != 1 in ToSquare).");
                }
                fmpz_clear(alpha);
            }
            fmpq_set_si(mat_flat[(std::size_t)(ti * L + chosen)], 1, 1);
            for (long j = 0; j < L; ++j) {
                if (j == chosen) continue;
                // Coefficient of the bilinear monomial l_chosen * l_j in d.
                // (Previous logic took `coeff_of(j, 1)` and then its
                // constant term — which silently dropped the contribution
                // from cross-loop bilinears like `-2*l_chosen*l_j` in
                // `(l_chosen - l_j)^2` props.)
                Mpoly cross = d.coeff_of(
                    std::vector<long>{chosen, j},
                    std::vector<unsigned long>{1, 1});
                if (cross.is_zero()) continue;
                fmpz_t c;
                fmpz_init(c);
                fmpz_mpoly_get_coeff_fmpz_ui(c, cross.raw(), z.data(),
                                              fc.ctx->raw());
                fmpq_t qval, two;
                fmpq_init(qval); fmpq_init(two);
                fmpq_set_fmpz(qval, c);
                fmpq_set_si(two, 2, 1);
                fmpq_div(qval, qval, two);
                fmpq_set(mat_flat[(std::size_t)(ti * L + j)], qval);
                fmpq_clear(qval); fmpq_clear(two);
                fmpz_clear(c);
            }
        }

        qft::DenseQMatrix M((std::size_t)n_tbl, (std::size_t)(L + n_tbl));
        for (long ti = 0; ti < n_tbl; ++ti) {
            for (long j = 0; j < L; ++j) {
                fmpq_set(M.at((std::size_t)ti, (std::size_t)j),
                          mat_flat[(std::size_t)(ti * L + j)]);
            }
            for (long j = 0; j < n_tbl; ++j) {
                fmpq_set_si(M.at((std::size_t)ti,
                                  (std::size_t)(L + j)),
                              (j == ti) ? -1 : 0, 1);
            }
        }
        for (auto* p : mat_flat) { fmpq_clear(p); flint_free(p); }

        M.rref();

        std::vector<std::pair<long, long>> pivots;
        for (long ti = 0; ti < n_tbl; ++ti) {
            long pivot_col = -1;
            for (long j = 0; j < L; ++j) {
                if (fmpq_is_one(M.at((std::size_t)ti, (std::size_t)j))) {
                    pivot_col = j;
                    break;
                }
            }
            if (pivot_col >= 0) pivots.emplace_back(ti, pivot_col);
        }

        if ((long)pivots.size() == n_tbl) {
            std::vector<long> row_pivot_col(n_tbl, -1);
            for (auto& [ti, jc] : pivots) row_pivot_col[ti] = jc;

            std::vector<std::pair<long, Mpoly>> redef;
            redef.reserve(pivots.size());
            for (auto& [ti, p_col] : pivots) {
                Mpoly rhs = Mpoly::zero(fc.ctx);
                for (long j = 0; j < n_tbl; ++j) {
                    const fmpq* mij = M.at((std::size_t)ti,
                                              (std::size_t)(L + j));
                    if (fmpq_is_zero(mij)) continue;
                    long lp = row_pivot_col[j];
                    if (lp < 0) continue;
                    if (!fmpz_is_one(fmpq_denref(mij))) {
                        throw std::runtime_error(
                            "factorize_family: non-integer "
                            "redef coefficient (augmented)");
                    }
                    fmpz_t neg_num;
                    fmpz_init(neg_num);
                    fmpz_neg(neg_num, fmpq_numref(mij));
                    Mpoly term = Mpoly::variable(fc.ctx, lp);
                    fmpz_mpoly_scalar_mul_fmpz(term.raw(), term.raw(),
                                                neg_num, fc.ctx->raw());
                    rhs += term;
                    fmpz_clear(neg_num);
                }
                for (long jc = 0; jc < L; ++jc) {
                    if (jc == p_col) continue;
                    const fmpq* mij = M.at((std::size_t)ti,
                                              (std::size_t)jc);
                    if (fmpq_is_zero(mij)) continue;
                    if (!fmpz_is_one(fmpq_denref(mij))) {
                        throw std::runtime_error(
                            "factorize_family: non-integer "
                            "redef coefficient (left-block)");
                    }
                    fmpz_t neg_num;
                    fmpz_init(neg_num);
                    fmpz_neg(neg_num, fmpq_numref(mij));
                    Mpoly term = Mpoly::variable(fc.ctx, jc);
                    fmpz_mpoly_scalar_mul_fmpz(term.raw(), term.raw(),
                                                neg_num, fc.ctx->raw());
                    rhs += term;
                    fmpz_clear(neg_num);
                }
                redef.emplace_back(p_col, std::move(rhs));
            }

            for (auto& [ti, jc] : pivots) {
                fcomp.loops.push_back(fc.loops[(std::size_t)jc]);
            }

            std::vector<long> var_pos;
            var_pos.reserve(comp.var.size());
            for (long v : comp.var) {
                auto it = var_to_pos.find(v);
                if (it != var_to_pos.end()) var_pos.push_back(it->second);
            }
            std::sort(var_pos.begin(), var_pos.end());
            fcomp.original_positions = var_pos;
            fcomp.propagators.reserve(var_pos.size());
            for (long pos : var_pos) {
                Mpoly d_orig = denominators[(std::size_t)pos].clone();
                Mpoly d_new = apply_loop_redef(fc, d_orig, redef);
                fcomp.propagators.push_back(std::move(d_new));
            }
        } else {
            // Defensive no-redef fallback — fires when the
            // loop-redefinition matrix M (constructed at lines
            // 137-186) has fewer pivots than `n_tbl` after RREF,
            // i.e. when the heuristic to-be-loop selection from
            // U-poly's first monomial doesn't span a complete
            // loop basis.  In that case we cannot apply the
            // standard loop redefinition; fall through to the
            // trivial "use original loops + original propagators"
            // structure.  The fallback is mathematically safe — it
            // leaves the family unchanged and downstream
            // FactorizeFamily consumers handle the unfactorized
            // output the same way they would a single-component
            // family.
            for (long jc = 0; jc < L; ++jc) {
                fcomp.loops.push_back(fc.loops[(std::size_t)jc]);
            }
            std::vector<long> var_pos;
            var_pos.reserve(comp.var.size());
            for (long v : comp.var) {
                auto it = var_to_pos.find(v);
                if (it != var_to_pos.end()) var_pos.push_back(it->second);
            }
            std::sort(var_pos.begin(), var_pos.end());
            fcomp.original_positions = var_pos;
            fcomp.propagators.reserve(var_pos.size());
            for (long pos : var_pos) {
                fcomp.propagators.push_back(
                    denominators[(std::size_t)pos].clone());
            }
        }

        fcomp.patterns.reserve(patts.size());
        for (const auto& patt : patts) {
            std::vector<long> p_new;
            for (long pos : fcomp.original_positions) {
                if ((long)patt.size() <= pos) {
                    throw std::runtime_error(
                        "factorize_family: pattern too short");
                }
                p_new.push_back(patt[(std::size_t)pos]);
            }
            fcomp.patterns.push_back(std::move(p_new));
        }

        out.push_back(std::move(fcomp));
    }

    return out;
}

}  // namespace amflow::pipeline
