// SPDX-License-Identifier: MIT
// qft::family_config — implementation.
//

#include "amflow/qft/family_config.hpp"

#include <algorithm>
#include <regex>
#include <set>
#include <sstream>
#include <stdexcept>
#include <utility>

#include <flint/flint.h>
#include <flint/fmpq.h>
#include <flint/fmpz.h>
#include <flint/fmpz_mpoly.h>

namespace amflow::qft {

using algebra::Mfrac;
using algebra::Mpoly;
using algebra::MpolyContext;

namespace {

std::set<std::string> scrape_identifiers(const std::string& s) {
    std::set<std::string> out;
    static const std::regex pat(R"([A-Za-z_][A-Za-z_0-9]*)");
    auto begin = std::sregex_iterator(s.begin(), s.end(), pat);
    auto end   = std::sregex_iterator{};
    for (auto it = begin; it != end; ++it) out.insert(it->str());
    return out;
}

Mpoly compose_mpoly(const Mpoly& src,
                    const std::shared_ptr<MpolyContext>& dst_ctx,
                    const std::vector<Mpoly>& gen_map) {
    auto src_ctx = src.ctx();
    if (static_cast<long>(gen_map.size()) != src_ctx->n_vars()) {
        throw std::invalid_argument("compose_mpoly: gen_map size mismatch");
    }
    for (const auto& g : gen_map) {
        if (g.ctx().get() != dst_ctx.get()) {
            throw std::invalid_argument("compose_mpoly: gen ctx mismatch");
        }
    }
    std::vector<fmpz_mpoly_struct*> ptrs(gen_map.size());
    for (std::size_t i = 0; i < gen_map.size(); ++i) {
        ptrs[i] = const_cast<fmpz_mpoly_struct*>(gen_map[i].raw());
    }
    Mpoly out(dst_ctx);
    int ok = fmpz_mpoly_compose_fmpz_mpoly(
        out.raw(), src.raw(), ptrs.data(),
        src_ctx->raw(), dst_ctx->raw());
    if (ok == 0) {
        throw std::runtime_error("compose_mpoly: FLINT compose failed");
    }
    return out;
}

std::string leg_monomial_key(const std::string& a, const std::string& b) {
    if (a <= b) return a + "*" + b;
    return b + "*" + a;
}

}  // namespace

FamilyConfig FamilyConfig::build(
    std::string family,
    std::vector<std::string> loops,
    std::vector<std::string> legs,
    std::vector<std::pair<std::string, std::string>> conservation,
    std::vector<std::pair<std::string, std::string>> replacement,
    std::vector<std::string> propagator_strings,
    std::vector<int> cut,
    std::vector<int> prescription) {

    FamilyConfig fc;
    fc.family   = std::move(family);
    fc.loops    = std::move(loops);
    fc.legs     = std::move(legs);

    auto check_distinct = [](const std::vector<std::string>& v,
                              const std::string& kind) {
        std::set<std::string> seen;
        for (const auto& s : v) {
            if (s.empty()) {
                throw std::invalid_argument(
                    kind + ": empty name not allowed");
            }
            if (!seen.insert(s).second) {
                throw std::invalid_argument(
                    kind + ": duplicate name '" + s + "'");
            }
        }
    };
    check_distinct(fc.loops, "FamilyConfig::loops");
    check_distinct(fc.legs,  "FamilyConfig::legs");

    std::set<std::string> leg_set(fc.legs.begin(), fc.legs.end());
    for (const auto& [lhs, _] : conservation) {
        if (!leg_set.count(lhs)) {
            throw std::invalid_argument(
                "conservation: lhs '" + lhs + "' is not a leg");
        }
    }

    std::set<std::string> known;
    for (const auto& s : fc.loops) known.insert(s);
    for (const auto& s : fc.legs)  known.insert(s);
    std::set<std::string> inv_set;
    auto scrape_into = [&](const std::string& s) {
        for (const auto& tok : scrape_identifiers(s)) {
            if (known.count(tok)) continue;
            inv_set.insert(tok);
        }
    };
    for (const auto& [lhs, rhs] : replacement)  { scrape_into(lhs); scrape_into(rhs); }
    for (const auto& s : propagator_strings)    { scrape_into(s); }
    for (const auto& [lhs, rhs] : conservation) { scrape_into(rhs); }
    fc.invariants.assign(inv_set.begin(), inv_set.end());

    std::set<std::string> conserved_keys;
    for (const auto& [lhs, _] : conservation) conserved_keys.insert(lhs);
    for (const auto& l : fc.legs) {
        if (!conserved_keys.count(l)) fc.reduced_legs.push_back(l);
    }

    fc.build_context();

    std::vector<std::string> full_vars = fc.loops;
    for (const auto& s : fc.legs)       full_vars.push_back(s);
    for (const auto& s : fc.invariants) full_vars.push_back(s);
    auto full_ctx = std::make_shared<MpolyContext>(std::move(full_vars));

    std::map<std::string, Mpoly> cons_full;
    for (const auto& [lhs, rhs_str] : conservation) {
        Mpoly rhs = Mpoly::from_string(full_ctx, rhs_str);
        cons_full.emplace(lhs, std::move(rhs));
    }

    auto make_sub_map = [&]() {
        std::vector<Mpoly> gen_map;
        gen_map.reserve(full_ctx->n_vars());
        for (long i = 0; i < full_ctx->n_vars(); ++i) {
            const std::string& name = full_ctx->var_name(i);
            auto it = cons_full.find(name);
            if (it != cons_full.end()) {
                gen_map.push_back(Mpoly::zero(fc.ctx));
            } else {
                long fi = fc.ctx->var_index(name);
                if (fi < 0) {
                    throw std::logic_error(
                        "FamilyConfig: variable mismatch: '" + name + "'");
                }
                gen_map.push_back(Mpoly::variable(fc.ctx, fi));
            }
        }
        return gen_map;
    };

    auto gen_map = make_sub_map();

    for (const auto& [lhs, rhs_full] : cons_full) {
        long full_idx = full_ctx->var_index(lhs);
        gen_map[full_idx] = compose_mpoly(rhs_full, fc.ctx, gen_map);
    }

    for (const auto& [lhs, rhs_full] : cons_full) {
        Mpoly rhs_in_fc = compose_mpoly(rhs_full, fc.ctx, gen_map);
        fc.conservation.emplace(lhs, std::move(rhs_in_fc));
    }

    fc.propagators_after_conservation.reserve(propagator_strings.size());
    fc.propagator_denominators.reserve(propagator_strings.size());
    for (const auto& s : propagator_strings) {
        fmpz_t scale;
        fmpz_init(scale);
        Mpoly p_full = Mpoly::from_string(full_ctx, s, scale);
        fc.propagators_after_conservation.push_back(
            compose_mpoly(p_full, fc.ctx, gen_map));
        char* scale_raw = fmpz_get_str(nullptr, 10, scale);
        fc.propagator_denominators.emplace_back(scale_raw);
        flint_free(scale_raw);
        fmpz_clear(scale);
    }

    std::vector<std::pair<Mpoly, Mfrac>> rep_in_fc;
    rep_in_fc.reserve(replacement.size());
    for (const auto& [lhs_str, rhs_str] : replacement) {
        Mpoly lhs_full = Mpoly::from_string(full_ctx, lhs_str);
        Mpoly lhs_fc   = compose_mpoly(lhs_full, fc.ctx, gen_map);

        fmpz_t scale;
        fmpz_init(scale);
        Mpoly rhs_full_p = Mpoly::from_string(full_ctx, rhs_str, scale);
        Mpoly rhs_fc_p   = compose_mpoly(rhs_full_p, fc.ctx, gen_map);
        Mfrac rhs_fc     = Mfrac::from_mpoly(std::move(rhs_fc_p));
        if (!fmpz_is_one(scale)) {
            fmpq_t qscale;
            fmpq_init(qscale);
            fmpz_one(fmpq_numref(qscale));
            fmpz_set(fmpq_denref(qscale), scale);
            Mfrac mul = Mfrac::from_fmpq(fc.ctx, qscale);
            rhs_fc *= mul;
            fmpq_clear(qscale);
        }
        fmpz_clear(scale);

        rep_in_fc.emplace_back(std::move(lhs_fc), std::move(rhs_fc));
    }

    long n_loop = static_cast<long>(fc.loops.size());
    long n_red  = static_cast<long>(fc.reduced_legs.size());
    long N      = fc.ctx->n_vars();
    std::vector<bool> resolved(rep_in_fc.size(), false);
    std::size_t remaining = rep_in_fc.size();
    while (remaining != 0) {
        bool progress = false;

        for (std::size_t idx = 0; idx < rep_in_fc.size(); ++idx) {
            if (resolved[idx]) continue;

            auto& lhs_fc = rep_in_fc[idx].first;
            auto rhs_remaining = rep_in_fc[idx].second.clone();

            const long len = fmpz_mpoly_length(lhs_fc.raw(), fc.ctx->raw());
            if (len == 0) {
                throw std::invalid_argument(
                    "replacement: lhs reduces to zero after conservation");
            }

            struct UnknownTerm {
                std::string key;
                Mfrac coeff;
            };
            std::vector<UnknownTerm> unknowns;
            unknowns.reserve(static_cast<std::size_t>(len));

            for (long term = 0; term < len; ++term) {
                std::vector<unsigned long> exp(N);
                fmpz_mpoly_get_term_exp_ui(
                    exp.data(), lhs_fc.raw(), term, fc.ctx->raw());

                long total_deg = 0;
                for (long i = 0; i < N; ++i) total_deg += static_cast<long>(exp[i]);
                if (total_deg != 2) {
                    throw std::invalid_argument(
                        "replacement: lhs must be quadratic in legs; got total degree "
                        + std::to_string(total_deg));
                }
                for (long i = 0; i < n_loop; ++i) {
                    if (exp[i] != 0) {
                        throw std::invalid_argument(
                            "replacement: lhs must not depend on loop momenta");
                    }
                }
                for (long i = n_loop + n_red; i < N; ++i) {
                    if (exp[i] != 0) {
                        throw std::invalid_argument(
                            "replacement: lhs must not depend on invariants");
                    }
                }

                std::vector<long> leg_indices_with_exp;
                long total_red_deg = 0;
                for (long i = 0; i < n_red; ++i) {
                    const long e = static_cast<long>(exp[n_loop + i]);
                    total_red_deg += e;
                    for (long k = 0; k < e; ++k) {
                        leg_indices_with_exp.push_back(i);
                    }
                }
                if (total_red_deg != 2 || leg_indices_with_exp.size() != 2) {
                    throw std::invalid_argument(
                        "replacement: lhs reduces to a non-leg-leg monomial");
                }

                std::string key = leg_monomial_key(
                    fc.reduced_legs[leg_indices_with_exp[0]],
                    fc.reduced_legs[leg_indices_with_exp[1]]);

                fmpz_t coeff_z;
                fmpz_init(coeff_z);
                fmpz_mpoly_get_term_coeff_fmpz(
                    coeff_z, lhs_fc.raw(), term, fc.ctx->raw());
                fmpq_t coeff_q;
                fmpq_init(coeff_q);
                fmpz_set(fmpq_numref(coeff_q), coeff_z);
                fmpz_one(fmpq_denref(coeff_q));
                fmpq_canonicalise(coeff_q);
                Mfrac coeff = Mfrac::from_fmpq(fc.ctx, coeff_q);
                fmpq_clear(coeff_q);
                fmpz_clear(coeff_z);

                auto known = fc.reduced_replacement.find(key);
                if (known != fc.reduced_replacement.end()) {
                    Mfrac contrib = known->second.clone();
                    contrib *= coeff;
                    rhs_remaining -= contrib;
                } else {
                    unknowns.push_back({std::move(key), std::move(coeff)});
                }
            }

            if (unknowns.empty()) {
                if (!rhs_remaining.is_zero()) {
                    throw std::invalid_argument(
                        "replacement: lhs '" + lhs_fc.to_string()
                        + "' conflicts with already resolved rules");
                }
                resolved[idx] = true;
                --remaining;
                progress = true;
                continue;
            }

            if (unknowns.size() == 1u) {
                rhs_remaining /= unknowns[0].coeff;
                fc.reduced_replacement.emplace(
                    std::move(unknowns[0].key), std::move(rhs_remaining));
                resolved[idx] = true;
                --remaining;
                progress = true;
            }
        }

        if (progress) continue;

        for (std::size_t idx = 0; idx < rep_in_fc.size(); ++idx) {
            if (resolved[idx]) continue;
            throw std::invalid_argument(
                "replacement: lhs '" + rep_in_fc[idx].first.to_string()
                + "' leaves multiple unknown leg-leg monomials; please add "
                  "enough replacement relations to isolate each scalar product");
        }
    }

    fc.build_sp_list();

    if (cut.empty()) {
        fc.cut.assign(fc.propagators_after_conservation.size(), 0);
    } else {
        if (cut.size() != fc.propagators_after_conservation.size()) {
            throw std::invalid_argument(
                "FamilyConfig::build: cut length must match #propagators");
        }
        for (int v : cut) {
            if (v != 0 && v != 1) {
                throw std::invalid_argument(
                    "FamilyConfig::build: cut entries must be 0 or 1");
            }
        }
        fc.cut = std::move(cut);
    }
    if (prescription.empty()) {
        fc.prescription.assign(fc.loops.size(), 1);
    } else {
        if (prescription.size() != fc.loops.size()) {
            throw std::invalid_argument(
                "FamilyConfig::build: prescription length must match #loops");
        }
        for (int v : prescription) {
            if (v != -1 && v != 0 && v != 1) {
                throw std::invalid_argument(
                    "FamilyConfig::build: prescription entries must be in {-1, 0, 1}");
            }
        }
        fc.prescription = std::move(prescription);
    }

    return fc;
}

void FamilyConfig::build_context() {
    std::vector<std::string> all_vars;
    all_vars.reserve(loops.size() + reduced_legs.size() + invariants.size());
    for (const auto& s : loops)        all_vars.push_back(s);
    for (const auto& s : reduced_legs) all_vars.push_back(s);
    for (const auto& s : invariants)   all_vars.push_back(s);
    if (all_vars.empty()) {
        throw std::invalid_argument(
            "FamilyConfig: at least one variable required (loops/legs/invariants)");
    }
    ctx = std::make_shared<MpolyContext>(std::move(all_vars));
}

void FamilyConfig::build_sp_list() {
    long n_loop = static_cast<long>(loops.size());
    long n_red  = static_cast<long>(reduced_legs.size());
    sp_list.clear();
    for (long i = 0; i < n_loop; ++i) {
        for (long j = i; j < n_loop; ++j) {
            Mpoly term = Mpoly::variable(ctx, i) * Mpoly::variable(ctx, j);
            sp_list.push_back(std::move(term));
        }
    }
    for (long i = 0; i < n_loop; ++i) {
        for (long j = 0; j < n_red; ++j) {
            Mpoly term = Mpoly::variable(ctx, i) * Mpoly::variable(ctx, n_loop + j);
            sp_list.push_back(std::move(term));
        }
    }
}

int FamilyConfig::prescription_of_loop(long loop_idx) const noexcept {
    if (loop_idx < 0 || (std::size_t)loop_idx >= prescription.size()) {
        return 1;
    }
    return prescription[(std::size_t)loop_idx];
}

std::optional<int> FamilyConfig::prescription_of_prop(const Mpoly& prop) const {
    if (prop.ctx().get() != ctx.get()) {
        throw std::invalid_argument(
            "FamilyConfig::prescription_of_prop: ctx mismatch");
    }

    long n_loop = (long)loops.size();
    long len = fmpz_mpoly_length(prop.raw(), ctx->raw());
    std::vector<unsigned long> exp((std::size_t)ctx->n_vars());
    std::vector<bool> seen_loop((std::size_t)n_loop, false);
    for (long t = 0; t < len; ++t) {
        fmpz_mpoly_get_term_exp_ui(exp.data(), prop.raw(), t, ctx->raw());
        for (long k = 0; k < n_loop; ++k) {
            if (exp[(std::size_t)k] > 0) {
                seen_loop[(std::size_t)k] = true;
            }
        }
    }

    std::vector<int> pres_list;
    for (long k = 0; k < n_loop; ++k) {
        if (seen_loop[(std::size_t)k]) {
            pres_list.push_back(prescription_of_loop(k));
        }
    }

    if (pres_list.empty()) return 0;
    bool all_zero = true;
    for (int p : pres_list) if (p != 0) { all_zero = false; break; }
    if (all_zero) return 0;

    bool any_pos = false, any_neg = false;
    for (int p : pres_list) {
        if (p > 0) any_pos = true;
        else if (p < 0) any_neg = true;
    }
    if (any_pos && !any_neg) return 1;
    if (!any_pos && any_neg) return -1;
    return std::nullopt;
}

Mpoly FamilyConfig::apply_conservation(const Mpoly& p) const {
    if (p.ctx().get() != ctx.get()) {
        throw std::invalid_argument(
            "FamilyConfig::apply_conservation: ctx mismatch");
    }
    return p.clone();
}

Mfrac FamilyConfig::apply_replacement(const Mpoly& p) const {
    return apply_replacement(Mfrac::from_mpoly(p.clone()));
}

Mfrac FamilyConfig::apply_replacement(const Mfrac& p) const {
    if (reduced_replacement.empty()) return p.clone();

    long n_loop = static_cast<long>(loops.size());
    long n_red  = static_cast<long>(reduced_legs.size());

    auto walk_mpoly = [&](const Mpoly& src) -> Mfrac {
        Mfrac acc = Mfrac::zero(ctx);
        long len = fmpz_mpoly_length(src.raw(), ctx->raw());
        std::vector<unsigned long> exp(ctx->n_vars());
        fmpz_t coeff;
        fmpz_init(coeff);
        for (long t = 0; t < len; ++t) {
            fmpz_mpoly_get_term_coeff_fmpz(coeff, src.raw(), t, ctx->raw());
            fmpz_mpoly_get_term_exp_ui(exp.data(), src.raw(), t, ctx->raw());

            std::vector<long> leg_powers;
            for (long i = 0; i < n_red; ++i) {
                long e = static_cast<long>(exp[n_loop + i]);
                for (long k = 0; k < e; ++k) leg_powers.push_back(i);
            }

            Mfrac leg_replacement = Mfrac::one(ctx);
            std::size_t paired = 0;
            for (; paired + 1 < leg_powers.size(); paired += 2) {
                std::string key = leg_monomial_key(
                    reduced_legs[leg_powers[paired]],
                    reduced_legs[leg_powers[paired + 1]]);
                auto it = reduced_replacement.find(key);
                if (it == reduced_replacement.end()) {
                    Mpoly mono = Mpoly::variable(
                        ctx, n_loop + leg_powers[paired]);
                    mono *= Mpoly::variable(
                        ctx, n_loop + leg_powers[paired + 1]);
                    leg_replacement *= Mfrac::from_mpoly(std::move(mono));
                } else {
                    leg_replacement *= it->second;
                }
            }
            for (; paired < leg_powers.size(); ++paired) {
                Mpoly mono = Mpoly::variable(ctx, n_loop + leg_powers[paired]);
                leg_replacement *= Mfrac::from_mpoly(std::move(mono));
            }

            std::vector<unsigned long> rest_exp = exp;
            for (long i = 0; i < n_red; ++i) rest_exp[n_loop + i] = 0;
            Mpoly rest_mono = Mpoly::zero(ctx);
            fmpz_mpoly_set_coeff_fmpz_ui(rest_mono.raw(), coeff,
                                         rest_exp.data(), ctx->raw());

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
            "FamilyConfig::apply_replacement: denominator vanished after substitution");
    }
    Mfrac quot(new_num.numerator() * new_den.denominator(),
               new_den.numerator() * new_num.denominator());
    return quot;
}

}  // namespace amflow::qft
