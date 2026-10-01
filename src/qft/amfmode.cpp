// SPDX-License-Identifier: MIT
// qft::amfmode — implementation.
//

#include "amflow/qft/amfmode.hpp"

#include <algorithm>
#include <functional>
#include <map>
#include <set>
#include <stdexcept>

#include <flint/fmpz.h>
#include <flint/fmpz_mpoly.h>

namespace amflow::qft {

using algebra::Mfrac;
using algebra::Mpoly;

const char* amf_mode_name(AMFMode m) noexcept {
    switch (m) {
        case AMFMode::Prescription: return "Prescription";
        case AMFMode::Mass:         return "Mass";
        case AMFMode::Propagator:   return "Propagator";
        case AMFMode::Branch:       return "Branch";
        case AMFMode::Loop:         return "Loop";
        case AMFMode::All:          return "All";
        case AMFMode::Bubble:       return "Bubble";
    }
    return "?";
}

AMFMode parse_amf_mode(const std::string& name) {
    if (name == "Prescription") return AMFMode::Prescription;
    if (name == "Mass")         return AMFMode::Mass;
    if (name == "Propagator")   return AMFMode::Propagator;
    if (name == "Branch")       return AMFMode::Branch;
    if (name == "Loop")         return AMFMode::Loop;
    if (name == "All")          return AMFMode::All;
    if (name == "Bubble")       return AMFMode::Bubble;
    throw std::invalid_argument(
        "parse_amf_mode: unknown AMFMode '" + name + "'.  Built-in modes are "
        "\"Prescription\", \"Mass\", \"Propagator\", \"Branch\", \"Loop\", "
        "\"All\", \"Bubble\".");
}

std::vector<TopSectorComponentInfo>
analyze_top_sector(const FamilyConfig& fc,
                   const std::vector<std::size_t>& topposi) {

    if (topposi.empty()) return {};

    std::vector<Mpoly> denominators;
    denominators.reserve(topposi.size());
    for (std::size_t k = 0; k < topposi.size(); ++k) {
        std::size_t i = topposi[k];
        if (i >= fc.propagators_after_conservation.size()) {
            throw std::out_of_range(
                "analyze_top_sector: propagator index out of range");
        }
        denominators.push_back(
            fc.propagators_after_conservation[i].clone());
    }

    std::vector<int> pres_full;
    pres_full.reserve(topposi.size());
    for (std::size_t k = 0; k < topposi.size(); ++k) {
        auto p = fc.prescription_of_prop(denominators[k]);
        if (!p.has_value()) {
            throw std::runtime_error(
                "analyze_top_sector: propagator at global index "
                + std::to_string(topposi[k])
                + " has conflicting prescription on its loops.");
        }
        pres_full.push_back(p.value());
    }

    std::vector<int> cut_full;
    cut_full.reserve(topposi.size());
    for (std::size_t k = 0; k < topposi.size(); ++k) {
        std::size_t i = topposi[k];
        if (i < fc.cut.size()) {
            cut_full.push_back(fc.cut[i]);
        } else {
            cut_full.push_back(0);
        }
    }

    auto comps = analyze_topology(fc, denominators);

    std::vector<TopSectorComponentInfo> out;
    out.reserve(comps.size());

    for (auto& ci : comps) {
        TopSectorComponentInfo ts;
        ts.u0      = ci.u0.clone();
        ts.loopnum = ci.loopnum;
        ts.vacQ    = ci.vacQ;

        long ctx_n = ci.u0.ctx()->n_vars();
        long n_x   = (long)topposi.size();
        long first_x = ctx_n - n_x;

        ts.var.reserve(ci.var.size());
        ts.prop_index.reserve(ci.var.size());
        ts.cut.reserve(ci.var.size());
        ts.pres.reserve(ci.var.size());
        for (long v : ci.var) {
            long local = v - first_x;
            if (local < 0 || local >= n_x) {
                throw std::runtime_error(
                    "analyze_top_sector: var index out of range");
            }
            ts.var.push_back(v);
            std::size_t global_pi = topposi[(std::size_t)local];
            ts.prop_index.push_back(global_pi);
            ts.cut.push_back(cut_full[(std::size_t)local]);
            ts.pres.push_back(pres_full[(std::size_t)local]);
        }

        ts.mass.reserve(ci.mass.size());
        for (auto& m : ci.mass) ts.mass.push_back(m.clone());

        out.push_back(std::move(ts));
    }

    return out;
}

bool vacuum_q(const TopSectorComponentInfo& info) {
    if (!info.vacQ) return false;
    for (int c : info.cut) {
        if (c != 0) return false;
    }
    return true;
}

bool single_mass_q(const TopSectorComponentInfo& info) {
    if (!vacuum_q(info)) return false;
    int n_one  = 0;
    int n_zero = 0;
    for (const auto& m : info.mass) {
        if (m.is_zero()) {
            ++n_zero;
        } else if (m.is_one()) {
            ++n_one;
        }
    }
    return n_one == 1 && n_zero == (int)info.mass.size() - 1;
}

bool phase_volume_q(const TopSectorComponentInfo& info) {
    for (int c : info.cut) {
        if (c != 1) return false;
    }
    return (long)info.var.size() == info.loopnum + 1;
}

bool ending_q(const TopSectorComponentInfo& info) {
    return single_mass_q(info) || phase_volume_q(info);
}

namespace {

std::vector<std::vector<std::size_t>>
put_first(std::vector<std::size_t>&& flat) {
    if (flat.empty()) return {};
    return {std::move(flat)};
}

std::vector<std::size_t>
local_indices_with_cut(const TopSectorComponentInfo& info, int cut_val) {
    std::vector<std::size_t> out;
    for (std::size_t k = 0; k < info.cut.size(); ++k) {
        if (info.cut[k] == cut_val) out.push_back(k);
    }
    return out;
}

// Detect 1-loop bubble sub-loops across the analyzed top sector: a loop
// momentum l_j that appears (quadratically — the AMFlow propagator
// convention is `(Σ q)² − m`, so l_j ∈ prop ⇔ coeff of l_j² ≠ 0) in
// EXACTLY two of the active top-sector propagators.  After integrating
// out l_j the result is a 1-loop bubble B0(P²; m_a, m_b) whose external
// invariant P² may depend on the *remaining* loop momenta — that is the
// case where U does NOT factor and `analyze_topology` therefore returns
// a single connected component.  (When P is purely external the bubble's
// U-factor `x_a + x_b` splits off and the existing vacuum-component skip
// already handles it; when l_j hits one prop it is a tadpole and U
// factors as `x_a · U_rest`.)
//
// Returns the global propagator indices to EXCLUDE from η-injection so
// the bubble stays at fixed mass and Kira's η-DE seed depth is set by
// the (L-1)-loop remainder only.  Never returns the entire top sector
// (η must land somewhere); never flags a bubble that is its own
// component (a 1-loop component has nothing left to η-inject and is
// already at the recursion floor).
std::set<std::size_t>
bubble_subloop_props(const FamilyConfig& fc,
                     const std::vector<TopSectorComponentInfo>& info_list) {
    std::set<std::size_t> top;
    for (const auto& ci : info_list) {
        for (auto pi : ci.prop_index) top.insert(pi);
    }
    std::set<std::size_t> out;
    const long L = (long)fc.n_loops();
    for (long j = 0; j < L; ++j) {
        std::vector<std::size_t> hits;
        for (std::size_t pi : top) {
            const Mpoly& d = fc.propagators_after_conservation[pi];
            if (!d.coeff_of(j, 2).is_zero()) hits.push_back(pi);
        }
        if (hits.size() != 2) continue;
        // Is this bubble its own already-factored component?  If so
        // `analyze_topology` produced a 1-loop component for it and the
        // vacuum/ending logic handles it — don't re-flag.
        bool own_component = false;
        for (const auto& ci : info_list) {
            if (ci.loopnum == 1 && ci.prop_index.size() == 2
                    && std::find(ci.prop_index.begin(), ci.prop_index.end(),
                                  hits[0]) != ci.prop_index.end()
                    && std::find(ci.prop_index.begin(), ci.prop_index.end(),
                                  hits[1]) != ci.prop_index.end()) {
                own_component = true;
                break;
            }
        }
        if (own_component) continue;
        for (auto pi : hits) out.insert(pi);
    }
    // Never strip everything — η has to go somewhere.
    if (out.size() >= top.size()) return {};
    return out;
}

std::vector<long>
same_branch_with(const TopSectorComponentInfo& info, long v) {
    Mpoly co = info.u0.numerator().coeff_of(v, 1);
    auto ctx = info.u0.ctx();
    long len = fmpz_mpoly_length(co.raw(), ctx->raw());
    std::vector<unsigned long> exp((std::size_t)ctx->n_vars());
    std::set<long> appearing;
    for (long t = 0; t < len; ++t) {
        fmpz_mpoly_get_term_exp_ui(exp.data(), co.raw(), t, ctx->raw());
        for (long w : info.var) {
            if (exp[(std::size_t)w] > 0) {
                appearing.insert(w);
            }
        }
    }
    std::vector<long> result;
    for (long w : info.var) {
        if (!appearing.count(w)) result.push_back(w);
    }
    std::sort(result.begin(), result.end());
    return result;
}

}  // namespace

std::vector<std::vector<std::size_t>>
all_possible_position(const FamilyConfig& fc,
                      const TopSectorComponentInfo& info,
                      AMFMode mode) {
    switch (mode) {
        case AMFMode::Prescription: {
            std::vector<std::size_t> flat;
            for (std::size_t k = 0; k < info.var.size(); ++k) {
                if (info.pres[k] == -1 && info.cut[k] == 0) {
                    flat.push_back(info.prop_index[k]);
                }
            }
            return put_first(std::move(flat));
        }

        case AMFMode::Propagator: {
            auto locals = local_indices_with_cut(info, 0);

            std::vector<std::pair<std::size_t, long>> with_len;
            with_len.reserve(locals.size());
            for (std::size_t k : locals) {
                long bl = (long)same_branch_with(info, info.var[k]).size();
                with_len.emplace_back(k, bl);
            }
            std::stable_sort(with_len.begin(), with_len.end(),
                              [](const auto& a, const auto& b) {
                                  return a.second > b.second;
                              });
            std::vector<std::vector<std::size_t>> out;
            out.reserve(with_len.size());
            for (auto& [k, _] : with_len) {
                out.push_back({info.prop_index[k]});
            }
            return out;
        }

        case AMFMode::Branch: {
            std::set<std::vector<long>> branches_set;
            std::vector<std::vector<long>> branches;
            for (long v : info.var) {
                auto b = same_branch_with(info, v);
                if (branches_set.insert(b).second) {
                    branches.push_back(std::move(b));
                }
            }

            std::set<long> cut0_vars;
            for (std::size_t k = 0; k < info.var.size(); ++k) {
                if (info.cut[k] == 0) cut0_vars.insert(info.var[k]);
            }
            std::vector<std::vector<long>> branches_filtered;
            for (auto& b : branches) {
                std::vector<long> filt;
                for (long v : b) {
                    if (cut0_vars.count(v)) filt.push_back(v);
                }
                if (!filt.empty()) branches_filtered.push_back(std::move(filt));
            }

            std::stable_sort(branches_filtered.begin(),
                              branches_filtered.end(),
                              [](const auto& a, const auto& b) {
                                  return a.size() <= b.size()
                                      ? a.size() < b.size()
                                      : false;
                              });

            std::map<long, std::size_t> var_to_local;
            for (std::size_t k = 0; k < info.var.size(); ++k) {
                var_to_local[info.var[k]] = k;
            }
            std::vector<std::vector<std::size_t>> out;
            for (auto& b : branches_filtered) {
                std::vector<std::size_t> g;
                for (long v : b) {
                    g.push_back(info.prop_index[var_to_local[v]]);
                }
                out.push_back(std::move(g));
            }
            return out;
        }

        case AMFMode::Loop: {
            long n_var = (long)info.var.size();
            long ln    = info.loopnum;
            if (ln <= 0 || n_var <= 0) return {};

            Mpoly u_num = info.u0.numerator();
            auto ctx = u_num.ctx();
            long nterms = fmpz_mpoly_length(u_num.raw(), ctx->raw());
            if (nterms == 0) return {};

            struct MonoExp { std::map<long, unsigned long> exp; };
            std::vector<MonoExp> monos((std::size_t)nterms);
            std::vector<unsigned long> exp((std::size_t)ctx->n_vars());
            for (long t = 0; t < nterms; ++t) {
                fmpz_mpoly_get_term_exp_ui(exp.data(), u_num.raw(),
                                            t, ctx->raw());
                for (long v : info.var) {
                    if (exp[(std::size_t)v] > 0) {
                        monos[(std::size_t)t].exp[v] = exp[(std::size_t)v];
                    }
                }
            }

            std::vector<std::vector<long>> subs;
            std::vector<long> tmp;
            std::function<void(std::size_t, std::size_t)> enumerate =
                [&](std::size_t start, std::size_t left) {
                    if (left == 0) { subs.push_back(tmp); return; }
                    for (std::size_t i = start;
                         i + left <= info.var.size(); ++i) {
                        tmp.push_back(info.var[i]);
                        enumerate(i + 1, left - 1);
                        tmp.pop_back();
                    }
                };
            if ((std::size_t)(ln - 1) > info.var.size()) return {};
            enumerate(0, (std::size_t)(ln - 1));

            std::set<std::vector<long>> loop_groups_set;
            std::vector<std::vector<long>> loop_groups;
            for (auto& sub : subs) {
                std::set<long> remaining;
                for (auto& m : monos) {
                    bool ok = true;
                    for (long v : sub) {
                        auto it = m.exp.find(v);
                        if (it == m.exp.end() || it->second < 1) {
                            ok = false; break;
                        }
                    }
                    if (!ok) continue;
                    for (auto& [v, e] : m.exp) {
                        unsigned long sub_e = 0;
                        for (long s : sub) if (s == v) sub_e++;
                        if (e > sub_e) {
                            remaining.insert(v);
                        }
                    }
                }
                if (remaining.empty()) continue;
                std::vector<long> g(remaining.begin(), remaining.end());
                std::sort(g.begin(), g.end());
                if (loop_groups_set.insert(g).second) {
                    loop_groups.push_back(std::move(g));
                }
            }

            std::set<long> cut0_vars;
            for (std::size_t k = 0; k < info.var.size(); ++k) {
                if (info.cut[k] == 0) cut0_vars.insert(info.var[k]);
            }
            std::vector<std::vector<long>> filtered;
            for (auto& g : loop_groups) {
                std::vector<long> ff;
                for (long v : g) {
                    if (cut0_vars.count(v)) ff.push_back(v);
                }
                if (!ff.empty()) filtered.push_back(std::move(ff));
            }

            std::stable_sort(filtered.begin(), filtered.end(),
                              [](const auto& a, const auto& b) {
                                  return a.size() < b.size();
                              });

            std::map<long, std::size_t> var_to_local;
            for (std::size_t k = 0; k < info.var.size(); ++k) {
                var_to_local[info.var[k]] = k;
            }
            std::vector<std::vector<std::size_t>> out;
            for (auto& g : filtered) {
                std::vector<std::size_t> gg;
                for (long v : g) {
                    gg.push_back(info.prop_index[var_to_local[v]]);
                }
                out.push_back(std::move(gg));
            }
            return out;
        }

        case AMFMode::Mass: {
            auto locals = local_indices_with_cut(info, 0);

            std::vector<std::pair<std::string, std::vector<std::size_t>>>
                by_mass;
            std::map<std::string, std::size_t> mass_group_index;
            for (std::size_t k : locals) {
                if (info.mass[k].is_zero()) continue;
                std::string key = info.mass[k].to_string();
                auto it = mass_group_index.find(key);
                if (it == mass_group_index.end()) {
                    mass_group_index.emplace(key, by_mass.size());
                    by_mass.emplace_back(std::move(key),
                                         std::vector<std::size_t>{k});
                } else {
                    by_mass[it->second].second.push_back(k);
                }
            }

            std::set<std::string> rr_vars;
            for (const auto& [_, rhs] : fc.reduced_replacement) {
                Mpoly num = rhs.numerator();
                Mpoly den = rhs.denominator();
                auto ctx = num.ctx();
                std::vector<unsigned long> exp((std::size_t)ctx->n_vars());
                auto walk = [&](const Mpoly& p) {
                    long len = fmpz_mpoly_length(p.raw(), ctx->raw());
                    for (long t = 0; t < len; ++t) {
                        fmpz_mpoly_get_term_exp_ui(exp.data(), p.raw(),
                                                    t, ctx->raw());
                        for (long v = 0; v < ctx->n_vars(); ++v) {
                            if (exp[(std::size_t)v] > 0) {
                                rr_vars.insert(ctx->var_name(v));
                            }
                        }
                    }
                };
                walk(num);
                walk(den);
            }

            std::vector<std::pair<std::string, std::vector<std::size_t>>> groups;
            for (auto& [mstr, locs] : by_mass) {
                const Mfrac& m = info.mass[locs.front()];
                Mpoly num = m.numerator();
                Mpoly den = m.denominator();
                auto ctx = num.ctx();
                std::vector<unsigned long> exp((std::size_t)ctx->n_vars());
                bool depend = false;
                auto walk = [&](const Mpoly& p) {
                    if (depend) return;
                    long len = fmpz_mpoly_length(p.raw(), ctx->raw());
                    for (long t = 0; t < len && !depend; ++t) {
                        fmpz_mpoly_get_term_exp_ui(exp.data(), p.raw(),
                                                    t, ctx->raw());
                        for (long v = 0; v < ctx->n_vars(); ++v) {
                            if (exp[(std::size_t)v] > 0) {
                                if (rr_vars.count(ctx->var_name(v))) {
                                    depend = true;
                                    break;
                                }
                            }
                        }
                    }
                };
                walk(num);
                walk(den);
                if (!depend) {
                    groups.emplace_back(mstr, locs);
                }
            }

            std::stable_sort(groups.begin(), groups.end(),
                              [](const auto& a, const auto& b) {
                                  return a.second.size() <= b.second.size()
                                      ? a.second.size() < b.second.size()
                                      : false;
                              });

            std::vector<std::vector<std::size_t>> out;
            out.reserve(groups.size());
            for (auto& [_, locs] : groups) {
                std::vector<std::size_t> gg;
                for (std::size_t k : locs) {
                    gg.push_back(info.prop_index[k]);
                }
                out.push_back(std::move(gg));
            }
            return out;
        }

        case AMFMode::Bubble:
            // Handled one level up in `amf_candidate` (needs the full
            // component list to detect connected bubble sub-loops).
            // Per-component this degrades to Prescription so any direct
            // caller still gets a sane answer.
            return all_possible_position(fc, info, AMFMode::Prescription);

        case AMFMode::All: {
            auto props = all_possible_position(fc, info, AMFMode::Propagator);
            std::vector<std::size_t> flat;
            for (auto& g : props) {
                for (std::size_t v : g) flat.push_back(v);
            }
            std::set<std::size_t> seen;
            std::vector<std::size_t> uniq;
            for (std::size_t v : flat) {
                if (seen.insert(v).second) uniq.push_back(v);
            }
            return put_first(std::move(uniq));
        }
    }
    return {};
}

std::vector<std::size_t>
amf_candidate_component(const FamilyConfig& fc,
                        const TopSectorComponentInfo& info,
                        AMFMode mode) {
    if (!vacuum_q(info)) {
        auto groups = all_possible_position(fc, info, mode);
        if (groups.empty()) return {};
        return groups.front();
    }
    if (!single_mass_q(info)) {
        auto groups = all_possible_position(fc, info, AMFMode::Branch);
        if (groups.empty()) return {};
        return groups.front();
    }
    return {};
}

std::vector<std::size_t>
amf_candidate(const FamilyConfig& fc,
              const std::vector<TopSectorComponentInfo>& info_list,
              AMFMode mode) {
    // Mirrors upstream `AMFCandidate` (AMFlow.m:614-617).  The Prescription
    // and All modes UNION candidates across components but explicitly
    // restrict the iteration to non-vacuum components:
    //
    //   If[mode==="Prescription" || mode==="All",
    //     candidates = AMFCandidateComponent[#, mode]&/@Select[info, !VacuumQ[#]&];
    //     candidates = PutFirst[Join@@candidates],
    //     candidates = AMFCandidateComponent[#, mode]&/@info]
    //
    // Without the vacuum filter, a vacuum-but-not-single-mass component
    // (e.g. a 1-loop tadpole with a symbolic mass) falls through to the
    // Branch fallback inside `amf_candidate_component` and contributes a
    // non-empty candidate, polluting the Prescription/All result.  For
    // banana_4L_mixed's first boundary sub-family (4 disconnected
    // 1-loop tadpoles with masses {mAsq, 1, mAsq, 1}), this produced
    // pos={0, 2} where MMA's filtered iteration produced pos={} and the
    // Mass mode (which iterates all components) was then tried and
    // returned the correct pos={0} per first non-empty component.
    if (mode == AMFMode::Bubble) {
        // C++-port extension.  Run the upstream default fallback chain
        // {Prescription, Mass, Propagator} to get the η-candidate set
        // exactly as stock AMFlow would, then DROP any propagator that
        // belongs to a detected 2-prop bubble sub-loop (see
        // `bubble_subloop_props` for the criterion and rationale).
        // When no bubble is detected — or when dropping the bubble
        // props would empty the candidate set — this is byte-identical
        // to the stock chain, so `amf_modes: ["Bubble"]` is a safe
        // superset of the upstream default.
        static const AMFMode kDefaultChain[] = {
            AMFMode::Prescription, AMFMode::Mass, AMFMode::Propagator
        };
        std::vector<std::size_t> base;
        for (AMFMode m : kDefaultChain) {
            base = amf_candidate(fc, info_list, m);
            if (!base.empty()) break;
        }
        if (base.empty()) return base;
        auto bubble = bubble_subloop_props(fc, info_list);
        if (bubble.empty()) return base;
        std::vector<std::size_t> filtered;
        filtered.reserve(base.size());
        for (std::size_t x : base) {
            if (!bubble.count(x)) filtered.push_back(x);
        }
        return filtered.empty() ? base : filtered;
    }

    if (mode == AMFMode::Prescription || mode == AMFMode::All) {
        std::vector<std::size_t> flat;
        std::set<std::size_t> seen;
        for (const auto& ci : info_list) {
            if (vacuum_q(ci)) continue;
            auto cand = amf_candidate_component(fc, ci, mode);
            for (std::size_t x : cand) {
                if (seen.insert(x).second) flat.push_back(x);
            }
        }
        return flat;
    }

    for (const auto& ci : info_list) {
        auto cand = amf_candidate_component(fc, ci, mode);
        if (!cand.empty()) return cand;
    }
    return {};
}

std::vector<std::size_t>
amf_position(const FamilyConfig& fc,
             const std::vector<std::size_t>& topposi,
             AMFMode mode) {
    if (topposi.empty()) return {};
    auto info = analyze_top_sector(fc, topposi);
    return amf_candidate(fc, info, mode);
}

std::vector<std::size_t>
amf_position(const FamilyConfig& fc,
             const std::vector<std::size_t>& topposi,
             const std::vector<AMFMode>& modes) {
    if (topposi.empty()) return {};
    auto info = analyze_top_sector(fc, topposi);
    for (AMFMode mode : modes) {
        auto cand = amf_candidate(fc, info, mode);
        if (!cand.empty()) return cand;
    }
    return {};
}

std::vector<int>
amf_eta_c(const FamilyConfig& fc,
          const std::vector<std::size_t>& posi) {
    std::vector<int> out(fc.propagators_after_conservation.size(), 0);
    for (std::size_t i : posi) {
        if (i >= out.size()) {
            throw std::out_of_range(
                "amf_eta_c: propagator index out of range");
        }
        out[i] = -1;
    }
    return out;
}

}  // namespace amflow::qft
