// SPDX-License-Identifier: MIT
// qft::findregion — FindAllRegion + ZeroRegionQ + RegionPower.
//

#include "amflow/qft/region.hpp"

#include <algorithm>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <map>
#include <set>
#include <stdexcept>

#include <flint/fmpq.h>
#include <flint/fmpz.h>
#include <flint/fmpz_mpoly.h>
#include <flint/fmpz_mpoly_q.h>
#include <flint/flint.h>

#include "amflow/numeric/log.hpp"
#include "amflow/qft/dense_q.hpp"
#include "amflow/qft/jintegral.hpp"
#include "amflow/qft/topology.hpp"

namespace amflow::qft {

using algebra::Mfrac;
using algebra::Mpoly;
using algebra::MpolyContext;

PowersContext make_powers_context(const FamilyConfig& fc) {
    std::vector<std::string> names;
    names.reserve((std::size_t)fc.ctx->n_vars() + 1);
    for (long i = 0; i < fc.ctx->n_vars(); ++i) {
        names.push_back(fc.ctx->var_name(i));
    }
    names.emplace_back("__amf_eps");
    auto p_ctx = std::make_shared<MpolyContext>(std::move(names));

    PowersContext pctx;
    pctx.ctx     = p_ctx;
    pctx.eps_var = fc.ctx->n_vars();
    pctx.lift_gens.reserve((std::size_t)fc.ctx->n_vars());
    for (long i = 0; i < fc.ctx->n_vars(); ++i) {
        pctx.lift_gens.push_back(Mpoly::variable(p_ctx, i));
    }
    return pctx;
}

namespace {

long max_degree_in_var(const Mpoly& p, long var) {
    auto ctx = p.ctx();
    long len = fmpz_mpoly_length(p.raw(), ctx->raw());
    std::vector<unsigned long> exp((std::size_t)ctx->n_vars());
    long d = 0;
    for (long t = 0; t < len; ++t) {
        fmpz_mpoly_get_term_exp_ui(exp.data(), p.raw(), t, ctx->raw());
        if ((long)exp[(std::size_t)var] > d) d = (long)exp[(std::size_t)var];
    }
    return d;
}

bool depends_on_var(const Mpoly& p, long var) {
    return max_degree_in_var(p, var) > 0;
}

// Eikonal region engine gate.  When
// AMFLOW_EIKONAL_REGIONS=1, find_all_region enumerates per-loop scale bits
// {0,1,2} (sigma in {0, 1/2, 1}; bit 2 == l ~ eta) instead of {0,1}, so the
// l ~ eta eikonal regions of families with eta-carrying linear (worldline)
// propagators are represented.  Default (env unset/0): behavior is
// bit-identical to the previous {0,1} enumeration.
// Math: the half-eta boundary-defect derivation, summarized in
// the guard comment below.
bool eikonal_regions_enabled() {
    const char* e = std::getenv("AMFLOW_EIKONAL_REGIONS");
    return e && *e && std::string(e) != "0";
}

Mfrac merge_eta_into_half_eta(const FamilyConfig& fc,
                              const RegionContext& rctx,
                              const Mfrac& raw) {
    long fc_eta_var = fc.ctx->var_index("eta");
    if (fc_eta_var < 0) return raw.clone();

    auto sub_eta = [&](const Mpoly& p) -> Mpoly {
        Mpoly out(rctx.ctx);
        long len = fmpz_mpoly_length(p.raw(), rctx.ctx->raw());
        std::vector<unsigned long> exp((std::size_t)rctx.ctx->n_vars());
        fmpz_t coeff;
        fmpz_init(coeff);
        for (long t = 0; t < len; ++t) {
            fmpz_mpoly_get_term_exp_ui(exp.data(), p.raw(), t, rctx.ctx->raw());
            fmpz_mpoly_get_term_coeff_fmpz(coeff, p.raw(), t, rctx.ctx->raw());
            long ek = (long)exp[(std::size_t)fc_eta_var];
            std::vector<unsigned long> e2 = exp;
            e2[(std::size_t)fc_eta_var] = 0;
            e2[(std::size_t)rctx.half_eta_var] += (unsigned long)(2 * ek);
            fmpz_mpoly_set_coeff_fmpz_ui(out.raw(), coeff, e2.data(), rctx.ctx->raw());
        }
        fmpz_clear(coeff);
        return out;
    };

    return Mfrac(sub_eta(raw.numerator()), sub_eta(raw.denominator()));
}

}  // namespace

std::vector<Region>
find_all_region(const FamilyConfig& fc,
                const RegionContext& rctx,
                const std::vector<std::size_t>& topposi) {
    long L = (long)fc.n_loops();
    std::size_t N = fc.propagators_after_conservation.size();

    // half_eta boundary-defect guard.  The
    // enumeration below scales each loop by he = sqrt(eta) only (scale bits
    // {0,1} in region_rule).  For families with a propagator LINEAR in the
    // loop momenta (uncut worldline props 2u.k), the eta->inf asymptotics
    // ALSO contains l ~ eta regions (2l.u - eta is homogeneous there); they
    // are not representable here and are silently dropped, so the flow is
    // started from incomplete boundary data.  Verified in review and
    // proved in closed form (artifact O = -i(pi/2)Gamma(eps)4^eps per odd
    // 1L factor = truth minus the exact resummed l~eta region; half-eta
    // boundary-defect derivation).
    // Consequence: odd-linear-count and coupled eikonal-sector output is NOT
    // a Feynman value; quadratic-only and gated even-linear content is safe.
    // This guard does not change any computed value: it warns once per
    // family, and AMFLOW_HALF_ETA_STRICT=1 upgrades the warning to an error.
    //
    // Region-engine demotion: with AMFLOW_EIKONAL_REGIONS=1 the
    // l ~ eta regions ARE enumerated below, so this blanket warning is
    // skipped; the guard instead fires at the drop sites (two-velocity /
    // multi-loop eikonal components, i.e. the stage-2 ending classes that
    // are still unrepresentable) further down.
    if (!eikonal_regions_enabled()) {
        bool any_uncut_linear = false;
        std::vector<unsigned long> gexp((std::size_t)fc.ctx->n_vars());
        for (std::size_t i = 0; i < N && !any_uncut_linear; ++i) {
            if (i < fc.cut.size() && fc.cut[i] == 1) continue;
            const Mpoly& p = fc.propagators_after_conservation[i];
            long len = fmpz_mpoly_length(p.raw(), fc.ctx->raw());
            long maxdeg = 0;
            for (long t = 0; t < len; ++t) {
                fmpz_mpoly_get_term_exp_ui(gexp.data(), p.raw(), t,
                                            fc.ctx->raw());
                long d = 0;
                for (long v = 0; v < L; ++v) {
                    d += (long)gexp[(std::size_t)v];
                }
                if (d > maxdeg) maxdeg = d;
            }
            if (maxdeg == 1) any_uncut_linear = true;
        }
        if (any_uncut_linear) {
            static std::set<std::string> warned;
            if (warned.insert(fc.family).second) {
                std::cerr << "[half_eta guard] family '" << fc.family
                          << "': uncut linear (worldline) propagator in an "
                             "eta->inf region analysis.  The l~eta regions "
                             "are NOT enumerated (only l~sqrt(eta)); "
                             "odd-linear-count and coupled eikonal-sector "
                             "results are not Feynman values (artifact "
                             "-i(pi/2)Gamma(eps)4^eps "
                             "per odd 1L factor).  Enable "
                             "AMFLOW_EIKONAL_REGIONS=1 to enumerate them.\n";
            }
            const char* strict = std::getenv("AMFLOW_HALF_ETA_STRICT");
            if (strict && std::string(strict) == "1") {
                throw std::runtime_error(
                    "find_all_region: AMFLOW_HALF_ETA_STRICT=1 and family '"
                    + fc.family + "' has an uncut linear propagator: the "
                    "eta->inf boundary is missing its l~eta regions "
                    "(set AMFLOW_EIKONAL_REGIONS=1).");
            }
        }
    }

    auto branch0 = branch_momenta(fc, fc.propagators_after_conservation);

    std::vector<Mpoly> dedup;
    std::vector<std::string> dedup_keys;
    auto add_unique = [&](const Mpoly& b) {
        std::string k = b.to_string();
        for (const auto& dk : dedup_keys) {
            if (dk == k) return;
        }
        dedup_keys.push_back(k);
        dedup.push_back(b.clone());
    };
    for (std::size_t k = 0; k < topposi.size(); ++k) {
        add_unique(branch0[topposi[k]]);
    }

    std::set<std::size_t> cutposi;
    for (std::size_t i = 0; i < N; ++i) {
        if (i < fc.cut.size() && fc.cut[i] == 1) {
            std::string k = branch0[i].to_string();
            for (std::size_t d = 0; d < dedup_keys.size(); ++d) {
                if (dedup_keys[d] == k) cutposi.insert(d);
            }
        }
    }

    if (L <= 0 || dedup.empty()) {
        return {};
    }
    std::vector<std::pair<LoopTransform, std::vector<std::size_t>>> trans_list;
    {
        std::vector<std::size_t> idx_tuple((std::size_t)L, 0);
        std::function<void(long)> rec = [&](long depth) {
            if (depth == L) {
                std::vector<Mpoly> cand;
                cand.reserve((std::size_t)L);
                for (long i = 0; i < L; ++i) {
                    cand.push_back(dedup[idx_tuple[(std::size_t)i]].clone());
                }
                LoopTransform tr = branch_to_loop(fc, rctx, cand);
                if (tr.ok) {
                    trans_list.emplace_back(std::move(tr), idx_tuple);
                }
                return;
            }
            for (std::size_t i = 0; i < dedup.size(); ++i) {
                idx_tuple[(std::size_t)depth] = i;
                rec(depth + 1);
            }
        };
        rec(0);
    }

    struct RegionEntry {
        LoopTransform     tr;
        std::vector<int>  scale;
        std::vector<int>  bscale;
    };

    std::vector<RegionEntry> all_entries;

    const bool eik_on = eikonal_regions_enabled();
    const long fc_eta_var = fc.ctx->var_index("eta");
    const long n_red = (long)fc.n_red_legs();

    // Per-monomial content of a (loop-transformed) propagator, used by the
    // eikonal engine's prune and drop checks.  Only built when the engine
    // is enabled.
    struct MonoInfo {
        std::vector<int> ldeg;   // per-loop degree
        long             tot;    // total loop degree
        long             eta;    // exponent of the fc 'eta' variable
        std::vector<long> legs;  // leg variables present
    };
    struct PropInfo {
        std::vector<MonoInfo> monos;
        bool is_linear   = false;  // max total loop degree == 1
        bool carries_eta = false;
        std::vector<char> loops_present;
    };

    auto& dedup_ref = dedup;
    for (auto& [tr, _] : trans_list) {
        // Scale-independent per-transform data (identity rule, transformed
        // branch momenta and, for the eikonal engine, transformed
        // propagator monomial info).
        std::vector<int> identity_scale((std::size_t)L, 0);
        std::vector<Mfrac> identity_rule = region_rule(rctx, tr,
                                                        identity_scale);
        std::vector<Mfrac> transformed;
        transformed.reserve(dedup_ref.size());
        for (auto& d : dedup_ref) {
            transformed.push_back(
                apply_region_rule(fc, rctx, d, identity_rule));
        }

        std::vector<PropInfo> prop_infos;
        if (eik_on) {
            std::vector<unsigned long> pex((std::size_t)rctx.ctx->n_vars());
            for (std::size_t kk = 0; kk < topposi.size(); ++kk) {
                std::size_t i = topposi[kk];
                if (i < fc.cut.size() && fc.cut[i] == 1) continue;
                Mfrac t = apply_region_rule(
                    fc, rctx, fc.propagators_after_conservation[i],
                    identity_rule);
                PropInfo pi;
                pi.loops_present.assign((std::size_t)L, 0);
                const Mpoly& num = t.numerator();
                long len = fmpz_mpoly_length(num.raw(), rctx.ctx->raw());
                long max_tot = 0;
                for (long tt = 0; tt < len; ++tt) {
                    fmpz_mpoly_get_term_exp_ui(pex.data(), num.raw(), tt,
                                                rctx.ctx->raw());
                    MonoInfo mi;
                    mi.ldeg.assign((std::size_t)L, 0);
                    mi.tot = 0;
                    for (long v = 0; v < L; ++v) {
                        int d = (int)pex[(std::size_t)v];
                        mi.ldeg[(std::size_t)v] = d;
                        mi.tot += d;
                        if (d > 0) pi.loops_present[(std::size_t)v] = 1;
                    }
                    for (long v = L; v < L + n_red; ++v) {
                        if (pex[(std::size_t)v] > 0) mi.legs.push_back(v);
                    }
                    mi.eta = (fc_eta_var >= 0)
                        ? (long)pex[(std::size_t)fc_eta_var] : 0;
                    if (mi.eta > 0) pi.carries_eta = true;
                    if (mi.tot > max_tot) max_tot = mi.tot;
                    pi.monos.push_back(std::move(mi));
                }
                pi.is_linear = (max_tot == 1);
                prop_infos.push_back(std::move(pi));
            }
        }

        // Drop-site guard bookkeeping (one warning per family+reason).
        auto drop_warn = [&](const std::string& reason,
                             const std::vector<int>& sc) {
            static std::set<std::string> warned;
            std::string sig;
            for (int s : sc) sig += std::to_string(s);
            if (warned.insert(fc.family + "|" + reason).second) {
                std::cerr << "[eikonal_regions guard] family '" << fc.family
                          << "': l~eta region (scale " << sig
                          << ") DROPPED: " << reason
                          << ".  The dropped content is a stage-2 ending "
                             "class (FOLLOWUP_BUILD_SPEC.md); results that "
                             "depend on it are not Feynman values.\n";
            }
            const char* strict = std::getenv("AMFLOW_HALF_ETA_STRICT");
            if (strict && std::string(strict) == "1") {
                throw std::runtime_error(
                    "find_all_region: AMFLOW_HALF_ETA_STRICT=1 and an "
                    "l~eta region of family '" + fc.family
                    + "' was dropped (" + reason + ")");
            }
        };

        std::vector<int> scale((std::size_t)L, 0);
        std::function<void(long)> rec_scale = [&](long depth) {
            if (depth == L) {
            bool any_sig1 = false;
            for (int s : scale) if (s == 2) any_sig1 = true;
            if (any_sig1) {
                // Prune: every l ~ eta loop must couple to an uncut,
                // eta-carrying LINEAR propagator; otherwise the region is
                // scaleless (massless quadratics + nothing eikonal) --
                // HALF_ETA_DERIVATION.md section 2/6 case analysis.
                for (long k = 0; k < L; ++k) {
                    if (scale[(std::size_t)k] != 2) continue;
                    bool ok_k = false;
                    for (const auto& pi : prop_infos) {
                        if (!pi.carries_eta || !pi.is_linear) continue;
                        if (pi.loops_present[(std::size_t)k]) {
                            ok_k = true;
                            break;
                        }
                    }
                    if (!ok_k) return;
                }
                // Leading-form component analysis of the sigma=1 loops:
                // group them by co-appearance in the leading (max
                // he-degree) monomials, and collect the velocities of the
                // eikonal-with-residual-mass leading forms.  Components
                // with two velocities (cusp class V1L/V2L(y)) or more
                // than one loop (the two-loop static-quark class of
                // Broadhurst and Grozin, Phys. Lett. B267 (1991) 105)
                // have no stage-1 ending: drop the region LOUDLY (guard
                // demotion).
                std::vector<long> uf((std::size_t)L);
                for (long k = 0; k < L; ++k) uf[(std::size_t)k] = k;
                std::function<long(long)> find = [&](long a) -> long {
                    while (uf[(std::size_t)a] != a) {
                        uf[(std::size_t)a] = uf[(std::size_t)uf[(std::size_t)a]];
                        a = uf[(std::size_t)a];
                    }
                    return a;
                };
                auto unite = [&](long a, long b) {
                    a = find(a); b = find(b);
                    if (a != b) uf[(std::size_t)a] = b;
                };
                std::vector<std::set<long>> loop_vels((std::size_t)L);
                for (const auto& pi : prop_infos) {
                    long maxhe = -1;
                    std::vector<long> hes;
                    hes.reserve(pi.monos.size());
                    for (const auto& mi : pi.monos) {
                        long he = 2 * mi.eta;
                        for (long v = 0; v < L; ++v) {
                            he += (long)scale[(std::size_t)v]
                                  * mi.ldeg[(std::size_t)v];
                        }
                        hes.push_back(he);
                        if (he > maxhe) maxhe = he;
                    }
                    std::set<long> s2loops;
                    bool lead_const = false;
                    std::set<long> lead_legs;
                    for (std::size_t m = 0; m < pi.monos.size(); ++m) {
                        if (hes[m] != maxhe) continue;
                        const auto& mi = pi.monos[m];
                        for (long v = 0; v < L; ++v) {
                            if (mi.ldeg[(std::size_t)v] > 0
                                    && scale[(std::size_t)v] == 2) {
                                s2loops.insert(v);
                            }
                        }
                        if (mi.tot == 0) lead_const = true;
                        if (mi.tot == 1) {
                            for (long lg : mi.legs) lead_legs.insert(lg);
                        }
                    }
                    long first = -1;
                    for (long v : s2loops) {
                        if (first < 0) first = v;
                        else unite(first, v);
                    }
                    if (lead_const && !lead_legs.empty()) {
                        for (long v : s2loops) {
                            loop_vels[(std::size_t)v].insert(
                                lead_legs.begin(), lead_legs.end());
                        }
                    }
                }
                std::map<long, std::vector<long>> comps;
                for (long k = 0; k < L; ++k) {
                    if (scale[(std::size_t)k] == 2) {
                        comps[find(k)].push_back(k);
                    }
                }
                for (const auto& [root, members] : comps) {
                    (void)root;
                    // Ending engine: 1-loop components with up to
                    // two velocities (V1L(y)) and 2-loop components (V2L
                    // tiers B/C) now HAVE endings — keep them.  Unsupported
                    // master shapes inside a kept region still throw
                    // loudly in eikonal_ending_value (no silent physics).
                    if (members.size() > 2) {
                        drop_warn("eikonal component with more than two "
                                  "loops (stage-3 endings not wired)",
                                  scale);
                        return;
                    }
                    std::set<long> vels;
                    for (long k : members) {
                        vels.insert(loop_vels[(std::size_t)k].begin(),
                                    loop_vels[(std::size_t)k].end());
                    }
                    if (vels.size() > 2) {
                        drop_warn("eikonal component with more than two "
                                  "velocities (endings not wired)",
                                  scale);
                        return;
                    }
                }
            }
            std::vector<int> bs = branch_scale(fc, transformed, scale);

            bool ok = true;
            for (std::size_t cp : cutposi) {
                if (cp >= bs.size()) continue;
                if (bs[cp] != 0) { ok = false; break; }
            }
            if (!ok) return;

            RegionEntry e;
            e.tr     = LoopTransform();
            e.tr.ok  = true;
            for (auto& mfm : tr.map) e.tr.map.push_back(mfm.clone());
            e.scale  = scale;
            e.bscale = std::move(bs);
            all_entries.push_back(std::move(e));
                return;
            }
            if (eik_on) {
                for (int bit : {0, 1, 2}) {
                    scale[(std::size_t)depth] = bit;
                    rec_scale(depth + 1);
                }
            } else {
                for (int bit : {0, 1}) {
                    scale[(std::size_t)depth] = bit;
                    rec_scale(depth + 1);
                }
            }
        };
        rec_scale(0);
    }

    std::vector<std::pair<std::vector<int>, std::vector<RegionEntry>>> by_bscale;
    for (auto& e : all_entries) {
        auto it = std::find_if(
            by_bscale.begin(), by_bscale.end(),
            [&](const auto& item) { return item.first == e.bscale; });
        if (it == by_bscale.end()) {
            by_bscale.emplace_back();
            by_bscale.back().first = e.bscale;
            it = std::prev(by_bscale.end());
        }
        it->second.push_back(std::move(e));
    }

    std::vector<int> pres0;
    pres0.reserve(dedup.size());
    for (auto& d : dedup) {
        auto p = fc.prescription_of_prop(d);
        if (!p.has_value()) {
            throw std::runtime_error(
                "find_all_region: deduped branch has conflicting prescription");
        }
        pres0.push_back(p.value());
    }

    auto check_pres_match = [&](const RegionEntry& e) -> bool {
        std::vector<int> identity_scale((std::size_t)L, 0);
        std::vector<Mfrac> identity_rule = region_rule(rctx, e.tr,
                                                        identity_scale);
        for (std::size_t i = 0; i < dedup.size(); ++i) {
            Mfrac t = apply_region_rule(fc, rctx, dedup[i], identity_rule);
            std::vector<int> pres_list;
            long n_loop = (long)fc.n_loops();
            for (long k = 0; k < n_loop; ++k) {
                if (depends_on_var(t.numerator(), k)
                        || depends_on_var(t.denominator(), k)) {
                    pres_list.push_back(fc.prescription_of_loop(k));
                }
            }
            std::optional<int> p_new;
            if (pres_list.empty()) {
                p_new = 0;
            } else {
                bool all_zero = true;
                for (int pp : pres_list) if (pp != 0) { all_zero = false; break; }
                if (all_zero) p_new = 0;
                else {
                    bool any_pos = false, any_neg = false;
                    for (int pp : pres_list) {
                        if (pp > 0) any_pos = true;
                        else if (pp < 0) any_neg = true;
                    }
                    if (any_pos && !any_neg) p_new = 1;
                    else if (!any_pos && any_neg) p_new = -1;
                }
            }
            if (!p_new.has_value()) return false;
            if (p_new.value() != pres0[i]) return false;
        }
        return true;
    };

    std::vector<Region> result;
    for (auto& [bs, entries] : by_bscale) {
        std::vector<RegionEntry> kept;
        for (auto& e : entries) {
            if (check_pres_match(e)) {
                kept.push_back(std::move(e));
            }
        }
        if (kept.empty()) {
            throw std::runtime_error(
                "find_all_region: some regions are prohibited by prescriptions");
        }
        Region r;
        r.scale = kept.front().scale;
        r.transform.ok = true;
        for (auto& mfm : kept.front().tr.map) r.transform.map.push_back(mfm.clone());
        result.push_back(std::move(r));
    }
    return result;
}

namespace {

Mpoly coeff_half_eta_mpoly(const Mfrac& f, long half_eta_var, long net_he) {
    auto ctx = f.ctx();

    long den_he = 0;
    {
        const Mpoly& den = f.denominator();
        long dlen = fmpz_mpoly_length(den.raw(), ctx->raw());
        std::vector<unsigned long> dexp((std::size_t)ctx->n_vars());
        bool first = true;
        for (long t = 0; t < dlen; ++t) {
            fmpz_mpoly_get_term_exp_ui(dexp.data(), den.raw(), t, ctx->raw());
            long d = (long)dexp[(std::size_t)half_eta_var];
            if (first) {
                den_he = d;
                first = false;
            } else if (d != den_he) {
                throw std::runtime_error(
                    "coeff_half_eta_mpoly: denominator is not a monomial in half_eta");
            }
        }
    }

    long target_num_he = net_he + den_he;
    Mpoly out = Mpoly::zero(ctx);
    const Mpoly& num = f.numerator();
    long nlen = fmpz_mpoly_length(num.raw(), ctx->raw());
    std::vector<unsigned long> nexp((std::size_t)ctx->n_vars());
    fmpz_t coeff;
    fmpz_init(coeff);
    for (long t = 0; t < nlen; ++t) {
        fmpz_mpoly_get_term_exp_ui(nexp.data(), num.raw(), t, ctx->raw());
        if ((long)nexp[(std::size_t)half_eta_var] != target_num_he) continue;
        fmpz_mpoly_get_term_coeff_fmpz(coeff, num.raw(), t, ctx->raw());
        std::vector<unsigned long> e2 = nexp;
        e2[(std::size_t)half_eta_var] = 0;
        fmpz_mpoly_set_coeff_fmpz_ui(out.raw(), coeff, e2.data(), ctx->raw());
    }
    fmpz_clear(coeff);
    return out;
}

Mpoly project_to_fc(const FamilyConfig& fc,
                     const RegionContext& rctx,
                     const Mpoly& p) {
    if (depends_on_var(p, rctx.eta_var) ||
        depends_on_var(p, rctx.half_eta_var)) {
        throw std::runtime_error(
            "project_to_fc: polynomial unexpectedly contains eta/half_eta");
    }
    Mpoly p1 = p.clone();
    fmpz_t zero;
    fmpz_init(zero);
    Mpoly tmp(rctx.ctx);
    fmpz_mpoly_evaluate_one_fmpz(tmp.raw(), p1.raw(),
                                  rctx.half_eta_var, zero, rctx.ctx->raw());
    p1 = std::move(tmp);
    Mpoly tmp2(rctx.ctx);
    fmpz_mpoly_evaluate_one_fmpz(tmp2.raw(), p1.raw(),
                                  rctx.eta_var, zero, rctx.ctx->raw());
    p1 = std::move(tmp2);
    fmpz_clear(zero);

    Mpoly out(fc.ctx);
    long len = fmpz_mpoly_length(p1.raw(), rctx.ctx->raw());
    std::vector<unsigned long> rexp((std::size_t)rctx.ctx->n_vars());
    std::vector<unsigned long> fexp((std::size_t)fc.ctx->n_vars(), 0);
    fmpz_t coeff;
    fmpz_init(coeff);
    for (long t = 0; t < len; ++t) {
        fmpz_mpoly_get_term_exp_ui(rexp.data(), p1.raw(),
                                    t, rctx.ctx->raw());
        fmpz_mpoly_get_term_coeff_fmpz(coeff, p1.raw(),
                                        t, rctx.ctx->raw());
        for (long i = 0; i < fc.ctx->n_vars(); ++i) {
            fexp[(std::size_t)i] = rexp[(std::size_t)i];
        }
        fmpz_mpoly_set_coeff_fmpz_ui(out.raw(), coeff,
                                      fexp.data(), fc.ctx->raw());
    }
    fmpz_clear(coeff);
    return out;
}

bool vanishes_at_loop_and_leg_zero(const FamilyConfig& fc, const Mpoly& p) {
    Mpoly q = p.clone();
    fmpz_t zero;
    fmpz_init(zero);
    long n_zero_vars = (long)fc.n_loops() + (long)fc.n_red_legs();
    for (long v = 0; v < n_zero_vars; ++v) {
        Mpoly tmp(fc.ctx);
        fmpz_mpoly_evaluate_one_fmpz(tmp.raw(), q.raw(), v, zero, fc.ctx->raw());
        q = std::move(tmp);
    }
    fmpz_clear(zero);
    return q.is_zero();
}

}  // namespace

bool zero_region_q(const FamilyConfig& fc,
                   const RegionContext& rctx,
                   const Region& region,
                   const std::vector<std::size_t>& topposi) {
    auto rule = region_rule(rctx, region.transform, region.scale);

    std::vector<Mpoly> expde;
    expde.reserve(topposi.size());
    for (std::size_t k = 0; k < topposi.size(); ++k) {
        std::size_t i = topposi[k];
        const Mpoly& prop = fc.propagators_after_conservation[i];
        Mfrac t = apply_region_rule(fc, rctx, prop, rule);
        t = merge_eta_into_half_eta(fc, rctx, t);
        // Linear (worldline) propagator fix: strip the ACTUAL max half_eta
        // degree (1 for linear/worldline props, 2 for quadratic), not a
        // hardcoded he^2.  Otherwise the leading coeff of a linear prop is
        // zero and the region is spuriously declared scaleless.
        long max_he = 0;
        {
            const Mpoly& tnum = t.numerator();
            long tlen = fmpz_mpoly_length(tnum.raw(), rctx.ctx->raw());
            std::vector<unsigned long> texp((std::size_t)rctx.ctx->n_vars());
            for (long tt = 0; tt < tlen; ++tt) {
                fmpz_mpoly_get_term_exp_ui(texp.data(), tnum.raw(), tt,
                                            rctx.ctx->raw());
                long e = (long)texp[(std::size_t)rctx.half_eta_var];
                if (e > max_he) max_he = e;
            }
        }
        Mfrac stripped = t.clone();
        if (max_he > 0) {
            Mpoly he = Mpoly::variable(rctx.ctx, rctx.half_eta_var);
            Mpoly he_pow = he.clone();
            for (long kk = 1; kk < max_he; ++kk) he_pow = he_pow * he;
            stripped = Mfrac(stripped.numerator(), stripped.denominator() * he_pow);
        }
        Mpoly lead = coeff_half_eta_mpoly(stripped, rctx.half_eta_var, 0);
        Mpoly proj = project_to_fc(fc, rctx, lead);
        expde.push_back(std::move(proj));
    }

    std::stable_sort(expde.begin(), expde.end(),
                     [&](const Mpoly& a, const Mpoly& b) {
                         bool ak = vanishes_at_loop_and_leg_zero(fc, a);
                         bool bk = vanishes_at_loop_and_leg_zero(fc, b);
                         return ak < bk;
                     });

    long M = (long)fc.sp_list.size();
    long rows = (long)expde.size();
    if (rows == 0) {
        return false;
    }
    std::vector<fmpq*> mat_flat((std::size_t)(rows * M), nullptr);
    for (long i = 0; i < rows * M; ++i) {
        mat_flat[(std::size_t)i] = (fmpq*)flint_malloc(sizeof(fmpq));
        fmpq_init(mat_flat[(std::size_t)i]);
    }
    {
        std::vector<fmpq*> row((std::size_t)M, nullptr);
        for (long i = 0; i < M; ++i) {
            row[(std::size_t)i] = (fmpq*)flint_malloc(sizeof(fmpq));
            fmpq_init(row[(std::size_t)i]);
        }
        for (long i = 0; i < rows; ++i) {
            coeff_over_splist(expde[(std::size_t)i], fc.sp_list, row);
            for (long k = 0; k < M; ++k) {
                fmpq_set(mat_flat[(std::size_t)(i * M + k)],
                          row[(std::size_t)k]);
            }
        }
        for (auto* p : row) { fmpq_clear(p); flint_free(p); }
    }

    auto group = maximal_group_rows((std::size_t)rows, (std::size_t)M, mat_flat);
    for (auto* p : mat_flat) { fmpq_clear(p); flint_free(p); }
    std::vector<Mpoly> pruned;
    pruned.reserve(group.size());
    for (std::size_t k = 0; k < group.size(); ++k) {
        pruned.push_back(expde[group[k]].clone());
    }

    if (pruned.empty()) {
        return true;
    }

    return zero_sector_q(fc, pruned);
}

std::vector<Mfrac>
region_power(const FamilyConfig& fc,
             const RegionContext& rctx,
             const PowersContext& pctx,
             const Region& region,
             const std::vector<JIntegral>& integrals) {
    std::size_t N = fc.propagators_after_conservation.size();

    auto rule = region_rule(rctx, region.transform, region.scale);
    // Linear (worldline) propagator fix: track the ACTUAL max half_eta
    // degree per propagator (1 for linear/worldline, 2 for quadratic) so
    // the eta-power below counts -nu_i * max_he_i / 2 instead of the
    // hardcoded -nu_i (which is correct only for quadratic).
    std::vector<long> factor_he(N, 0);
    for (std::size_t i = 0; i < N; ++i) {
        Mfrac t = apply_region_rule(fc, rctx,
            fc.propagators_after_conservation[i], rule);
        t = merge_eta_into_half_eta(fc, rctx, t);
        Mpoly num = t.numerator();
        auto ctx = num.ctx();
        long len = fmpz_mpoly_length(num.raw(), ctx->raw());
        if (len == 0) { factor_he[i] = 0; continue; }
        std::vector<unsigned long> exp((std::size_t)ctx->n_vars());
        long max_he = 0;
        for (long t2 = 0; t2 < len; ++t2) {
            fmpz_mpoly_get_term_exp_ui(exp.data(), num.raw(),
                                        t2, ctx->raw());
            long e = (long)exp[(std::size_t)rctx.half_eta_var];
            if (e > max_he) max_he = e;
        }
        factor_he[i] = max_he;
    }

    long sum_scale = 0;
    for (int s : region.scale) sum_scale += s;

    long eps = pctx.eps_var;
    Mpoly two_minus_eps = Mpoly::constant(pctx.ctx, 2);
    {
        Mpoly e = Mpoly::variable(pctx.ctx, eps);
        two_minus_eps -= e;
    }

    std::vector<Mfrac> out;
    out.reserve(integrals.size());
    for (const auto& integ : integrals) {
        if ((long)integ.indices().size() != (long)N) {
            throw std::invalid_argument(
                "region_power: integral indices length != n_propagators");
        }
        // sum over half-eta units: each prop contributes nu_i * max_he_i
        // half-powers; final eta-power = (2-eps)*sum_scale - sum_he_units/2.
        long sum_he_units = 0;
        for (std::size_t i = 0; i < N; ++i) {
            sum_he_units += integ.indices()[i] * factor_he[i];
        }
        Mpoly val_num(pctx.ctx);
        {
            fmpz_t s;
            fmpz_init_set_si(s, 2 * sum_scale);
            fmpz_mpoly_scalar_mul_fmpz(val_num.raw(),
                                        two_minus_eps.raw(),
                                        s, pctx.ctx->raw());
            fmpz_clear(s);
        }
        Mpoly minus_sum = Mpoly::constant(pctx.ctx, -sum_he_units);
        val_num += minus_sum;
        out.push_back(Mfrac(std::move(val_num),
                              Mpoly::constant(pctx.ctx, 2)));
    }
    return out;
}

}  // namespace amflow::qft
