// SPDX-License-Identifier: MIT
// Layer 16: AMFSystem pipeline -- implementation.
//
// See amflow/pipeline/amfsystem.hpp for the conceptual model.

#include "amflow/pipeline/amfsystem.hpp"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <unistd.h>
#include <exception>
#include <filesystem>
#include <functional>
#include <fstream>
#include <iostream>
#include <map>
#include <mutex>
#include <set>
#include <sstream>
#include <stdexcept>
#include <thread>
#include <utility>
#include <vector>

#include <flint/acb.h>
#include <flint/acb_calc.h>
#include <flint/acb_hypgeom.h>
#include <flint/arb.h>
#include <flint/fmpq.h>
#include <flint/fmpz.h>
#include <flint/fmpz_mpoly.h>
#include <flint/fmpz_mpoly_q.h>

#include "amflow/qft/amfmode.hpp"
#include "amflow/ibp/load_adaptive.hpp"
#include "amflow/ibp/reduce.hpp"
#include "amflow/qft/boundary.hpp"
#include "amflow/pipeline/factorize.hpp"
#include "amflow/qft/family_config.hpp"
#include "amflow/qft/family_uf.hpp"
#include "amflow/ode/inf.hpp"
#include "amflow/numeric/log.hpp"
#include "amflow/numeric/options.hpp"
#include "amflow/numeric/rational.hpp"
#include "amflow/qft/region.hpp"
#include "amflow/ode/amflow.hpp"
#include "amflow/qft/topology.hpp"
#include "amflow/qft/vacuum.hpp"
#include "amflow/ibp/ibp_cache.hpp"
#include "amflow/pipeline/system_checkpoint.hpp"

namespace amflow::pipeline {

// stage-2 endings (defined near the end of this file): subsector
// corners + one-dot variants of the given targets as a Kira
// preferred-masters list, ONLY for eikonal-vacuum families with the region
// engine on (else empty).
std::vector<qft::JIntegral>
build_eikonal_corner_preferred(const qft::FamilyConfig& fc,
                               const std::vector<qft::JIntegral>& targets);

namespace fs = std::filesystem;

namespace {
// Kira keeps per-directory auxiliary state (numeric substitution count).
// Boundary reduction reuses similar path stems; append a unique tag so
// each BlackBoxReduce run gets a fresh working directory.
std::atomic<long long> g_boundary_reduce_unique{0};

// Resolve the effective ε-parallel worker count from the option struct
// and the `AMFLOW_EPS_PARALLEL` env override.  Clamped to >= 1.
//
// Upgrade F2: `AMFLOW_EPS_PARALLEL=auto` picks the worker count from
// current spare capacity, n = clamp(floor((nproc − loadavg1)/FARM_N),
// 1, ∞) — the per-call-site `min(n, n_eps)` already caps to the ε-list
// length.  See amflow/ibp/load_adaptive.hpp.
long resolve_eps_parallel(const AMFSystemOptions& opts) {
    long n = opts.eps_parallel;
    if (const char* e = std::getenv("AMFLOW_EPS_PARALLEL")) {
        if (std::string(e) == "auto") {
            n = ibp::load_adaptive_threads(/*cap=*/0);
            std::cerr << "[amflow][load-adaptive] eps_parallel=auto -> " << n
                      << " (nproc=" << ibp::load_adaptive_nproc()
                      << " load1=" << ibp::load_adaptive_loadavg1()
                      << " farm_n=" << ibp::load_adaptive_farm_n() << ")\n";
            return n > 1 ? n : 1;
        }
        try { n = std::stol(e); } catch (...) {}
    }
    return n > 1 ? n : 1;
}

// Resolve the effective node-tree (Lever 2) worker count from the option
// struct and the `AMFLOW_NODE_PARALLEL` env override.  Clamped to >= 1.
// Same contract as resolve_eps_parallel, including `auto` via the
// load-adaptive scheduler; the call site caps to the sibling count.
long resolve_node_parallel(const AMFSystemOptions& opts) {
    long n = opts.node_parallel;
    if (const char* e = std::getenv("AMFLOW_NODE_PARALLEL")) {
        if (std::string(e) == "auto") {
            n = ibp::load_adaptive_threads(/*cap=*/0);
            std::cerr << "[amflow][load-adaptive] node_parallel=auto -> " << n
                      << " (nproc=" << ibp::load_adaptive_nproc()
                      << " load1=" << ibp::load_adaptive_loadavg1()
                      << " farm_n=" << ibp::load_adaptive_farm_n() << ")\n";
            return n > 1 ? n : 1;
        }
        try { n = std::stol(e); } catch (...) {}
    }
    return n > 1 ? n : 1;
}
}  // namespace

// ===========================================================================
//  EndingScheme name
// ===========================================================================

const char* ending_scheme_name(EndingScheme s) noexcept {
    switch (s) {
        case EndingScheme::Tradition:  return "Tradition";
        case EndingScheme::Cutkosky:   return "Cutkosky";
        case EndingScheme::SingleMass: return "SingleMass";
        case EndingScheme::Trivial:    return "Trivial";
    }
    return "?";
}

// ===========================================================================
//  Helpers
// ===========================================================================

namespace {

std::string mpoly_str(const algebra::Mpoly& p) { return p.to_string(); }

std::string mfrac_str(const algebra::Mfrac& f) {
    algebra::Mpoly num = f.numerator();
    algebra::Mpoly den = f.denominator();
    if (den.is_one()) return num.to_string();
    return "(" + num.to_string() + ")/(" + den.to_string() + ")";
}

std::vector<std::pair<std::string, std::string>>
extract_conservation_strings(const qft::FamilyConfig& fc) {
    std::vector<std::pair<std::string, std::string>> out;
    for (const auto& [key, mpoly] : fc.conservation) {
        out.emplace_back(key, mpoly.to_string());
    }
    return out;
}

std::vector<std::pair<std::string, std::string>>
extract_replacement_strings(const qft::FamilyConfig& fc) {
    std::vector<std::pair<std::string, std::string>> out;
    for (const auto& [key, val] : fc.reduced_replacement) {
        out.emplace_back(key, mfrac_str(val));
    }
    return out;
}

numeric::AcbValue midpoint_only(const numeric::AcbValue& src) {
    // [bug U1] When AMFLOW_XORDER_TAIL_BOUND=1, the eta-propagated master
    // values carry a *meaningful* Arb radius (the x_order truncation tail),
    // so do NOT strip it.  Stripping here once reported a value
    // as "+/- 4e-341" that was in fact good to only ~32 digits.
    static const bool keep_radius = []() -> bool {
        const char* v = std::getenv("AMFLOW_XORDER_TAIL_BOUND");
        return v != nullptr && *v != '\0' && *v != '0';
    }();
    if (keep_radius) return src.clone();
    numeric::AcbValue out;
    arb_set_arf(acb_realref(out.raw()), arb_midref(acb_realref(src.raw())));
    arb_set_arf(acb_imagref(out.raw()), arb_midref(acb_imagref(src.raw())));
    return out;
}

std::vector<std::string>
extract_propagator_strings(const qft::FamilyConfig& fc) {
    std::vector<std::string> out;
    for (const auto& p : fc.propagators_after_conservation) {
        out.push_back(p.to_string());
    }
    return out;
}

std::string jintegral_key(const qft::JIntegral& j) {
    std::ostringstream oss;
    oss << j.family();
    for (long v : j.indices()) oss << "|" << v;
    return oss.str();
}

// Build the eta-injected family.  Mathematica:
//   AMFlowInfo["Propagator"] = ReducedPropagator + $Eta * etac
// We add an `eta` invariant to fc.ctx (assuming it isn't already there)
// and rewrite each propagator string.  If etac is all-zero, returns a
// fresh copy of fc unchanged.
std::unique_ptr<qft::FamilyConfig>
build_injected_fc(const qft::FamilyConfig& parent, const std::vector<int>& etac) {
    if (etac.size() != parent.propagators_after_conservation.size()) {
        throw std::invalid_argument(
            "build_injected_fc: etac size mismatch");
    }
    bool any = false;
    for (int v : etac) if (v != 0) { any = true; break; }
    if (!any) {
        return std::make_unique<qft::FamilyConfig>(qft::FamilyConfig::build(
            parent.family, parent.loops, parent.legs,
            extract_conservation_strings(parent),
            extract_replacement_strings(parent),
            extract_propagator_strings(parent),
            parent.cut, parent.prescription));
    }
    std::vector<std::string> prop_strs;
    prop_strs.reserve(parent.propagators_after_conservation.size());
    for (std::size_t i = 0; i < parent.propagators_after_conservation.size();
            ++i) {
        std::ostringstream oss;
        oss << "(" << mpoly_str(parent.propagators_after_conservation[i]) << ")";
        if (etac[i] != 0) {
            if (etac[i] > 0) oss << " + ";
            else             oss << " - ";
            long ae = std::abs(etac[i]);
            if (ae != 1) oss << ae << "*";
            oss << "eta";
        }
        prop_strs.push_back(oss.str());
    }
    return std::make_unique<qft::FamilyConfig>(qft::FamilyConfig::build(
        parent.family,
        parent.loops, parent.legs,
        extract_conservation_strings(parent),
        extract_replacement_strings(parent),
        prop_strs,
        parent.cut, parent.prescription));
}

// Project an algebra::Mfrac on a "wide" ctx down to dst_ctx by walking
// monomials.  Variables in src that aren't in dst:
//   * "__amf_*" / "__feyn_*" / "__zsq_*" -> drop (set to 0)
//   * everything else -> error.
algebra::Mfrac project_mfrac_by_name(const algebra::Mfrac& src,
                              const std::shared_ptr<algebra::MpolyContext>& dst_ctx) {
    auto src_ctx = src.ctx();
    if (src_ctx.get() == dst_ctx.get()) return src.clone();

    std::map<std::string, long> name_to_dst;
    for (long i = 0; i < dst_ctx->n_vars(); ++i) {
        name_to_dst[dst_ctx->var_name(i)] = i;
    }
    std::vector<long> src_to_dst((std::size_t)src_ctx->n_vars(), -2);
    for (long i = 0; i < src_ctx->n_vars(); ++i) {
        const std::string& name = src_ctx->var_name(i);
        auto it = name_to_dst.find(name);
        if (it != name_to_dst.end()) {
            src_to_dst[(std::size_t)i] = it->second;
        } else if (name.rfind("__amf_", 0) == 0
                || name.rfind("__feyn_", 0) == 0
                || name.rfind("__zsq_", 0) == 0) {
            src_to_dst[(std::size_t)i] = -1;   // drop
        }
    }
    auto walk = [&](const algebra::Mpoly& p) -> algebra::Mpoly {
        algebra::Mpoly out(dst_ctx);
        long len = fmpz_mpoly_length(p.raw(), src_ctx->raw());
        std::vector<unsigned long> exp((std::size_t)src_ctx->n_vars());
        std::vector<unsigned long> dexp((std::size_t)dst_ctx->n_vars(), 0);
        fmpz_t coeff;
        fmpz_init(coeff);
        for (long t = 0; t < len; ++t) {
            fmpz_mpoly_get_term_exp_ui(exp.data(), p.raw(), t, src_ctx->raw());
            fmpz_mpoly_get_term_coeff_fmpz(coeff, p.raw(), t, src_ctx->raw());
            std::fill(dexp.begin(), dexp.end(), 0);
            bool drop = false;
            for (long i = 0; i < src_ctx->n_vars(); ++i) {
                if (exp[(std::size_t)i] == 0) continue;
                long pi = src_to_dst[(std::size_t)i];
                if (pi == -1) { drop = true; break; }
                if (pi == -2) {
                    fmpz_clear(coeff);
                    throw std::runtime_error(
                        "project_mfrac_by_name: residual var '"
                        + src_ctx->var_name(i) + "' not in dst");
                }
                dexp[(std::size_t)pi] = exp[(std::size_t)i];
            }
            if (drop) continue;
            fmpz_mpoly_set_coeff_fmpz_ui(out.raw(), coeff,
                                          dexp.data(), dst_ctx->raw());
        }
        fmpz_clear(coeff);
        return out;
    };
    algebra::Mpoly num = walk(src.numerator());
    algebra::Mpoly den = walk(src.denominator());
    if (den.is_zero()) {
        throw std::runtime_error("project_mfrac_by_name: denominator vanished");
    }
    return algebra::Mfrac(std::move(num), std::move(den));
}

// FmpqHolder: small RAII wrapper.
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
    fmpq* raw()             { return q_; }
    const fmpq* raw() const { return q_; }
private:
    fmpq_t q_;
};

// Test if a variable v has positive exponent in any term of `p`.
bool mpoly_has_var(const algebra::Mpoly& p, long v) {
    auto ctx = p.ctx();
    long len = fmpz_mpoly_length(p.raw(), ctx->raw());
    std::vector<unsigned long> exp((std::size_t)ctx->n_vars());
    for (long t = 0; t < len; ++t) {
        fmpz_mpoly_get_term_exp_ui(exp.data(), p.raw(), t, ctx->raw());
        if (exp[(std::size_t)v] > 0) return true;
    }
    return false;
}

// Substitute every variable in `src` *except* those flagged `keep_names`
// to the rational value supplied via numeric_q (by name).  For "d",
// derive d = 4 - 2*eps when not given directly.  Internal "__amf_*"
// vars not in numeric_q default to 0.  Variables that don't actually
// appear in `src` (positive exponent in num or den) are skipped.
algebra::Mfrac mfrac_substitute(const algebra::Mfrac& src,
                          const std::map<std::string, FmpqHolder>& numeric_q,
                          const std::set<std::string>& keep_names) {
    auto ctx = src.ctx();
    algebra::Mfrac cur = src.clone();
    for (long v = 0; v < ctx->n_vars(); ++v) {
        const std::string& name = ctx->var_name(v);
        if (keep_names.count(name)) continue;
        // Skip variables that don't appear at all (so we don't demand
        // a numeric value for them).  Note: `cur` has already been
        // partially substituted for earlier `v'`s; recompute.
        if (!mpoly_has_var(cur.numerator(), v)
                && !mpoly_has_var(cur.denominator(), v)) {
            continue;
        }
        auto it = numeric_q.find(name);
        if (it == numeric_q.end()) {
            if (name == "d") {
                auto eit = numeric_q.find("eps");
                if (eit == numeric_q.end()) eit = numeric_q.find("__amf_eps");
                if (eit == numeric_q.end()) {
                    throw std::runtime_error(
                        "mfrac_substitute: missing 'd' and 'eps'");
                }
                fmpq_t d_val, two_eps;
                fmpq_init(d_val); fmpq_init(two_eps);
                fmpq_set_si(d_val, 4, 1);
                fmpq_set(two_eps, eit->second.raw());
                fmpq_mul_si(two_eps, two_eps, 2);
                fmpq_sub(d_val, d_val, two_eps);
                fmpq_clear(two_eps);
                cur = cur.substitute(v, d_val);
                fmpq_clear(d_val);
                continue;
            }
            if (name.rfind("__amf_", 0) == 0
                    || name.rfind("__feyn_", 0) == 0
                    || name.rfind("__zsq_", 0) == 0) {
                fmpq_t zero; fmpq_init(zero);
                cur = cur.substitute(v, zero);
                fmpq_clear(zero);
                continue;
            }
            throw std::runtime_error(
                "mfrac_substitute: no value for '" + name + "'");
        }
        cur = cur.substitute(v, it->second.raw());
    }
    return cur;
}

// Convert an algebra::Mfrac in a single variable (after all other substitutions)
// to a numeric::RationalFunction over Q.  `eta_var` < 0 means the result is a
// constant.
numeric::RationalFunction mfrac_to_eta_rational(const algebra::Mfrac& src, long eta_var) {
    auto ctx = src.ctx();
    long n_vars = ctx->n_vars();
    auto walk = [&](const algebra::Mpoly& p) -> numeric::FmpqPoly {
        numeric::FmpqPoly out;
        long len = fmpz_mpoly_length(p.raw(), ctx->raw());
        std::vector<unsigned long> exp((std::size_t)n_vars);
        fmpz_t coeff;
        fmpz_init(coeff);
        for (long t = 0; t < len; ++t) {
            fmpz_mpoly_get_term_exp_ui(exp.data(), p.raw(), t, ctx->raw());
            fmpz_mpoly_get_term_coeff_fmpz(coeff, p.raw(), t, ctx->raw());
            for (long v = 0; v < n_vars; ++v) {
                if (v == eta_var) continue;
                if (exp[(std::size_t)v] != 0) {
                    fmpz_clear(coeff);
                    throw std::runtime_error(
                        "mfrac_to_eta_rational: residual symbolic var '"
                        + ctx->var_name(v) + "'");
                }
            }
            long k = (eta_var >= 0) ? (long)exp[(std::size_t)eta_var] : 0;
            fmpq_t qc, qexist;
            fmpq_init(qc); fmpq_init(qexist);
            fmpq_set_fmpz(qc, coeff);
            out.coeff(k, qexist);
            fmpq_add(qc, qc, qexist);
            out.set_coeff_fmpq(k, qc);
            fmpq_clear(qc); fmpq_clear(qexist);
        }
        fmpz_clear(coeff);
        return out;
    };
    numeric::FmpqPoly num_p = walk(src.numerator());
    numeric::FmpqPoly den_p = walk(src.denominator());
    return numeric::RationalFunction(std::move(num_p), std::move(den_p));
}

// Approximate an acb value as an fmpq (RationalizePre digits).
//
// PRECISION-MISMATCH FIX (see src/ode/inf.cpp:acb_real_to_fmpq_local): when
// `rationalize_pre` > working_prec * log10(2), the rationalization captures
// the binary-representation noise in the input as a "real" rational, which
// propagates as a 1e-25 residual through the AMFSystem build and amplifies
// into 1e7-1e21 errors on masters with all-zero BC.  Cap rationalize_digits
// at what working_prec can resolve.
void acb_to_fmpq_rational(fmpq_t out, const acb_t v, long prec) {
    int digits = numeric::rationalize_pre();
    if (digits <= 0) digits = 20;

    long working_digits = static_cast<long>(prec * 0.301029995663981195L);
    if (working_digits >= 5) working_digits -= 5;
    if (working_digits < 1) working_digits = 1;
    if (digits > working_digits) {
        digits = static_cast<int>(working_digits);
    }

    const long bits = std::max<long>(prec, numeric::decimal_digits_to_bits(digits) + 32);
    arb_t scale, scaled, xb;
    arb_init(scale);
    arb_init(scaled);
    arb_init(xb);
    arb_set_si(scale, 10);
    arb_pow_ui(scale, scale, static_cast<unsigned long>(digits), bits);
    arb_set(xb, acb_realref(v));
    arb_mul(scaled, xb, scale, bits);

    fmpz_t numer, pow10;
    fmpz_init(numer);
    fmpz_init(pow10);
    arf_get_fmpz(numer, arb_midref(scaled), ARF_RND_NEAR);
    fmpz_set_si(pow10, 10);
    fmpz_pow_ui(pow10, pow10, static_cast<unsigned long>(digits));
    fmpq_set_fmpz_frac(out, numer, pow10);
    fmpq_canonicalise(out);

    fmpz_clear(numer);
    fmpz_clear(pow10);
    arb_clear(scale);
    arb_clear(scaled);
    arb_clear(xb);
}

}  // namespace

// ===========================================================================
//  Per-eps rational frame helpers (exposed for unit tests)
// ===========================================================================

namespace detail {

bool numd_exact_eps(fmpq_t out) {
    const std::string numd = ibp::diffeq_numd_env();
    if (numd.empty()) return false;
    fmpq_t d_q;
    fmpq_init(d_q);
    const bool parsed = (fmpq_set_str(d_q, numd.c_str(), 10) == 0);
    if (parsed) {
        fmpq_canonicalise(d_q);
        fmpq_t four;
        fmpq_init(four);
        fmpq_set_si(four, 4, 1);
        fmpq_sub(d_q, four, d_q);        // 4 - d
        fmpq_div_2exp(out, d_q, 1);      // (4 - d)/2 = eps
        fmpq_clear(four);
    }
    fmpq_clear(d_q);
    return parsed;
}

void eps_sample_to_fmpq(fmpq_t out, const numeric::AcbValue& eps, long prec) {
    // Frame consistency with a numeric-d differential equation.  When
    // AMFLOW_DIFFEQ_NUMD=P/Q is set, the DE tables arrive with d = P/Q
    // baked in exactly, so the DE sits at eps = (4 - P/Q)/2 exactly and
    // the d -> 4 - 2*eps substitution of `mfrac_substitute` is a no-op on
    // them.  The
    // region exponents and boundary coefficients, however, are evaluated
    // at whatever rational is stored here.  Rationalizing the sample to
    // `rationalize_pre` decimal digits (a truncation, not an exact
    // recovery: 1/113 does not survive it) would put them in a slightly
    // different frame; the difference b*(eps_q - eps) then enters the
    // Taylor system's diagonal and, once it exceeds the sparse chop
    // threshold, is pivoted on and the solve returns a wrong value or an
    // exact zero.  Use the exact value from the env instead — the startup
    // assert has already verified that the single eps sample equals it.
    if (numd_exact_eps(out)) return;
    acb_to_fmpq_rational(out, eps.raw(), prec);
}

}  // namespace detail

// ===========================================================================
//  Counter for system IDs
// ===========================================================================

namespace {
std::atomic<long> g_amf_system_id_counter{0};
}

// ===========================================================================
//  AMFSystem ctor / dtor / move
// ===========================================================================

AMFSystem::AMFSystem(qft::FamilyConfig fc,
                       std::vector<qft::JIntegral> preferred,
                       std::vector<int>       etac,
                       EndingScheme           ending_used,
                       const AMFSystemOptions& opts)
    : preferred_(std::move(preferred)),
      etac_(std::move(etac)),
      ending_(ending_used),
      opts_(opts),
      system_id_(g_amf_system_id_counter.fetch_add(1)) {
    bool all_zero = std::all_of(etac_.begin(), etac_.end(),
                                  [](int v) { return v == 0; });
    is_ending_ = all_zero;
    fc_ = std::make_unique<qft::FamilyConfig>(std::move(fc));
}

AMFSystem::AMFSystem(AMFSystem&&) noexcept            = default;
AMFSystem& AMFSystem::operator=(AMFSystem&&) noexcept = default;
AMFSystem::~AMFSystem() = default;

// ===========================================================================
//  EndingQ
// ===========================================================================

// Forward decl; implementation below uses topology analysis.
namespace {
bool single_mass_ending_q_impl(const qft::FamilyConfig& fc,
                                  const std::vector<qft::JIntegral>& preferred,
                                  const AMFSystemOptions& opts);
bool cutkosky_ending_q_impl(const qft::FamilyConfig& fc,
                              const std::vector<qft::JIntegral>& preferred);
}

bool ending_q(const qft::FamilyConfig& fc,
              const std::vector<qft::JIntegral>& preferred,
              EndingScheme scheme,
              const AMFSystemOptions& opts) {
    if (preferred.empty()) return true;
    auto top = qft::get_top_position(preferred);
    if (top.empty()) return true;

    switch (scheme) {
        case EndingScheme::Tradition: {
            auto pos = qft::amf_position(fc, top, opts.amf_modes);
            return pos.empty();
        }
        case EndingScheme::SingleMass:
            return single_mass_ending_q_impl(fc, preferred, opts);
        case EndingScheme::Cutkosky:
            return cutkosky_ending_q_impl(fc, preferred);
        case EndingScheme::Trivial:
            // Mirrors upstream AMFlow.m:1031:
            //   AMFSystemEndingQ[preferred_, "Trivial"] := False;
            // Trivial is the always-applicable fallback; it never
            // declares "ending here" because it is the *handler* that
            // takes over when no other scheme can.
            return false;
    }
    return false;
}

namespace {

// Apply numeric substitution (from opts.bb.numeric_values) to an algebra::Mfrac
// living on a context that prefixes fc.ctx (so its variables map by
// position to fc.ctx variables).  Variables not in numeric_values are
// left symbolic.  Returns a new algebra::Mfrac on the same ctx.
algebra::Mfrac apply_numeric(const algebra::Mfrac& src,
                     const std::map<std::string, FmpqHolder>& numeric_q) {
    auto ctx = src.ctx();
    algebra::Mfrac cur = src.clone();
    for (long v = 0; v < ctx->n_vars(); ++v) {
        const std::string& name = ctx->var_name(v);
        auto it = numeric_q.find(name);
        if (it == numeric_q.end()) continue;
        cur = cur.substitute(v, it->second.raw());
    }
    return cur;
}

// Build numeric_q from opts.bb.numeric_values.
std::map<std::string, FmpqHolder>
build_numeric_q(const ibp::ReduceOptions& bb) {
    std::map<std::string, FmpqHolder> out;
    auto parse = [&](const std::string& s, fmpq_t v) -> bool {
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
            // Decimal.
            std::size_t dot = s.find('.');
            if (dot == std::string::npos) return false;
            std::string sint = s.substr(0, dot) + s.substr(dot + 1);
            long denom = 1;
            for (std::size_t i = dot + 1; i < s.size(); ++i) denom *= 10;
            long num = std::stol(sint);
            fmpq_set_si(v, num, denom);
            return true;
        } catch (...) { return false; }
    };
    for (const auto& [name, val_str] : bb.numeric_values) {
        FmpqHolder h;
        if (parse(val_str, h.raw())) out.emplace(name, std::move(h));
    }
    return out;
}

std::map<std::string, FmpqHolder>
build_family_numeric_q(const ibp::ReduceOptions& bb) {
    auto out = build_numeric_q(bb);
    out.erase("eps");
    out.erase("__amf_eps");
    out.erase("d");
    out.erase("eta");
    out.erase("__amf_eta");
    return out;
}

std::map<std::string, FmpqHolder>
build_numeric_q_for_eps(const ibp::ReduceOptions& bb,
                        const numeric::AcbValue& eps,
                        long prec) {
    std::map<std::string, FmpqHolder> out;
    {
        FmpqHolder h;
        detail::eps_sample_to_fmpq(h.raw(), eps, prec);
        out.emplace("eps", h);
        out.emplace("__amf_eps", h);
    }

    auto parse = [&](const std::string& s, fmpq_t v) -> bool {
        std::size_t slash = s.find('/');
        try {
            if (slash != std::string::npos) {
                long p = std::stol(s.substr(0, slash));
                long q = std::stol(s.substr(slash + 1));
                fmpq_set_si(v, p, q);
                return true;
            }
            std::size_t pos = 0;
            long p = std::stol(s, &pos);
            if (pos == s.size()) {
                fmpq_set_si(v, p, 1);
                return true;
            }
            std::stod(s, &pos);
            if (pos != s.size()) return false;
            std::size_t dot = s.find('.');
            if (dot == std::string::npos) return false;
            std::string sint = s.substr(0, dot) + s.substr(dot + 1);
            long denom = 1;
            for (std::size_t i = dot + 1; i < s.size(); ++i) denom *= 10;
            long num = std::stol(sint);
            fmpq_set_si(v, num, denom);
            return true;
        } catch (...) {
            return false;
        }
    };

    for (const auto& [name, val_str] : bb.numeric_values) {
        if (out.count(name)) continue;
        FmpqHolder h;
        if (!parse(val_str, h.raw())) {
            throw std::runtime_error(
                "build_numeric_q_for_eps: cannot parse numeric '"
                + val_str + "' for '" + name + "'");
        }
        out.emplace(name, std::move(h));
    }
    return out;
}

void evaluate_const_mfrac(const algebra::Mfrac& m,
                          const std::map<std::string, FmpqHolder>& numeric_q,
                          long prec,
                          numeric::AcbValue& out) {
    std::set<std::string> empty_keep;
    algebra::Mfrac sub = mfrac_substitute(m, numeric_q, empty_keep);
    algebra::Mpoly num = sub.numerator();
    algebra::Mpoly den = sub.denominator();
    fmpz_t nint, dint;
    fmpz_init(nint);
    fmpz_init(dint);
    std::vector<unsigned long> z((std::size_t)sub.ctx()->n_vars(), 0);
    fmpz_mpoly_get_coeff_fmpz_ui(nint, num.raw(), z.data(), sub.ctx()->raw());
    fmpz_mpoly_get_coeff_fmpz_ui(dint, den.raw(), z.data(), sub.ctx()->raw());
    if (std::getenv("AMF_DBG_X00") && fmpz_is_zero(dint)) {
        // calcx00 audit: an exact
        // zero denominator here is the precision-independent 0/0 -> NaN
        // source that resurfaces downstream as the calcx00 abort.
        char* ns = fmpz_get_str(nullptr, 10, nint);
        fprintf(stderr,
                "[mfrac0] DENOM constant term is ZERO at substitution "
                "(num.const=%s, num.nterms=%ld den.nterms=%ld) -> 0/0/inf\n",
                ns, (long)fmpz_mpoly_length(num.raw(), sub.ctx()->raw()),
                (long)fmpz_mpoly_length(den.raw(), sub.ctx()->raw()));
        flint_free(ns);
    }
    fmpq_t q;
    fmpq_init(q);
    fmpq_set_fmpz_frac(q, nint, dint);
    acb_set_fmpq(out.raw(), q, prec);
    fmpq_clear(q);
    fmpz_clear(nint);
    fmpz_clear(dint);
}

// ===========================================================================
//  EikonalVacuum ending engine (region engine)
// ===========================================================================
//
// Gated by AMFLOW_EIKONAL_REGIONS=1 (same gate as the sigma in {0,1/2,1}
// region enumeration in qft::find_all_region).  Boundary sub-families
// generated by the l ~ eta regions contain eikonal propagators with unit
// residual mass (2 l.u - 1).  These families are SELF-SIMILAR under eta
// re-injection (the l ~ eta region of 2l.u - 1 - eta reproduces the family),
// so they MUST be ended directly, never recursed.  Closed forms (verified
// two independent routes, the half-eta fix, routes V7a/V7b; this is the
// standard one-loop HQET integral, see e.g. A. G. Grozin, Heavy Quark
// Effective Theory (Springer Tracts Mod. Phys. 201, 2004), and for linear
// propagators in AMF see Z.-F. Liu, Y.-Q. Ma, arXiv:2201.11636):
//
//   T(a,b) = Int d^Dk/(i pi^{D/2}) 1/[(-k^2-i0)^a (1-2k.u-i0)^b]
//          = Gamma(D/2-a) Gamma(2a+b-D) / (Gamma(a) Gamma(b)),  D = 4-2eps,
//
// with the family-basis translation (props as constructed by
// boundary_integrals: cq*L^2 and ce*2L.u + c0, |cq|=|c0|=1, |ce|=2):
//
//   J[a on quad, b on eikonal] = cq^-a c0^-b (-1)^a T(a,b).
//
// Disconnected components multiply; single-mass tadpole components
// (L^2 - 1)^a contribute (-1)^a Gamma(a-2+eps)/Gamma(a); constant
// propagators c contribute c^-index.  Anything outside these shapes
// throws LOUDLY (no silent wrong physics).
namespace {

bool eikonal_regions_on() {
    const char* e = std::getenv("AMFLOW_EIKONAL_REGIONS");
    return e && *e && std::string(e) != "0";
}

// True iff some propagator is linear in the loop momenta AND carries a
// nonzero constant monomial (unit residual mass): the eikonal-vacuum
// fingerprint.  Such propagators only occur in l~eta boundary families.
bool has_eikonal_mass_prop(const qft::FamilyConfig& fc) {
    const long L = (long)fc.n_loops();
    const long n_vars = fc.ctx->n_vars();
    std::vector<unsigned long> exp((std::size_t)n_vars);
    for (const auto& p : fc.propagators_after_conservation) {
        long len = fmpz_mpoly_length(p.raw(), fc.ctx->raw());
        long maxdeg = 0;
        bool has_const = false;
        for (long t = 0; t < len; ++t) {
            fmpz_mpoly_get_term_exp_ui(exp.data(), p.raw(), t,
                                        fc.ctx->raw());
            long d = 0;
            bool any_var = false;
            for (long v = 0; v < n_vars; ++v) {
                if (exp[(std::size_t)v] > 0) any_var = true;
                if (v < L) d += (long)exp[(std::size_t)v];
            }
            if (d > maxdeg) maxdeg = d;
            if (!any_var) has_const = true;
        }
        if (maxdeg == 1 && has_const) return true;
    }
    return false;
}

// ===========================================================================
//  Stage-2 endings: two-velocity (cusp)
//  eikonal vacua V1L(y) and the wired V2L(y) tiers B/C.
// ===========================================================================
//
// Verified closed forms (stage-2 endings; every formula checked >= 30
// digits against independent quadrature routes at 3+ (y,eps) points):
//
//  V1L(a,b,c;y) = Int d^Dk/(i pi^{D/2}) /[(-k^2)^a (1-2k.u1)^b (1-2k.u2)^c]
//               = T(a,b+c) * 3F2(2-eps-a, b, c;
//                                (b+c)/2, (b+c+1)/2; (1-y)/2),   u1.u2 = y,
//  generic complex a (the Tier-B bubble shift needs a = int - 2 + eps).
//  Orientations: props (1 - r1 2k.u1)(1 - r2 2k.u2) give y_eff = r1 r2 y.
//
//  Tier B (2-loop component, one loop purely quadratic):
//      = Gb(a2,a3) * V1L(a1+a2+a3-D/2, b, c; y),
//      Gb(x,w) = G(D/2-x)G(D/2-w)G(x+w-D/2) / (G(x)G(w)G(D-x-w)).
//
//  Tier C (2-loop component, inner loop = {crossed quad^a3, one eikonal^b2},
//  no own quadratic):
//      = T(a3,b2) * T(a1,b1+bh)
//        * 3F2(2-eps-a1, b1, bh; (b1+bh)/2, (b1+bh+1)/2; z),
//      bh = 2 a3 + b2 - D,   z = (1 + sigma_eff y)/2,
//      sigma_eff = r_outer * sigma_crossed * r_inner.
//  z > 1 (the radiative sheet, Im != 0) is the z + i0 boundary value,
//  evaluated through the classical 1/z continuation (NIST 16.8.8) — the
//  parameter differences carry eps, so the formula is non-degenerate at the
//  flow's eps samples.  All i0 bookkeeping is anchored by the verified
//  spacelike-velocity identity T_sp(a,B) = T(a,B) e^{-i pi (a-D/2)}.
//  Unsupported shapes (three velocities in a slot, nested cusp slots,
//  Tier D = both loops carrying own quads AND eikonals, >2-loop components,
//  z <= -1 or z == 1 kinematics) THROW loudly.

// T(a,B) = G(D/2-a) G(2a+B-D) rgamma(a) rgamma(B), generic complex a, B.
void eik_T_acb(acb_t res, const acb_t a, const acb_t B,
               const numeric::AcbValue& eps, long prec) {
    acb_t Dh, arg, g, rg;
    acb_init(Dh); acb_init(arg); acb_init(g); acb_init(rg);
    acb_set_si(Dh, 2);
    acb_sub(Dh, Dh, eps.raw(), prec);          // D/2
    acb_sub(arg, Dh, a, prec);
    acb_gamma(res, arg, prec);                 // G(D/2-a)
    acb_mul_2exp_si(arg, Dh, 1);               // D
    acb_neg(arg, arg);
    acb_addmul_si(arg, a, 2, prec);
    acb_add(arg, arg, B, prec);                // 2a+B-D
    acb_gamma(g, arg, prec);
    acb_mul(res, res, g, prec);
    acb_rgamma(rg, a, prec);
    acb_mul(res, res, rg, prec);
    acb_rgamma(rg, B, prec);
    acb_mul(res, res, rg, prec);
    acb_clear(Dh); acb_clear(arg); acb_clear(g); acb_clear(rg);
}

// Massless bubble Gb(x,w) (integer indices).
void eik_bubble(acb_t res, long x, long w,
                const numeric::AcbValue& eps, long prec) {
    acb_t Dh, arg, g, rg;
    acb_init(Dh); acb_init(arg); acb_init(g); acb_init(rg);
    acb_set_si(Dh, 2);
    acb_sub(Dh, Dh, eps.raw(), prec);          // D/2
    acb_sub_si(arg, Dh, x, prec);
    acb_gamma(res, arg, prec);                 // G(D/2-x)
    acb_sub_si(arg, Dh, w, prec);
    acb_gamma(g, arg, prec);                   // G(D/2-w)
    acb_mul(res, res, g, prec);
    acb_neg(arg, Dh);
    acb_add_si(arg, arg, x + w, prec);         // x+w-D/2
    acb_gamma(g, arg, prec);
    acb_mul(res, res, g, prec);
    acb_set_si(arg, x);
    acb_rgamma(rg, arg, prec);
    acb_mul(res, res, rg, prec);
    acb_set_si(arg, w);
    acb_rgamma(rg, arg, prec);
    acb_mul(res, res, rg, prec);
    acb_mul_2exp_si(arg, Dh, 1);               // D
    acb_sub_si(arg, arg, x + w, prec);         // D-x-w
    acb_rgamma(rg, arg, prec);
    acb_mul(res, res, rg, prec);
    acb_clear(Dh); acb_clear(arg); acb_clear(g); acb_clear(rg);
}

// 3F2(al[0..2]; rh[0..1]; z) for real rational z > 1 on the z + i0 side:
// classical three-term 1/z continuation (NIST 16.8.8),
//   3F2(..;z) = sum_j C_j (-z)^{-al_j}
//               * 3F2(al_j, al_j-rh_0+1, al_j-rh_1+1;
//                     al_j-al_k+1, al_j-al_l+1; 1/z),
//   C_j = G(rh_0)G(rh_1)G(al_k-al_j)G(al_l-al_j)
//         / (G(al_k)G(al_l)G(rh_0-al_j)G(rh_1-al_j)),
//   (-z)^{-al_j} = z^{-al_j} e^{+i pi al_j}   (z + i0 side).
// Caller must guarantee pairwise non-integer al differences.
void eik_hyp3f2_cont_above(acb_t res, acb_srcptr al, acb_srcptr rh,
                           const acb_t z, long prec) {
    acb_t sum, C, t, u, zi, F;
    acb_init(sum); acb_init(C); acb_init(t); acb_init(u);
    acb_init(zi); acb_init(F);
    acb_ptr A = _acb_vec_init(3);
    acb_ptr Bv = _acb_vec_init(2);
    acb_one(t);
    acb_div(zi, t, z, prec);                   // 1/z
    acb_zero(sum);
    for (int jj = 0; jj < 3; ++jj) {
        const int kk = (jj + 1) % 3, ll = (jj + 2) % 3;
        // C_j
        acb_gamma(C, rh + 0, prec);
        acb_gamma(t, rh + 1, prec);
        acb_mul(C, C, t, prec);
        acb_sub(u, al + kk, al + jj, prec);
        acb_gamma(t, u, prec);
        acb_mul(C, C, t, prec);
        acb_sub(u, al + ll, al + jj, prec);
        acb_gamma(t, u, prec);
        acb_mul(C, C, t, prec);
        acb_rgamma(t, al + kk, prec);
        acb_mul(C, C, t, prec);
        acb_rgamma(t, al + ll, prec);
        acb_mul(C, C, t, prec);
        acb_sub(u, rh + 0, al + jj, prec);
        acb_rgamma(t, u, prec);
        acb_mul(C, C, t, prec);
        acb_sub(u, rh + 1, al + jj, prec);
        acb_rgamma(t, u, prec);
        acb_mul(C, C, t, prec);
        // (-z)^{-al_j} = z^{-al_j} e^{+i pi al_j}
        acb_neg(u, al + jj);
        acb_pow(t, z, u, prec);
        acb_mul(C, C, t, prec);
        acb_exp_pi_i(t, al + jj, prec);
        acb_mul(C, C, t, prec);
        // inner 3F2 at 1/z
        acb_set(A + 0, al + jj);
        acb_sub(A + 1, al + jj, rh + 0, prec);
        acb_add_si(A + 1, A + 1, 1, prec);
        acb_sub(A + 2, al + jj, rh + 1, prec);
        acb_add_si(A + 2, A + 2, 1, prec);
        acb_sub(Bv + 0, al + jj, al + kk, prec);
        acb_add_si(Bv + 0, Bv + 0, 1, prec);
        acb_sub(Bv + 1, al + jj, al + ll, prec);
        acb_add_si(Bv + 1, Bv + 1, 1, prec);
        acb_hypgeom_pfq(F, A, 3, Bv, 2, zi, 0, prec);
        acb_addmul(sum, C, F, prec);
    }
    acb_set(res, sum);
    _acb_vec_clear(A, 3);
    _acb_vec_clear(Bv, 2);
    acb_clear(sum); acb_clear(C); acb_clear(t); acb_clear(u);
    acb_clear(zi); acb_clear(F);
}

// 2F1(al, be; ga; z + i0) for rational z in (1,2): the 1-z connection
// formula (non-degenerate: ga-al-be carries eps),
//   = C1 2F1(al,be; al+be-ga+1; 1-z)
//   + C2 |1-z|^{ga-al-be} e^{-i pi (ga-al-be)} 2F1(ga-al,ga-be; ga-al-be+1; 1-z)
// Verified vs branch-explicit lambda-contour quadrature
// (independent mixed-case verification, >= 30d, 15 points).
void eik_hyp2f1_above(acb_t res, const acb_t al, const acb_t be,
                      const acb_t ga, const acb_t z, long prec) {
    acb_t u, nu, C1, C2, t, g, F, acc;
    acb_init(u); acb_init(nu); acb_init(C1); acb_init(C2);
    acb_init(t); acb_init(g); acb_init(F); acb_init(acc);
    acb_one(u);
    acb_sub(u, u, z, prec);                    // u = 1-z in (-1,0)
    acb_sub(nu, ga, al, prec);
    acb_sub(nu, nu, be, prec);                 // nu = ga-al-be
    // C1 = G(ga) G(nu) rgamma(ga-al) rgamma(ga-be)
    acb_gamma(C1, ga, prec);
    acb_gamma(t, nu, prec);
    acb_mul(C1, C1, t, prec);
    acb_sub(t, ga, al, prec);
    acb_rgamma(g, t, prec);
    acb_mul(C1, C1, g, prec);
    acb_sub(t, ga, be, prec);
    acb_rgamma(g, t, prec);
    acb_mul(C1, C1, g, prec);
    // C2 = G(ga) G(-nu) rgamma(al) rgamma(be)
    acb_gamma(C2, ga, prec);
    acb_neg(t, nu);
    acb_gamma(t, t, prec);
    acb_mul(C2, C2, t, prec);
    acb_rgamma(g, al, prec);
    acb_mul(C2, C2, g, prec);
    acb_rgamma(g, be, prec);
    acb_mul(C2, C2, g, prec);
    // acc = C1 * 2F1(al, be; al+be-ga+1; u)
    acb_add(t, al, be, prec);
    acb_sub(t, t, ga, prec);
    acb_add_si(t, t, 1, prec);
    acb_hypgeom_2f1(F, al, be, t, u, 0, prec);
    acb_mul(acc, C1, F, prec);
    // += C2 * |u|^nu e^{-i pi nu} * 2F1(ga-al, ga-be; nu+1; u)
    acb_neg(t, u);                             // |u| = -u > 0
    acb_pow(t, t, nu, prec);
    acb_mul(C2, C2, t, prec);
    acb_neg(g, nu);
    acb_exp_pi_i(g, g, prec);                  // e^{-i pi nu}
    acb_mul(C2, C2, g, prec);
    acb_sub(t, ga, al, prec);
    acb_sub(g, ga, be, prec);
    acb_add_si(F, nu, 1, prec);
    {
        acb_t F2;
        acb_init(F2);
        acb_hypgeom_2f1(F2, t, g, F, u, 0, prec);
        acb_addmul(acc, C2, F2, prec);
        acb_clear(F2);
    }
    acb_set(res, acc);
    acb_clear(u); acb_clear(nu); acb_clear(C1); acb_clear(C2);
    acb_clear(t); acb_clear(g); acb_clear(F); acb_clear(acc);
}

// ===========================================================================
//  Ending-class integrals: the two remaining production ending classes.
//  Mathematics: the ending-class derivations (nested cusp + leg-quad
//  vertex); every formula was checked to >= 30 digits against independent
//  quadrature routes, and is validated here against 40/45-digit reference
//  value tables.
//
//  Class 2 — eikonal x leg-quad 1L vertex:
//    C2(a,b,c) = Int d^Dk 1/[(-k^2)^a (-(k-q)^2)^b (1-2k.u)^c],
//    q.u = 0, q^2 = -1 (guarded exact):  single Barnes contour
//    (exponentially convergent, valid at ALL eps incl. the flow's tiny
//    real samples) + crossed-pole residue corrections (sign +1,
//    validated on 33/33 cases at 55 digits).
//
//  Class 1 — nested cusp: 2-loop component whose
//  inner loop is {crossed quad^a1, one eikonal^b} and whose outer loop
//  carries its own quad^a plus TWO eikonal slots.  The E-step
//  continuation E(w) = T(a1,b)(w-i0)^{D-2a1-b} converts the inner loop to
//  an opposite-orientation eikonal power bh = 2a1+b-D on the outer loop:
//    V1L3(a; c,ch,bh; y) = Int d^Dk 1/[(-k^2)^a (1-2k.u1)^c (1-2k.u2)^ch
//                                       ((1+2k.u2)-i0)^bh],
//  evaluated by the reduced 1-fold form (route R):
//    V1L3 = pref * 2^{1-ch-bh} * (J+ + J-),
//  with J+ ELEMENTARY (the 2F1 in I+ terminates: rgamma(1-ch) = 0 exactly
//  for integer ch; the surviving connection piece is a degree-(ch-1)
//  polynomial) and J- a 1-fold over 2F1s.  All quadrature is tanh-sinh in
//  the exact DISTANCE variable per segment (distance-frame rule), and
//  every singular edge d^p uses the exact one-term subtraction
//    Int_0^h d^p A(d) dd = A(0) h^{p+1}/(p+1) + Int_0^h d^p (A(d)-A(0)) dd
//  (Beta-continuation) — mandatory at the flow's eps ~ 1e-14 where the
//  edge exponents sit at -1+O(eps) and carry O(1/eps) mass that no purely
//  numerical rule can integrate.  The J- m->0 endpoint uses the exact 1-z
//  connection split (tau = 1-ch-bh never an integer), giving one m^p piece
//  plus one analytic piece.  Radii: tanh-sinh folds the last level
//  difference (x4) into the result ball; the Barnes contour uses rigorous
//  acb_calc_integrate with the truncation tail added to the radius.
//
//  Scope kept loud (unchanged throws): Tier D, >=3-loop components,
//  >2 velocities, y_eff <= 1 (incl. z <= -1), q^2 != -1 rescale,
//  numerator/ISP recursion, Re(s) <= -1 (y+i0 contour route; not
//  exercised by the flow's real-positive eps grid nor the tables).
// ===========================================================================

// True if the ball may contain a pole of Gamma (a non-positive integer).
bool eik_maybe_gamma_pole(const acb_t x) {
    if (arb_is_positive(acb_realref(x))) return false;
    if (!arb_contains_zero(acb_imagref(x))) return false;
    return arb_contains_int(acb_realref(x)) != 0;
}

struct EikC2Param {
    acb_t Dh;      // D/2
    acb_t sg;      // 2a+2b+c-D
    acb_t rgden;   // rgamma(a) rgamma(b) rgamma(c)
    long a, b;
};

// Barnes integrand for C2 (class 2).  Meromorphic: returns an
// indeterminate ball whenever a Gamma argument may sit on a pole, which
// both keeps the order==1 holomorphy contract honest and forces
// subdivision.
int eik_c2_integrand(acb_ptr out, const acb_t z, void* vp, slong order,
                     slong prec) {
    (void)order;
    const EikC2Param* P = static_cast<const EikC2Param*>(vp);
    acb_t t, g, acc;
    acb_init(t); acb_init(g); acb_init(acc);
    bool bad = false;
    // G(-z)
    acb_neg(t, z);
    if (eik_maybe_gamma_pole(t)) bad = true;
    else acb_gamma(acc, t, prec);
    // * G(D/2-a-b-z)
    if (!bad) {
        acb_sub(t, P->Dh, z, prec);
        acb_sub_si(t, t, P->a + P->b, prec);
        if (eik_maybe_gamma_pole(t)) bad = true;
        else { acb_gamma(g, t, prec); acb_mul(acc, acc, g, prec); }
    }
    // * G(sg+2z)
    if (!bad) {
        acb_mul_2exp_si(t, z, 1);
        acb_add(t, t, P->sg, prec);
        if (eik_maybe_gamma_pole(t)) bad = true;
        else { acb_gamma(g, t, prec); acb_mul(acc, acc, g, prec); }
    }
    // * G(a+z) * G(b+z)
    if (!bad) {
        acb_add_si(t, z, P->a, prec);
        if (eik_maybe_gamma_pole(t)) bad = true;
        else { acb_gamma(g, t, prec); acb_mul(acc, acc, g, prec); }
    }
    if (!bad) {
        acb_add_si(t, z, P->b, prec);
        if (eik_maybe_gamma_pole(t)) bad = true;
        else { acb_gamma(g, t, prec); acb_mul(acc, acc, g, prec); }
    }
    // * rgamma(a+b+2z) * rgden   (entire — no pole guard needed)
    if (!bad) {
        acb_mul_2exp_si(t, z, 1);
        acb_add_si(t, t, P->a + P->b, prec);
        acb_rgamma(g, t, prec);
        acb_mul(acc, acc, g, prec);
        acb_mul(acc, acc, P->rgden, prec);
        acb_set(out, acc);
    } else {
        acb_indeterminate(out);
    }
    acb_clear(t); acb_clear(g); acb_clear(acc);
    return 0;
}

// C2(a,b,c) via the Barnes route (a,b,c >= 1 integers; any eps):
//   C2 = (1/2 pi i) Int_{c0} dz G(-z)G(D/2-a-b-z)G(sg+2z)G(a+z)G(b+z)
//                                / [G(a)G(b)G(c)G(a+b+2z)]
//        + Sum_n [(-1)^n/n!] (residue corrections for family-2 poles
//          z_n = D/2-a-b+n with Re(z_n) < c0; sign +1, validated),
//   c0 = -0.9 min(a, b, Re(sg)/2, 1/2),  sg = 2a+2b+c-D.
void eik_class2_barnes(acb_t res, long a, long b, long c,
                       const numeric::AcbValue& eps, long prec) {
    const long wp = prec + 60;
    EikC2Param P;
    acb_init(P.Dh); acb_init(P.sg); acb_init(P.rgden);
    P.a = a; P.b = b;
    acb_set_si(P.Dh, 2);
    acb_sub(P.Dh, P.Dh, eps.raw(), wp);              // D/2
    acb_mul_2exp_si(P.sg, P.Dh, 1);                  // D
    acb_neg(P.sg, P.sg);
    acb_add_si(P.sg, P.sg, 2 * a + 2 * b + c, wp);   // sg
    {
        acb_t t, g;
        acb_init(t); acb_init(g);
        acb_set_si(t, a); acb_rgamma(P.rgden, t, wp);
        acb_set_si(t, b); acb_rgamma(g, t, wp);
        acb_mul(P.rgden, P.rgden, g, wp);
        acb_set_si(t, c); acb_rgamma(g, t, wp);
        acb_mul(P.rgden, P.rgden, g, wp);
        acb_clear(t); acb_clear(g);
    }
    // c0 (double is plenty: it only needs to separate the pole families)
    double sg_re = arf_get_d(arb_midref(acb_realref(P.sg)), ARF_RND_NEAR);
    double c0 = -0.9 * std::min(std::min((double)a, (double)b),
                                std::min(sg_re / 2.0, 0.5));
    // contour truncation: integrand decays like e^{-2 pi |Im z|}
    double T = 0.12 * (double)wp + 30.0;
    acb_t lo, hi, I;
    acb_init(lo); acb_init(hi); acb_init(I);
    acb_set_d_d(lo, c0, -T);
    acb_set_d_d(hi, c0, T);
    // absolute tolerance keyed to the integrand scale at z = c0
    mag_t tol, mtmp;
    mag_init(tol); mag_init(mtmp);
    {
        acb_t f0, z0;
        acb_init(f0); acb_init(z0);
        acb_set_d(z0, c0);
        eik_c2_integrand(f0, z0, &P, 0, wp);
        acb_get_mag(tol, f0);
        mag_set_ui(mtmp, 1);
        mag_max(tol, tol, mtmp);
        mag_mul_2exp_si(tol, tol, -(prec + 10));
        acb_clear(f0); acb_clear(z0);
    }
    int ok = acb_calc_integrate(I, eik_c2_integrand, &P, lo, hi,
                                wp, tol, nullptr, wp);
    if (ok != 0 || !acb_is_finite(I)) {
        acb_clear(P.Dh); acb_clear(P.sg); acb_clear(P.rgden);
        acb_clear(lo); acb_clear(hi); acb_clear(I);
        mag_clear(tol); mag_clear(mtmp);
        throw std::runtime_error(
            "eikonal_ending_value: class-2 Barnes contour did not converge "
            "for C2(" + std::to_string(a) + "," + std::to_string(b) + ","
            + std::to_string(c) + ") — refusing to return unverified digits");
    }
    // truncation tail: |f| at the endpoints, added to the radius
    {
        acb_t fe;
        acb_init(fe);
        eik_c2_integrand(fe, hi, &P, 0, wp);
        acb_get_mag(mtmp, fe);
        eik_c2_integrand(fe, lo, &P, 0, wp);
        mag_t m2;
        mag_init(m2);
        acb_get_mag(m2, fe);
        mag_max(mtmp, mtmp, m2);
        mag_mul_2exp_si(mtmp, mtmp, 2);
        acb_add_error_mag(I, mtmp);
        mag_clear(m2);
        acb_clear(fe);
    }
    // I / (2 pi i)
    {
        acb_t twopii;
        acb_init(twopii);
        acb_const_pi(twopii, wp);
        acb_mul_onei(twopii, twopii);
        acb_mul_2exp_si(twopii, twopii, 1);
        acb_div(I, I, twopii, wp);
        acb_clear(twopii);
    }
    // residue corrections: z_n = D/2 - a - b + n with Re(z_n) < c0
    {
        acb_t zn, term, t, g, corr;
        acb_init(zn); acb_init(term); acb_init(t); acb_init(g);
        acb_init(corr);
        acb_zero(corr);
        fmpz_t fact;
        fmpz_init(fact);
        fmpz_one(fact);
        for (long n = 0;; ++n) {
            acb_sub_si(zn, P.Dh, a + b - n, wp);      // z_n
            double zr = arf_get_d(arb_midref(acb_realref(zn)), ARF_RND_NEAR);
            if (zr >= c0) break;
            if (n > 0) fmpz_mul_si(fact, fact, n);
            // G(-z_n) G(sg+2 z_n) G(a+z_n) G(b+z_n) rgden rgamma(a+b+2z_n)
            acb_neg(t, zn);
            acb_gamma(term, t, wp);
            acb_mul_2exp_si(t, zn, 1);
            acb_add(t, t, P.sg, wp);
            acb_gamma(g, t, wp);
            acb_mul(term, term, g, wp);
            acb_add_si(t, zn, a, wp);
            acb_gamma(g, t, wp);
            acb_mul(term, term, g, wp);
            acb_add_si(t, zn, b, wp);
            acb_gamma(g, t, wp);
            acb_mul(term, term, g, wp);
            acb_mul_2exp_si(t, zn, 1);
            acb_add_si(t, t, a + b, wp);
            acb_rgamma(g, t, wp);
            acb_mul(term, term, g, wp);
            acb_mul(term, term, P.rgden, wp);
            // * (-1)^n / n!
            acb_div_fmpz(term, term, fact, wp);
            if (n % 2 != 0) acb_neg(term, term);
            acb_add(corr, corr, term, wp);
            if (n > 64) {
                fmpz_clear(fact);
                acb_clear(zn); acb_clear(term); acb_clear(t); acb_clear(g);
                acb_clear(corr);
                acb_clear(P.Dh); acb_clear(P.sg); acb_clear(P.rgden);
                acb_clear(lo); acb_clear(hi); acb_clear(I);
                mag_clear(tol); mag_clear(mtmp);
                throw std::runtime_error(
                    "eikonal_ending_value: class-2 Barnes residue loop "
                    "runaway — unsupported index/eps combination");
            }
        }
        acb_add(I, I, corr, wp);
        fmpz_clear(fact);
        acb_clear(zn); acb_clear(term); acb_clear(t); acb_clear(g);
        acb_clear(corr);
    }
    acb_set_round(res, I, prec);
    acb_clear(P.Dh); acb_clear(P.sg); acb_clear(P.rgden);
    acb_clear(lo); acb_clear(hi); acb_clear(I);
    mag_clear(tol); mag_clear(mtmp);
}

void eik_val_dbg(const char* what, const qft::FamilyConfig& fc,
                 const qft::JIntegral& j, const acb_t v,
                 const numeric::AcbValue& eps) {
    if (!std::getenv("AMFLOW_EIK_DEBUG")) return;
    char* rs = arb_get_str(acb_realref(v), 12, ARB_STR_NO_RADIUS);
    char* is = arb_get_str(acb_imagref(v), 12, ARB_STR_NO_RADIUS);
    char* es = arb_get_str(acb_realref(eps.raw()), 8, ARB_STR_NO_RADIUS);
    std::cerr << "[eik-val] " << what << " fam=" << fc.family
              << " master=" << jintegral_key(j) << " eps~" << es
              << " val=(" << rs << "," << is << ")\n";
    flint_free(rs); flint_free(is); flint_free(es);
}

bool eik_dbg9() {
    static const bool on = (std::getenv("AMF_DBG_EIK9") != nullptr);
    return on;
}

void eik_dbg9_print(const char* what, const acb_t v) {
    if (!eik_dbg9()) return;
    char* rs = arb_get_str(acb_realref(v), 25, 0);
    char* is = arb_get_str(acb_imagref(v), 25, 0);
    std::cerr << "[eik9] " << what << " = " << rs << " + " << is << " i\n";
    flint_free(rs); flint_free(is);
}

// ---------------------------------------------------------------------------
// tanh-sinh quadrature of Int_0^h F(d) dd in the exact DISTANCE variable d
// from the (possibly singular) endpoint d = 0.  Node distances are computed
// through q = e^{-pi sinh u}, d = h q/(1+q) (never 1 - tanh: no
// cancellation; distance-frame rule).  Level doubling; the last level
// difference (x4) is folded into the returned ball's radius — the same
// truncation-control scheme as the reference mpmath oracle; end-to-end
// digit validation is against the reference value tables.  F must be
// integrable with any
// algebraic endpoint behavior no worse than d^p, Re(p) > -1 + margin
// (singular edges are pre-subtracted by eik_int_edge below).
void eik_ts_edge(acb_t res,
                 const std::function<void(acb_t, const arb_t, long)>& F,
                 const arb_t h, long prec, const char* what) {
    const long wp = prec + 40;
    const double sU = 0.2206 * (double)(prec + 40);
    const double Ud = std::log(sU + std::sqrt(sU * sU + 1.0));
    arb_t pi, u, sh, chh, q, opq, d, w, tmp, tauk;
    arb_init(pi); arb_init(u); arb_init(sh); arb_init(chh);
    arb_init(q); arb_init(opq); arb_init(d); arb_init(w); arb_init(tmp);
    arb_init(tauk);
    arb_const_pi(pi, wp);
    acb_t S, Sprev, fv, node;
    acb_init(S); acb_init(Sprev); acb_init(fv); acb_init(node);
    mag_t md, ms, cutoff;
    mag_init(md); mag_init(ms); mag_init(cutoff);
    mag_set_ui_2exp_si(cutoff, 1, -(prec + 50));
    bool converged = false;
    long nodes_used = 0;
    const long kmin = 3, kmax = 13;
    for (long k = kmin; k <= kmax; ++k) {
        arb_one(tauk);
        arb_mul_2exp_si(tauk, tauk, -k);
        const long jmax = (long)std::ceil(Ud * std::pow(2.0, (double)k));
        acb_zero(S);
        for (long j = -jmax; j <= jmax; ++j) {
            arb_set_si(u, j);
            arb_mul_2exp_si(u, u, -k);
            arb_sinh_cosh(sh, chh, u, wp);
            arb_mul(tmp, pi, sh, wp);
            arb_neg(tmp, tmp);
            arb_exp(q, tmp, wp);                      // q = e^{-pi sinh u}
            arb_add_si(opq, q, 1, wp);                // 1+q
            // weight = h (pi/2) cosh(u) 2q/(1+q)^2 * tau_k
            arb_mul(w, pi, chh, wp);
            arb_mul(w, w, q, wp);
            arb_div(w, w, opq, wp);
            arb_div(w, w, opq, wp);
            arb_mul(w, w, h, wp);
            arb_mul(w, w, tauk, wp);
            // negligible-weight cutoff (both wings)
            arb_get_mag(md, w);
            if (mag_cmp(md, cutoff) < 0) continue;
            // d = h q/(1+q)
            arb_div(d, q, opq, wp);
            arb_mul(d, d, h, wp);
            F(fv, d, wp);
            if (!acb_is_finite(fv)) {
                arb_clear(pi); arb_clear(u); arb_clear(sh); arb_clear(chh);
                arb_clear(q); arb_clear(opq); arb_clear(d); arb_clear(w);
                arb_clear(tmp); arb_clear(tauk);
                acb_clear(S); acb_clear(Sprev); acb_clear(fv);
                acb_clear(node);
                mag_clear(md); mag_clear(ms); mag_clear(cutoff);
                throw std::runtime_error(
                    std::string("eikonal_ending_value: non-finite tanh-sinh "
                                "node in ") + what);
            }
            acb_mul_arb(node, fv, w, wp);
            acb_add(S, S, node, wp);
            if (++nodes_used > 2000000) {
                arb_clear(pi); arb_clear(u); arb_clear(sh); arb_clear(chh);
                arb_clear(q); arb_clear(opq); arb_clear(d); arb_clear(w);
                arb_clear(tmp); arb_clear(tauk);
                acb_clear(S); acb_clear(Sprev); acb_clear(fv);
                acb_clear(node);
                mag_clear(md); mag_clear(ms); mag_clear(cutoff);
                throw std::runtime_error(
                    std::string("eikonal_ending_value: tanh-sinh node budget "
                                "exceeded in ") + what);
            }
        }
        if (k > kmin) {
            // convergence is measured on the MIDPOINTS (the discretization
            // error); the accumulated node radii stay in S itself and are
            // returned honestly either way.
            acb_get_mid(fv, S);
            acb_get_mid(node, Sprev);
            acb_sub(fv, fv, node, wp);
            acb_get_mag(md, fv);
            acb_get_mag(ms, S);
            if (eik_dbg9()) {
                std::cerr << "[eik-ts-lvl] " << what << " k=" << k
                          << " |dS|~2^" << (double)mag_get_d_log2_approx(md)
                          << " |S|~2^" << (double)mag_get_d_log2_approx(ms)
                          << "\n";
            }
            mag_mul_2exp_si(ms, ms, -(prec + 2));
            mag_set_ui_2exp_si(cutoff, 1, -(prec + 50));
            // converged when |S_k - S_{k-1}| <= max(2^-(prec+2)|S|, floor,
            // 8 x the accumulated ball radius of S — below that the level
            // difference is inside the noise the radius already carries)
            mag_max(ms, ms, cutoff);
            {
                mag_t mrad;
                mag_init(mrad);
                mag_max(mrad, arb_radref(acb_realref(S)),
                        arb_radref(acb_imagref(S)));
                mag_mul_2exp_si(mrad, mrad, 3);
                mag_max(ms, ms, mrad);
                mag_clear(mrad);
            }
            if (mag_cmp(md, ms) <= 0) {
                converged = true;
                acb_set(res, S);
                mag_mul_2exp_si(md, md, 2);
                acb_add_error_mag(res, md);
                break;
            }
        }
        acb_set(Sprev, S);
    }
    if (!converged) {
        // honest: return the last level with the (large) level difference
        // in the radius — downstream digit checks then fail loudly.
        acb_sub(fv, S, Sprev, wp);
        acb_get_mag(md, fv);
        mag_mul_2exp_si(md, md, 2);
        acb_set(res, S);
        acb_add_error_mag(res, md);
        std::cerr << "[eik-ts] WARNING: tanh-sinh not converged in " << what
                  << " (radius carries the level difference)\n";
    }
    arb_clear(pi); arb_clear(u); arb_clear(sh); arb_clear(chh);
    arb_clear(q); arb_clear(opq); arb_clear(d); arb_clear(w);
    arb_clear(tmp); arb_clear(tauk);
    acb_clear(S); acb_clear(Sprev); acb_clear(fv); acb_clear(node);
    mag_clear(md); mag_clear(ms); mag_clear(cutoff);
}

// Int_0^h d^p A(d) dd with A analytic on [0,h] and generic complex p,
// Re(p) > -2 (p+1 bounded away from 0):
//   = A(0) h^{p+1}/(p+1)  +  Int_0^h d^p (A(d) - A(0)) dd,
// the first term the exact Beta-continuation carrying the edge mass (at
// the flow's eps the exponent sits at -1+O(eps): O(1/eps) mass), the
// remainder ~ d^{p+1} integrable by eik_ts_edge.  A is evaluated at
// elevated precision (the A(d)-A(0) cancellation costs up to prec bits
// at the deepest nodes).
void eik_int_edge(acb_t res, const acb_t p, const arb_t h,
                  const std::function<void(acb_t, const arb_t, long)>& A,
                  long prec, const char* what) {
    const long hp = 2 * prec + 80;
    acb_t A0, p1, hb, piece;
    acb_init(A0); acb_init(p1); acb_init(hb); acb_init(piece);
    {
        arb_t zero;
        arb_init(zero);
        A(A0, zero, hp);
        arb_clear(zero);
    }
    acb_add_si(p1, p, 1, hp);
    if (acb_contains_zero(p1)) {
        acb_clear(A0); acb_clear(p1); acb_clear(hb); acb_clear(piece);
        throw std::runtime_error(
            std::string("eikonal_ending_value: edge exponent pinch (p = -1) "
                        "in ") + what);
    }
    acb_set_arb(hb, h);
    acb_pow(piece, hb, p1, hp);
    acb_mul(piece, piece, A0, hp);
    acb_div(piece, piece, p1, hp);
    acb_t rem;
    acb_init(rem);
    // depth-adaptive precision: the A(d)-A(0) cancellation loses about
    // -log2(d/h) bits, so only the (few, deep) small-d nodes need the full
    // elevated precision — the bulk of the nodes run near prec.
    auto F = [&](acb_t out, const arb_t d, long cbp) {
        (void)cbp;
        if (arb_is_zero(d)) {
            // exact endpoint never sampled by tanh-sinh; belt and braces
            acb_zero(out);
            return;
        }
        long loss = 0;
        {
            arb_t rat;
            arb_init(rat);
            arb_div(rat, d, h, 53);
            const long e = arf_abs_bound_lt_2exp_si(arb_midref(rat));
            if (e < 0) loss = -e;
            arb_clear(rat);
        }
        const long np = std::min(hp, prec + 80 + loss + loss / 8);
        acb_t av, db, dp;
        acb_init(av); acb_init(db); acb_init(dp);
        A(av, d, np);
        acb_sub(av, av, A0, np);
        acb_set_arb(db, d);
        acb_pow(dp, db, p, np);
        acb_mul(out, av, dp, np);
        acb_clear(av); acb_clear(db); acb_clear(dp);
    };
    eik_ts_edge(rem, F, h, prec, what);
    acb_add(res, piece, rem, prec + 40);
    acb_clear(A0); acb_clear(p1); acb_clear(hb); acb_clear(piece);
    acb_clear(rem);
}

// tanh-sinh with LEVEL REUSE (two-velocity E-step; used ONLY by the new class-3
// paths — existing evaluators keep eik_ts_edge so their outputs stay
// byte-identical): S_k = S_{k-1}/2 + tau_k * sum over odd-j nodes.  Same
// node map, weights, cutoff, convergence and radius bookkeeping as
// eik_ts_edge; ~2x fewer integrand evaluations.
void eik_ts_edge3(acb_t res,
                  const std::function<void(acb_t, const arb_t, long)>& F,
                  const arb_t h, long prec, const char* what) {
    const long wp = prec + 40;
    const double sU = 0.2206 * (double)(prec + 40);
    const double Ud = std::log(sU + std::sqrt(sU * sU + 1.0));
    arb_t pi, u, sh, chh, q, opq, d, w, tmp, tauk;
    arb_init(pi); arb_init(u); arb_init(sh); arb_init(chh);
    arb_init(q); arb_init(opq); arb_init(d); arb_init(w); arb_init(tmp);
    arb_init(tauk);
    arb_const_pi(pi, wp);
    acb_t S, Sprev, fv, node;
    acb_init(S); acb_init(Sprev); acb_init(fv); acb_init(node);
    mag_t md, ms, cutoff;
    mag_init(md); mag_init(ms); mag_init(cutoff);
    mag_set_ui_2exp_si(cutoff, 1, -(prec + 50));
    bool converged = false;
    long nodes_used = 0;
    const long kmin = 3, kmax = 13;
    for (long k = kmin; k <= kmax; ++k) {
        arb_one(tauk);
        arb_mul_2exp_si(tauk, tauk, -k);
        const long jmax = (long)std::ceil(Ud * std::pow(2.0, (double)k));
        if (k == kmin) {
            acb_zero(S);
        } else {
            acb_mul_2exp_si(S, S, -1);          // reuse: halved weights
        }
        for (long j = -jmax; j <= jmax; ++j) {
            if (k > kmin && (j % 2 == 0)) continue;   // even = reused
            arb_set_si(u, j);
            arb_mul_2exp_si(u, u, -k);
            arb_sinh_cosh(sh, chh, u, wp);
            arb_mul(tmp, pi, sh, wp);
            arb_neg(tmp, tmp);
            arb_exp(q, tmp, wp);                      // q = e^{-pi sinh u}
            arb_add_si(opq, q, 1, wp);                // 1+q
            arb_mul(w, pi, chh, wp);
            arb_mul(w, w, q, wp);
            arb_div(w, w, opq, wp);
            arb_div(w, w, opq, wp);
            arb_mul(w, w, h, wp);
            arb_mul(w, w, tauk, wp);
            arb_get_mag(md, w);
            if (mag_cmp(md, cutoff) < 0) continue;
            arb_div(d, q, opq, wp);
            arb_mul(d, d, h, wp);
            F(fv, d, wp);
            if (!acb_is_finite(fv)) {
                arb_clear(pi); arb_clear(u); arb_clear(sh); arb_clear(chh);
                arb_clear(q); arb_clear(opq); arb_clear(d); arb_clear(w);
                arb_clear(tmp); arb_clear(tauk);
                acb_clear(S); acb_clear(Sprev); acb_clear(fv);
                acb_clear(node);
                mag_clear(md); mag_clear(ms); mag_clear(cutoff);
                throw std::runtime_error(
                    std::string("eikonal_ending_value: non-finite "
                                "tanh-sinh node in ") + what);
            }
            acb_mul_arb(node, fv, w, wp);
            acb_add(S, S, node, wp);
            if (++nodes_used > 2000000) {
                arb_clear(pi); arb_clear(u); arb_clear(sh); arb_clear(chh);
                arb_clear(q); arb_clear(opq); arb_clear(d); arb_clear(w);
                arb_clear(tmp); arb_clear(tauk);
                acb_clear(S); acb_clear(Sprev); acb_clear(fv);
                acb_clear(node);
                mag_clear(md); mag_clear(ms); mag_clear(cutoff);
                throw std::runtime_error(
                    std::string("eikonal_ending_value: tanh-sinh node "
                                "budget exceeded in ") + what);
            }
        }
        if (k > kmin) {
            acb_get_mid(fv, S);
            acb_get_mid(node, Sprev);
            acb_sub(fv, fv, node, wp);
            acb_get_mag(md, fv);
            acb_get_mag(ms, S);
            if (eik_dbg9()) {
                std::cerr << "[eik-ts3-lvl] " << what << " k=" << k
                          << " |dS|~2^" << (double)mag_get_d_log2_approx(md)
                          << " |S|~2^" << (double)mag_get_d_log2_approx(ms)
                          << "\n";
            }
            mag_mul_2exp_si(ms, ms, -(prec + 2));
            mag_set_ui_2exp_si(cutoff, 1, -(prec + 50));
            mag_max(ms, ms, cutoff);
            {
                mag_t mrad;
                mag_init(mrad);
                mag_max(mrad, arb_radref(acb_realref(S)),
                        arb_radref(acb_imagref(S)));
                // 2^6 margin (eik_ts_edge uses 2^3): with the fixed-
                // precision outer [F-T] nodes the level differences hover
                // at ~16x the accumulated-radius noise floor and a 8x
                // margin burns the segment to kmax (measured); the
                // returned radius carries the difference honestly either
                // way.
                mag_mul_2exp_si(mrad, mrad, 6);
                mag_max(ms, ms, mrad);
                mag_clear(mrad);
            }
            if (mag_cmp(md, ms) <= 0) {
                converged = true;
                acb_set(res, S);
                mag_mul_2exp_si(md, md, 2);
                acb_add_error_mag(res, md);
                break;
            }
        }
        acb_set(Sprev, S);
    }
    if (!converged) {
        acb_sub(fv, S, Sprev, wp);
        acb_get_mag(md, fv);
        mag_mul_2exp_si(md, md, 2);
        acb_set(res, S);
        acb_add_error_mag(res, md);
        std::cerr << "[eik-ts3] WARNING: tanh-sinh not converged in "
                  << what << " (radius carries the level difference)\n";
    }
    arb_clear(pi); arb_clear(u); arb_clear(sh); arb_clear(chh);
    arb_clear(q); arb_clear(opq); arb_clear(d); arb_clear(w);
    arb_clear(tmp); arb_clear(tauk);
    acb_clear(S); acb_clear(Sprev); acb_clear(fv); acb_clear(node);
    mag_clear(md); mag_clear(ms); mag_clear(cutoff);
}

// eik_int_edge on top of the level-reuse quadrature (TC3-only variant).
void eik_int_edge3(acb_t res, const acb_t p, const arb_t h,
                   const std::function<void(acb_t, const arb_t, long)>& A,
                   long prec, const char* what) {
    const long hp = 2 * prec + 80;
    acb_t A0, p1, hb, piece;
    acb_init(A0); acb_init(p1); acb_init(hb); acb_init(piece);
    {
        arb_t zero;
        arb_init(zero);
        A(A0, zero, hp);
        arb_clear(zero);
    }
    acb_add_si(p1, p, 1, hp);
    if (acb_contains_zero(p1)) {
        acb_clear(A0); acb_clear(p1); acb_clear(hb); acb_clear(piece);
        throw std::runtime_error(
            std::string("eikonal_ending_value: edge exponent pinch "
                        "(p = -1) in ") + what);
    }
    acb_set_arb(hb, h);
    acb_pow(piece, hb, p1, hp);
    acb_mul(piece, piece, A0, hp);
    acb_div(piece, piece, p1, hp);
    acb_t rem;
    acb_init(rem);
    auto F = [&](acb_t out, const arb_t d, long cbp) {
        (void)cbp;
        if (arb_is_zero(d)) {
            acb_zero(out);
            return;
        }
        long loss = 0;
        {
            arb_t rat;
            arb_init(rat);
            arb_div(rat, d, h, 53);
            const long e = arf_abs_bound_lt_2exp_si(arb_midref(rat));
            if (e < 0) loss = -e;
            arb_clear(rat);
        }
        const long np = std::min(hp, prec + 80 + loss + loss / 8);
        acb_t av, db, dp;
        acb_init(av); acb_init(db); acb_init(dp);
        A(av, d, np);
        acb_sub(av, av, A0, np);
        acb_set_arb(db, d);
        acb_pow(dp, db, p, np);
        acb_mul(out, av, dp, np);
        acb_clear(av); acb_clear(db); acb_clear(dp);
    };
    eik_ts_edge3(rem, F, h, prec, what);
    acb_add(res, piece, rem, prec + 40);
    acb_clear(A0); acb_clear(p1); acb_clear(hb); acb_clear(piece);
    acb_clear(rem);
}

// ---------------------------------------------------------------------------
// Class-3 (two-velocity E-step) evaluator TC3: the third-class TC3
// derivation (route R3 below), and the bare-eikonal vertex F(eps).
// ---------------------------------------------------------------------------

// F(eps) = -(sqrt(pi)/2) G(1/2+eps) G(1/2-eps)^2 / G(1-2eps): the finite
// 1-loop bare-eikonal vertex Fbare(1,1,1) (the half-eta fix, two
// independent derivations; validated to 75.4 digits).
void eik_bare_F(acb_t res, const numeric::AcbValue& eps, long prec) {
    acb_t t, g, acc;
    acb_init(t); acb_init(g); acb_init(acc);
    // G(1/2+eps)
    acb_one(t);
    acb_mul_2exp_si(t, t, -1);
    acb_add(t, t, eps.raw(), prec);
    acb_gamma(acc, t, prec);
    // * G(1/2-eps)^2
    acb_one(t);
    acb_mul_2exp_si(t, t, -1);
    acb_sub(t, t, eps.raw(), prec);
    acb_gamma(g, t, prec);
    acb_mul(acc, acc, g, prec);
    acb_mul(acc, acc, g, prec);
    // * rgamma(1-2eps)
    acb_one(t);
    acb_submul_si(t, eps.raw(), 2, prec);
    acb_rgamma(g, t, prec);
    acb_mul(acc, acc, g, prec);
    // * -sqrt(pi)/2
    {
        arb_t sp;
        arb_init(sp);
        arb_const_sqrt_pi(sp, prec);
        acb_mul_arb(acc, acc, sp, prec);
        arb_clear(sp);
    }
    acb_mul_2exp_si(acc, acc, -1);
    acb_neg(res, acc);
    acb_clear(t); acb_clear(g); acb_clear(acc);
}

// power-substitution order for an algebraic endpoint exponent p (derived
// k-substitution): d = h w^k turns d^p mass into w^{k(1+p)-1}; k ~ 3/(1+p).
long eik_ksub(double p) {
    if (p <= -1.0) {
        throw std::runtime_error(
            "eikonal_ending_value: non-integrable endpoint exponent in "
            "class-3 quadrature (internal)");
    }
    if (p >= 0.5) return 1;
    long k = (long)std::floor(3.0 / (1.0 + p)) + 1;
    return std::min(k, (long)90);
}

// Int_0^h F(d) dd with an algebraic d^p endpoint at d = 0 handled by the
// EXACT power substitution d = h w^kpow (kpow from eik_ksub; kpow == 1 is
// plain eik_ts_edge).  d is computed from the node w by exact powering
// (distance-frame rule: no cancellation anywhere).
void eik_ts_pow(acb_t res,
                const std::function<void(acb_t, const arb_t, long)>& F,
                const arb_t h, long kpow, long prec, const char* what) {
    if (kpow <= 1) {
        eik_ts_edge3(res, F, h, prec, what);
        return;
    }
    arb_t one;
    arb_init(one);
    arb_one(one);
    auto Fw = [&](acb_t out, const arb_t w, long ep) {
        arb_t d, jac;
        arb_init(d); arb_init(jac);
        arb_pow_ui(jac, w, (ulong)(kpow - 1), ep);   // w^{k-1}
        arb_mul(d, jac, w, ep);                      // w^k
        arb_mul(d, d, h, ep);                        // d = h w^k
        arb_mul(jac, jac, h, ep);
        arb_mul_si(jac, jac, kpow, ep);              // k h w^{k-1}
        F(out, d, ep);
        acb_mul_arb(out, out, jac, ep);
        arb_clear(d); arb_clear(jac);
    };
    eik_ts_edge3(res, Fw, one, prec, what);
    arb_clear(one);
}

// lam = 1 pinch coefficient of the class-3 lam-fold:
//   I(s2) = Int_-inf^inf ((w^2 - 1/4) - i0)^{s2} dw
//         = 2^{-2s2-1} G(s2+1) [ G(-s2-1/2)/G(1/2)
//                               + e^{-i pi s2} G(1/2)/G(s2+3/2) ].
// Validated 8/8 against an independent oracle (incl. s2 -> -1, the
// flow-eps analogue).
void eik_tc3_pinch_I(acb_t res, const acb_t s2, long prec) {
    acb_t t, g, br, acc;
    acb_init(t); acb_init(g); acb_init(br); acb_init(acc);
    arb_t sp;
    arb_init(sp);
    arb_const_sqrt_pi(sp, prec);
    // G(-s2-1/2)/G(1/2)
    acb_neg(t, s2);
    acb_one(g);
    acb_mul_2exp_si(g, g, -1);
    acb_sub(t, t, g, prec);                       // -s2-1/2
    acb_gamma(br, t, prec);
    acb_div_arb(br, br, sp, prec);
    // + e^{-i pi s2} G(1/2)/G(s2+3/2)
    acb_set_si(t, 3);
    acb_mul_2exp_si(t, t, -1);
    acb_add(t, t, s2, prec);                      // s2+3/2
    acb_rgamma(g, t, prec);
    acb_mul_arb(g, g, sp, prec);
    acb_neg(t, s2);
    acb_exp_pi_i(acc, t, prec);                   // e^{-i pi s2}
    acb_mul(g, g, acc, prec);
    acb_add(br, br, g, prec);
    // * 2^{-2s2-1} G(s2+1)
    acb_add_si(t, s2, 1, prec);
    acb_gamma(acc, t, prec);
    acb_mul(br, br, acc, prec);
    acb_mul_2exp_si(t, s2, 1);
    acb_add_si(t, t, 1, prec);
    acb_neg(t, t);                                // -2s2-1
    acb_set_si(g, 2);
    acb_pow(g, g, t, prec);
    acb_mul(res, br, g, prec);
    acb_clear(t); acb_clear(g); acb_clear(br); acb_clear(acc);
    arb_clear(sp);
}

// Root-edge integral Int_0^h B(d) d^{s2} dd with B(d) containing a
// (dm +- d)^{s2} factor (expansion radius dm):
//   (i) exact 2-term subtraction on [0, hs], hs = min(h, dm/2):
//       B0 hs^{s2+1}/(s2+1) + B0 Lg hs^{s2+2}/(s2+2)
//       + Int_0^hs [B(d) - B0 (1 + Lg d)] d^{s2} dd   (B0, Lg EXACT);
//   (ii) tail [hs, h] in the LOG variable d = hs e^{v ln(h/hs)} (spans
//        O(|log dm|) decades near the pinch; an affine segment cannot
//        carry the d^{2s2} mass — measured nan/9d).
// NEVER subtract with h > dm/2 (1e117-scale cancellations, measured).
void eik_tc3_edge2(acb_t res, const acb_t s2, const arb_t h, const arb_t dm,
                   const std::function<void(acb_t, const arb_t, long)>& B,
                   const acb_t B0, const acb_t Lg, long prec,
                   const char* what) {
    const long wp = prec + 40;
    const long hp = 2 * prec + 200;
    arb_t hs, tmp;
    arb_init(hs); arb_init(tmp);
    arb_mul_2exp_si(hs, dm, -1);
    arb_min(hs, hs, h, wp);
    acb_t p1, p2, hb, piece, term;
    acb_init(p1); acb_init(p2); acb_init(hb); acb_init(piece);
    acb_init(term);
    acb_add_si(p1, s2, 1, hp);
    acb_add_si(p2, s2, 2, hp);
    if (acb_contains_zero(p1) || acb_contains_zero(p2)) {
        acb_clear(p1); acb_clear(p2); acb_clear(hb); acb_clear(piece);
        acb_clear(term);
        arb_clear(hs); arb_clear(tmp);
        throw std::runtime_error(
            std::string("eikonal_ending_value: class-3 root-edge exponent "
                        "pinch in ") + what);
    }
    // B0 hs^{s2+1}/(s2+1) + B0 Lg hs^{s2+2}/(s2+2)
    acb_set_arb(hb, hs);
    acb_pow(piece, hb, p1, hp);
    acb_mul(piece, piece, B0, hp);
    acb_div(piece, piece, p1, hp);
    acb_pow(term, hb, p2, hp);
    acb_mul(term, term, B0, hp);
    acb_mul(term, term, Lg, hp);
    acb_div(term, term, p2, hp);
    acb_add(piece, piece, term, hp);
    // remainder on [0, hs].  NOTE precision neutrality: the bracket
    // cancellation at deep nodes loses ~2 log2(hs/d) bits of the node
    // VALUE, but the tanh-sinh weight (~d) shrinks the node's
    // contribution by the same factor — fixed precision is error-neutral
    // (balls stay honest either way; measured 15x cost of the escalated
    // variant at 653 bits).
    acb_t rem;
    acb_init(rem);
    auto Fr = [&](acb_t out, const arb_t d, long cbp) {
        (void)cbp;
        if (arb_is_zero(d)) {
            acb_zero(out);
            return;
        }
        // depth-adaptive: the bracket cancellation needs ~2 log2(hs/d)
        // extra bits or the level differences stall at the cancellation
        // floor and the segment burns to kmax (measured).  This is CHEAP
        // (B is a few pows — unlike the outer [F-T], where the analogous
        // escalation would re-run the whole nested t-fold at high
        // precision; there the ts weight makes fixed precision neutral).
        long loss = 0;
        {
            arb_t rat;
            arb_init(rat);
            arb_div(rat, d, hs, 53);
            const long ee = arf_abs_bound_lt_2exp_si(arb_midref(rat));
            if (ee < 0) loss = -ee;
            arb_clear(rat);
        }
        const long np = std::min(hp, prec + 80 + 2 * loss);
        acb_t bv, sub, db, dp;
        acb_init(bv); acb_init(sub); acb_init(db); acb_init(dp);
        B(bv, d, np);
        acb_mul_arb(sub, Lg, d, np);
        acb_add_si(sub, sub, 1, np);
        acb_mul(sub, sub, B0, np);
        acb_sub(bv, bv, sub, np);
        acb_set_arb(db, d);
        acb_pow(dp, db, s2, np);
        acb_mul(out, bv, dp, np);
        acb_clear(bv); acb_clear(sub); acb_clear(db); acb_clear(dp);
    };
    eik_ts_edge3(rem, Fr, hs, prec, what);
    acb_add(res, piece, rem, wp);
    // tail [hs, h] in the LOG variable
    if (arb_lt(hs, h)) {
        arb_t lr, one;
        arb_init(lr); arb_init(one);
        arb_div(lr, h, hs, wp);
        arb_log(lr, lr, wp);
        arb_one(one);
        auto Ft = [&](acb_t out, const arb_t v, long ep) {
            arb_t d;
            arb_init(d);
            arb_mul(d, v, lr, ep);
            arb_exp(d, d, ep);
            arb_mul(d, d, hs, ep);                 // d = hs e^{v lr}
            acb_t bv, db, dp, p1v;
            acb_init(bv); acb_init(db); acb_init(dp); acb_init(p1v);
            B(bv, d, ep);
            acb_set_arb(db, d);
            acb_add_si(p1v, s2, 1, ep);
            acb_pow(dp, db, p1v, ep);              // d^{s2+1}
            acb_mul(out, bv, dp, ep);
            acb_mul_arb(out, out, lr, ep);
            acb_clear(bv); acb_clear(db); acb_clear(dp); acb_clear(p1v);
            arb_clear(d);
        };
        acb_t tail;
        acb_init(tail);
        eik_ts_edge3(tail, Ft, one, prec, what);
        acb_add(res, res, tail, wp);
        acb_clear(tail);
        arb_clear(lr); arb_clear(one);
    }
    acb_clear(p1); acb_clear(p2); acb_clear(hb); acb_clear(piece);
    acb_clear(term); acb_clear(rem);
    arb_clear(hs); arb_clear(tmp);
}

// TC3(a1,a2,b,c,e; y) by route R3:
//   TC3 = pref3 * Int_0^1 dlam lam^{b-1}(1-lam)^{c-1}(u_lam^2)^{s1}
//         * Int_0^1 dt t^{e-1}(1-t)^{bh3-1}(V_t^2 - i0)^{s2},
// pure power kernel, radiative window t in (taum, taup) with phase
// e^{-i pi s2}, exact product geometry.  Numerics per C.1 (2-term root
// edges + LOG tail; K=2 Beta-continuation at t=1) plus the EXACT lam=1
// pinch subtraction (coefficient I(s2), gate_pinchI/gate_proto_pinch —
// the earlier spec's plain power-substituted lam-fold cannot carry the
// (1-lam)^{c-1+2s2+1} window mass at the flow's eps ~ 1e-14).
// y = 1 (gate use only): lam-fold collapses to B(b,c), double-root
// t-kernel; requires Re(2 s2) > -1.
void eik_tc3_R3(acb_t res, long a1, long a2, long b, long c, long e,
                const fmpq_t yq, const numeric::AcbValue& eps, long prec,
                const std::string& fam) {
    if (b < 1 || c < 1 || e < 1 || a1 < 1 || a2 < 1) {
        throw std::runtime_error(
            "eikonal_ending_value: class-3 TC3 of family '" + fam
            + "' with a non-positive index — outside the R3 domain "
            "(two-velocity E-step scope)");
    }
    if (fmpq_cmp_si(yq, 1) < 0) {
        throw std::runtime_error(
            "eikonal_ending_value: class-3 TC3 of family '" + fam
            + "' at y < 1 — not wired (two-velocity E-step scope)");
    }
    const bool y1 = (fmpq_cmp_si(yq, 1) == 0);
    // evaluate at the boundary-value budget (a measured lesson, iii)
    if (prec > 700) prec = prec / 2 + 128;
    const long wp = prec + 60;
    const long hp = 2 * prec + 200;
    acb_t Dh, s1, s2, bh3, g3, eip, pref, t1, t2, rg;
    acb_init(Dh); acb_init(s1); acb_init(s2); acb_init(bh3);
    acb_init(g3); acb_init(eip); acb_init(pref); acb_init(t1);
    acb_init(t2); acb_init(rg);
    acb_set_si(Dh, 2);
    acb_sub(Dh, Dh, eps.raw(), hp);                  // D/2 = 2 - eps
    acb_set_si(s1, a1);
    acb_sub(s1, s1, Dh, hp);                         // s1 = a1 - D/2
    acb_set_si(s2, a2);
    acb_sub(s2, s2, Dh, hp);                         // s2 = a2 - D/2
    acb_set_si(bh3, 2 * a1 + b + c);
    acb_submul_si(bh3, Dh, 2, hp);                   // bh3 = 2a1+b+c-D
    acb_set_si(g3, 2 * a1 + 2 * a2 + b + c + e);
    acb_submul_si(g3, Dh, 4, hp);                    // g3 = ... - 2D
    const double s2re = arf_get_d(arb_midref(acb_realref(s2)),
                                  ARF_RND_NEAR);
    const double bh3re = arf_get_d(arb_midref(acb_realref(bh3)),
                                   ARF_RND_NEAR);
    // guards (strict, no epsilon margins: flow eps ~ 1e-14; a measured
    // lesson, i)
    if (s2re <= -1.0) {
        acb_clear(Dh); acb_clear(s1); acb_clear(s2); acb_clear(bh3);
        acb_clear(g3); acb_clear(eip); acb_clear(pref); acb_clear(t1);
        acb_clear(t2); acb_clear(rg);
        throw std::runtime_error(
            "eikonal_ending_value: class-3 TC3 of family '" + fam
            + "' at Re(s2) <= -1 — the y+i delta contour does NOT "
            "regulate the lam=1 pinch (a measured negative result, "
            "gate_GD.log: diverges ~ delta^-0.9); no continuation wired "
            "(two-velocity E-step scope; cannot occur on the flow's real-positive "
            "eps grid)");
    }
    if ((double)c + 2.0 * s2re <= -1.0) {
        acb_clear(Dh); acb_clear(s1); acb_clear(s2); acb_clear(bh3);
        acb_clear(g3); acb_clear(eip); acb_clear(pref); acb_clear(t1);
        acb_clear(t2); acb_clear(rg);
        throw std::runtime_error(
            "eikonal_ending_value: class-3 TC3 of family '" + fam
            + "' at Re(c + 2 s2) <= -1 (lam=1 pinch non-integrable) — "
            "not wired (two-velocity E-step scope; cannot occur on the AMFlow "
            "grid)");
    }
    if (y1 && 2.0 * s2re <= -1.0) {
        acb_clear(Dh); acb_clear(s1); acb_clear(s2); acb_clear(bh3);
        acb_clear(g3); acb_clear(eip); acb_clear(pref); acb_clear(t1);
        acb_clear(t2); acb_clear(rg);
        throw std::runtime_error(
            "eikonal_ending_value: class-3 TC3 of family '" + fam
            + "' at y = 1 with Re(2 s2) <= -1 (double-root kernel) — "
            "not wired (two-velocity E-step scope)");
    }
    // pref3 = G(D/2-a1) G(D/2-a2) G(g3) / (G(a1)G(a2)G(b)G(c)G(e))
    acb_sub_si(t1, Dh, a1, hp);
    acb_neg(t1, t1);
    acb_neg(t1, t1);
    acb_gamma(pref, t1, hp);
    acb_sub_si(t1, Dh, a2, hp);
    acb_gamma(t2, t1, hp);
    acb_mul(pref, pref, t2, hp);
    acb_gamma(t2, g3, hp);
    acb_mul(pref, pref, t2, hp);
    {
        const long dens[5] = {a1, a2, b, c, e};
        for (int i = 0; i < 5; ++i) {
            acb_set_si(t1, dens[i]);
            acb_rgamma(rg, t1, hp);
            acb_mul(pref, pref, rg, hp);
        }
    }
    acb_neg(t1, s2);
    acb_exp_pi_i(eip, t1, hp);                       // e^{-i pi s2}
    // real kinematics
    arb_t y, sq;
    arb_init(y); arb_init(sq);
    arb_set_fmpq(y, yq, hp);
    arb_mul(sq, y, y, hp);
    arb_sub_si(sq, sq, 1, hp);
    arb_sqrt(sq, sq, hp);                            // sqrt(y^2-1), 0 at y=1
    const bool tsub = (bh3re <= 0.05);
    const bool edge2 = (s2re <= -0.5) && !y1;
    const bool psub = (s2re <= -0.5) && !y1;         // lam=1 pinch subtr.
    // Ap-free constants
    acb_t Ipin;
    acb_init(Ipin);
    if (psub) eik_tc3_pinch_I(Ipin, s2, hp);

    // ---- inner t-fold Jt at exact (lam, oml), precision ep ----
    auto Jt_eval = [&](acb_t Jt, const arb_t lam, const arb_t oml,
                       long ep) {
        const long ep2 = ep + 40;
        arb_t w, rho2, P, mup, mum, taum, taup, omtaum, omtaup, dm;
        arb_t h1, h2, h3, ta, tb;
        arb_init(w); arb_init(rho2); arb_init(P); arb_init(mup);
        arb_init(mum); arb_init(taum); arb_init(taup); arb_init(omtaum);
        arb_init(omtaup); arb_init(dm); arb_init(h1); arb_init(h2);
        arb_init(h3); arb_init(ta); arb_init(tb);
        arb_mul(w, oml, y, ep2);
        arb_add(w, w, lam, ep2);                     // w = lam + oml*y
        arb_mul(rho2, lam, oml, ep2);
        arb_mul_2exp_si(rho2, rho2, 1);
        arb_sub_si(ta, y, 1, ep2);
        arb_mul(rho2, rho2, ta, ep2);
        arb_add_si(rho2, rho2, 1, ep2);              // 1+2 lam oml (y-1)
        arb_mul_2exp_si(P, w, 1);
        arb_add(P, P, rho2, ep2);
        arb_add_si(P, P, 1, ep2);                    // P = 1+2w+rho2
        arb_mul(ta, oml, sq, ep2);                   // oml*sqrt(y^2-1)
        arb_add(mup, w, ta, ep2);
        arb_sub(mum, w, ta, ep2);
        arb_add_si(tb, mup, 1, ep2);                 // 1+mup
        arb_inv(omtaup, tb, ep2);                    // 1-taup = 1/(1+mup)
        arb_div(taup, mup, tb, ep2);
        arb_add_si(tb, mum, 1, ep2);                 // 1+mum
        arb_inv(omtaum, tb, ep2);                    // 1-taum = 1/(1+mum)
        arb_div(taum, mum, tb, ep2);
        // dm = taup-taum = 2 oml sqrt(y^2-1)/P (exact product form)
        arb_mul_2exp_si(dm, ta, 1);
        arb_div(dm, dm, P, ep2);
        const bool dm0 = arb_is_zero(dm);            // y = 1 only
        arb_mul_2exp_si(h1, taum, -1);               // taum/2
        arb_mul_2exp_si(h2, dm, -1);                 // dm/2
        arb_mul_2exp_si(h3, omtaup, -1);             // (1-taup)/2
        acb_t Ps2, seg, B0, Lg, tt;
        acb_init(Ps2); acb_init(seg); acb_init(B0); acb_init(Lg);
        acb_init(tt);
        acb_set_arb(tt, P);
        acb_pow(Ps2, tt, s2, ep2);                   // P^{s2}
        acb_zero(Jt);
        // --- S1a: plain [0, taum/2] (t from 0; analytic: e >= 1) ---
        auto F1a = [&](acb_t out, const arb_t t, long p_) {
            arb_t d1, d2;
            arb_init(d1); arb_init(d2);
            arb_sub(d1, taum, t, p_);
            arb_sub(d2, taup, t, p_);
            acb_t k1v, k2v, k3v;
            acb_init(k1v); acb_init(k2v); acb_init(k3v);
            acb_set_arb(k1v, t);
            acb_pow_si(k1v, k1v, e - 1, p_);         // t^{e-1}
            arb_t omt;
            arb_init(omt);
            arb_sub_si(omt, t, 1, p_);
            arb_neg(omt, omt);
            acb_set_arb(k2v, omt);
            acb_sub_si(tt, bh3, 1, p_);
            acb_pow(k2v, k2v, tt, p_);               // (1-t)^{bh3-1}
            arb_mul(d1, d1, d2, p_);
            arb_mul(d1, d1, P, p_);
            acb_set_arb(k3v, d1);
            acb_pow(k3v, k3v, s2, p_);               // (P(taum-t)(taup-t))^s2
            acb_mul(out, k1v, k2v, p_);
            acb_mul(out, out, k3v, p_);
            acb_clear(k1v); acb_clear(k2v); acb_clear(k3v);
            arb_clear(d1); arb_clear(d2); arb_clear(omt);
        };
        eik_ts_edge3(seg, F1a, h1, ep, "TC3 t [0,taum/2]");
        acb_add(Jt, Jt, seg, ep2);
        // generic root-edge kernel: t = t0 +- d, sign structure via
        // (P (dm +- d))^{s2}; sgn_t: t = t0 + sgn_t*d; sgn_dm: (dm + sgn_d*d)
        auto Bmk = [&](const arb_t t0, int sgn_t, int sgn_d, bool ph)
            -> std::function<void(acb_t, const arb_t, long)> {
            return [&, sgn_t, sgn_d, ph, t0](acb_t out, const arb_t d,
                                             long p_) {
                arb_t t, omt, dd;
                arb_init(t); arb_init(omt); arb_init(dd);
                if (sgn_t > 0) arb_add(t, t0, d, p_);
                else arb_sub(t, t0, d, p_);
                arb_sub_si(omt, t, 1, p_);
                arb_neg(omt, omt);                   // 1-t
                if (sgn_d > 0) arb_add(dd, dm, d, p_);
                else arb_sub(dd, dm, d, p_);
                arb_mul(dd, dd, P, p_);              // P(dm +- d)
                acb_t k1v, k2v, k3v, ttl;
                acb_init(k1v); acb_init(k2v); acb_init(k3v);
                acb_init(ttl);
                acb_set_arb(k1v, t);
                acb_pow_si(k1v, k1v, e - 1, p_);
                acb_set_arb(k2v, omt);
                acb_sub_si(ttl, bh3, 1, p_);
                acb_pow(k2v, k2v, ttl, p_);
                acb_set_arb(k3v, dd);
                acb_pow(k3v, k3v, s2, p_);
                acb_mul(out, k1v, k2v, p_);
                acb_mul(out, out, k3v, p_);
                if (ph) acb_mul(out, out, eip, p_);
                acb_clear(k1v); acb_clear(k2v); acb_clear(k3v);
                acb_clear(ttl);
                arb_clear(t); arb_clear(omt); arb_clear(dd);
            };
        };
        // exact B0 and log-derivative Lg at a root t0 (tau- or tau+):
        //   B0 = t0^{e-1} (1-t0)^{bh3-1} (P dm)^{s2} [* eip]
        //   Lg = sgn_t [ (e-1)/t0 - (bh3-1)/(1-t0) ] + sgn_d * s2/dm
        auto B0Lg = [&](const arb_t t0, const arb_t omt0, int sgn_t,
                        int sgn_d, bool ph) {
            arb_t pd;
            arb_init(pd);
            arb_mul(pd, P, dm, ep2);
            acb_t k1v, k2v, k3v;
            acb_init(k1v); acb_init(k2v); acb_init(k3v);
            acb_set_arb(k1v, t0);
            acb_pow_si(k1v, k1v, e - 1, ep2);
            acb_set_arb(k2v, omt0);
            acb_sub_si(tt, bh3, 1, ep2);
            acb_pow(k2v, k2v, tt, ep2);
            acb_set_arb(k3v, pd);
            acb_pow(k3v, k3v, s2, ep2);
            acb_mul(B0, k1v, k2v, ep2);
            acb_mul(B0, B0, k3v, ep2);
            if (ph) acb_mul(B0, B0, eip, ep2);
            // Lg
            acb_set_si(Lg, e - 1);
            acb_div_arb(Lg, Lg, t0, ep2);            // (e-1)/t0
            acb_sub_si(tt, bh3, 1, ep2);
            acb_div_arb(tt, tt, omt0, ep2);          // (bh3-1)/(1-t0)
            acb_sub(Lg, Lg, tt, ep2);
            if (sgn_t < 0) acb_neg(Lg, Lg);
            acb_div_arb(tt, s2, dm, ep2);
            if (sgn_d < 0) acb_neg(tt, tt);
            acb_add(Lg, Lg, tt, ep2);
            acb_clear(k1v); acb_clear(k2v); acb_clear(k3v);
            arb_clear(pd);
        };
        // --- S1b: edge at taum from the left: t = taum - d ---
        if (!dm0) {
            if (edge2) {
                B0Lg(taum, omtaum, -1, +1, false);
                eik_tc3_edge2(seg, s2, h1, dm, Bmk(taum, -1, +1, false),
                              B0, Lg, ep, "TC3 edge taum-");
            } else {
                auto Bf = Bmk(taum, -1, +1, false);
                auto Ff = [&](acb_t out, const arb_t d, long p_) {
                    Bf(out, d, p_);
                    acb_t db, dp;
                    acb_init(db); acb_init(dp);
                    acb_set_arb(db, d);
                    acb_pow(dp, db, s2, p_);
                    acb_mul(out, out, dp, p_);
                    acb_clear(db); acb_clear(dp);
                };
                eik_ts_pow(seg, Ff, h1, eik_ksub(s2re), ep,
                           "TC3 edge taum- (plain)");
            }
            acb_add(Jt, Jt, seg, ep2);
            // --- S2a/S2b: window edges (phase e^{-i pi s2}) ---
            if (edge2) {
                B0Lg(taum, omtaum, +1, -1, true);
                eik_tc3_edge2(seg, s2, h2, dm, Bmk(taum, +1, -1, true),
                              B0, Lg, ep, "TC3 window taum+");
                acb_add(Jt, Jt, seg, ep2);
                B0Lg(taup, omtaup, -1, -1, true);
                eik_tc3_edge2(seg, s2, h2, dm, Bmk(taup, -1, -1, true),
                              B0, Lg, ep, "TC3 window taup-");
                acb_add(Jt, Jt, seg, ep2);
            } else {
                auto Bf1 = Bmk(taum, +1, -1, true);
                auto Ff1 = [&](acb_t out, const arb_t d, long p_) {
                    Bf1(out, d, p_);
                    acb_t db, dp;
                    acb_init(db); acb_init(dp);
                    acb_set_arb(db, d);
                    acb_pow(dp, db, s2, p_);
                    acb_mul(out, out, dp, p_);
                    acb_clear(db); acb_clear(dp);
                };
                eik_ts_pow(seg, Ff1, h2, eik_ksub(s2re), ep,
                           "TC3 window taum+ (plain)");
                acb_add(Jt, Jt, seg, ep2);
                auto Bf2 = Bmk(taup, -1, -1, true);
                auto Ff2 = [&](acb_t out, const arb_t d, long p_) {
                    Bf2(out, d, p_);
                    acb_t db, dp;
                    acb_init(db); acb_init(dp);
                    acb_set_arb(db, d);
                    acb_pow(dp, db, s2, p_);
                    acb_mul(out, out, dp, p_);
                    acb_clear(db); acb_clear(dp);
                };
                eik_ts_pow(seg, Ff2, h2, eik_ksub(s2re), ep,
                           "TC3 window taup- (plain)");
                acb_add(Jt, Jt, seg, ep2);
            }
            // --- S3a: edge at taup from the right: t = taup + d ---
            if (edge2) {
                B0Lg(taup, omtaup, +1, +1, false);
                eik_tc3_edge2(seg, s2, h3, dm, Bmk(taup, +1, +1, false),
                              B0, Lg, ep, "TC3 edge taup+");
            } else {
                auto Bf = Bmk(taup, +1, +1, false);
                auto Ff = [&](acb_t out, const arb_t d, long p_) {
                    Bf(out, d, p_);
                    acb_t db, dp;
                    acb_init(db); acb_init(dp);
                    acb_set_arb(db, d);
                    acb_pow(dp, db, s2, p_);
                    acb_mul(out, out, dp, p_);
                    acb_clear(db); acb_clear(dp);
                };
                eik_ts_pow(seg, Ff, h3, eik_ksub(s2re), ep,
                           "TC3 edge taup+ (plain)");
            }
            acb_add(Jt, Jt, seg, ep2);
        } else {
            // y = 1: double root at taum = taup = 1/2; kernel d^{2 s2}
            // on both sides (window absent; Re(2 s2) > -1 guarded)
            auto Fd = [&](int sgn_t) {
                return [&, sgn_t](acb_t out, const arb_t d, long p_) {
                    arb_t t, omt;
                    arb_init(t); arb_init(omt);
                    if (sgn_t > 0) arb_add(t, taum, d, p_);
                    else arb_sub(t, taum, d, p_);
                    arb_sub_si(omt, t, 1, p_);
                    arb_neg(omt, omt);
                    acb_t k1v, k2v, k3v, s22;
                    acb_init(k1v); acb_init(k2v); acb_init(k3v);
                    acb_init(s22);
                    acb_set_arb(k1v, t);
                    acb_pow_si(k1v, k1v, e - 1, p_);
                    acb_set_arb(k2v, omt);
                    acb_sub_si(s22, bh3, 1, p_);
                    acb_pow(k2v, k2v, s22, p_);
                    acb_mul_2exp_si(s22, s2, 1);     // 2 s2
                    acb_t db, dp;
                    acb_init(db); acb_init(dp);
                    acb_set_arb(db, d);
                    acb_pow(dp, db, s22, p_);        // d^{2s2}
                    acb_mul(out, k1v, k2v, p_);
                    acb_mul(out, out, dp, p_);
                    acb_mul(out, out, Ps2, p_);      // P^{s2}
                    acb_clear(k1v); acb_clear(k2v); acb_clear(k3v);
                    acb_clear(s22); acb_clear(db); acb_clear(dp);
                    arb_clear(t); arb_clear(omt);
                };
            };
            const long k2s = eik_ksub(2.0 * s2re);
            eik_ts_pow(seg, Fd(-1), h1, k2s, ep, "TC3 y1 root-");
            acb_add(Jt, Jt, seg, ep2);
            eik_ts_pow(seg, Fd(+1), h3, k2s, ep, "TC3 y1 root+");
            acb_add(Jt, Jt, seg, ep2);
        }
        // --- S3b: t = 1 edge, weight d^{bh3-1}; A(0) = 1 EXACTLY ---
        auto Aa2 = [&](acb_t out, const arb_t d, long p_) {
            arb_t omd, f1, f2;
            arb_init(omd); arb_init(f1); arb_init(f2);
            arb_sub_si(omd, d, 1, p_);
            arb_neg(omd, omd);                       // 1-d
            arb_sub(f1, omtaup, d, p_);              // (1-taup)-d
            arb_sub(f2, omtaum, d, p_);              // (1-taum)-d
            arb_mul(f1, f1, f2, p_);
            arb_mul(f1, f1, P, p_);
            acb_t k1v, k3v;
            acb_init(k1v); acb_init(k3v);
            acb_set_arb(k1v, omd);
            acb_pow_si(k1v, k1v, e - 1, p_);
            acb_set_arb(k3v, f1);
            acb_pow(k3v, k3v, s2, p_);
            acb_mul(out, k1v, k3v, p_);
            acb_clear(k1v); acb_clear(k3v);
            arb_clear(omd); arb_clear(f1); arb_clear(f2);
        };
        if (tsub) {
            // K=2 Beta-continuation; A(1)=1, Ap = -(e-1) - s2 (2+2w)
            acb_t Ap, hb, pb1, pb2, piece;
            acb_init(Ap); acb_init(hb); acb_init(pb1); acb_init(pb2);
            acb_init(piece);
            arb_t twp;
            arb_init(twp);
            arb_mul_2exp_si(twp, w, 1);
            arb_add_si(twp, twp, 2, ep2);            // 2+2w
            acb_mul_arb(Ap, s2, twp, ep2);
            acb_neg(Ap, Ap);
            acb_sub_si(Ap, Ap, e - 1, ep2);
            acb_set_arb(hb, h3);
            acb_pow(pb1, hb, bh3, ep2);
            acb_div(pb1, pb1, bh3, ep2);             // h^{bh3}/bh3
            acb_add_si(tt, bh3, 1, ep2);
            acb_pow(pb2, hb, tt, ep2);
            acb_div(pb2, pb2, tt, ep2);
            acb_mul(pb2, pb2, Ap, ep2);
            acb_add(piece, pb1, pb2, ep2);
            acb_add(Jt, Jt, piece, ep2);
            auto Frem = [&](acb_t out, const arb_t d, long p_) {
                if (arb_is_zero(d)) {
                    acb_zero(out);
                    return;
                }
                // depth-adaptive (cheap pows; see eik_tc3_edge2 note)
                long loss = 0;
                {
                    arb_t rat;
                    arb_init(rat);
                    arb_div(rat, d, h3, 53);
                    const long ee = arf_abs_bound_lt_2exp_si(
                        arb_midref(rat));
                    if (ee < 0) loss = -ee;
                    arb_clear(rat);
                }
                const long np = std::min(2 * ep + 200, p_ + 80 + 2 * loss);
                acb_t av, sub, db, dp, bm1;
                acb_init(av); acb_init(sub); acb_init(db); acb_init(dp);
                acb_init(bm1);
                Aa2(av, d, np);
                acb_mul_arb(sub, Ap, d, np);
                acb_add_si(sub, sub, 1, np);
                acb_sub(av, av, sub, np);
                acb_set_arb(db, d);
                acb_sub_si(bm1, bh3, 1, np);
                acb_pow(dp, db, bm1, np);
                acb_mul(out, av, dp, np);
                acb_clear(av); acb_clear(sub); acb_clear(db);
                acb_clear(dp); acb_clear(bm1);
            };
            eik_ts_edge3(seg, Frem, h3, ep, "TC3 t=1 K2 remainder");
            acb_add(Jt, Jt, seg, ep2);
            acb_clear(Ap); acb_clear(hb); acb_clear(pb1); acb_clear(pb2);
            acb_clear(piece);
            arb_clear(twp);
        } else {
            auto Ff = [&](acb_t out, const arb_t d, long p_) {
                Aa2(out, d, p_);
                acb_t db, dp, bm1;
                acb_init(db); acb_init(dp); acb_init(bm1);
                acb_set_arb(db, d);
                acb_sub_si(bm1, bh3, 1, p_);
                acb_pow(dp, db, bm1, p_);
                acb_mul(out, out, dp, p_);
                acb_clear(db); acb_clear(dp); acb_clear(bm1);
            };
            eik_ts_pow(seg, Ff, h3, eik_ksub(bh3re - 1.0), ep,
                       "TC3 t=1 (plain)");
            acb_add(Jt, Jt, seg, ep2);
        }
        acb_clear(Ps2); acb_clear(seg); acb_clear(B0); acb_clear(Lg);
        acb_clear(tt);
        arb_clear(w); arb_clear(rho2); arb_clear(P); arb_clear(mup);
        arb_clear(mum); arb_clear(taum); arb_clear(taup);
        arb_clear(omtaum); arb_clear(omtaup); arb_clear(dm);
        arb_clear(h1); arb_clear(h2); arb_clear(h3); arb_clear(ta);
        arb_clear(tb);
    };

    acb_t J, seg;
    acb_init(J); acb_init(seg);
    acb_zero(J);
    arb_t half;
    arb_init(half);
    arb_one(half);
    arb_mul_2exp_si(half, half, -1);

    if (y1) {
        // lam-fold collapses: rho2 == 1, Jt lam-independent:
        // J = B(b,c) * Jt
        arb_t lam1, oml1;
        arb_init(lam1); arb_init(oml1);
        arb_mul_2exp_si(lam1, half, 0);              // 1/2 (any lam)
        arb_set(oml1, half);
        Jt_eval(seg, lam1, oml1, wp);
        acb_set_si(t1, b);
        acb_set_si(t2, c);
        acb_t bb, cc, bcs;
        acb_init(bb); acb_init(cc); acb_init(bcs);
        acb_gamma(bb, t1, hp);
        acb_gamma(cc, t2, hp);
        acb_set_si(bcs, b + c);
        acb_rgamma(bcs, bcs, hp);
        acb_mul(bb, bb, cc, hp);
        acb_mul(bb, bb, bcs, hp);                    // B(b,c)
        acb_mul(J, seg, bb, hp);
        acb_clear(bb); acb_clear(cc); acb_clear(bcs);
        arb_clear(lam1); arb_clear(oml1);
    } else {
        // outer weight at (lam, oml), times Jt
        auto Fout = [&](acb_t out, const arb_t lam, const arb_t oml,
                        long p_) {
            arb_t rho2, ta;
            arb_init(rho2); arb_init(ta);
            arb_mul(rho2, lam, oml, p_);
            arb_mul_2exp_si(rho2, rho2, 1);
            arb_sub_si(ta, y, 1, p_);
            arb_mul(rho2, rho2, ta, p_);
            arb_add_si(rho2, rho2, 1, p_);
            acb_t k1v, k2v, k3v, jt;
            acb_init(k1v); acb_init(k2v); acb_init(k3v); acb_init(jt);
            acb_set_arb(k1v, lam);
            acb_pow_si(k1v, k1v, b - 1, p_);
            acb_set_arb(k2v, oml);
            acb_pow_si(k2v, k2v, c - 1, p_);
            acb_set_arb(k3v, rho2);
            acb_pow(k3v, k3v, s1, p_);
            Jt_eval(jt, lam, oml, p_);
            acb_mul(out, k1v, k2v, p_);
            acb_mul(out, out, k3v, p_);
            acb_mul(out, out, jt, p_);
            acb_clear(k1v); acb_clear(k2v); acb_clear(k3v); acb_clear(jt);
            arb_clear(rho2); arb_clear(ta);
        };
        // ---- lam = 0 half: [0, 1/2], analytic (b >= 1 integer) ----
        auto F0 = [&](acb_t out, const arb_t d, long p_) {
            arb_t oml;
            arb_init(oml);
            arb_sub_si(oml, d, 1, p_);
            arb_neg(oml, oml);                       // 1-d
            Fout(out, d, oml, p_);
            arb_clear(oml);
        };
        eik_ts_edge3(seg, F0, half, prec, "TC3 lam [0,1/2]");
        acb_add(J, J, seg, wp);
        // ---- lam = 1 half: d = 1-lam ----
        if (psub) {
            // exact pinch subtraction: T(d) = A_T(d) d^{c+2s2},
            //   A_T = (1-d)^{b-1} rho2^{s1} I(s2) P^{s2} phi(taubar)
            //         (2 sq / P)^{2s2+1},
            //   taubar = (w+rho2)/P, 1-taubar = (1+w)/P (exact forms)
            auto AT = [&](acb_t out, const arb_t d, long p_) {
                arb_t omd, rho2, ta, w, P, taub, omtaub;
                arb_init(omd); arb_init(rho2); arb_init(ta); arb_init(w);
                arb_init(P); arb_init(taub); arb_init(omtaub);
                arb_sub_si(omd, d, 1, p_);
                arb_neg(omd, omd);                   // 1-d = lam
                arb_mul(rho2, omd, d, p_);
                arb_mul_2exp_si(rho2, rho2, 1);
                arb_sub_si(ta, y, 1, p_);
                arb_mul(rho2, rho2, ta, p_);
                arb_add_si(rho2, rho2, 1, p_);
                arb_mul(w, d, y, p_);
                arb_add(w, w, omd, p_);              // w = lam + d*y
                arb_mul_2exp_si(P, w, 1);
                arb_add(P, P, rho2, p_);
                arb_add_si(P, P, 1, p_);
                arb_add(taub, w, rho2, p_);
                arb_div(taub, taub, P, p_);          // (w+rho2)/P
                arb_add_si(omtaub, w, 1, p_);
                arb_div(omtaub, omtaub, P, p_);      // (1+w)/P
                acb_t k1v, k2v, k3v, k4v, k5v, ttl;
                acb_init(k1v); acb_init(k2v); acb_init(k3v);
                acb_init(k4v); acb_init(k5v); acb_init(ttl);
                acb_set_arb(k1v, omd);
                acb_pow_si(k1v, k1v, b - 1, p_);     // (1-d)^{b-1}
                acb_set_arb(k2v, rho2);
                acb_pow(k2v, k2v, s1, p_);           // rho2^{s1}
                acb_set_arb(k3v, P);
                acb_pow(k3v, k3v, s2, p_);           // P^{s2}
                acb_set_arb(k4v, taub);
                acb_pow_si(k4v, k4v, e - 1, p_);     // taub^{e-1}
                acb_set_arb(k5v, omtaub);
                acb_sub_si(ttl, bh3, 1, p_);
                acb_pow(k5v, k5v, ttl, p_);          // (1-taub)^{bh3-1}
                acb_mul(k4v, k4v, k5v, p_);          // phi(taubar)
                // (2 sq / P)^{2s2+1}
                arb_mul_2exp_si(ta, sq, 1);
                arb_div(ta, ta, P, p_);
                acb_set_arb(k5v, ta);
                acb_mul_2exp_si(ttl, s2, 1);
                acb_add_si(ttl, ttl, 1, p_);
                acb_pow(k5v, k5v, ttl, p_);
                acb_mul(out, k1v, k2v, p_);
                acb_mul(out, out, k3v, p_);
                acb_mul(out, out, k4v, p_);
                acb_mul(out, out, k5v, p_);
                acb_mul(out, out, Ipin, p_);
                acb_clear(k1v); acb_clear(k2v); acb_clear(k3v);
                acb_clear(k4v); acb_clear(k5v); acb_clear(ttl);
                arb_clear(omd); arb_clear(rho2); arb_clear(ta);
                arb_clear(w); arb_clear(P); arb_clear(taub);
                arb_clear(omtaub);
            };
            acb_t pexp;
            acb_init(pexp);
            acb_mul_2exp_si(pexp, s2, 1);
            acb_add_si(pexp, pexp, c, hp);           // p = c + 2 s2
            // [F - T] with depth-adaptive precision (the T subtraction
            // cancels the d^{c+2s2} window mass: ~|2s2+1| log2(h/d) bits
            // at deep nodes).  COST NOTE: this re-runs the nested t-fold
            // at up to hp bits on the deep outer nodes — measured
            // 92 s/entry at 320 bits but ~15+ min at the flow's 664 bits;
            // a fixed-precision variant (weight-neutrality argument)
            // stalls the inner segments at a noise floor and burns them
            // to kmax (worse).  Speeding this up is the named follow-up.
            const double cancel = -(2.0 * s2re + 1.0);
            auto FmT = [&](acb_t out, const arb_t d, long p_) {
                if (arb_is_zero(d)) {
                    acb_zero(out);
                    return;
                }
                long loss = 0;
                {
                    arb_t rat;
                    arb_init(rat);
                    arb_div(rat, d, half, 53);
                    const long ee = arf_abs_bound_lt_2exp_si(
                        arb_midref(rat));
                    if (ee < 0) loss = (long)(cancel * (double)(-ee)) + 8;
                    arb_clear(rat);
                }
                const long np = std::min(hp, p_ + 80 + loss);
                arb_t oml;
                arb_init(oml);
                arb_sub_si(oml, d, 1, np);
                arb_neg(oml, oml);                   // lam = 1-d
                acb_t fv, tv, db, dp;
                acb_init(fv); acb_init(tv); acb_init(db); acb_init(dp);
                Fout(fv, oml, d, np);                // (lam, oml=d)
                AT(tv, d, np);
                acb_set_arb(db, d);
                acb_pow(dp, db, pexp, np);
                acb_mul(tv, tv, dp, np);
                acb_sub(out, fv, tv, np);
                acb_clear(fv); acb_clear(tv); acb_clear(db); acb_clear(dp);
                arb_clear(oml);
            };
            eik_ts_edge3(seg, FmT, half, prec, "TC3 lam=1 [F-T]");
            acb_add(J, J, seg, wp);
            // + Int_0^{1/2} A_T(d) d^{c+2s2} dd (K=1 Beta-continuation)
            eik_int_edge3(seg, pexp, half, AT, prec, "TC3 lam=1 pinch");
            acb_add(J, J, seg, wp);
            acb_clear(pexp);
        } else {
            const double pm = (double)(c - 1)
                + std::min(0.0, 1.0 + 2.0 * s2re);
            auto F1 = [&](acb_t out, const arb_t d, long p_) {
                arb_t oml;
                arb_init(oml);
                arb_sub_si(oml, d, 1, p_);
                arb_neg(oml, oml);
                Fout(out, oml, d, p_);
                arb_clear(oml);
            };
            eik_ts_pow(seg, F1, half, eik_ksub(pm), prec,
                       "TC3 lam [1/2,1]");
            acb_add(J, J, seg, wp);
        }
    }
    acb_mul(res, pref, J, wp);
    acb_clear(Dh); acb_clear(s1); acb_clear(s2); acb_clear(bh3);
    acb_clear(g3); acb_clear(eip); acb_clear(pref); acb_clear(t1);
    acb_clear(t2); acb_clear(rg); acb_clear(Ipin);
    acb_clear(J); acb_clear(seg);
    arb_clear(y); arb_clear(sq); arb_clear(half);
}

// V1L3(a; c, ch, bh; y) by route R (reduced 1-fold form).
// a >= 1, c >= 1, ch >= 1 integers; bh = 2eps + n generic (never integer);
// y rational > 1.  Valid for Re(s) > -1 (s = a - D/2) and Re(ch+bh) > -1
// — both hold on the flow's real-positive eps grid and the g1 tables;
// outside that, throws loudly (y+i0 contour route not wired tonight).
void eik_v1l3_R(acb_t res, long a, long c, long ch, const acb_t bh,
                const fmpq_t yq, const numeric::AcbValue& eps, long prec,
                const std::string& fam) {
    if (fmpq_cmp_si(yq, 1) <= 0) {
        throw std::runtime_error(
            "eikonal_ending_value: nested-cusp V1L3 of family '" + fam
            + "' at y_eff <= 1 — pinch/continuation not wired (ending-class "
            "scope)");
    }
    // Evaluation precision: the flow passes working_prec_bits() (= 2x the
    // single_pre fit budget).  The quadrature cost is superlinear in
    // precision, and the fit only consumes ~single_pre digits of the
    // boundary value, so evaluate at prec/2 + 128 bits when the caller
    // asks for more — the returned ball's (honest) radius reflects it.
    if (prec > 700) prec = prec / 2 + 128;
    const long wp = prec + 60;
    const long hp = 2 * prec + 80;
    acb_t Dh, s, g, p, tau, eip, t1, t2;
    acb_init(Dh); acb_init(s); acb_init(g); acb_init(p); acb_init(tau);
    acb_init(eip); acb_init(t1); acb_init(t2);
    acb_set_si(Dh, 2);
    acb_sub(Dh, Dh, eps.raw(), hp);                   // D/2
    acb_set_si(s, a);
    acb_sub(s, s, Dh, hp);                            // s = a - D/2
    acb_mul_2exp_si(g, Dh, 1);
    acb_neg(g, g);
    acb_add_si(g, g, 2 * a + c + ch, hp);
    acb_add(g, g, bh, hp);                            // g = 2a+c+ch+bh-D
    acb_add_si(p, bh, ch - 1, hp);                    // p = ch+bh-1
    acb_neg(tau, p);                                  // tau = 1-ch-bh
    {
        double sre = arf_get_d(arb_midref(acb_realref(s)), ARF_RND_NEAR);
        double pre = arf_get_d(arb_midref(acb_realref(p)), ARF_RND_NEAR);
        // NB no epsilon margin here: the flow's eps samples are ~1e-14, so
        // s = -1 + eps and p = -1 + 2eps sit within any sloppy margin.  The
        // real-axis route is valid for Re(s) > -1 strictly (edge pinches
        // are separately guarded by the exact ball check in eik_int_edge),
        // and the one-term subtraction covers Re(p) > -2.
        if (sre <= -1.0) {
            acb_clear(Dh); acb_clear(s); acb_clear(g); acb_clear(p);
            acb_clear(tau); acb_clear(eip); acb_clear(t1); acb_clear(t2);
            throw std::runtime_error(
                "eikonal_ending_value: nested-cusp V1L3 of family '" + fam
                + "' at Re(s) <= -1 needs the y+i0 contour continuation — "
                "not wired (ending-class scope; cannot occur on the flow's "
                "real-positive eps grid)");
        }
        if (pre <= -1.999) {
            acb_clear(Dh); acb_clear(s); acb_clear(g); acb_clear(p);
            acb_clear(tau); acb_clear(eip); acb_clear(t1); acb_clear(t2);
            throw std::runtime_error(
                "eikonal_ending_value: nested-cusp V1L3 of family '" + fam
                + "' at Re(ch+bh) <= -1 (beyond the one-term edge "
                "subtraction) — not wired (ending-class scope)");
        }
    }
    acb_neg(t1, s);
    acb_exp_pi_i(eip, t1, hp);                        // e^{-i pi s}
    // kinematics: x = y + sqrt(y^2-1), xi = 1/x, w0 = x - xi  (all real)
    arb_t y, x, xi, w0, h1, h2, h3, tmp;
    arb_init(y); arb_init(x); arb_init(xi); arb_init(w0);
    arb_init(h1); arb_init(h2); arb_init(h3); arb_init(tmp);
    arb_set_fmpq(y, yq, hp);
    arb_mul(tmp, y, y, hp);
    arb_sub_si(tmp, tmp, 1, hp);
    arb_sqrt(tmp, tmp, hp);
    arb_add(x, y, tmp, hp);
    arb_inv(xi, x, hp);
    arb_sub(w0, x, xi, hp);                           // > 0 (y > 1)
    arb_mul_2exp_si(h1, xi, -1);
    arb_mul_2exp_si(h3, w0, -1);
    arb_min(h1, h1, h3, hp);                          // h1 = min(xi/2, w0/2)
    arb_set(h2, h1);                                  // h2 = min(w0/2, xi/2)
    // constants of the connection split (tau = 1-ch-bh, never integer)
    //   Bm  = G(ch) G(1-bh+g-ch) rg(1-bh+g)
    //   CAm = G(1-bh+g) G(tau)  rg(1-bh) rg(1-bh+g-ch)
    //   CBm = G(1-bh+g) G(-tau) rg(g) rg(ch)
    //   Cplus = G(1-ch+g-bh) G(ch+bh-1) rg(g)      [ = B(bh,.)*CB+ ]
    acb_t Bm, CAm, CBm, Cplus, Cm1, Cm2;
    acb_init(Bm); acb_init(CAm); acb_init(CBm); acb_init(Cplus);
    acb_init(Cm1); acb_init(Cm2);
    {
        acb_t A1v, A2v, A3v, rg;
        acb_init(A1v); acb_init(A2v); acb_init(A3v); acb_init(rg);
        // A1v = 1-bh+g ; A2v = 1-bh+g-ch ; A3v = 1-ch+g-bh
        acb_sub(A1v, g, bh, hp);
        acb_add_si(A1v, A1v, 1, hp);
        acb_add_si(A2v, A1v, -ch, hp);
        acb_set(A3v, A2v);                            // same quantity
        acb_set_si(t1, ch);
        acb_gamma(Bm, t1, hp);
        acb_gamma(t2, A2v, hp);
        acb_mul(Bm, Bm, t2, hp);
        acb_rgamma(rg, A1v, hp);
        acb_mul(Bm, Bm, rg, hp);
        acb_gamma(CAm, A1v, hp);
        acb_gamma(t2, tau, hp);
        acb_mul(CAm, CAm, t2, hp);
        acb_neg(t1, bh);
        acb_add_si(t1, t1, 1, hp);                    // 1-bh
        acb_rgamma(rg, t1, hp);
        acb_mul(CAm, CAm, rg, hp);
        acb_rgamma(rg, A2v, hp);
        acb_mul(CAm, CAm, rg, hp);
        acb_gamma(CBm, A1v, hp);
        acb_gamma(t2, p, hp);                         // G(-tau) = G(p)
        acb_mul(CBm, CBm, t2, hp);
        acb_rgamma(rg, g, hp);
        acb_mul(CBm, CBm, rg, hp);
        acb_set_si(t1, ch);
        acb_rgamma(rg, t1, hp);
        acb_mul(CBm, CBm, rg, hp);
        acb_gamma(Cplus, A3v, hp);
        acb_gamma(t2, p, hp);
        acb_mul(Cplus, Cplus, t2, hp);
        acb_rgamma(rg, g, hp);
        acb_mul(Cplus, Cplus, rg, hp);
        // Cm1 = Bm*CAm (piece-1 constant), Cm2 = Bm*CBm (piece-2 constant)
        acb_mul(Cm1, Bm, CAm, hp);
        acb_mul(Cm2, Bm, CBm, hp);
        acb_clear(A1v); acb_clear(A2v); acb_clear(A3v); acb_clear(rg);
    }
    // Fb+ polynomial coefficients: 2F1(1-ch, 1-ch+g-bh; 1+tau; w1) =
    //   sum_{k=0}^{ch-1} (1-ch)_k (1-ch+g-bh)_k / ((1+tau)_k k!) w1^k
    std::vector<acb_struct> fbp((std::size_t)ch);
    {
        acb_t coef, num1, num2, den1, kk;
        acb_init(coef); acb_init(num1); acb_init(num2); acb_init(den1);
        acb_init(kk);
        acb_one(coef);
        acb_set_si(num1, 1 - ch);
        acb_sub(num2, g, bh, hp);
        acb_add_si(num2, num2, 1 - ch, hp);
        acb_add_si(den1, tau, 1, hp);
        for (long k = 0; k < ch; ++k) {
            acb_init(&fbp[(std::size_t)k]);
            acb_set(&fbp[(std::size_t)k], coef);
            // coef *= (1-ch+k)(1-ch+g-bh+k) / ((1+tau+k)(k+1))
            acb_add_si(t1, num1, k, hp);
            acb_mul(coef, coef, t1, hp);
            acb_add_si(t1, num2, k, hp);
            acb_mul(coef, coef, t1, hp);
            acb_add_si(t1, den1, k, hp);
            acb_div(coef, coef, t1, hp);
            acb_div_si(coef, coef, k + 1, hp);
        }
        acb_clear(coef); acb_clear(num1); acb_clear(num2); acb_clear(den1);
        acb_clear(kk);
    }
    auto fbp_eval = [&](acb_t out, const acb_t w1, long ep) {
        acb_zero(out);
        for (long k = ch - 1; k >= 0; --k) {
            acb_mul(out, out, w1, ep);
            acb_add(out, out, &fbp[(std::size_t)k], ep);
        }
    };
    // I-(m) via the direct 2F1 (used for m >= ~xi/2; w = (1-m)/(1+m)):
    //   I- = Bm (2m)^p (1+m)^{-g} 2F1(g, ch; 1-bh+g; w)
    auto Iminus = [&](acb_t out, const acb_t m, long ep) {
        acb_t w, om, f, tt;
        acb_init(w); acb_init(om); acb_init(f); acb_init(tt);
        acb_add_si(om, m, 1, ep);                     // 1+m
        acb_sub_si(w, m, 1, ep);
        acb_neg(w, w);
        acb_div(w, w, om, ep);                        // (1-m)/(1+m)
        acb_sub(tt, g, bh, ep);
        acb_add_si(tt, tt, 1, ep);                    // 1-bh+g
        {
            acb_t aa, bb;
            acb_init(aa); acb_init(bb);
            acb_set(aa, g);
            acb_set_si(bb, ch);
            acb_hypgeom_2f1(f, aa, bb, tt, w, 0, ep);
            acb_clear(aa); acb_clear(bb);
        }
        acb_mul(out, f, Bm, ep);
        acb_mul_2exp_si(tt, m, 1);                    // 2m
        acb_pow(tt, tt, p, ep);
        acb_mul(out, out, tt, ep);
        acb_neg(tt, g);
        acb_pow(tt, om, tt, ep);
        acb_mul(out, out, tt, ep);
        acb_clear(w); acb_clear(om); acb_clear(f); acb_clear(tt);
    };
    // ---- J+ (elementary integrand) ---------------------------------------
    //   f+(m) = ((m+x)(m+xi))^s (1+m)^{-(g+tau)} Fb+(2m/(1+m))
    acb_t gt;
    acb_init(gt);
    acb_add(gt, g, tau, hp);                          // g+tau
    auto fplus = [&](acb_t out, const acb_t m, long ep) {
        acb_t f1, f2, om, w1, pv;
        acb_init(f1); acb_init(f2); acb_init(om); acb_init(w1);
        acb_init(pv);
        acb_add_arb(f1, m, x, ep);
        acb_add_arb(f2, m, xi, ep);
        acb_mul(f1, f1, f2, ep);
        acb_pow(f1, f1, s, ep);
        acb_add_si(om, m, 1, ep);
        acb_neg(pv, gt);
        acb_pow(f2, om, pv, ep);
        acb_mul(f1, f1, f2, ep);
        acb_mul_2exp_si(w1, m, 1);
        acb_div(w1, w1, om, ep);
        fbp_eval(f2, w1, ep);
        acb_mul(out, f1, f2, ep);
        acb_clear(f1); acb_clear(f2); acb_clear(om); acb_clear(w1);
        acb_clear(pv);
    };
    acb_t Jp, Jm, seg;
    acb_init(Jp); acb_init(Jm); acb_init(seg);
    acb_zero(Jp);
    {
        // [0,1]
        arb_t one;
        arb_init(one);
        arb_one(one);
        auto F01 = [&](acb_t out, const arb_t d, long ep) {
            acb_t m;
            acb_init(m);
            acb_set_arb(m, d);
            fplus(out, m, ep);
            acb_clear(m);
        };
        eik_ts_edge(seg, F01, one, prec, "V1L3 J+ [0,1]");
        eik_dbg9_print("J+ [0,1]", seg);
        acb_add(Jp, Jp, seg, wp);
        // [1, x+4]: m = 1+d
        arb_t len;
        arb_init(len);
        arb_add_si(len, x, 3, hp);                    // (x+4) - 1
        auto F1x = [&](acb_t out, const arb_t d, long ep) {
            acb_t m;
            acb_init(m);
            acb_set_arb(m, d);
            acb_add_si(m, m, 1, ep);
            fplus(out, m, ep);
            acb_clear(m);
        };
        eik_ts_edge(seg, F1x, len, prec, "V1L3 J+ [1,x+4]");
        eik_dbg9_print("J+ [1,x+4]", seg);
        acb_add(Jp, Jp, seg, wp);
        // tail via t = 1/m on (0, 1/(x+4)]:
        //   t^{c-1} ((1-x t)(1-xi t))^s (1+t)^{-(g+tau)} Fb+(2/(1+t))
        arb_t tl;
        arb_init(tl);
        arb_add_si(tl, x, 4, hp);
        arb_inv(tl, tl, hp);
        auto Ft = [&](acb_t out, const arb_t td, long ep) {
            acb_t tb, f1, f2, ot, w1, pv;
            acb_init(tb); acb_init(f1); acb_init(f2); acb_init(ot);
            acb_init(w1); acb_init(pv);
            acb_set_arb(tb, td);
            acb_mul_arb(f1, tb, x, ep);
            acb_add_si(f1, f1, 1, ep);                // 1 + x t
            acb_mul_arb(f2, tb, xi, ep);
            acb_add_si(f2, f2, 1, ep);                // 1 + xi t
            acb_mul(f1, f1, f2, ep);
            acb_pow(f1, f1, s, ep);
            acb_add_si(ot, tb, 1, ep);                // 1+t
            acb_neg(pv, gt);
            acb_pow(f2, ot, pv, ep);
            acb_mul(f1, f1, f2, ep);
            acb_set_si(w1, 2);
            acb_div(w1, w1, ot, ep);
            fbp_eval(f2, w1, ep);
            acb_mul(f1, f1, f2, ep);
            if (c > 1) {
                acb_pow_si(f2, tb, c - 1, ep);
                acb_mul(f1, f1, f2, ep);
            }
            acb_set(out, f1);
            acb_clear(tb); acb_clear(f1); acb_clear(f2); acb_clear(ot);
            acb_clear(w1); acb_clear(pv);
        };
        eik_ts_edge(seg, Ft, tl, prec, "V1L3 J+ tail");
        eik_dbg9_print("J+ tail", seg);
        acb_add(Jp, Jp, seg, wp);
        arb_clear(one); arb_clear(len); arb_clear(tl);
    }
    acb_mul(Jp, Jp, Cplus, wp);
    // ---- J- ---------------------------------------------------------------
    acb_zero(Jm);
    // K-(m) on [0, xi): ((x-m)(xi-m))^s, safe direct subtractions for
    // m <= xi - h1 (distances >= h1 ~ O(xi)).
    // seg0 [0, m0], m0 = xi/2: exact 1-z connection split.
    //   piece1: m^p * A1(m), A1 = 2^p Cm1 K-(m) (1+m)^{-g} Fa(w1)
    //   piece2:        A2(m),  A2 = Cm2 K-(m) (1+m)^{-(g+tau)} Fb-(w1)
    //   Fa  = 2F1(g, ch; ch+bh; w1),  Fb- = 2F1(1-bh, 1-bh+g-ch; 1+tau; w1)
    {
        arb_t m0;
        arb_init(m0);
        arb_mul_2exp_si(m0, xi, -1);
        auto Km_lo = [&](acb_t out, const acb_t m, long ep) {
            acb_t f1, f2;
            acb_init(f1); acb_init(f2);
            acb_neg(f1, m);
            acb_add_arb(f1, f1, x, ep);               // x - m
            acb_neg(f2, m);
            acb_add_arb(f2, f2, xi, ep);              // xi - m
            acb_mul(f1, f1, f2, ep);
            acb_pow(out, f1, s, ep);
            acb_clear(f1); acb_clear(f2);
        };
        auto A1 = [&](acb_t out, const arb_t d, long ep) {
            acb_t m, om, w1, f, kv, aa, bb, cc;
            acb_init(m); acb_init(om); acb_init(w1); acb_init(f);
            acb_init(kv); acb_init(aa); acb_init(bb); acb_init(cc);
            acb_set_arb(m, d);
            Km_lo(kv, m, ep);
            acb_add_si(om, m, 1, ep);
            acb_mul_2exp_si(w1, m, 1);
            acb_div(w1, w1, om, ep);
            acb_set(aa, g);
            acb_set_si(bb, ch);
            acb_add_si(cc, bh, ch, ep);               // ch+bh = 1-tau
            acb_hypgeom_2f1(f, aa, bb, cc, w1, 0, ep);
            acb_mul(f, f, kv, ep);
            acb_neg(kv, g);
            acb_pow(kv, om, kv, ep);
            acb_mul(f, f, kv, ep);
            acb_mul(f, f, Cm1, ep);
            acb_set_si(kv, 2);
            acb_pow(kv, kv, p, ep);
            acb_mul(out, f, kv, ep);
            acb_clear(m); acb_clear(om); acb_clear(w1); acb_clear(f);
            acb_clear(kv); acb_clear(aa); acb_clear(bb); acb_clear(cc);
        };
        eik_int_edge(seg, p, m0, A1, prec, "V1L3 J- [0,m0] m^p piece");
        eik_dbg9_print("J- [0,m0] m^p", seg);
        acb_add(Jm, Jm, seg, wp);
        auto A2 = [&](acb_t out, const arb_t d, long ep) {
            acb_t m, om, w1, f, kv, aa, bb, cc;
            acb_init(m); acb_init(om); acb_init(w1); acb_init(f);
            acb_init(kv); acb_init(aa); acb_init(bb); acb_init(cc);
            acb_set_arb(m, d);
            Km_lo(kv, m, ep);
            acb_add_si(om, m, 1, ep);
            acb_mul_2exp_si(w1, m, 1);
            acb_div(w1, w1, om, ep);
            acb_neg(aa, bh);
            acb_add_si(aa, aa, 1, ep);                // 1-bh
            acb_sub(bb, g, bh, ep);
            acb_add_si(bb, bb, 1 - ch, ep);           // 1-bh+g-ch
            acb_add_si(cc, tau, 1, ep);               // 1+tau
            acb_hypgeom_2f1(f, aa, bb, cc, w1, 0, ep);
            acb_mul(f, f, kv, ep);
            acb_neg(kv, gt);
            acb_pow(kv, om, kv, ep);
            acb_mul(f, f, kv, ep);
            acb_mul(out, f, Cm2, ep);
            acb_clear(m); acb_clear(om); acb_clear(w1); acb_clear(f);
            acb_clear(kv); acb_clear(aa); acb_clear(bb); acb_clear(cc);
        };
        eik_ts_edge(seg, A2, m0, prec, "V1L3 J- [0,m0] analytic piece");
        eik_dbg9_print("J- [0,m0] analytic", seg);
        acb_add(Jm, Jm, seg, wp);
        // seg1 [m0, xi-h1] (may be empty)
        arb_t lo1, len1;
        arb_init(lo1); arb_init(len1);
        arb_set(lo1, m0);
        arb_sub(len1, xi, h1, hp);
        arb_sub(len1, len1, m0, hp);
        if (arf_sgn(arb_midref(len1)) > 0) {
            auto Fm1 = [&](acb_t out, const arb_t d, long ep) {
                acb_t m, iv;
                acb_init(m); acb_init(iv);
                acb_set_arb(m, d);
                acb_add_arb(m, m, lo1, ep);
                Km_lo(out, m, ep);
                Iminus(iv, m, ep);
                acb_mul(out, out, iv, ep);
                acb_clear(m); acb_clear(iv);
            };
            eik_ts_edge(seg, Fm1, len1, prec, "V1L3 J- [m0,xi-h1]");
        eik_dbg9_print("J- [m0,xi-h1]", seg);
            acb_add(Jm, Jm, seg, wp);
        }
        // seg2 [xi-h1, xi]: d from xi, K = (d (w0+d))^s
        auto A3 = [&](acb_t out, const arb_t d, long ep) {
            acb_t m, f;
            acb_init(m); acb_init(f);
            acb_set_arb(m, d);
            acb_neg(m, m);
            acb_add_arb(m, m, xi, ep);                // m = xi - d
            acb_set_arb(f, d);
            acb_add_arb(f, f, w0, ep);                // w0 + d
            acb_pow(f, f, s, ep);
            Iminus(out, m, ep);
            acb_mul(out, out, f, ep);
            acb_clear(m); acb_clear(f);
        };
        eik_int_edge(seg, s, h1, A3, prec, "V1L3 J- edge xi-");
        eik_dbg9_print("J- edge xi-", seg);
        acb_add(Jm, Jm, seg, wp);
        // seg3 [xi, xi+h2]: inside the window, K = eip (d (w0-d))^s
        auto A4 = [&](acb_t out, const arb_t d, long ep) {
            acb_t m, f;
            acb_init(m); acb_init(f);
            acb_set_arb(m, d);
            acb_add_arb(m, m, xi, ep);                // m = xi + d
            acb_set_arb(f, d);
            acb_neg(f, f);
            acb_add_arb(f, f, w0, ep);                // w0 - d
            acb_pow(f, f, s, ep);
            acb_mul(f, f, eip, ep);
            Iminus(out, m, ep);
            acb_mul(out, out, f, ep);
            acb_clear(m); acb_clear(f);
        };
        eik_int_edge(seg, s, h2, A4, prec, "V1L3 J- edge xi+");
        eik_dbg9_print("J- edge xi+", seg);
        acb_add(Jm, Jm, seg, wp);
        // seg4 [xi+h2, x-h3] (may be empty): K = eip ((m-xi)(x-m))^s with
        // exact distances (m-xi) = h2+d, (x-m) = (len4-d)+h3
        arb_t lo4, len4;
        arb_init(lo4); arb_init(len4);
        arb_add(lo4, xi, h2, hp);
        arb_sub(len4, x, h3, hp);
        arb_sub(len4, len4, lo4, hp);
        if (arf_sgn(arb_midref(len4)) > 0) {
            auto Fm4 = [&](acb_t out, const arb_t d, long ep) {
                acb_t m, f1, f2;
                acb_init(m); acb_init(f1); acb_init(f2);
                acb_set_arb(m, d);
                acb_add_arb(m, m, lo4, ep);           // m = lo4 + d
                acb_set_arb(f1, d);
                acb_add_arb(f1, f1, h2, ep);          // m - xi
                acb_set_arb(f2, d);
                acb_neg(f2, f2);
                acb_add_arb(f2, f2, len4, ep);
                acb_add_arb(f2, f2, h3, ep);          // x - m
                acb_mul(f1, f1, f2, ep);
                acb_pow(f1, f1, s, ep);
                acb_mul(f1, f1, eip, ep);
                Iminus(out, m, ep);
                acb_mul(out, out, f1, ep);
                acb_clear(m); acb_clear(f1); acb_clear(f2);
            };
            eik_ts_edge(seg, Fm4, len4, prec, "V1L3 J- [xi+h2,x-h3]");
        eik_dbg9_print("J- mid", seg);
            acb_add(Jm, Jm, seg, wp);
        }
        // seg5 [x-h3, x]: d from x inside the window
        auto A5 = [&](acb_t out, const arb_t d, long ep) {
            acb_t m, f;
            acb_init(m); acb_init(f);
            acb_set_arb(m, d);
            acb_neg(m, m);
            acb_add_arb(m, m, x, ep);                 // m = x - d
            acb_set_arb(f, d);
            acb_neg(f, f);
            acb_add_arb(f, f, w0, ep);                // w0 - d
            acb_pow(f, f, s, ep);
            acb_mul(f, f, eip, ep);
            Iminus(out, m, ep);
            acb_mul(out, out, f, ep);
            acb_clear(m); acb_clear(f);
        };
        eik_int_edge(seg, s, h3, A5, prec, "V1L3 J- edge x-");
        eik_dbg9_print("J- edge x-", seg);
        acb_add(Jm, Jm, seg, wp);
        // seg6 [x, x+4]: d from x outside the window
        arb_t four;
        arb_init(four);
        arb_set_si(four, 4);
        auto A6 = [&](acb_t out, const arb_t d, long ep) {
            acb_t m, f;
            acb_init(m); acb_init(f);
            acb_set_arb(m, d);
            acb_add_arb(m, m, x, ep);                 // m = x + d
            acb_set_arb(f, d);
            acb_add_arb(f, f, w0, ep);                // w0 + d
            acb_pow(f, f, s, ep);
            Iminus(out, m, ep);
            acb_mul(out, out, f, ep);
            acb_clear(m); acb_clear(f);
        };
        eik_int_edge(seg, s, four, A6, prec, "V1L3 J- edge x+");
        eik_dbg9_print("J- edge x+", seg);
        acb_add(Jm, Jm, seg, wp);
        // seg7 tail via t = 1/m on (0, 1/(x+4)]:
        //   2^p t^{c-1} ((1-x t)(1-xi t))^s (1+t)^{-g} Bm
        //       * 2F1(g, ch; 1-bh+g; -(1-t)/(1+t))
        arb_t tl;
        arb_init(tl);
        arb_add_si(tl, x, 4, hp);
        arb_inv(tl, tl, hp);
        auto Ftm = [&](acb_t out, const arb_t td, long ep) {
            acb_t tb, f1, f2, ot, w, cc;
            acb_init(tb); acb_init(f1); acb_init(f2); acb_init(ot);
            acb_init(w); acb_init(cc);
            acb_set_arb(tb, td);
            acb_mul_arb(f1, tb, x, ep);
            acb_neg(f1, f1);
            acb_add_si(f1, f1, 1, ep);                // 1 - x t
            acb_mul_arb(f2, tb, xi, ep);
            acb_neg(f2, f2);
            acb_add_si(f2, f2, 1, ep);                // 1 - xi t
            acb_mul(f1, f1, f2, ep);
            acb_pow(f1, f1, s, ep);
            acb_add_si(ot, tb, 1, ep);                // 1+t
            acb_sub_si(w, tb, 1, ep);
            acb_div(w, w, ot, ep);                    // (t-1)/(1+t) = -(1-t)/(1+t)
            acb_sub(cc, g, bh, ep);
            acb_add_si(cc, cc, 1, ep);
            {
                acb_t aa, bb, f;
                acb_init(aa); acb_init(bb); acb_init(f);
                acb_set(aa, g);
                acb_set_si(bb, ch);
                acb_hypgeom_2f1(f, aa, bb, cc, w, 0, ep);
                acb_mul(f1, f1, f, ep);
                acb_clear(aa); acb_clear(bb); acb_clear(f);
            }
            acb_neg(f2, g);
            acb_pow(f2, ot, f2, ep);
            acb_mul(f1, f1, f2, ep);
            acb_mul(f1, f1, Bm, ep);
            acb_set_si(f2, 2);
            acb_pow(f2, f2, p, ep);
            acb_mul(f1, f1, f2, ep);
            if (c > 1) {
                acb_pow_si(f2, tb, c - 1, ep);
                acb_mul(f1, f1, f2, ep);
            }
            acb_set(out, f1);
            acb_clear(tb); acb_clear(f1); acb_clear(f2); acb_clear(ot);
            acb_clear(w); acb_clear(cc);
        };
        eik_ts_edge(seg, Ftm, tl, prec, "V1L3 J- tail");
        eik_dbg9_print("J- tail", seg);
        acb_add(Jm, Jm, seg, wp);
        arb_clear(m0); arb_clear(lo1); arb_clear(len1);
        arb_clear(lo4); arb_clear(len4); arb_clear(four); arb_clear(tl);
    }
    // ---- assemble: pref * 2^{1-ch-bh} * (J+ + J-) -------------------------
    {
        acb_t pref, rg, arg;
        acb_init(pref); acb_init(rg); acb_init(arg);
        acb_sub_si(arg, Dh, a, hp);
        acb_neg(arg, arg);
        acb_neg(arg, arg);
        acb_gamma(pref, arg, hp);                     // G(D/2-a)
        acb_gamma(rg, g, hp);
        acb_mul(pref, pref, rg, hp);
        acb_set_si(arg, a);
        acb_rgamma(rg, arg, hp);
        acb_mul(pref, pref, rg, hp);
        acb_set_si(arg, c);
        acb_rgamma(rg, arg, hp);
        acb_mul(pref, pref, rg, hp);
        acb_set_si(arg, ch);
        acb_rgamma(rg, arg, hp);
        acb_mul(pref, pref, rg, hp);
        acb_rgamma(rg, bh, hp);
        acb_mul(pref, pref, rg, hp);
        // * 2^{1-ch-bh} = 2^{tau}
        acb_set_si(arg, 2);
        acb_pow(arg, arg, tau, hp);
        acb_mul(pref, pref, arg, hp);
        acb_add(rg, Jp, Jm, hp);
        acb_mul(pref, pref, rg, hp);
        acb_set_round(res, pref, prec);
        acb_clear(pref); acb_clear(rg); acb_clear(arg);
    }
    for (long k = 0; k < ch; ++k) acb_clear(&fbp[(std::size_t)k]);
    acb_clear(Jp); acb_clear(Jm); acb_clear(seg);
    acb_clear(Bm); acb_clear(CAm); acb_clear(CBm); acb_clear(Cplus);
    acb_clear(Cm1); acb_clear(Cm2);
    acb_clear(gt);
    arb_clear(y); arb_clear(x); arb_clear(xi); arb_clear(w0);
    arb_clear(h1); arb_clear(h2); arb_clear(h3); arb_clear(tmp);
    acb_clear(Dh); acb_clear(s); acb_clear(g); acb_clear(p);
    acb_clear(tau); acb_clear(eip); acb_clear(t1); acb_clear(t2);
}

// Numeric rational scalar product of two reduced legs (family replacement
// rules + blackbox numeric_values).  va/vb are context VARIABLE indices.
void eik_leg_dot(const qft::FamilyConfig& fc, const AMFSystemOptions& opts,
                 const numeric::AcbValue& eps, long prec,
                 long va, long vb, fmpq_t res) {
    const long L = (long)fc.n_loops();
    const std::string& na = fc.reduced_legs[(std::size_t)(va - L)];
    const std::string& nb = fc.reduced_legs[(std::size_t)(vb - L)];
    const std::string key = (na <= nb) ? (na + "*" + nb) : (nb + "*" + na);
    auto it = fc.reduced_replacement.find(key);
    if (it == fc.reduced_replacement.end()) {
        throw std::runtime_error(
            "eikonal_ending_value: family '" + fc.family
            + "' has no replacement rule for the velocity product '" + key
            + "' needed by a two-velocity (cusp) ending");
    }
    auto numeric_q = build_numeric_q_for_eps(opts.bb, eps, prec);
    std::set<std::string> empty_keep;
    algebra::Mfrac sub = mfrac_substitute(it->second, numeric_q, empty_keep);
    algebra::Mpoly num = sub.numerator();
    algebra::Mpoly den = sub.denominator();
    if (!fmpz_mpoly_is_fmpz(num.raw(), sub.ctx()->raw())
            || !fmpz_mpoly_is_fmpz(den.raw(), sub.ctx()->raw())) {
        throw std::runtime_error(
            "eikonal_ending_value: velocity product '" + key
            + "' does not evaluate to a rational number under "
            "numeric_values (got '" + sub.to_string() + "')");
    }
    fmpz_t nz, dz;
    fmpz_init(nz); fmpz_init(dz);
    std::vector<unsigned long> zexp((std::size_t)sub.ctx()->n_vars(), 0);
    fmpz_mpoly_get_coeff_fmpz_ui(nz, num.raw(), zexp.data(), sub.ctx()->raw());
    fmpz_mpoly_get_coeff_fmpz_ui(dz, den.raw(), zexp.data(), sub.ctx()->raw());
    if (fmpz_is_zero(dz)) {
        fmpz_clear(nz); fmpz_clear(dz);
        throw std::runtime_error(
            "eikonal_ending_value: zero denominator for '" + key + "'");
    }
    fmpq_set_fmpz_frac(res, nz, dz);
    fmpz_clear(nz); fmpz_clear(dz);
}

// T(a,b+c) * 3F2(2-eps-a, b, c; (b+c)/2, (b+c+1)/2; z), z rational.
// Dispatch: -1 < z < 1 defining series; z > 1 the +i0 continuation:
//   * `cont_ok` (c carries eps, non-degenerate): NIST 16.8.8 1/z formula;
//   * integer (b,c) passed as (bi,ci) >= 1: parameter-degenerate, but the
//     3F2 COLLAPSES to a 2F1 for (1,1)/(2,1)/(1,2) [-> 2F1(2-eps-a,1;3/2;z)]
//     and (2,2) [-> 2F1(2-eps-a,2;5/2;z)] via Pochhammer cancellation
//     (verified independently, >= 30d); evaluated on the z+i0 side with the
//     1-z connection formula (needs z < 2, i.e. |y| < 3).
// z == 1 (pinch), z <= -1, and unreduced degenerate shapes throw.
void eik_pair_value(acb_t res, const acb_t a, const acb_t b, const acb_t c,
                    const fmpq_t zq, bool cont_ok, long bi, long ci,
                    const numeric::AcbValue& eps, long prec,
                    const std::string& fam) {
    acb_ptr al = _acb_vec_init(3);
    acb_ptr rh = _acb_vec_init(2);
    acb_t z, bc, F, Tv;
    acb_init(z); acb_init(bc); acb_init(F); acb_init(Tv);
    acb_set_fmpq(z, zq, prec);
    acb_add(bc, b, c, prec);
    // al = {2-eps-a, b, c}; rh = {(b+c)/2, (b+c+1)/2}
    acb_set_si(al + 0, 2);
    acb_sub(al + 0, al + 0, eps.raw(), prec);
    acb_sub(al + 0, al + 0, a, prec);
    acb_set(al + 1, b);
    acb_set(al + 2, c);
    acb_mul_2exp_si(rh + 0, bc, -1);
    acb_add_si(rh + 1, bc, 1, prec);
    acb_mul_2exp_si(rh + 1, rh + 1, -1);
    const int cmp_hi = fmpq_cmp_si(zq, 1);
    const int cmp_lo = fmpq_cmp_si(zq, -1);
    if (cmp_hi < 0 && cmp_lo > 0) {
        acb_hypgeom_pfq(F, al, 3, rh, 2, z, 0, prec);
    } else if (cmp_hi > 0) {
        if (cont_ok) {
            eik_hyp3f2_cont_above(F, al, rh, z, prec);
        } else {
            // integer (b,c): 3F2 -> 2F1 reductions
            const bool r32 = (bi == 1 && ci == 1) || (bi == 2 && ci == 1)
                          || (bi == 1 && ci == 2);
            const bool r52 = (bi == 2 && ci == 2);
            if ((!r32 && !r52) || fmpq_cmp_si(zq, 2) >= 0) {
                _acb_vec_clear(al, 3); _acb_vec_clear(rh, 2);
                acb_clear(z); acb_clear(bc); acb_clear(F); acb_clear(Tv);
                throw std::runtime_error(
                    "eikonal_ending_value: mixed-orientation two-velocity "
                    "ending of family '" + fam + "' with (b,c)=("
                    + std::to_string(bi) + "," + std::to_string(ci)
                    + ") or z >= 2 — 2F1 reduction not wired (stage-2 "
                    "ending scope)");
            }
            acb_t be, ga;
            acb_init(be); acb_init(ga);
            if (r32) {
                acb_one(be);
                acb_set_si(ga, 3);
            } else {
                acb_set_si(be, 2);
                acb_set_si(ga, 5);
            }
            acb_mul_2exp_si(ga, ga, -1);       // 3/2 or 5/2
            eik_hyp2f1_above(F, al + 0, be, ga, z, prec);
            acb_clear(be); acb_clear(ga);
        }
    } else if (cmp_hi == 0) {
        // z == 1: same-velocity mixed-orientation / BG-class crossed
        // coupling.  ((1-2lam)^2 - i0)^s == |1-2lam|^{2s} (positive base),
        // closed form via Euler integrals at argument -1
        // (verified verify_z1_bg*.log, >= 30d):
        //   F(1) = G(b+c) G(2s+1) 2^{1-b-c}
        //          * [ rg(c) rg(2s+1+b) 2F1(1-c, 2s+1; 2s+1+b; -1)
        //            + rg(b) rg(2s+1+c) 2F1(1-b, 2s+1; 2s+1+c; -1) ],
        //   s = a-2+eps.
        acb_t s2, t1, t2, g, F2, one_m;
        acb_init(s2); acb_init(t1); acb_init(t2); acb_init(g);
        acb_init(F2); acb_init(one_m);
        acb_add_si(s2, a, -2, prec);
        acb_add(s2, s2, eps.raw(), prec);      // s
        acb_mul_2exp_si(s2, s2, 1);
        acb_add_si(s2, s2, 1, prec);           // 2s+1
        acb_set_si(one_m, -1);
        // term 1: rg(c) rg(2s+1+b) 2F1(1-c, 2s+1; 2s+1+b; -1)
        acb_neg(t1, c);
        acb_add_si(t1, t1, 1, prec);           // 1-c
        acb_add(t2, s2, b, prec);              // 2s+1+b
        acb_hypgeom_2f1(F2, t1, s2, t2, one_m, 0, prec);
        acb_rgamma(g, c, prec);
        acb_mul(F2, F2, g, prec);
        acb_rgamma(g, t2, prec);
        acb_mul(F2, F2, g, prec);
        acb_set(F, F2);
        // term 2: rg(b) rg(2s+1+c) 2F1(1-b, 2s+1; 2s+1+c; -1)
        acb_neg(t1, b);
        acb_add_si(t1, t1, 1, prec);           // 1-b
        acb_add(t2, s2, c, prec);              // 2s+1+c
        acb_hypgeom_2f1(F2, t1, s2, t2, one_m, 0, prec);
        acb_rgamma(g, b, prec);
        acb_mul(F2, F2, g, prec);
        acb_rgamma(g, t2, prec);
        acb_mul(F2, F2, g, prec);
        acb_add(F, F, F2, prec);
        // * G(b+c) G(2s+1) 2^{1-b-c}
        acb_gamma(g, bc, prec);
        acb_mul(F, F, g, prec);
        acb_gamma(g, s2, prec);
        acb_mul(F, F, g, prec);
        acb_neg(t1, bc);
        acb_add_si(t1, t1, 1, prec);           // 1-b-c
        acb_set_si(t2, 2);
        acb_pow(t2, t2, t1, prec);
        acb_mul(F, F, t2, prec);
        acb_clear(s2); acb_clear(t1); acb_clear(t2); acb_clear(g);
        acb_clear(F2); acb_clear(one_m);
    } else {
        char* zs_raw = fmpq_get_str(nullptr, 10, zq);
        std::string zs(zs_raw);
        flint_free(zs_raw);
        _acb_vec_clear(al, 3); _acb_vec_clear(rh, 2);
        acb_clear(z); acb_clear(bc); acb_clear(F); acb_clear(Tv);
        throw std::runtime_error(
            "eikonal_ending_value: two-velocity ending of family '"
            + fam + "' at unsupported kinematics (3F2 argument z = "
            + zs + " <= -1 needs the |y| >= 3 continuation — not wired, "
            "stage-2 ending scope)");
    }
    eik_T_acb(Tv, a, bc, eps, prec);
    acb_mul(res, Tv, F, prec);
    _acb_vec_clear(al, 3); _acb_vec_clear(rh, 2);
    acb_clear(z); acb_clear(bc); acb_clear(F); acb_clear(Tv);
}

// Compute the value of an ending master on an eikonal-vacuum family as a
// product over disconnected components.  Stage-1 shapes (one velocity per
// component) keep the stage-1 formulas; stage-2 adds V1L(y) (1-loop,
// two velocities) and the V2L tiers B and C (2-loop components).  Anything
// else throws LOUDLY (no silent wrong physics).  `props` is the propagator
// list in a chosen loop-momentum routing (see the unimodular-retry wrapper
// below).
void eikonal_ending_value_props(const qft::FamilyConfig& fc,
                                const std::vector<algebra::Mpoly>& props,
                                const qft::JIntegral& j,
                                const AMFSystemOptions& opts,
                                const numeric::AcbValue& eps,
                                numeric::AcbValue& out,
                                long prec) {
    const auto& idx = j.indices();
    if (idx.size() != props.size()) {
        throw std::runtime_error(
            "eikonal_ending_value: index/propagator length mismatch for '"
            + jintegral_key(j) + "'");
    }
    const long L = (long)fc.n_loops();
    const long n_red = (long)fc.n_red_legs();
    const long n_vars = fc.ctx->n_vars();

    struct PClass {
        long deg = 0;          // max total loop degree
        long loop = -1;        // single loop var (if exactly one)
        std::set<long> loops;
        std::set<long> legs;   // legs in the loop-linear monomials
        bool has_const = false;
        long c0 = 0;           // constant monomial coefficient
        long cq = 0;           // coefficient of the pure loop^2 monomial(s)
        bool cq_set = false;   // (stage-2) a loop^2 monomial was seen
        bool cq_conflict = false;  // two loop^2 monomials with different coeff
        long cx = 0;           // (stage-2) coefficient of a l_i*l_j
                               // cross monomial (crossed quadratic)
        long ce = 0;           // coefficient of the loop-leg monomial
        long lq = 0;           // (ending-class engine) coefficient of a pure leg^2
        long lq_leg = -1;      //   monomial (the q^2 of a (k±q)^2 quad)
        bool lq_multi = false; //   second/malformed leg-quad monomial
        bool bad = false;      // unrecognized monomial content
    };

    std::vector<PClass> cls(props.size());
    std::vector<unsigned long> exp((std::size_t)n_vars);
    fmpz_t coeff;
    fmpz_init(coeff);
    for (std::size_t i = 0; i < props.size(); ++i) {
        if (idx[i] == 0) continue;
        const auto& p = props[i];
        PClass& c = cls[i];
        long len = fmpz_mpoly_length(p.raw(), fc.ctx->raw());
        for (long t = 0; t < len; ++t) {
            fmpz_mpoly_get_term_exp_ui(exp.data(), p.raw(), t, fc.ctx->raw());
            fmpz_mpoly_get_term_coeff_fmpz(coeff, p.raw(), t, fc.ctx->raw());
            long d = 0;
            bool any_leg = false, any_other = false;
            for (long v = 0; v < n_vars; ++v) {
                if (exp[(std::size_t)v] == 0) continue;
                if (v < L) d += (long)exp[(std::size_t)v];
                else if (v < L + n_red) any_leg = true;
                else any_other = true;
            }
            if (d > c.deg) c.deg = d;
            for (long v = 0; v < L; ++v) {
                if (exp[(std::size_t)v] > 0) c.loops.insert(v);
            }
            if (d == 0 && !any_leg && !any_other) {
                c.has_const = true;
                c.c0 = fmpz_get_si(coeff);
            } else if (d == 2 && !any_leg && !any_other) {
                // pure loop-quadratic monomial (l^2 or l1.l2)
                bool sq = false;
                for (long v = 0; v < L; ++v) {
                    if (exp[(std::size_t)v] == 2) sq = true;
                }
                if (sq) {
                    long cv = fmpz_get_si(coeff);
                    if (c.cq_set && c.cq != cv) c.cq_conflict = true;
                    c.cq = cv;
                    c.cq_set = true;
                } else {
                    // cross monomial l_i * l_j (stage-2: crossed
                    // quadratic of a coupled 2-loop eikonal component)
                    c.cx = fmpz_get_si(coeff);
                }
            } else if (d == 1 && any_leg && !any_other) {
                c.ce = fmpz_get_si(coeff);
                for (long v = L; v < L + n_red; ++v) {
                    if (exp[(std::size_t)v] > 0) c.legs.insert(v);
                }
            } else if (d == 1 && !any_leg && !any_other) {
                // bare loop-linear monomial without a leg: unsupported
                c.bad = true;
            } else if (d == 0) {
                // leg/invariant-only monomial.  A pure leg^2 monomial is
                // the q^2 of a (k±q)^2 leg quadratic (ending class 2);
                // anything else stays unsupported.
                long lv = -1;
                bool ok_lq = !any_other;
                for (long v = L; ok_lq && v < L + n_red; ++v) {
                    if (exp[(std::size_t)v] == 0) continue;
                    if (exp[(std::size_t)v] == 2 && lv < 0) lv = v;
                    else ok_lq = false;
                }
                if (ok_lq && lv >= 0) {
                    if (c.lq_leg >= 0) c.lq_multi = true;
                    else {
                        c.lq = fmpz_get_si(coeff);
                        c.lq_leg = lv;
                    }
                } else {
                    c.bad = true;
                }
            } else {
                c.bad = true;
            }
        }
        if (c.loops.size() == 1) c.loop = *c.loops.begin();
        if (idx[i] < 0) {
            fmpz_clear(coeff);
            throw std::runtime_error(
                "eikonal_ending_value: master '" + jintegral_key(j)
                + "' carries a NUMERATOR index on propagator "
                + std::to_string(i) + " — u-projection numerator recursion "
                "is not implemented (stage-1 scope); this should not occur "
                "for Kira corner masters");
        }
        if (c.bad) {
            fmpz_clear(coeff);
            throw std::runtime_error(
                "eikonal_ending_value: propagator " + std::to_string(i)
                + " ('" + props[i].to_string() + "') of family '"
                + fc.family + "' has unsupported monomial content for an "
                "eikonal-vacuum ending");
        }
    }
    fmpz_clear(coeff);

    // Union-find over loops via positive-index props.
    std::vector<long> uf((std::size_t)L);
    for (long k = 0; k < L; ++k) uf[(std::size_t)k] = k;
    std::function<long(long)> find = [&](long a) -> long {
        while (uf[(std::size_t)a] != a) {
            uf[(std::size_t)a] = uf[(std::size_t)uf[(std::size_t)a]];
            a = uf[(std::size_t)a];
        }
        return a;
    };
    for (std::size_t i = 0; i < props.size(); ++i) {
        if (idx[i] <= 0) continue;
        long first = -1;
        for (long v : cls[i].loops) {
            if (first < 0) { first = v; continue; }
            long a = find(first), b = find(v);
            if (a != b) uf[(std::size_t)a] = b;
        }
    }

    // If a loop of the family is touched by NO positive-index propagator,
    // the integrand does not depend on it: scaleless, value 0.
    {
        std::vector<char> covered((std::size_t)L, 0);
        for (std::size_t i = 0; i < props.size(); ++i) {
            if (idx[i] <= 0) continue;
            for (long v : cls[i].loops) covered[(std::size_t)v] = 1;
        }
        for (long k = 0; k < L; ++k) {
            if (!covered[(std::size_t)k]) {
                acb_zero(out.raw());
                return;
            }
        }
    }

    // Per-component accumulation (stage 1 + stage 2).
    struct EikSlot {
        long leg = -1;         // velocity leg variable index
        int  r = 0;            // orientation: prop = c0 * (1 - r * 2 l.u)
        long b = 0;            // index sum
        bool bare = false;     // (two-velocity E-step / bare-eikonal fix) pure +-2 l.u,
                               // no residual constant
    };
    struct LoopData {
        long a = 0;            // own massless-quad index sum
        long am = 0;           // massive-quad index sum
        long mq_cq = 1, mq_c0 = 0;   // massive quad normalization
        long bq = 0;           // (ending-class engine) leg-quad (k±q)^2 index sum
        long q_leg = -1;       //   its leg variable
        std::vector<EikSlot> eiks;
    };
    struct Comp {
        std::map<long, LoopData> ld;   // loop var -> data
        long ax = 0;           // crossed-quad (k_i + sigma k_j)^2 index sum
        int  ax_sigma = 0;     // sign of the cross monomial
        bool ax_conflict = false;
    };
    std::map<long, Comp> comps;
    numeric::AcbValue const_factor;    // constants + |c|>1 normalizations
    acb_one(const_factor.raw());
    long norm_sign = 1;                // accumulated +-1 normalization

    auto mul_const_pow = [&](long base, long idxp) {
        if (base == 1) return;
        if (base == -1) {
            if (idxp % 2 != 0) norm_sign = -norm_sign;
            return;
        }
        acb_t bb;
        acb_init(bb);
        acb_set_si(bb, base);
        acb_pow_si(bb, bb, -(slong)idxp, prec);
        acb_mul(const_factor.raw(), const_factor.raw(), bb, prec);
        acb_clear(bb);
    };

    for (std::size_t i = 0; i < props.size(); ++i) {
        if (idx[i] <= 0) continue;
        const PClass& c = cls[i];
        if (c.deg == 0) {
            // constant propagator: factor c0^-idx
            if (c.c0 == 0) {
                throw std::runtime_error(
                    "eikonal_ending_value: zero constant propagator");
            }
            mul_const_pow(c.c0, idx[i]);
            continue;
        }
        if (c.deg == 1) {
            // eikonal propagator
            if (!c.has_const) {
                // BARE eikonal (pure +-2 l.u, no
                // residual constant): NOT an automatic zero — the
                // scaleless argument is valid only for a VACUUM component
                // (an automatic-zero early return here would silently
                // delete one loop's eikonal-region content — stock
                // upstream returns O*F instead of F^2, shown at 66d).  Classify as
                // a bare EikSlot; the component valuation decides
                // (0 / Fbare / loud throw).
                if (c.loops.size() != 1) {
                    throw std::runtime_error(
                        "eikonal_ending_value: bare eikonal propagator '"
                        + props[i].to_string() + "' spans more than one "
                        "loop in family '" + fc.family + "' — not wired "
                        "(two-velocity E-step scope)");
                }
                if (c.ce != 2 && c.ce != -2) {
                    throw std::runtime_error(
                        "eikonal_ending_value: bare eikonal propagator '"
                        + props[i].to_string() + "' does not have a +-2 "
                        "l.u monomial — not wired (two-velocity E-step scope)");
                }
                if (c.legs.size() != 1) {
                    throw std::runtime_error(
                        "eikonal_ending_value: bare eikonal propagator '"
                        + props[i].to_string() + "' has "
                        + std::to_string(c.legs.size())
                        + " velocity legs — not wired (two-velocity E-step scope)");
                }
                const int rb = (c.ce == 2) ? +1 : -1;
                LoopData& ldd = comps[find(c.loop)].ld[c.loop];
                bool merged = false;
                for (auto& es : ldd.eiks) {
                    if (es.bare && es.leg == *c.legs.begin()
                            && es.r == rb) {
                        es.b += idx[i];
                        merged = true;
                        break;
                    }
                }
                if (!merged) {
                    ldd.eiks.push_back(
                        EikSlot{*c.legs.begin(), rb, (long)idx[i], true});
                }
                continue;
            }
            if (c.loops.size() != 1) {
                throw std::runtime_error(
                    "eikonal_ending_value: eikonal propagator '"
                    + props[i].to_string() + "' spans more than one loop "
                    "in family '" + fc.family + "' — not wired (stage-3 "
                    "shape)");
            }
            if (c.ce != 2 && c.ce != -2) {
                throw std::runtime_error(
                    "eikonal_ending_value: eikonal propagator '"
                    + props[i].to_string() + "' does not have a +-2 l.u "
                    "monomial");
            }
            if (c.c0 != 1 && c.c0 != -1) {
                throw std::runtime_error(
                    "eikonal_ending_value: eikonal residual mass is not "
                    "+-1 in '" + props[i].to_string() + "'");
            }
            if (c.legs.size() != 1) {
                throw std::runtime_error(
                    "eikonal_ending_value: eikonal propagator '"
                    + props[i].to_string() + "' has "
                    + std::to_string(c.legs.size()) + " velocity legs");
            }
            // c0 (1 - r 2 l.u):  ce/c0 = -2r
            const int r = (c.ce * c.c0 == -2) ? +1 : -1;
            mul_const_pow(c.c0, idx[i]);
            LoopData& ldd = comps[find(c.loop)].ld[c.loop];
            bool merged = false;
            for (auto& e : ldd.eiks) {
                if (e.leg == *c.legs.begin() && e.r == r) {
                    e.b += idx[i];
                    merged = true;
                    break;
                }
            }
            if (!merged) {
                ldd.eiks.push_back(
                    EikSlot{*c.legs.begin(), r, (long)idx[i]});
            }
            continue;
        }
        // deg == 2: quadratic
        if (!c.cq_set || c.cq_conflict || (c.cq != 1 && c.cq != -1)) {
            throw std::runtime_error(
                "eikonal_ending_value: quadratic propagator '"
                + props[i].to_string() + "' has non-unit or inconsistent "
                "leading coefficients");
        }
        if (c.loops.size() == 1) {
            LoopData& ldd = comps[find(c.loop)].ld[c.loop];
            if (c.cx != 0) {
                throw std::runtime_error(
                    "eikonal_ending_value: single-loop quadratic with a "
                    "cross monomial in '" + props[i].to_string() + "'");
            }
            if (c.lq_leg >= 0 || !c.legs.empty()) {
                // (ending class 2) leg quadratic: must be a COMPLETE
                // square cq*(k ± q)^2 — exactly one leg, ce = ±2 cq on that
                // same leg, leg-quad coefficient == cq, no constant.
                if (c.lq_multi || c.has_const || c.lq_leg < 0
                        || c.legs.size() != 1
                        || *c.legs.begin() != c.lq_leg
                        || c.lq != c.cq
                        || (c.ce != 2 * c.cq && c.ce != -2 * c.cq)) {
                    throw std::runtime_error(
                        "eikonal_ending_value: quadratic propagator '"
                        + props[i].to_string() + "' of family '" + fc.family
                        + "' carries external-leg monomials but is not a "
                        "complete (k±q)^2 square — not wired (ending-class "
                        "scope)");
                }
                if (ldd.bq > 0) {
                    throw std::runtime_error(
                        "eikonal_ending_value: two distinct leg quadratics "
                        "in one component of family '" + fc.family
                        + "' — not wired (ending-class scope)");
                }
                ldd.bq += idx[i];
                ldd.q_leg = c.lq_leg;
                mul_const_pow(c.cq, idx[i]);
                // (l^2)^-b = (-1)^b (-l^2)^-b
                if (idx[i] % 2 != 0) norm_sign = -norm_sign;
            } else if (c.has_const) {
                if (ldd.am > 0) {
                    throw std::runtime_error(
                        "eikonal_ending_value: two distinct massive "
                        "quadratics in one component");
                }
                ldd.am += idx[i];
                ldd.mq_cq = c.cq;
                ldd.mq_c0 = c.c0;
                // normalization applied at valuation (stage-1 identical)
            } else {
                ldd.a += idx[i];
                mul_const_pow(c.cq, idx[i]);
                // (l^2)^-a = (-1)^a (-l^2)^-a
                if (idx[i] % 2 != 0) norm_sign = -norm_sign;
            }
        } else if (c.loops.size() == 2) {
            // crossed quadratic cq*(k_i + sigma k_j)^2
            if (c.lq_leg >= 0 || !c.legs.empty()) {
                throw std::runtime_error(
                    "eikonal_ending_value: two-loop quadratic '"
                    + props[i].to_string() + "' carries external-leg "
                    "monomials — not wired (ending-class scope)");
            }
            if (c.has_const) {
                throw std::runtime_error(
                    "eikonal_ending_value: massive crossed quadratic '"
                    + props[i].to_string() + "' — not wired");
            }
            if (c.cx != 2 * c.cq && c.cx != -2 * c.cq) {
                throw std::runtime_error(
                    "eikonal_ending_value: two-loop quadratic '"
                    + props[i].to_string() + "' is not of the form "
                    "(k_i +- k_j)^2");
            }
            const int sigma = (c.cx == 2 * c.cq) ? +1 : -1;
            Comp& cm = comps[find(*c.loops.begin())];
            for (long v : c.loops) cm.ld[v];   // register both loops
            if (cm.ax > 0 && cm.ax_sigma != sigma) cm.ax_conflict = true;
            cm.ax += idx[i];
            cm.ax_sigma = sigma;
            mul_const_pow(c.cq, idx[i]);
            if (idx[i] % 2 != 0) norm_sign = -norm_sign;
        } else {
            throw std::runtime_error(
                "eikonal_ending_value: quadratic propagator '"
                + props[i].to_string() + "' spans more than two loops — "
                "not wired (stage-3 shape)");
        }
    }

    // Velocity-normalization guard for legs entering two-velocity forms.
    auto check_unit_velocity = [&](long legv) {
        fmpq_t nq;
        fmpq_init(nq);
        eik_leg_dot(fc, opts, eps, prec, legv, legv, nq);
        const bool ok = fmpq_is_one(nq);
        fmpq_clear(nq);
        if (!ok) {
            throw std::runtime_error(
                "eikonal_ending_value: velocity leg '"
                + fc.reduced_legs[(std::size_t)(legv - L)]
                + "' of family '" + fc.family + "' is not unit-normalized "
                "(u^2 != 1) — two-velocity endings assume unit velocities");
        }
    };

    // One-loop eikonal-component value: T(a,b) (one slot) or V1L (two
    // slots).  `a_acb` may be generic complex (Tier-B bubble shift).
    auto value_1loop_eik = [&](const acb_t a_acb,
                               const std::vector<EikSlot>& eiks,
                               acb_t res) {
        if (eiks.size() == 1) {
            acb_t bb;
            acb_init(bb);
            acb_set_si(bb, eiks[0].b);
            eik_T_acb(res, a_acb, bb, eps, prec);
            acb_clear(bb);
            return;
        }
        if (eiks.size() != 2) {
            throw std::runtime_error(
                "eikonal_ending_value: more than two eikonal slots on one "
                "loop in family '" + fc.family + "' — not wired");
        }
        const EikSlot& e1 = eiks[0];
        const EikSlot& e2 = eiks[1];
        fmpq_t yq, zq;
        fmpq_init(yq);
        fmpq_init(zq);
        if (e1.leg == e2.leg) {
            fmpq_set_si(yq, 1, 1);
        } else {
            check_unit_velocity(e1.leg);
            check_unit_velocity(e2.leg);
            eik_leg_dot(fc, opts, eps, prec, e1.leg, e2.leg, yq);
        }
        // y_eff = r1 r2 y ; z = (1 - y_eff)/2
        if (e1.r * e2.r < 0) fmpq_neg(yq, yq);
        fmpq_set_si(zq, 1, 1);
        fmpq_sub(zq, zq, yq);
        fmpq_div_2exp(zq, zq, 1);
        acb_t bb, cc;
        acb_init(bb); acb_init(cc);
        acb_set_si(bb, e1.b);
        acb_set_si(cc, e2.b);
        // integer b,c: the z>1 continuation is parameter-degenerate;
        // eik_pair_value dispatches the wired 2F1 reductions there
        eik_pair_value(res, a_acb, bb, cc, zq, /*cont_ok=*/false,
                       e1.b, e2.b, eps, prec, fc.family);
        acb_clear(bb); acb_clear(cc);
        fmpq_clear(yq); fmpq_clear(zq);
    };

    // Assemble: product over components.
    acb_set(out.raw(), const_factor.raw());
    acb_t D_half, val, g1, rg, arg;
    acb_init(D_half); acb_init(val); acb_init(g1);
    acb_init(rg); acb_init(arg);
    // D/2 = 2 - eps
    acb_set_si(D_half, 2);
    acb_sub(D_half, D_half, eps.raw(), prec);

    for (auto& [root, cm] : comps) {
        (void)root;
        if (cm.ax_conflict) {
            throw std::runtime_error(
                "eikonal_ending_value: component with both (k_i+k_j)^2 and "
                "(k_i-k_j)^2 crossed quadratics in family '" + fc.family
                + "' — not wired");
        }
        if (cm.ld.size() == 1) {
            LoopData& ldd = cm.ld.begin()->second;
            long nbare = 0;
            for (const auto& es : ldd.eiks) {
                if (es.bare) ++nbare;
            }
            if (nbare > 0) {
                // ---- two-velocity E-step fix: component-aware bare-eik ----
                if ((long)ldd.eiks.size() != nbare) {
                    throw std::runtime_error(
                        "eikonal_ending_value: bare and residual-mass "
                        "eikonal slots mixed on one loop of family '"
                        + fc.family + "' — not wired (two-velocity E-step scope); "
                        "master " + jintegral_key(j));
                }
                if (ldd.am > 0) {
                    throw std::runtime_error(
                        "eikonal_ending_value: bare eikonal with a massive "
                        "quadratic in one component of family '" + fc.family
                        + "' — not wired (two-velocity E-step scope); master "
                        + jintegral_key(j));
                }
                if (ldd.bq == 0) {
                    // NO leg content: the genuine scaleless case (vacuum
                    // component) — component factor 0
                    if (std::getenv("AMFLOW_EIK_DEBUG")) {
                        std::cerr << "[eik-val] bare-scaleless component "
                                  << "(loop " << cm.ld.begin()->first
                                  << ") family=" << fc.family
                                  << " master=" << jintegral_key(j)
                                  << "\n";
                    }
                    acb_zero(out.raw());
                    break;
                }
                // bare vertex Fbare(a, b, c): own massless quad a,
                // complete (k±q)^2 leg quad b, bare eik c — same guards
                // as class 2 (q^2 = -1, q.u = 0, u^2 = 1)
                if (nbare != 1) {
                    throw std::runtime_error(
                        "eikonal_ending_value: more than one bare eikonal "
                        "slot on a leg-quad component of family '"
                        + fc.family + "' — not wired (two-velocity E-step scope); "
                        "master " + jintegral_key(j));
                }
                {
                    fmpq_t q2;
                    fmpq_init(q2);
                    eik_leg_dot(fc, opts, eps, prec, ldd.q_leg, ldd.q_leg,
                                q2);
                    const bool ok = (fmpq_cmp_si(q2, -1) == 0);
                    fmpq_clear(q2);
                    if (!ok) {
                        throw std::runtime_error(
                            "eikonal_ending_value: bare-eik leg quadratic "
                            "of family '" + fc.family + "' with q^2 != -1 "
                            "— rescale not wired (two-velocity E-step scope); "
                            "master " + jintegral_key(j));
                    }
                }
                {
                    fmpq_t qu;
                    fmpq_init(qu);
                    eik_leg_dot(fc, opts, eps, prec, ldd.q_leg,
                                ldd.eiks[0].leg, qu);
                    const bool ok = fmpq_is_zero(qu);
                    fmpq_clear(qu);
                    if (!ok) {
                        throw std::runtime_error(
                            "eikonal_ending_value: bare-eik component of "
                            "family '" + fc.family + "' with q.u != 0 — "
                            "not derived (two-velocity E-step scope); master "
                            + jintegral_key(j));
                    }
                }
                check_unit_velocity(ldd.eiks[0].leg);
                if (ldd.a == 1 && ldd.bq == 1 && ldd.eiks[0].b == 1) {
                    // STAGE 1: Fbare(1,1,1) = F(eps), the
                    // double-derived vertex, validated to 75 digits
                    eik_bare_F(val, eps, prec);
                    eik_val_dbg("bare-F", fc, j, val, eps);
                    acb_mul(out.raw(), out.raw(), val, prec);
                } else {
                    throw std::runtime_error(
                        "eikonal_ending_value: bare-eikonal vertex Fbare("
                        + std::to_string(ldd.a) + ","
                        + std::to_string(ldd.bq) + ","
                        + std::to_string(ldd.eiks[0].b) + ") of family '"
                        + fc.family + "' — only (1,1,1) is implemented; "
                        "the general closed form is not implemented "
                        "(two-velocity E-step scope); "
                        "master " + jintegral_key(j));
                }
                continue;
            }
            if (ldd.bq > 0) {
                // ---- ending class 2: q-quad slot ----
                //   (-k^2)^-a (-(k±q)^2)^-bq (1-2k.u)^-c, q.u = 0, q^2 = -1
                // (orientation of the eikonal and the ±q sign are k -> -k /
                //  q -> -q symmetries of the value).
                if (ldd.am > 0) {
                    throw std::runtime_error(
                        "eikonal_ending_value: leg quadratic and massive "
                        "quadratic share a component in family '" + fc.family
                        + "' — not wired (ending-class scope); master "
                        + jintegral_key(j));
                }
                if (ldd.eiks.size() > 1) {
                    throw std::runtime_error(
                        "eikonal_ending_value: leg quadratic with more than "
                        "one eikonal slot in family '" + fc.family
                        + "' — not derived (ending-class scope); master "
                        + jintegral_key(j));
                }
                {
                    fmpq_t q2;
                    fmpq_init(q2);
                    eik_leg_dot(fc, opts, eps, prec, ldd.q_leg, ldd.q_leg,
                                q2);
                    const bool ok = (fmpq_cmp_si(q2, -1) == 0);
                    fmpq_clear(q2);
                    if (!ok) {
                        throw std::runtime_error(
                            "eikonal_ending_value: leg quadratic of family '"
                            + fc.family + "' with q^2 != -1 — the q2 rescale "
                            "is not wired (ending-class scope, "
                            "spec A); master "
                            + jintegral_key(j));
                    }
                }
                long cidx = 0;
                if (!ldd.eiks.empty()) {
                    cidx = ldd.eiks[0].b;
                    fmpq_t qu;
                    fmpq_init(qu);
                    eik_leg_dot(fc, opts, eps, prec, ldd.q_leg,
                                ldd.eiks[0].leg, qu);
                    const bool ok = fmpq_is_zero(qu);
                    fmpq_clear(qu);
                    if (!ok) {
                        throw std::runtime_error(
                            "eikonal_ending_value: leg quadratic of family '"
                            + fc.family + "' with q.u != 0 — not derived "
                            "(ending-class scope); master " + jintegral_key(j));
                    }
                    check_unit_velocity(ldd.eiks[0].leg);
                }
                if (cidx == 0 && ldd.a == 0) {
                    // (-(k±q)^2)^-b alone: shift k -> k∓q, scaleless
                    acb_zero(out.raw());
                    break;
                }
                if (cidx == 0) {
                    // massless bubble at q^2 = -1 (unit |q^2|: no rescale)
                    eik_bubble(val, ldd.a, ldd.bq, eps, prec);
                    eik_val_dbg("class2-Gb", fc, j, val, eps);
                } else if (ldd.a == 0) {
                    // shift k -> k±q (eikonal invariant, q.u = 0): T(bq,c)
                    acb_t aa, bb;
                    acb_init(aa); acb_init(bb);
                    acb_set_si(aa, ldd.bq);
                    acb_set_si(bb, cidx);
                    eik_T_acb(val, aa, bb, eps, prec);
                    acb_clear(aa); acb_clear(bb);
                    eik_val_dbg("class2-T", fc, j, val, eps);
                } else {
                    eik_class2_barnes(val, ldd.a, ldd.bq, cidx, eps, prec);
                    eik_val_dbg("class2-C2", fc, j, val, eps);
                }
                acb_mul(out.raw(), out.raw(), val, prec);
            } else if (!ldd.eiks.empty()) {
                if (ldd.am > 0) {
                    throw std::runtime_error(
                        "eikonal_ending_value: massive quadratic and "
                        "eikonal propagator share a component — "
                        "unsupported shape");
                }
                acb_set_si(arg, ldd.a);
                value_1loop_eik(arg, ldd.eiks, val);
                acb_mul(out.raw(), out.raw(), val, prec);
            } else if (ldd.am > 0) {
                if (ldd.a > 0) {
                    throw std::runtime_error(
                        "eikonal_ending_value: massless and massive "
                        "quadratics share a component — unsupported master "
                        "shape (should have been reduced)");
                }
                // m^2 = -c0/cq must be +1
                if (ldd.mq_c0 * ldd.mq_cq != -1) {
                    throw std::runtime_error(
                        "eikonal_ending_value: massive quadratic with "
                        "m^2 != 1");
                }
                const long a = ldd.am;
                // cq^-a (-1)^a Gamma(a-2+eps) rgamma(a)
                acb_set_si(arg, a - 2);
                acb_add(arg, arg, eps.raw(), prec);
                acb_gamma(g1, arg, prec);
                acb_set_si(arg, a);
                acb_rgamma(rg, arg, prec);
                acb_mul(val, g1, rg, prec);
                long sign = 1;
                if (ldd.mq_cq == -1 && (a % 2 != 0)) sign = -sign;
                if (a % 2 != 0) sign = -sign;
                if (sign < 0) acb_neg(val, val);
                acb_mul(out.raw(), out.raw(), val, prec);
            } else {
                // massless-only component: scaleless
                acb_zero(out.raw());
                break;
            }
        } else if (cm.ld.size() == 2) {
            if (cm.ax == 0) {
                throw std::runtime_error(
                    "eikonal_ending_value: internal error — two-loop "
                    "component without a crossed quadratic");
            }
            auto it0 = cm.ld.begin();
            auto it1 = std::next(it0);
            LoopData* A = &it0->second;
            LoopData* B = &it1->second;
            for (const LoopData* ldp : {(const LoopData*)A,
                                        (const LoopData*)B}) {
                for (const auto& es : ldp->eiks) {
                    if (es.bare) {
                        throw std::runtime_error(
                            "eikonal_ending_value: bare eikonal inside a "
                            "coupled 2-loop component of family '"
                            + fc.family + "' — not wired (two-velocity E-step "
                            "scope); master "
                            + jintegral_key(j));
                    }
                }
            }
            if (A->am > 0 || B->am > 0) {
                throw std::runtime_error(
                    "eikonal_ending_value: massive quadratic inside a "
                    "coupled 2-loop eikonal component — not wired");
            }
            if (A->bq > 0 || B->bq > 0) {
                throw std::runtime_error(
                    "eikonal_ending_value: leg quadratic inside a coupled "
                    "2-loop eikonal component in family '" + fc.family
                    + "' — not wired (ending-class scope); master "
                    + jintegral_key(j));
            }
            if (A->eiks.empty() || B->eiks.empty()) {
                // ---- Tier B: one loop purely quadratic ----
                if (A->eiks.empty() && B->eiks.empty()) {
                    // massless 2-loop vacuum: scaleless
                    acb_zero(out.raw());
                    break;
                }
                if (A->eiks.empty()) std::swap(A, B);  // A carries eikonals
                if (B->a == 0) {
                    // inner loop = crossed quad only: massless tadpole
                    acb_zero(out.raw());
                    break;
                }
                // Gb(B.a, ax) * V1L(A.a + B.a + ax - D/2, ...)
                eik_bubble(g1, B->a, cm.ax, eps, prec);
                acb_set_si(arg, A->a + B->a + cm.ax);
                acb_sub(arg, arg, D_half, prec);
                value_1loop_eik(arg, A->eiks, val);
                acb_mul(val, val, g1, prec);
                acb_mul(out.raw(), out.raw(), val, prec);
            } else {
                // ---- Tier C: inner loop = crossed quad + one eikonal ----
                LoopData* inner = nullptr;
                LoopData* outer = nullptr;
                if (A->a == 0 && B->a == 0) {
                    // after the inner integration the outer loop has no
                    // quadratic left: eikonal-only, scaleless
                    acb_zero(out.raw());
                    break;
                } else if (A->a == 0) {
                    inner = A; outer = B;
                } else if (B->a == 0) {
                    inner = B; outer = A;
                } else {
                    throw std::runtime_error(
                        "eikonal_ending_value: coupled 2-loop eikonal "
                        "component with own quadratics AND eikonals on "
                        "both loops (V2L Tier D) in family '" + fc.family
                        + "' — not wired (stage-2 ending scope); master "
                        + jintegral_key(j));
                }
                if (inner->eiks.size() == 2) {
                    // ---- two-velocity E-step class 3: two-velocity E-step ----
                    // canonical (after the GL(2,Z) wrapper):
                    //   TC3(a1,a2,b,c,e) = IntInt Dk1 Dk2 /
                    //     [(-(k1+k2)^2)^{a1} (-k2^2)^{a2} (1-2k1.u1)^{b}
                    //      (1-2k1.u2)^{c} (1-2k2.u1)^{e}],
                    // assembled value = TC3 (all
                    // Gammas inside pref3; NO separate T-factor).
                    const EikSlot& i0 = inner->eiks[0];
                    const EikSlot& i1 = inner->eiks[1];
                    if (i0.leg == i1.leg) {
                        throw std::runtime_error(
                            "eikonal_ending_value: two-velocity E-step "
                            "loop of family '" + fc.family + "' with a "
                            "SAME-velocity slot pair — not in the class-3 "
                            "inventory (two-velocity E-step scope); master "
                            + jintegral_key(j));
                    }
                    if (i0.r != i1.r) {
                        throw std::runtime_error(
                            "eikonal_ending_value: two-velocity E-step "
                            "loop of family '" + fc.family + "' with "
                            "OPPOSITE slot orientations (u_lam^2 can "
                            "vanish; lam-mix not valid) — not wired "
                            "(two-velocity E-step scope); master "
                            + jintegral_key(j));
                    }
                    if (outer->eiks.size() != 1) {
                        throw std::runtime_error(
                            "eikonal_ending_value: class-3 component of "
                            "family '" + fc.family + "' with "
                            + std::to_string(outer->eiks.size())
                            + " eikonal slots on the quad-carrying loop — "
                            "outside the implemented class list "
                            "(two-velocity scope); master " + jintegral_key(j));
                    }
                    const EikSlot& eo = outer->eiks[0];
                    // outer velocity must match ONE E-step velocity
                    // (u1 := the matched one; mirror swap b <-> c)
                    const EikSlot* sb = nullptr;   // b-slot (u1, matched)
                    const EikSlot* sc = nullptr;   // c-slot (u2)
                    if (eo.leg == i0.leg) { sb = &i0; sc = &i1; }
                    else if (eo.leg == i1.leg) { sb = &i1; sc = &i0; }
                    else {
                        throw std::runtime_error(
                            "eikonal_ending_value: class-3 component of "
                            "family '" + fc.family + "' with THREE "
                            "distinct velocities — outside the implemented class list "
                            "(two-velocity scope); master "
                            + jintegral_key(j));
                    }
                    // orientation: after k1 -> r_in k1, k2 -> r_e k2 the
                    // crossed quad is (k1 + r_in r_e sigma k2)^2; the
                    // canonical class needs +1 (other signs: let the
                    // GL(2,Z) wrapper retry / stay loud)
                    if (i0.r * eo.r * cm.ax_sigma != 1) {
                        throw std::runtime_error(
                            "eikonal_ending_value: class-3 component of "
                            "family '" + fc.family + "' in a non-canonical "
                            "crossed-quad orientation — not wired "
                            "(two-velocity E-step scope); master "
                            + jintegral_key(j));
                    }
                    check_unit_velocity(i0.leg);
                    check_unit_velocity(i1.leg);
                    fmpq_t yq3;
                    fmpq_init(yq3);
                    eik_leg_dot(fc, opts, eps, prec, i0.leg, i1.leg, yq3);
                    if (fmpq_cmp_si(yq3, 1) <= 0) {
                        fmpq_clear(yq3);
                        throw std::runtime_error(
                            "eikonal_ending_value: class-3 component of "
                            "family '" + fc.family + "' at y <= 1 — "
                            "pinch/continuation not wired (two-velocity "
                            "E-step scope); master " + jintegral_key(j));
                    }
                    eik_tc3_R3(val, cm.ax, outer->a, sb->b, sc->b, eo.b,
                               yq3, eps, prec, fc.family);
                    fmpq_clear(yq3);
                    eik_val_dbg("class3-R3", fc, j, val, eps);
                    acb_mul(out.raw(), out.raw(), val, prec);
                    continue;
                }
                if (inner->eiks.size() != 1) {
                    throw std::runtime_error(
                        "eikonal_ending_value: coupled 2-loop eikonal "
                        "component with more than two eikonal slots on the "
                        "E-step (quad-free) loop in family '" + fc.family
                        + "' — not wired (two-velocity E-step scope); master "
                        + jintegral_key(j));
                }
                if (outer->eiks.size() > 2) {
                    throw std::runtime_error(
                        "eikonal_ending_value: coupled 2-loop eikonal "
                        "component with more than two eikonal slots on the "
                        "outer loop in family '" + fc.family
                        + "' — not wired (ending-class scope); master "
                        + jintegral_key(j));
                }
                if (outer->eiks.size() == 2) {
                    // ---- ending class 1: nested cusp ----
                    // E-step continuation on the inner
                    // loop:  Int d^Dk_A (-(k_A + sigma k_B)^2)^{-ax}
                    //          (1 - r_i 2k_A.v)^{-b}
                    //   = T(ax, b) * ((1 + rho 2k_B.v) - i0)^{-bh},
                    //   rho = r_i * sigma,  bh = 2 ax + b - D (never an
                    //   integer).  In the (1 - r 2k.u) slot basis the new
                    //   outer-loop slot has orientation r_bh = -rho.
                    const EikSlot& ein = inner->eiks[0];
                    const int rho = ein.r * cm.ax_sigma;
                    const int r_bh = -rho;
                    acb_t bh, TE, aouter, bb1, cc1;
                    acb_init(bh); acb_init(TE); acb_init(aouter);
                    acb_init(bb1); acb_init(cc1);
                    acb_set_si(bh, 2 * cm.ax + ein.b);
                    acb_submul_si(bh, D_half, 2, prec);
                    {
                        acb_t axq, b2;
                        acb_init(axq); acb_init(b2);
                        acb_set_si(axq, cm.ax);
                        acb_set_si(b2, ein.b);
                        eik_T_acb(TE, axq, b2, eps, prec);
                        acb_clear(axq); acb_clear(b2);
                    }
                    acb_set_si(aouter, outer->a);
                    const EikSlot& s0 = outer->eiks[0];
                    const EikSlot& s1 = outer->eiks[1];
                    if (s0.leg == s1.leg) {
                        // outer carries a same-velocity mixed pair itself
                        if (ein.leg != s0.leg) {
                            acb_clear(bh); acb_clear(TE); acb_clear(aouter);
                            acb_clear(bb1); acb_clear(cc1);
                            throw std::runtime_error(
                                "eikonal_ending_value: nested cusp with an "
                                "integer same-velocity pair plus a distinct "
                                "E-step velocity in family '" + fc.family
                                + "' — tau-degenerate route R, not wired "
                                "(ending-class scope); master "
                                + jintegral_key(j));
                        }
                        // all three slots on ONE velocity: merge bh into
                        // the matching-orientation slot -> the z = 1
                        // same-velocity pair (BG-kernel class,
                        // generic non-integer index)
                        check_unit_velocity(s0.leg);
                        acb_set_si(bb1, s0.b);
                        acb_set_si(cc1, s1.b);
                        if (r_bh == s0.r) acb_add(bb1, bb1, bh, prec);
                        else acb_add(cc1, cc1, bh, prec);
                        fmpq_t zq1;
                        fmpq_init(zq1);
                        fmpq_set_si(zq1, 1, 1);
                        eik_pair_value(val, aouter, bb1, cc1, zq1,
                                       /*cont_ok=*/true, -1, -1, eps, prec,
                                       fc.family);
                        fmpq_clear(zq1);
                    } else {
                        const EikSlot* M = nullptr;   // matches E-step leg
                        const EikSlot* Uo = nullptr;  // the odd velocity
                        if (ein.leg == s0.leg) { M = &s0; Uo = &s1; }
                        else if (ein.leg == s1.leg) { M = &s1; Uo = &s0; }
                        else {
                            acb_clear(bh); acb_clear(TE); acb_clear(aouter);
                            acb_clear(bb1); acb_clear(cc1);
                            throw std::runtime_error(
                                "eikonal_ending_value: nested cusp with "
                                "three distinct velocities in family '"
                                + fc.family + "' — not wired (ending-class "
                                "scope); master " + jintegral_key(j));
                        }
                        check_unit_velocity(M->leg);
                        check_unit_velocity(Uo->leg);
                        fmpq_t yq2;
                        fmpq_init(yq2);
                        eik_leg_dot(fc, opts, eps, prec, Uo->leg, M->leg,
                                    yq2);
                        if (Uo->r * M->r < 0) fmpq_neg(yq2, yq2); // y_eff
                        if (r_bh == M->r) {
                            // bh joins the matched slot with the SAME
                            // orientation: two-velocity V1L with a
                            // non-integer second index (the
                            // Tier-C 3F2 class), z = (1 - y_eff)/2
                            fmpq_t zq2;
                            fmpq_init(zq2);
                            fmpq_set_si(zq2, 1, 1);
                            fmpq_sub(zq2, zq2, yq2);
                            fmpq_div_2exp(zq2, zq2, 1);
                            acb_set_si(bb1, Uo->b);
                            acb_add_si(cc1, bh, M->b, prec);
                            eik_pair_value(val, aouter, bb1, cc1, zq2,
                                           /*cont_ok=*/true, -1, -1, eps,
                                           prec, fc.family);
                            fmpq_clear(zq2);
                        } else {
                            // OPPOSITE orientation on the matched leg: the
                            // V1L3 nested-cusp class (route R), y_eff > 1
                            eik_v1l3_R(val, outer->a, Uo->b, M->b, bh, yq2,
                                       eps, prec, fc.family);
                        }
                        fmpq_clear(yq2);
                    }
                    acb_mul(val, val, TE, prec);
                    eik_val_dbg("class1-nested", fc, j, val, eps);
                    acb_mul(out.raw(), out.raw(), val, prec);
                    acb_clear(bh); acb_clear(TE); acb_clear(aouter);
                    acb_clear(bb1); acb_clear(cc1);
                    continue;
                }
                const EikSlot& ei = inner->eiks[0];
                const EikSlot& eo = outer->eiks[0];
                const int sig_eff = eo.r * cm.ax_sigma * ei.r;
                fmpq_t yq, zq;
                fmpq_init(yq);
                fmpq_init(zq);
                if (eo.leg == ei.leg) {
                    fmpq_set_si(yq, 1, 1);
                    check_unit_velocity(eo.leg);
                } else {
                    check_unit_velocity(eo.leg);
                    check_unit_velocity(ei.leg);
                    eik_leg_dot(fc, opts, eps, prec, eo.leg, ei.leg, yq);
                }
                // z = (1 + sig_eff * y)/2
                if (sig_eff < 0) fmpq_neg(yq, yq);
                fmpq_set_si(zq, 1, 1);
                fmpq_add(zq, zq, yq);
                fmpq_div_2exp(zq, zq, 1);
                // bh = 2 ax + b2 - D
                acb_t bh, a1, b1, axq, b2;
                acb_init(bh); acb_init(a1); acb_init(b1);
                acb_init(axq); acb_init(b2);
                acb_set_si(bh, 2 * cm.ax + ei.b);
                acb_submul_si(bh, D_half, 2, prec);
                acb_set_si(a1, outer->a);
                acb_set_si(b1, eo.b);
                acb_set_si(axq, cm.ax);
                acb_set_si(b2, ei.b);
                // T(ax, b2) * T(a1, b1+bh) * 3F2(...; z)
                eik_T_acb(g1, axq, b2, eps, prec);
                eik_pair_value(val, a1, b1, bh, zq, /*cont_ok=*/true,
                               -1, -1, eps, prec, fc.family);
                acb_mul(val, val, g1, prec);
                acb_mul(out.raw(), out.raw(), val, prec);
                acb_clear(bh); acb_clear(a1); acb_clear(b1);
                acb_clear(axq); acb_clear(b2);
                fmpq_clear(yq); fmpq_clear(zq);
            }
        } else {
            throw std::runtime_error(
                "eikonal_ending_value: eikonal component with more than "
                "two loops in family '" + fc.family + "' — not wired "
                "(stage-3 scope); master " + jintegral_key(j));
        }
    }
    if (norm_sign < 0) acb_neg(out.raw(), out.raw());
    acb_clear(D_half); acb_clear(val); acb_clear(g1);
    acb_clear(rg); acb_clear(arg);
}

// Unimodular-retry wrapper: some boundary families carry difference-momentum
// eikonals (e.g. 2(k1-k2).u - 1) or crossed quads in a routing the shape
// classifier does not recognize.  A unimodular loop-momentum redefinition
// leaves the integral value invariant (|det| = 1, measure preserved), so on
// a shape throw we retry the valuation under the GL(2,Z) transforms with
// entries in {-1,0,1}.  If none classifies, the ORIGINAL error is rethrown
// (still loud).
void eikonal_ending_value(const qft::FamilyConfig& fc,
                          const qft::JIntegral& j,
                          const AMFSystemOptions& opts,
                          const numeric::AcbValue& eps,
                          numeric::AcbValue& out,
                          long prec) {
    try {
        eikonal_ending_value_props(fc, fc.propagators_after_conservation,
                                   j, opts, eps, out, prec);
        eik_val_dbg("ending", fc, j, out.raw(), eps);
        return;
    } catch (const std::runtime_error& e0) {
        const long L = (long)fc.n_loops();
        if (L != 2) throw;
        const long n_vars = fc.ctx->n_vars();
        // substitution polynomials: identity for non-loop vars
        std::vector<algebra::Mpoly> gens;
        gens.reserve((std::size_t)n_vars);
        for (long v = 0; v < n_vars; ++v) {
            gens.push_back(algebra::Mpoly::variable(fc.ctx, v));
        }
        for (int m00 = -1; m00 <= 1; ++m00)
        for (int m01 = -1; m01 <= 1; ++m01)
        for (int m10 = -1; m10 <= 1; ++m10)
        for (int m11 = -1; m11 <= 1; ++m11) {
            const int det = m00*m11 - m01*m10;
            if (det != 1 && det != -1) continue;
            if (m00 == 1 && m01 == 0 && m10 == 0 && m11 == 1) continue;
            // k1 -> m00 k1 + m01 k2 ; k2 -> m10 k1 + m11 k2
            algebra::Mpoly s0 = algebra::Mpoly::zero(fc.ctx);
            algebra::Mpoly s1 = algebra::Mpoly::zero(fc.ctx);
            if (m00 != 0) s0 += algebra::Mpoly::constant(fc.ctx, m00)*gens[0];
            if (m01 != 0) s0 += algebra::Mpoly::constant(fc.ctx, m01)*gens[1];
            if (m10 != 0) s1 += algebra::Mpoly::constant(fc.ctx, m10)*gens[0];
            if (m11 != 0) s1 += algebra::Mpoly::constant(fc.ctx, m11)*gens[1];
            std::vector<const fmpz_mpoly_struct*> C;
            C.reserve((std::size_t)n_vars);
            C.push_back(s0.raw());
            C.push_back(s1.raw());
            for (long v = 2; v < n_vars; ++v) C.push_back(gens[(std::size_t)v].raw());
            std::vector<algebra::Mpoly> tprops;
            tprops.reserve(fc.propagators_after_conservation.size());
            bool comp_ok = true;
            for (const auto& p : fc.propagators_after_conservation) {
                algebra::Mpoly tp = algebra::Mpoly::zero(fc.ctx);
                if (!fmpz_mpoly_compose_fmpz_mpoly(
                        tp.raw(), p.raw(),
                        const_cast<fmpz_mpoly_struct* const*>(C.data()),
                        fc.ctx->raw(), fc.ctx->raw())) {
                    comp_ok = false;
                    break;
                }
                tprops.push_back(std::move(tp));
            }
            if (!comp_ok) continue;
            try {
                eikonal_ending_value_props(fc, tprops, j, opts, eps, out,
                                           prec);
                eik_val_dbg("ending-gl2z", fc, j, out.raw(), eps);
                return;
            } catch (const std::runtime_error& et) {
                if (std::getenv("AMFLOW_EIK_DEBUG")) {
                    std::cerr << "[eik-debug] M=(" << m00 << "," << m01
                              << ";" << m10 << "," << m11 << "): "
                              << et.what() << "\n";
                }
            }
        }
        if (std::getenv("AMFLOW_EIK_DEBUG")) {
            std::cerr << "[eik-debug] family '" << fc.family
                      << "' master " << jintegral_key(j) << " props:\n";
            const auto& idx = j.indices();
            for (std::size_t i = 0;
                 i < fc.propagators_after_conservation.size(); ++i) {
                std::cerr << "  [" << i << "] idx="
                          << (i < idx.size() ? idx[i] : 0) << "  "
                          << fc.propagators_after_conservation[i].to_string()
                          << "\n";
            }
        }
        throw;   // rethrow the identity-routing error (loud)
    }
}

}  // namespace

// ---------------------------------------------------------------------------
// Ending-class validation hooks: direct evaluation of the two ending classes
// against the 40/45-digit reference tables.  Exported
// with C linkage for a standalone validation driver (links libamflow.a);
// no pipeline path calls these — CLI behavior is unchanged.
// ---------------------------------------------------------------------------
extern "C" {

// C2(a,b,c) at eps = eps_re + i eps_im (decimal strings).  Writes
// "re_str im_str" (arb_get_str at `digits`, radius-aware) into buf.
void amflow_eik_gate_c2(long a, long b, long c,
                        const char* eps_re, const char* eps_im,
                        long prec, long digits, char* buf, long buflen) {
    numeric::AcbValue eps;
    arb_set_str(acb_realref(eps.raw()), eps_re, prec);
    arb_set_str(acb_imagref(eps.raw()), eps_im, prec);
    acb_t v;
    acb_init(v);
    eik_class2_barnes(v, a, b, c, eps, prec);
    char* rs = arb_get_str(acb_realref(v), digits, 0);
    char* is = arb_get_str(acb_imagref(v), digits, 0);
    std::snprintf(buf, (std::size_t)buflen, "%s %s", rs, is);
    flint_free(rs);
    flint_free(is);
    acb_clear(v);
}

// Assembled class-1 master E = T(a1,b) * V1L3(a2; c, ch, bh; y),
// bh = 2 a1 + b - D, y = ynum/yden.  c == 0 routes through the SAME code
// the engine uses (the z=1 pair value with a non-integer second index —
// the previously ungated path the spec flags); c >= 1 through route R.
void amflow_eik_gate_e(long a1, long a2, long b, long c, long ch,
                       long ynum, long yden,
                       const char* eps_re, const char* eps_im,
                       long prec, long digits, char* buf, long buflen) {
    numeric::AcbValue eps;
    arb_set_str(acb_realref(eps.raw()), eps_re, prec);
    arb_set_str(acb_imagref(eps.raw()), eps_im, prec);
    acb_t v, bh, Dh, T1, aa, bb;
    acb_init(v); acb_init(bh); acb_init(Dh); acb_init(T1);
    acb_init(aa); acb_init(bb);
    acb_set_si(Dh, 2);
    acb_sub(Dh, Dh, eps.raw(), prec);
    acb_set_si(bh, 2 * a1 + b);
    acb_submul_si(bh, Dh, 2, prec);              // bh = 2a1+b-D
    if (c == 0) {
        // engine path: T(a2, ch+bh) * F(z=1) == V1L3_M closed continuation
        fmpq_t zq;
        fmpq_init(zq);
        fmpq_set_si(zq, 1, 1);
        acb_t a2b, chb;
        acb_init(a2b); acb_init(chb);
        acb_set_si(a2b, a2);
        acb_set_si(chb, ch);
        eik_pair_value(v, a2b, chb, bh, zq, /*cont_ok=*/true, -1, -1,
                       eps, prec, "g1_gate");
        acb_clear(a2b); acb_clear(chb);
        fmpq_clear(zq);
    } else {
        fmpq_t yq;
        fmpq_init(yq);
        fmpq_set_si(yq, ynum, yden);
        eik_v1l3_R(v, a2, c, ch, bh, yq, eps, prec, "g1_gate");
        fmpq_clear(yq);
    }
    acb_set_si(aa, a1);
    acb_set_si(bb, b);
    eik_T_acb(T1, aa, bb, eps, prec);
    acb_mul(v, v, T1, prec);
    char* rs = arb_get_str(acb_realref(v), digits, 0);
    char* is = arb_get_str(acb_imagref(v), digits, 0);
    std::snprintf(buf, (std::size_t)buflen, "%s %s", rs, is);
    flint_free(rs);
    flint_free(is);
    acb_clear(v); acb_clear(bh); acb_clear(Dh); acb_clear(T1);
    acb_clear(aa); acb_clear(bb);
}

// Two-velocity E-step gate hooks — class 3 + bare-eikonal fix.

// TC3 gate: mode 0 = assembled TC3 via route R3 (the engine path);
// mode 1 = y=1 Dirichlet-type closed identity
//          T(a1,b+c) G(D/2-a2)G(g3)/(G(a2)G(e)G(bh3)) * bg_closed(e,bh3,s2)
//          (bg_closed = the z=1 BG kernel);
// mode 2 = e=0 ENGINE composite (Tier-B path: Gb(a2,a1) *
//          pair_value(a1+a2-D/2; b,c; z=(1-y)/2)) — e ignored;
// mode 3 = e=0 closed 1-fold (C.4) — e ignored.
void amflow_eik_gate_tc3(long a1, long a2, long b, long c, long e,
                         long ynum, long yden, long mode,
                         const char* eps_re, const char* eps_im,
                         long prec, long digits, char* buf, long buflen) {
    numeric::AcbValue eps;
    arb_set_str(acb_realref(eps.raw()), eps_re, prec);
    arb_set_str(acb_imagref(eps.raw()), eps_im, prec);
    fmpq_t yq;
    fmpq_init(yq);
    fmpq_set_si(yq, ynum, yden);
    acb_t v, Dh, t1, t2, g;
    acb_init(v); acb_init(Dh); acb_init(t1); acb_init(t2); acb_init(g);
    acb_set_si(Dh, 2);
    acb_sub(Dh, Dh, eps.raw(), prec);
    try {
        if (mode == 0) {
            eik_tc3_R3(v, a1, a2, b, c, e, yq, eps, prec, "g1_gate_tc3");
        } else if (mode == 1) {
            // closed y=1 form
            acb_t s2, bh3, g3, bg, T1, bb, cc;
            acb_init(s2); acb_init(bh3); acb_init(g3); acb_init(bg);
            acb_init(T1); acb_init(bb); acb_init(cc);
            acb_set_si(s2, a2);
            acb_sub(s2, s2, Dh, prec);
            acb_set_si(bh3, 2 * a1 + b + c);
            acb_submul_si(bh3, Dh, 2, prec);
            acb_set_si(g3, 2 * a1 + 2 * a2 + b + c + e);
            acb_submul_si(g3, Dh, 4, prec);
            // bg_closed(p1,p2,s) = 2^{1-p1-p2} [B(2s+1,p1) 2F1(1-p2,2s+1;
            //   2s+1+p1;-1) + B(2s+1,p2) 2F1(1-p1,2s+1;2s+1+p2;-1)]
            acb_t twos1, p1v, p2v, f, term, beta, mone;
            acb_init(twos1); acb_init(p1v); acb_init(p2v); acb_init(f);
            acb_init(term); acb_init(beta); acb_init(mone);
            acb_set_si(mone, -1);
            acb_mul_2exp_si(twos1, s2, 1);
            acb_add_si(twos1, twos1, 1, prec);        // 2s+1
            acb_set_si(p1v, e);
            acb_set(p2v, bh3);
            acb_zero(bg);
            for (int piece = 0; piece < 2; ++piece) {
                const acb_struct* pa = (piece == 0) ? p1v : p2v;
                const acb_struct* pb = (piece == 0) ? p2v : p1v;
                // B(2s+1, pa)
                acb_gamma(beta, twos1, prec);
                acb_gamma(t1, pa, prec);
                acb_mul(beta, beta, t1, prec);
                acb_add(t1, twos1, pa, prec);
                acb_rgamma(t2, t1, prec);
                acb_mul(beta, beta, t2, prec);
                // 2F1(1-pb, 2s+1; 2s+1+pa; -1)
                acb_neg(t2, pb);
                acb_add_si(t2, t2, 1, prec);
                acb_hypgeom_2f1(f, t2, twos1, t1, mone, 0, prec);
                acb_mul(term, beta, f, prec);
                acb_add(bg, bg, term, prec);
            }
            acb_set_si(t1, 1 - e);
            acb_sub(t1, t1, bh3, prec);               // 1-p1-p2
            acb_set_si(t2, 2);
            acb_pow(t2, t2, t1, prec);
            acb_mul(bg, bg, t2, prec);
            // T(a1, b+c)
            acb_set_si(t1, a1);
            acb_set_si(t2, b + c);
            eik_T_acb(T1, t1, t2, eps, prec);
            // G(D/2-a2) G(g3) rg(a2) rg(e) rg(bh3)
            acb_sub_si(t1, Dh, a2, prec);
            acb_gamma(v, t1, prec);
            acb_gamma(t2, g3, prec);
            acb_mul(v, v, t2, prec);
            acb_set_si(t1, a2);
            acb_rgamma(t2, t1, prec);
            acb_mul(v, v, t2, prec);
            acb_set_si(t1, e);
            acb_rgamma(t2, t1, prec);
            acb_mul(v, v, t2, prec);
            acb_rgamma(t2, bh3, prec);
            acb_mul(v, v, t2, prec);
            acb_mul(v, v, T1, prec);
            acb_mul(v, v, bg, prec);
            acb_clear(twos1); acb_clear(p1v); acb_clear(p2v); acb_clear(f);
            acb_clear(term); acb_clear(beta); acb_clear(mone);
            acb_clear(s2); acb_clear(bh3); acb_clear(g3); acb_clear(bg);
            acb_clear(T1); acb_clear(bb); acb_clear(cc);
        } else if (mode == 2) {
            // e=0 ENGINE composite: Gb(a2, a1) * pair_value(a1+a2-D/2;
            // b, c; z = (1-y)/2) — the Tier-B path the classifier routes
            // e=0 shapes through
            acb_t arg, bb, cc, Gb;
            acb_init(arg); acb_init(bb); acb_init(cc); acb_init(Gb);
            eik_bubble(Gb, a2, a1, eps, prec);
            acb_set_si(arg, a1 + a2);
            acb_sub(arg, arg, Dh, prec);
            acb_set_si(bb, b);
            acb_set_si(cc, c);
            fmpq_t zq;
            fmpq_init(zq);
            fmpq_set_si(zq, 1, 1);
            fmpq_sub(zq, zq, yq);
            fmpq_div_2exp(zq, zq, 1);
            eik_pair_value(v, arg, bb, cc, zq, /*cont_ok=*/false, b, c,
                           eps, prec, "g1_gate_tc3_e0");
            fmpq_clear(zq);
            acb_mul(v, v, Gb, prec);
            acb_clear(arg); acb_clear(bb); acb_clear(cc); acb_clear(Gb);
        } else if (mode == 3) {
            // e=0 closed 1-fold (C.4): pref * Int_0^1 lam^{b-1}(1-lam)^{c-1}
            // (u_lam^2)^{s1+s2} dlam
            acb_t s12, pref;
            acb_init(s12); acb_init(pref);
            acb_set_si(s12, a1 + a2);
            acb_submul_si(s12, Dh, 2, prec);          // s1+s2 = a1+a2-D
            acb_sub_si(t1, Dh, a1, prec);
            acb_gamma(pref, t1, prec);
            acb_sub_si(t1, Dh, a2, prec);
            acb_gamma(t2, t1, prec);
            acb_mul(pref, pref, t2, prec);
            // G(2a2+bh3-D) = G(2(a1+a2)+b+c-2D)
            acb_set_si(t1, 2 * (a1 + a2) + b + c);
            acb_submul_si(t1, Dh, 4, prec);
            acb_gamma(t2, t1, prec);
            acb_mul(pref, pref, t2, prec);
            const long dens[4] = {a1, a2, b, c};
            for (int i = 0; i < 4; ++i) {
                acb_set_si(t1, dens[i]);
                acb_rgamma(t2, t1, prec);
                acb_mul(pref, pref, t2, prec);
            }
            arb_t yv, half;
            arb_init(yv); arb_init(half);
            arb_set_fmpq(yv, yq, prec + 60);
            arb_one(half);
            arb_mul_2exp_si(half, half, -1);
            auto K = [&](acb_t outv, const arb_t lam, const arb_t oml,
                         long p_) {
                arb_t rho2, ta;
                arb_init(rho2); arb_init(ta);
                arb_mul(rho2, lam, oml, p_);
                arb_mul_2exp_si(rho2, rho2, 1);
                arb_sub_si(ta, yv, 1, p_);
                arb_mul(rho2, rho2, ta, p_);
                arb_add_si(rho2, rho2, 1, p_);
                acb_t k1v, k2v, k3v;
                acb_init(k1v); acb_init(k2v); acb_init(k3v);
                acb_set_arb(k1v, lam);
                acb_pow_si(k1v, k1v, b - 1, p_);
                acb_set_arb(k2v, oml);
                acb_pow_si(k2v, k2v, c - 1, p_);
                acb_set_arb(k3v, rho2);
                acb_pow(k3v, k3v, s12, p_);
                acb_mul(outv, k1v, k2v, p_);
                acb_mul(outv, outv, k3v, p_);
                acb_clear(k1v); acb_clear(k2v); acb_clear(k3v);
                arb_clear(rho2); arb_clear(ta);
            };
            acb_t J, seg;
            acb_init(J); acb_init(seg);
            auto F0 = [&](acb_t outv, const arb_t d, long p_) {
                arb_t oml;
                arb_init(oml);
                arb_sub_si(oml, d, 1, p_);
                arb_neg(oml, oml);
                K(outv, d, oml, p_);
                arb_clear(oml);
            };
            auto F1 = [&](acb_t outv, const arb_t d, long p_) {
                arb_t oml;
                arb_init(oml);
                arb_sub_si(oml, d, 1, p_);
                arb_neg(oml, oml);
                K(outv, oml, d, p_);
                arb_clear(oml);
            };
            eik_ts_edge3(J, F0, half, prec, "tc3 e0 lam[0,1/2]");
            eik_ts_edge3(seg, F1, half, prec, "tc3 e0 lam[1/2,1]");
            acb_add(J, J, seg, prec);
            acb_mul(v, pref, J, prec);
            acb_clear(J); acb_clear(seg);
            acb_clear(s12); acb_clear(pref);
            arb_clear(yv); arb_clear(half);
        } else {
            std::snprintf(buf, (std::size_t)buflen, "ERROR bad mode");
            fmpq_clear(yq);
            acb_clear(v); acb_clear(Dh); acb_clear(t1); acb_clear(t2);
            acb_clear(g);
            return;
        }
        char* rs = arb_get_str(acb_realref(v), digits, 0);
        char* is = arb_get_str(acb_imagref(v), digits, 0);
        std::snprintf(buf, (std::size_t)buflen, "%s %s", rs, is);
        flint_free(rs);
        flint_free(is);
    } catch (const std::exception& ex) {
        std::snprintf(buf, (std::size_t)buflen, "THROW: %s", ex.what());
    }
    fmpq_clear(yq);
    acb_clear(v); acb_clear(Dh); acb_clear(t1); acb_clear(t2); acb_clear(g);
}

// Bare-eikonal consistency checks:
// mode 0: synthetic ended family {2*k1*u1, k1^2}, master [1,1] — the
//         genuine scaleless case: expect "0 0" (+ the [eik-val]
//         bare-scaleless debug line under AMFLOW_EIK_DEBUG);
// mode 1: {2*k1*u1, k1^2, (k1-q)^2}, master [1,1,1] — expect F(eps);
// mode 2: same family, master [2,1,1] — expect THROW (stage-2 gap loud);
// mode 3: bare eik on a coupled 2-loop component — expect THROW;
// mode 4: F(eps) closed form directly (transcription oracle check);
// mode 5: {2*k1*u1, k1^2-1}, master [1,1] — bare + massive quad: THROW;
// mode 6: {2*k1*u1, 2*k1*u1-1, k1^2, (k1-q)^2}, master [1,1,1,1] — mixed
//         bare/residual-mass slots on one loop: THROW.
void amflow_eik_gate_bare(long mode, const char* eps_re, const char* eps_im,
                          long prec, long digits, char* buf, long buflen) {
    numeric::AcbValue eps;
    arb_set_str(acb_realref(eps.raw()), eps_re, prec);
    arb_set_str(acb_imagref(eps.raw()), eps_im, prec);
    try {
        if (mode == 4) {
            acb_t f;
            acb_init(f);
            eik_bare_F(f, eps, prec);
            char* rs = arb_get_str(acb_realref(f), digits, 0);
            char* is = arb_get_str(acb_imagref(f), digits, 0);
            std::snprintf(buf, (std::size_t)buflen, "%s %s", rs, is);
            flint_free(rs);
            flint_free(is);
            acb_clear(f);
            return;
        }
        std::vector<std::string> props;
        std::vector<long> idx;
        std::vector<std::string> loops = {"k1"};
        if (mode == 0) {
            props = {"2*k1*u1", "k1^2"};
            idx = {1, 1};
        } else if (mode == 1) {
            props = {"2*k1*u1", "k1^2", "(k1-q)^2"};
            idx = {1, 1, 1};
        } else if (mode == 2) {
            props = {"2*k1*u1", "k1^2", "(k1-q)^2"};
            idx = {2, 1, 1};
        } else if (mode == 3) {
            loops = {"k1", "k2"};
            props = {"2*k1*u1", "(k1+k2)^2", "k2^2", "2*k2*u1-1"};
            idx = {1, 1, 1, 1};
        } else if (mode == 5) {
            props = {"2*k1*u1", "k1^2-1"};
            idx = {1, 1};
        } else if (mode == 6) {
            props = {"2*k1*u1", "2*k1*u1-1", "k1^2", "(k1-q)^2"};
            idx = {1, 1, 1, 1};
        } else if (mode >= 10 && mode <= 12) {
            // FLOW-DISPATCH gate for class 3: a production reference family
            // shape VERBATIM (identity routing throws; the GL(2,Z)
            // wrapper must find k1 -> k1+k2 and the recognizer must map
            // the slots, incl. the b<->c mirror on the dotted variants):
            //   mode 10: master [1,1,1,1,1]  -> -TC3(1,1,1,1,1)
            //   mode 11: master [1,1,1,2,1]  -> +TC3(1,1,2,1,1) (b-dot)
            //   mode 12: master [1,2,1,1,1]  -> +TC3(1,1,1,2,1) (c-dot)
            // (sign = (-1)^{#odd-index c0=-1 eiks + #odd-index quads})
            loops = {"k1", "k2"};
            props = {"2*k2*u1-1", "2*k1*u2-2*k2*u2-1", "k1^2",
                     "2*k1*u1-2*k2*u1-1", "k2^2"};
            if (mode == 10) idx = {1, 1, 1, 1, 1};
            if (mode == 11) idx = {1, 1, 1, 2, 1};
            if (mode == 12) idx = {1, 2, 1, 1, 1};
            auto fc3 = qft::FamilyConfig::build(
                "synthc3", loops, {"u1", "u2"}, {},
                {{"u1^2", "1"}, {"u2^2", "1"}, {"u1*u2", "5/4"}},
                props);
            qft::JIntegral j3("synthc3", idx);
            AMFSystemOptions opts3;
            numeric::AcbValue out3;
            eikonal_ending_value(fc3, j3, opts3, eps, out3, prec);
            char* rs3 = arb_get_str(acb_realref(out3.raw()), digits, 0);
            char* is3 = arb_get_str(acb_imagref(out3.raw()), digits, 0);
            std::snprintf(buf, (std::size_t)buflen, "%s %s", rs3, is3);
            flint_free(rs3);
            flint_free(is3);
            return;
        } else {
            std::snprintf(buf, (std::size_t)buflen, "ERROR bad mode");
            return;
        }
        auto fc = qft::FamilyConfig::build(
            "synthbare", loops, {"u1", "q"}, {},
            {{"u1^2", "1"}, {"q^2", "-1"}, {"u1*q", "0"}},
            props);
        qft::JIntegral j("synthbare", idx);
        AMFSystemOptions opts;
        numeric::AcbValue outv;
        eikonal_ending_value(fc, j, opts, eps, outv, prec);
        char* rs = arb_get_str(acb_realref(outv.raw()), digits, 0);
        char* is = arb_get_str(acb_imagref(outv.raw()), digits, 0);
        std::snprintf(buf, (std::size_t)buflen, "%s %s", rs, is);
        flint_free(rs);
        flint_free(is);
    } catch (const std::exception& ex) {
        std::snprintf(buf, (std::size_t)buflen, "THROW: %s", ex.what());
    }
}

}  // extern "C"

bool try_builtin_ending_value(const qft::FamilyConfig& fc,
                              const qft::JIntegral& j,
                              const AMFSystemOptions& opts,
                              std::size_t eps_index,
                              const numeric::AcbValue& eps,
                              numeric::AcbValue& out) {
    const std::string key = jintegral_key(j);

    // Cut-ending vacuum: closed-form maximal cut of the
    // equal-mass 3-loop tetrahedron (K4 "Mercedes") vacuum that the eta->inf
    // ending under the maximal cut [1,1,1,1,1,1,0..] reduces to.  Computed
    // from the exact `eps` of each sample so it is correct on the whole grid.
    //
    // Baikov maximal cut (P. A. Baikov, hep-ph/9611449; cut prescription as
    // in Frellesvig-Papadopoulos, arXiv:1701.07356), measure
    // prod d^d l/(i pi^{d/2}):
    //   MaxCut[a] = K3 * 1/prod((a_i-1)!) * d^{a-1}_z [ U(z)^{-eps} ]|_{z=0}
    //   K3 = 2^{2-2eps} i / (pi^2 Gamma(1-eps) Gamma(3-2eps))   (i = 1/i^3)
    // NOTE (unresolved): the code below computes the power as 2^{2-eps}, not
    // the 2^{2-2eps} written here; which exponent is correct has not been
    // settled in this tree. The feature is experimental and its overall
    // normalization is user-supplied (CONV below), so the discrepancy is
    // recorded rather than silently "fixed"; see PATCHES.md.
    // U(z) is the loop Gram det in propagator variables (degree 3, 6 vars).
    // We multiply by an env-set convention factor CONV (default i) absorbing
    // the measure-i and any AMFlow family-cut normalization (not validated
    // in this tree).
    // Activated by setting env AMFLOW_MERCEDES_CUT to the family name;
    // applies to any top-sector master  fam|a1|a2|a3|a4|a5|a6|0|...|0  (a_i>=1).
    if (const char* env = std::getenv("AMFLOW_MERCEDES_CUT"); env && *env) {
        const auto& idx = j.indices();
        bool is_top = (j.family() == env) && (idx.size() == 15);
        if (is_top) {
            for (int p = 0; p < 6; ++p) if (idx[p] < 1) is_top = false;
            for (int p = 6; p < 15; ++p) if (idx[p] != 0) is_top = false;
        }
        if (is_top) {
            const long prec = numeric::working_prec_bits();
            int a[6];
            int amax = 1;
            for (int p = 0; p < 6; ++p) { a[p] = (int)idx[p]; if (a[p] > amax) amax = a[p]; }

            // ---- U(z) monomials: list of (exps[6], num, den) with coeff num/den ----
            // (generated from sympy; U has 47 monomials, coeffs +/-1/4 and 1/2 const)
            struct Mono { int e[6]; int num; int den; };
            static const Mono U[] = {
              {{2,0,1,0,0,0},-1,4},{{2,0,0,0,0,0},-1,4},{{1,1,1,0,0,0},1,4},
              {{1,1,0,1,0,0},1,4},{{1,1,0,0,0,1},-1,4},{{1,1,0,0,0,0},1,4},
              {{1,0,2,0,0,0},-1,4},{{1,0,1,1,0,0},1,4},{{1,0,1,0,1,0},1,4},
              {{1,0,1,0,0,1},1,4},{{1,0,0,1,1,0},-1,4},{{1,0,0,1,0,0},1,4},
              {{1,0,0,0,1,1},1,4},{{1,0,0,0,1,0},1,4},{{1,0,0,0,0,1},1,4},
              {{1,0,0,0,0,0},1,4},{{0,2,0,1,0,0},-1,4},{{0,2,0,0,0,0},-1,4},
              {{0,1,1,1,0,0},1,4},{{0,1,1,0,1,0},-1,4},{{0,1,1,0,0,0},1,4},
              {{0,1,0,2,0,0},-1,4},{{0,1,0,1,1,0},1,4},{{0,1,0,1,0,1},1,4},
              {{0,1,0,0,1,1},1,4},{{0,1,0,0,1,0},1,4},{{0,1,0,0,0,1},1,4},
              {{0,1,0,0,0,0},1,4},{{0,0,2,0,0,0},-1,4},{{0,0,1,1,0,1},-1,4},
              {{0,0,1,1,0,0},1,4},{{0,0,1,0,1,1},1,4},{{0,0,1,0,1,0},1,4},
              {{0,0,1,0,0,1},1,4},{{0,0,1,0,0,0},1,4},{{0,0,0,2,0,0},-1,4},
              {{0,0,0,1,1,1},1,4},{{0,0,0,1,1,0},1,4},{{0,0,0,1,0,1},1,4},
              {{0,0,0,1,0,0},1,4},{{0,0,0,0,2,1},-1,4},{{0,0,0,0,2,0},-1,4},
              {{0,0,0,0,1,2},-1,4},{{0,0,0,0,1,0},1,4},{{0,0,0,0,0,2},-1,4},
              {{0,0,0,0,0,1},1,4},{{0,0,0,0,0,0},1,2}
            };
            const int NU = (int)(sizeof(U)/sizeof(U[0]));

            // We need the multivariate Taylor coefficient of U(z)^{-eps} at z=0
            // for the monomial z^{a-1}.  Write U = U0 (1 + W), U0 = 1/2,
            // W = (U - U0)/U0.  Then U^{-eps} = U0^{-eps} sum_{k} C(-eps,k) W^k.
            // Max total degree needed = sum(a_i-1) <= 6*(amax-1).  Since W has no
            // constant term, only k <= deg matters.  We accumulate the target
            // coefficient as a polynomial-in-W expansion truncated at the needed
            // multidegree.  Implementation: dense map over exponent vectors.
            // Needed exponent vector:
            int need[6];
            int needtot = 0;
            for (int p = 0; p < 6; ++p) { need[p] = a[p]-1; needtot += need[p]; }

            // Represent a multivariate polynomial as map<key,acb>; key packs the
            // 6 exponents (each < 16) into one 32-bit int.
            auto packk = [](const int e[6]) {
                long k = 0; for (int p = 0; p < 6; ++p) k = (k<<4) | (e[p] & 0xF); return k;
            };
            auto le_need = [&](const int e[6]) {
                for (int p = 0; p < 6; ++p) if (e[p] > need[p]) return false; return true;
            };

            // W = (U - 1/2)/(1/2) = 2*U - 1  (constant term cancels).
            // Build W monomials (drop constant) truncated to <= need.
            std::map<long, fmpq*> Wm;   // exponent-key -> rational coeff (2*Umono)
            auto add_q = [&](std::map<long,fmpq*>& m, const int e[6], int num, int den) {
                long k = packk(e);
                auto it = m.find(k);
                if (it == m.end()) { fmpq* q = (fmpq*)flint_malloc(sizeof(fmpq)); fmpq_init(q);
                    fmpq_set_si(q, num, den); m[k] = q; }
                else { fmpq_t t; fmpq_init(t); fmpq_set_si(t, num, den);
                    fmpq_add(it->second, it->second, t); fmpq_clear(t); }
            };
            for (int m = 0; m < NU; ++m) {
                bool allzero = true; for (int p=0;p<6;++p) if (U[m].e[p]) allzero=false;
                if (allzero) continue;             // skip constant 1/2
                if (!le_need(U[m].e)) continue;
                add_q(Wm, U[m].e, 2*U[m].num, U[m].den);   // 2 * coeff
            }

            // Series U0^{-eps} * sum_{k>=0} binom(-eps,k) W^k, keep coeff of z^need.
            // Iterate: term_0 = 1 (coeff map: {0:1}); term_{k} = term_{k-1} * W
            // (truncated), and accumulate binom(-eps,k)*coeff(need) into acb result.
            // binom(-eps,k) = prod_{j=0}^{k-1} (-eps - j)/(k - j) ... = (-eps)(-eps-1)...(-eps-k+1)/k!
            std::map<long, fmpq*> cur;      // current W^k coefficients (rational)
            { int z6[6]={0,0,0,0,0,0}; fmpq* q=(fmpq*)flint_malloc(sizeof(fmpq)); fmpq_init(q);
              fmpq_set_si(q,1,1); cur[packk(z6)]=q; }

            acb_t Racc; acb_init(Racc); acb_zero(Racc);     // sum binom(-eps,k)*[z^need]W^k
            acb_t binom; acb_init(binom); acb_set_si(binom, 1);  // binom(-eps,0)=1
            long needkey = packk(need);

            for (int kk = 0; kk <= needtot; ++kk) {
                // contribution of W^kk: coeff of z^need in cur, times binom(-eps,kk)
                auto it = cur.find(needkey);
                if (it != cur.end()) {
                    acb_t cval; acb_init(cval);
                    acb_set_fmpq(cval, it->second, prec);
                    acb_addmul(Racc, cval, binom, prec);
                    acb_clear(cval);
                }
                if (kk == needtot) break;
                // update binom(-eps,kk+1) = binom * (-eps - kk)/(kk+1)
                { acb_t f; acb_init(f); acb_neg(f, eps.raw()); acb_sub_si(f, f, kk, prec);
                  acb_mul(binom, binom, f, prec); acb_div_si(binom, binom, kk+1, prec); acb_clear(f); }
                // cur = cur * W  (truncated to <= need)
                std::map<long, fmpq*> nxt;
                for (auto& [ck, cq] : cur) {
                    int ce[6]; long t=ck; for (int p=5;p>=0;--p){ ce[p]=t&0xF; t>>=4; }
                    for (auto& [wk, wq] : Wm) {
                        int we[6]; long s=wk; for (int p=5;p>=0;--p){ we[p]=s&0xF; s>>=4; }
                        int ne[6]; bool ok=true;
                        for (int p=0;p<6;++p){ ne[p]=ce[p]+we[p]; if (ne[p]>need[p]){ok=false;break;} }
                        if (!ok) continue;
                        long nk = packk(ne);
                        fmpq_t prod; fmpq_init(prod); fmpq_mul(prod, cq, wq);
                        auto nit = nxt.find(nk);
                        if (nit == nxt.end()) { fmpq* q=(fmpq*)flint_malloc(sizeof(fmpq)); fmpq_init(q);
                            fmpq_set(q, prod); nxt[nk]=q; }
                        else fmpq_add(nit->second, nit->second, prod);
                        fmpq_clear(prod);
                    }
                }
                for (auto& [k2,q2] : cur) { fmpq_clear(q2); flint_free(q2); }
                cur.swap(nxt);
            }
            for (auto& [k2,q2] : cur) { fmpq_clear(q2); flint_free(q2); }
            for (auto& [k2,q2] : Wm) { fmpq_clear(q2); flint_free(q2); }
            acb_clear(binom);

            // K3 = 2^{2-2eps} i / (pi^2 Gamma(1-eps) Gamma(3-2eps)); U0^{-eps}=2^{eps}.
            // MaxCut = K3 * 2^{eps} * Racc / prod((a_i-1)!)
            //        = 2^{2-eps} i /(pi^2 Gamma(1-eps)Gamma(3-2eps)) * Racc / prod(fact)
            acb_t pi2, g1, g2, twop, val;
            acb_init(pi2); acb_init(g1); acb_init(g2); acb_init(twop); acb_init(val);
            acb_const_pi(pi2, prec); acb_mul(pi2, pi2, pi2, prec);
            acb_neg(g1, eps.raw()); acb_add_si(g1, g1, 1, prec); acb_gamma(g1, g1, prec); // G(1-eps)
            acb_mul_si(g2, eps.raw(), -2, prec); acb_add_si(g2, g2, 3, prec); acb_gamma(g2, g2, prec); // G(3-2eps)
            { acb_t base, ex; acb_init(base); acb_init(ex); acb_set_si(base, 2);
              acb_neg(ex, eps.raw()); acb_add_si(ex, ex, 2, prec); acb_pow(twop, base, ex, prec); // 2^{2-eps} (header comment says 2^{2-2eps}; unresolved, see NOTE above)
              acb_clear(base); acb_clear(ex); }
            acb_mul(val, twop, Racc, prec);
            acb_t den; acb_init(den); acb_mul(den, pi2, g1, prec); acb_mul(den, den, g2, prec);
            acb_div(val, val, den, prec);
            // divide by prod (a_i-1)!
            long fac = 1; for (int p=0;p<6;++p){ for(int q=2;q<=a[p]-1;++q) fac*=q; }
            acb_div_si(val, val, fac, prec);
            // convention factor CONV (default i): absorbs measure-i + cut-norm
            double cre = 0.0, cim = 1.0;
            if (const char* r = std::getenv("AMFLOW_MERCEDES_CONV_RE")) cre = std::atof(r);
            if (const char* im = std::getenv("AMFLOW_MERCEDES_CONV_IM")) cim = std::atof(im);
            acb_t conv; acb_init(conv); acb_set_d_d(conv, cre, cim);
            acb_mul(out.raw(), val, conv, prec);
            if (std::getenv("AMFLOW_MERCEDES_DEBUG")) {
                char* es = arb_get_str(acb_realref(eps.raw()), 20, 0);
                char* vr = arb_get_str(acb_realref(val), 28, 0);
                char* vi = arb_get_str(acb_imagref(val), 28, 0);
                std::cerr << "[mercedes] key=" << key << " eps=" << es
                          << " val(preconv)= " << vr << " + " << vi << "i"
                          << " fac=" << fac << "\n";
                flint_free(es); flint_free(vr); flint_free(vi);
            }
            acb_clear(pi2); acb_clear(g1); acb_clear(g2); acb_clear(twop);
            acb_clear(val); acb_clear(den); acb_clear(conv); acb_clear(Racc);
            return true;
        }
    }

    auto eit = opts.explicit_boundary.find(key);
    if (eit != opts.explicit_boundary.end() && eit->second && !eit->second->empty()) {
        // A single-element vector is broadcast to every eps sample (an
        // eps-independent boundary constant); otherwise the vector must be
        // long enough to index by eps_index.
        if (eit->second->size() == 1) {
            out = (*eit->second)[0].clone();
            return true;
        }
        if (eps_index >= eit->second->size()) {
            throw std::runtime_error(
                "AMFSystem::solve: explicit_boundary too short for '" + key + "'");
        }
        out = (*eit->second)[eps_index].clone();
        return true;
    }

    bool any_pos = false;
    for (long v : j.indices()) {
        if (v > 0) {
            any_pos = true;
            break;
        }
    }
    if (!any_pos) {
        acb_one(out.raw());
        return true;
    }

    // EikonalVacuum ending (region engine): must run BEFORE the
    // count-keyed (loop_count, prop_count) vacuum dispatch below — that
    // table keys on index counts only and would misfire on eikonal
    // families (e.g. an l~eta mixed-region family with 2 loops and 3
    // positive indices is NOT the 2-loop massive sunset).  This is the
    // re-audit of the P0 (L,L) product-of-tadpoles special case demanded
    // by FOLLOWUP_BUILD_SPEC.md section 4: with the engine on, any family
    // carrying an eikonal-mass propagator takes this path exclusively.
    if (eikonal_regions_on() && has_eikonal_mass_prop(fc)) {
        eikonal_ending_value(fc, j, opts, eps, out,
                             numeric::working_prec_bits());
        return true;
    }

    long n_props = 0;
    bool all_one_or_zero = true;
    for (long v : j.indices()) {
        if (v == 1) {
            ++n_props;
        } else if (v != 0) {
            all_one_or_zero = false;
            break;
        }
    }

    // === Single-mass vacuum boundary (possibly FACTORIZED) ===
    //
    // The eta->infinity ending corner of a multi-loop graph is a single-mass
    // vacuum: every active line is massless except the one carrying the
    // auxiliary mass.  For a CONNECTED such vacuum the value is the closed
    // form `qft::vacuum(L, n)` (banana family n=L+1, plus the (1,1) tadpole
    // and the (3,5) corner).  But for fully-massless graphs the ending corner
    // frequently FACTORIZES across the loops into a product of smaller
    // single-mass vacua sharing no loop momentum — e.g. the massless planar
    // double box bottoms out at two disjoint massive tadpoles
    // {k^2-1}*{l^2-1}, a (2,2) vacuum that is NOT a banana and is not in any
    // (L,n) table, but equals Gamma(eps-1)^2 = vacuum(1,1)^2.
    //
    // We detect the factorization directly from the topology: two loops are
    // in the same component iff some active propagator depends on both.  Each
    // connected component is itself a single-mass vacuum `vacuum(L_c, n_c)`
    // (its loop labels are irrelevant — `qft::vacuum` only consumes counts),
    // and the integral is the product over components.  A single component
    // reproduces the old connected lookup exactly, so this strictly extends
    // the previous behaviour.
    if (all_one_or_zero && n_props > 0) {
        const long L = (long)fc.n_loops();
        // For each active propagator, the set of loops it depends on.
        std::vector<std::vector<long>> prop_loops;       // active props only
        const auto& idx = j.indices();
        for (std::size_t p = 0; p < idx.size(); ++p) {
            if (idx[p] != 1) continue;
            if (p >= fc.propagators_after_conservation.size()) break;
            const auto& prop = fc.propagators_after_conservation[p];
            std::vector<long> deps;
            for (long li = 0; li < L; ++li) {
                if (prop.degree(li) > 0) deps.push_back(li);
            }
            prop_loops.push_back(std::move(deps));
        }
        // Union-find over the L loops, linked by shared active propagators.
        std::vector<long> uf(L);
        for (long i = 0; i < L; ++i) uf[i] = i;
        std::function<long(long)> find = [&](long x) {
            while (uf[x] != x) { uf[x] = uf[uf[x]]; x = uf[x]; }
            return x;
        };
        auto unite = [&](long a, long b) { uf[find(a)] = find(b); };
        for (const auto& deps : prop_loops) {
            for (std::size_t q = 1; q < deps.size(); ++q) unite(deps[0], deps[q]);
        }
        // Group active props by their component root.  A prop with no loop
        // dependence (pure ISP=const after the mass shift) cannot occur for a
        // power-1 vacuum master, but guard anyway by skipping it.
        std::map<long, std::pair<long, long>> comp;  // root -> (L_c, n_c)
        bool ok = true;
        for (std::size_t a = 0; a < prop_loops.size(); ++a) {
            if (prop_loops[a].empty()) { ok = false; break; }
            long root = find(prop_loops[a][0]);
            comp[root].second += 1;  // n_c
        }
        // Count loops per component.
        for (long li = 0; li < L; ++li) comp[find(li)].first += 1;
        // Each component must be a recognised single-mass vacuum.
        if (ok) {
            for (const auto& [root, ln] : comp) {
                const long Lc = ln.first;
                const long nc = ln.second;
                // A component with loops but no active propagator (nc==0) is a
                // scaleless massless sub-vacuum -> the whole integral is 0.
                if (nc == 0) {
                    out.set_zero();
                    return true;
                }
                if (!qft::vacuum_known(Lc, nc)) { ok = false; break; }
            }
        }
        if (ok) {
            const long prec = numeric::working_prec_bits();
            numeric::AcbValue acc;
            acb_one(acc.raw());
            for (const auto& [root, ln] : comp) {
                numeric::AcbValue piece;
                qft::vacuum(ln.first, ln.second, eps.raw(), piece.raw(), prec);
                acb_mul(acc.raw(), acc.raw(), piece.raw(), prec);
            }
            out = std::move(acc);
            return true;
        }
    }

    return false;
}

ibp::ReduceResult&
get_ending_reduction(const qft::FamilyConfig& fc,
                     const qft::JIntegral& target,
                     long system_id,
                     const AMFSystemOptions& opts,
                     std::map<std::string, ibp::ReduceResult>& cache) {
    const std::string key = jintegral_key(target);
    auto it = cache.find(key);
    if (it != cache.end()) return it->second;

    ibp::ReduceOptions reduce_opts = opts.bb;
    if (!opts.cache_root.empty()) {
        reduce_opts.work_dir = opts.cache_root + "/system_"
                             + std::to_string(system_id)
                             + "_ending_reduce_" + key;
    } else if (reduce_opts.work_dir.empty()) {
        reduce_opts.work_dir = (fs::temp_directory_path() /
            ("amflow_amfsys_" + std::to_string(system_id)
             + "_ending_reduce_" + key)).string();
    } else {
        reduce_opts.work_dir += "/ending_reduce_" + key;
    }

    // stage-2 endings: on eikonal-vacuum families bias Kira toward
    // denominator masters (subsector corners + one-dot variants) — the
    // stage-2 ending values are wired for denominator indices only
    // (numerator u-projection is the spec'd stage-1/2 gap).  Gated to
    // eikonal families with the engine on: default paths untouched.
    std::vector<qft::JIntegral> ending_preferred =
        build_eikonal_corner_preferred(fc, {target});
    // AMFLOW_KIRA_FFSAVE=1: skip the pre-wipe; ibp::reduce's own
    // input-match logic decides resume-vs-fresh per embedded Kira dir
    // (see reduce.hpp).  Knob OFF (default): stock wipe.
    if (!ibp::kira_ffsave_enabled()) {
        fs::remove_all(reduce_opts.work_dir);
    }

    auto inserted = cache.emplace(
        key, ibp::reduce(fc, {target}, ending_preferred,
                              qft::get_top_sector({target}), reduce_opts));
    return inserted.first->second;
}

numeric::AcbValue solve_ending_master_value(
        const qft::FamilyConfig& fc,
        const qft::JIntegral& target,
        long system_id,
        const AMFSystemOptions& opts,
        std::size_t eps_index,
        const numeric::AcbValue& eps,
        const std::map<std::string, FmpqHolder>& numeric_q,
        std::map<std::string, ibp::ReduceResult>& reduction_cache,
        std::map<std::string, numeric::AcbValue>& value_cache,
        std::set<std::string>& active_keys,
        long prec) {
    const std::string key = jintegral_key(target);
    auto vit = value_cache.find(key);
    if (vit != value_cache.end()) return vit->second.clone();

    if (!active_keys.insert(key).second) {
        throw std::runtime_error(
            "AMFSystem::solve: cyclic ending reduction for '" + key + "'");
    }

    numeric::AcbValue out;
    if (try_builtin_ending_value(fc, target, opts, eps_index, eps, out)) {
        value_cache.emplace(key, midpoint_only(out));
        active_keys.erase(key);
        return out;
    }

    ibp::ReduceResult& red =
        get_ending_reduction(fc, target, system_id, opts, reduction_cache);
    if (red.rules.size() != 1 || red.rules[0].empty()) {
        active_keys.erase(key);
        throw std::runtime_error(
            "AMFSystem::solve: ending system master '" + key
            + "' has no Vacuum entry, no explicit_boundary value, "
            "and no usable reduction rule");
    }

    numeric::AcbValue rhs_sum;
    acb_zero(rhs_sum.raw());
    numeric::AcbValue self_coef;
    acb_zero(self_coef.raw());
    bool saw_non_self = false;
    for (const auto& term : red.rules[0]) {
        numeric::AcbValue coef;
        evaluate_const_mfrac(term.coef, numeric_q, prec, coef);
        if (jintegral_key(term.integ) == key) {
            acb_add(self_coef.raw(), self_coef.raw(), coef.raw(), prec);
            continue;
        }
        saw_non_self = true;
        numeric::AcbValue rhs = solve_ending_master_value(
            fc, term.integ, system_id, opts, eps_index, eps, numeric_q,
            reduction_cache, value_cache, active_keys, prec);
        numeric::AcbValue prod;
        acb_mul(prod.raw(), coef.raw(), rhs.raw(), prec);
        acb_add(rhs_sum.raw(), rhs_sum.raw(), prod.raw(), prec);
    }

    numeric::AcbValue denom;
    acb_one(denom.raw());
    acb_sub(denom.raw(), denom.raw(), self_coef.raw(), prec);

    if (!saw_non_self && acb_is_one(self_coef.raw())) {
        active_keys.erase(key);
        throw std::runtime_error(
            "AMFSystem::solve: ending system master '" + key
            + "' has no Vacuum entry, no explicit_boundary value, "
            "and reduction did not lower it");
    }

    // denom = 1 - self_coef.  When the requested ε makes this vanish (the
    // generic situation at ε=0 exactly for a SingleMass / vacuum-type leaf
    // whose ending IBP reduction becomes an identity, 1 - self_coef -> 0), the
    // division yields Inf/NaN.  The boundary constant of such a leaf genuinely
    // carries a 1/ε pole at ε=0; the assembled physical integral is finite
    // only AFTER subsector recombination.  Evaluating it at ε=0 *exactly* is
    // therefore ill-posed -- no amount of working precision can rationalize a
    // NaN.  Surface this as a clear, actionable error instead of letting the
    // NaN propagate downstream and resurface as the misleading calcx00
    // "rationalized boundary matching ... interval too wide" abort.
    acb_div(out.raw(), rhs_sum.raw(), denom.raw(), prec);
    // A non-finite quotient means denom = 1 - self_coef vanished: the ending
    // IBP reduction degenerated to an identity at this epsilon.  At eps=0
    // exactly this is the generic situation for a SingleMass / vacuum-type
    // leaf whose boundary constant carries a 1/eps pole (the assembled
    // physical integral is eps^0-finite only after subsector recombination),
    // so a direct eps=0 evaluation is ill-posed -- no working-precision
    // increase can rationalize a NaN.  Surface a clear, actionable error
    // instead of letting the NaN propagate and resurface downstream as the
    // misleading calcx00 "interval too wide" abort.
    if (!acb_is_finite(out.raw())) {
        active_keys.erase(key);
        char* dn = arb_get_str(acb_realref(denom.raw()), 12, 0);
        char* eps_s = arb_get_str(acb_realref(eps.raw()), 12, 0);
        std::ostringstream msg;
        msg << "AMFSystem::solve: ending master '" << key
            << "' evaluated to a non-finite value at this epsilon "
               "(1 - self_coef ~ " << dn << ", eps ~ " << eps_s
            << ").  The boundary constant is divergent here, typically a 1/eps "
               "pole at eps=0; a direct eps=0 evaluation is ill-posed.  Use a "
               "small nonzero-eps grid and extrapolate to eps=0 (e.g. a base-2 "
               "grid eps = 2^-4 .. 2^-18) instead of eps_samples=[0].";
        flint_free(dn);
        flint_free(eps_s);
        throw std::runtime_error(msg.str());
    }
    value_cache.emplace(key, midpoint_only(out));
    active_keys.erase(key);
    return out;
}

// Test whether a qft::TopSectorComponentInfo, after applying numeric
// substitution to mass entries, has exactly one mass=1 and all others
// mass=0.
bool vacuum_q_numeric(const qft::TopSectorComponentInfo& info,
                      const std::map<std::string, FmpqHolder>&) {
    return qft::vacuum_q(info);
}

// Substitute Numeric values into the mass list before applying the
// literal "exactly one mass=1, all others mass=0" test.  Upstream's
// `SingleMassQ` (`AMFlow.m`) tests literally and rejects any
// symbolic mass entry; this C++ enhancement recognises the case
// where the user-supplied `Numeric` map resolves the symbols to
// the literal `0` / `1` shape, allowing SingleMass scheme to fire
// on families like `{l^2 - msq}` with `Numeric = {msq -> 1}`.
//
// Equivalence: with empty `numeric_q`, this collapses to the
// literal upstream test (no substitution).  When `numeric_q`
// resolves the symbolic masses to literal 0/1, the substituted
// shape is mathematically a single-mass family and SingleMass
// scheme is a valid flow path — the final integral is the same
// as upstream's, just dispatched via SingleMass earlier instead
// of falling through to the next scheme.  Locked by
// `tests/test_amflow_amfsystem.cpp` `SingleMassEnhancement_*` and
// `tests/test_qft_amfmode.cpp` `SingleMassQ_*` (the latter
// confirms the qft layer's `single_mass_q` remains
// upstream-literal).
bool single_mass_q_numeric(const qft::TopSectorComponentInfo& info,
                              const std::map<std::string, FmpqHolder>& numeric_q) {
    if (!vacuum_q_numeric(info, numeric_q)) return false;
    int n_one = 0, n_zero = 0;
    for (const auto& m : info.mass) {
        algebra::Mfrac sub = apply_numeric(m, numeric_q);
        if (sub.is_zero()) ++n_zero;
        else if (sub.is_one()) ++n_one;
        else return false;   // some mass is neither 0 nor 1 -- not single-mass
    }
    return n_one == 1 && n_zero == (int)info.mass.size() - 1;
}

bool ending_q_numeric(const qft::TopSectorComponentInfo& info,
                       const std::map<std::string, FmpqHolder>& numeric_q) {
    return single_mass_q_numeric(info, numeric_q) || qft::phase_volume_q(info);
}

// Mirrors AMFSystemEndingQ[preferred_, "Cutkosky"] (.m line 996-1000):
//   ending = ! (all components are EndingQ and exactly one component is
//               PhaseVolumeQ)
bool cutkosky_ending_q_impl(const qft::FamilyConfig& fc,
                              const std::vector<qft::JIntegral>& preferred) {
    auto top = qft::get_top_position(preferred);
    if (top.empty()) return true;
    auto info = qft::analyze_top_sector(fc, top);
    if (info.empty()) return true;

    long phase_count = 0;
    for (const auto& comp : info) {
        if (!qft::ending_q(comp)) return true;
        if (qft::phase_volume_q(comp)) ++phase_count;
    }
    return phase_count != 1;
}

// Mirrors AMFSystemEndingQ[preferred_, "SingleMass"] (.m line 1002-1007):
//   ending = ! ( all components ending && cutcom==0 && Loop>0 && all idx>=0 )
// i.e. SingleMass is *applicable* (ending=False) only when all four
// conditions hold. In particular, a heterogeneous-mass vacuum is still
// ending in the Mathematica reference because EndingQ = SingleMassQ ||
// PhaseVolumeQ, so such a component fails condition 1 and never reaches
// the later FactorizeFamily setup branch.
bool single_mass_ending_q_impl(const qft::FamilyConfig& fc,
                                  const std::vector<qft::JIntegral>& preferred,
                                  const AMFSystemOptions& opts) {
    auto top = qft::get_top_position(preferred);
    if (top.empty()) return true;          // empty top -> already vacuum-like, ending
    auto info = qft::analyze_top_sector(fc, top);
    if (info.empty()) return true;         // no analyzed components -> SingleMass not applicable

    auto numeric_q = build_numeric_q(opts.bb);

    // Early exit for the same Mathematica condition-1 failure: a pure
    // vacuum component with masses other than {1,0,...,0} is VacuumQ but
    // not SingleMassQ, hence not EndingQ.
    for (const auto& comp : info) {
        if (vacuum_q_numeric(comp, numeric_q)
                && !single_mass_q_numeric(comp, numeric_q)) {
            return true;
        }
    }

    // condition 1: every component is "EndingQ" under numeric
    // (= single-mass-vacuum or phase-volume)
    for (const auto& comp : info) {
        if (!ending_q_numeric(comp, numeric_q)) return true;  // not ending -> SingleMass not applicable
    }
    // condition 2: no phase-volume (cut) component
    for (const auto& comp : info) {
        if (qft::phase_volume_q(comp)) return true;
    }
    // condition 3: Length[Loop] > 0
    if (fc.n_loops() == 0) return true;
    // condition 4: all preferred indices non-negative
    for (const auto& j : preferred) {
        for (long v : j.indices()) if (v < 0) return true;
    }
    // All four conditions satisfied -> SingleMass is applicable -> NOT ending.
    return false;
}

// ===========================================================================
//  AMFLOW_BOUNDARY_DUMP (env-gated, read-only) -- JSON dump helpers
// ===========================================================================
//
// When AMFLOW_BOUNDARY_DUMP=<dir> is set (non-empty), each non-ending
// AMFSystem writes
//   <dir>/system_<id>_boundary.json   symbolic eta=infinity boundary tables
//                                     (once, at the end of setup()), and
//   <dir>/system_<id>_bc_eps<k>.json  the numeric bc_sorted handed to
//                                     ode::amflow (once per eps sample).
// With the env var unset or empty every hook is a no-op (byte-identical
// behavior).  Write failures print a warning and continue; they never
// abort a solve.

const std::string& boundary_dump_dir() {
    static const std::string dir = [] {
        const char* e = std::getenv("AMFLOW_BOUNDARY_DUMP");
        return std::string(e ? e : "");
    }();
    return dir;
}

std::string json_escape(const std::string& s) {
    std::string out;
    out.reserve(s.size() + 8);
    for (char c : s) {
        switch (c) {
            case '"':  out += "\\\""; break;
            case '\\': out += "\\\\"; break;
            case '\n': out += "\\n";  break;
            case '\r': out += "\\r";  break;
            case '\t': out += "\\t";  break;
            default:
                if ((unsigned char)c < 0x20) {
                    char buf[8];
                    std::snprintf(buf, sizeof(buf), "\\u%04x",
                                  (unsigned)(unsigned char)c);
                    out += buf;
                } else {
                    out += c;
                }
        }
    }
    return out;
}

std::string bdump_arb_str40(const arb_t x) {
    char* s = arb_get_str(x, 40, 0);
    std::string out = s ? s : "?";
    if (s) flint_free(s);
    return out;
}

std::string bdump_fmpq_str(const fmpq_t q) {
    char* s = fmpq_get_str(nullptr, 10, q);
    std::string out = s ? s : "?";
    if (s) flint_free(s);
    return out;
}

// Try to read `m` as  a + b*eps  where the regulator variable is named
// "__amf_eps" (powers context) or "eps", the denominator is a rational
// constant, and no other variable appears.  On success writes the exact
// rationals a, b as decimal p/q strings.  The verbatim Mfrac string is
// always dumped alongside, so failure here only omits the (a, b) pair.
bool bdump_mfrac_linear_in_eps(const algebra::Mfrac& m,
                               std::string& a_out, std::string& b_out) {
    auto ctx = m.ctx();
    long eps_var = ctx->var_index("__amf_eps");
    if (eps_var < 0) eps_var = ctx->var_index("eps");
    algebra::Mpoly num = m.numerator();
    algebra::Mpoly den = m.denominator();
    long len_d = fmpz_mpoly_length(den.raw(), ctx->raw());
    if (len_d != 1) return false;
    std::vector<unsigned long> exp((std::size_t)ctx->n_vars());
    fmpz_mpoly_get_term_exp_ui(exp.data(), den.raw(), 0, ctx->raw());
    for (long v = 0; v < ctx->n_vars(); ++v) {
        if (exp[(std::size_t)v] != 0) return false;
    }
    fmpz_t dc, coeff, a_num, b_num;
    fmpz_init(dc); fmpz_init(coeff); fmpz_init(a_num); fmpz_init(b_num);
    fmpz_mpoly_get_term_coeff_fmpz(dc, den.raw(), 0, ctx->raw());
    bool ok = !fmpz_is_zero(dc);
    long len_n = fmpz_mpoly_length(num.raw(), ctx->raw());
    for (long t = 0; ok && t < len_n; ++t) {
        fmpz_mpoly_get_term_exp_ui(exp.data(), num.raw(), t, ctx->raw());
        fmpz_mpoly_get_term_coeff_fmpz(coeff, num.raw(), t, ctx->raw());
        unsigned long e_eps =
            (eps_var >= 0) ? exp[(std::size_t)eps_var] : 0ul;
        for (long v = 0; v < ctx->n_vars(); ++v) {
            if (v == eps_var) continue;
            if (exp[(std::size_t)v] != 0) { ok = false; break; }
        }
        if (!ok) break;
        if (e_eps == 0)      fmpz_add(a_num, a_num, coeff);
        else if (e_eps == 1) fmpz_add(b_num, b_num, coeff);
        else ok = false;
    }
    if (ok) {
        fmpq_t a_q, b_q;
        fmpq_init(a_q); fmpq_init(b_q);
        fmpq_set_fmpz_frac(a_q, a_num, dc);
        fmpq_set_fmpz_frac(b_q, b_num, dc);
        a_out = bdump_fmpq_str(a_q);
        b_out = bdump_fmpq_str(b_q);
        fmpq_clear(a_q); fmpq_clear(b_q);
    }
    fmpz_clear(dc); fmpz_clear(coeff); fmpz_clear(a_num); fmpz_clear(b_num);
    return ok;
}

// Emit  "mu": "<verbatim Mfrac>"[, "mu_a": "...", "mu_b": "..."]  fields.
void bdump_emit_mu(std::ostream& os, const algebra::Mfrac& m) {
    os << "\"mu\": \"" << json_escape(mfrac_str(m)) << "\"";
    std::string a, b;
    if (bdump_mfrac_linear_in_eps(m, a, b)) {
        os << ", \"mu_a\": \"" << a << "\", \"mu_b\": \"" << b << "\"";
    }
}

// Read-only classification of an ending-system master.  Mirrors the branch
// order of try_builtin_ending_value + the solve_ending_master_value
// fallback, but performs NO evaluation, NO reduction, and NO caching.
std::string classify_ending_master_json(const qft::FamilyConfig& fc,
                                        const qft::JIntegral& j,
                                        const AMFSystemOptions& opts) {
    const std::string key = jintegral_key(j);

    // 1. Mercedes maximal-cut closed form (env-gated).
    if (const char* env = std::getenv("AMFLOW_MERCEDES_CUT"); env && *env) {
        const auto& idx = j.indices();
        bool is_top = (j.family() == env) && (idx.size() == 15);
        if (is_top) {
            for (int p = 0; p < 6; ++p) if (idx[p] < 1) is_top = false;
            for (int p = 6; p < 15; ++p) if (idx[p] != 0) is_top = false;
        }
        if (is_top) return "{\"kind\": \"mercedes_cut\"}";
    }

    // 2. Explicit boundary table.
    auto eit = opts.explicit_boundary.find(key);
    if (eit != opts.explicit_boundary.end() && eit->second
            && !eit->second->empty()) {
        return "{\"kind\": \"explicit_boundary\"}";
    }

    // 3. All indices non-positive -> integral == 1.
    bool any_pos = false;
    for (long v : j.indices()) if (v > 0) { any_pos = true; break; }
    if (!any_pos) return "{\"kind\": \"one\"}";

    // 4. Power-one single-mass vacuum, possibly factorized across loops
    //    (same component analysis as try_builtin_ending_value).
    long n_props = 0;
    bool all_one_or_zero = true;
    for (long v : j.indices()) {
        if (v == 1) ++n_props;
        else if (v != 0) { all_one_or_zero = false; break; }
    }
    if (all_one_or_zero && n_props > 0) {
        const long L = (long)fc.n_loops();
        std::vector<std::vector<long>> prop_loops;
        const auto& idx = j.indices();
        for (std::size_t p = 0; p < idx.size(); ++p) {
            if (idx[p] != 1) continue;
            if (p >= fc.propagators_after_conservation.size()) break;
            const auto& prop = fc.propagators_after_conservation[p];
            std::vector<long> deps;
            for (long li = 0; li < L; ++li) {
                if (prop.degree(li) > 0) deps.push_back(li);
            }
            prop_loops.push_back(std::move(deps));
        }
        std::vector<long> uf((std::size_t)L);
        for (long i = 0; i < L; ++i) uf[(std::size_t)i] = i;
        std::function<long(long)> find = [&](long x) {
            while (uf[(std::size_t)x] != x) {
                uf[(std::size_t)x] = uf[(std::size_t)uf[(std::size_t)x]];
                x = uf[(std::size_t)x];
            }
            return x;
        };
        auto unite = [&](long a, long b) {
            uf[(std::size_t)find(a)] = find(b);
        };
        for (const auto& deps : prop_loops) {
            for (std::size_t q = 1; q < deps.size(); ++q) {
                unite(deps[0], deps[q]);
            }
        }
        std::map<long, std::pair<long, long>> comp;  // root -> (L_c, n_c)
        bool ok = true;
        for (std::size_t a = 0; a < prop_loops.size(); ++a) {
            if (prop_loops[a].empty()) { ok = false; break; }
            comp[find(prop_loops[a][0])].second += 1;
        }
        if (ok) {
            for (long li = 0; li < L; ++li) comp[find(li)].first += 1;
            bool zero = false;
            for (const auto& [root, ln] : comp) {
                if (ln.second == 0) { zero = true; break; }
                if (!qft::vacuum_known(ln.first, ln.second)) {
                    ok = false;
                    break;
                }
            }
            if (zero) return "{\"kind\": \"zero_scaleless\"}";
            if (ok) {
                std::ostringstream os;
                os << "{\"kind\": \"vacuum_product\", \"components\": [";
                bool first = true;
                for (const auto& [root, ln] : comp) {
                    if (!first) os << ", ";
                    first = false;
                    os << "[" << ln.first << ", " << ln.second << "]";
                }
                os << "]}";
                return os.str();
            }
        }
    }

    // 5. Fallback: value comes from a Kira ending reduction + recursion
    //    (solve_ending_master_value).
    return "{\"kind\": \"ending_reduction\"}";
}

}  // namespace

bool ending_q_all(const qft::FamilyConfig& fc,
                  const std::vector<qft::JIntegral>& preferred,
                  const AMFSystemOptions& opts) {
    for (auto s : opts.ending_schemes) {
        if (!ending_q(fc, preferred, s, opts)) return false;
    }
    return true;
}

// ---------------------------------------------------------------------------
//  GRAVITYFLOW_FIXEDPOINT guard helper
// ---------------------------------------------------------------------------
//
// Canonical propagator signature of a family: the propagator polynomials
// after momentum conservation, rendered to their canonical strings and
// sorted (order-insensitive).  Two families with equal signatures carry the
// byte-identical propagator content — the self-similarity criterion of the
// eta->inf fixed point (a boundary family reproducing its parent).

static std::string gravityflow_prop_signature(const qft::FamilyConfig& fc) {
    std::vector<std::string> ps;
    ps.reserve(fc.propagators_after_conservation.size());
    for (const auto& p : fc.propagators_after_conservation) {
        ps.push_back(p.to_string());
    }
    std::sort(ps.begin(), ps.end());
    std::string sig;
    for (const auto& s : ps) {
        sig += s;
        sig += ';';
    }
    return sig;
}

// ===========================================================================
//  Setup -- top level
// ===========================================================================
//
//   * If is_ending_: nothing to do.  solve() will look up Vacuum or
//     explicit_boundary.
//   * Else:
//       a. build fc_with_eta_ = inject(fc_, etac_);
//       b. build_diffeq() -- compute dM/dη on fc_with_eta_;
//       c. build_boundary() -- find regions, compute boundary
//          integrals, recursively setup sub-systems.

void AMFSystem::setup() {
    if (!is_ending_) {
        // We'll show after injection.
    }
    if (is_ending_) return;
    if (fc_with_eta_) return;
    // Self-similarity assert (region engine, FOLLOWUP_BUILD_SPEC
    // risk note): an EikonalVacuum-shaped family must NEVER reach eta
    // injection — the l~eta region regenerates the family (infinite
    // recursion with silently wrong boundary content).  Loud, not silent.
    if (eikonal_regions_on() && has_eikonal_mass_prop(*fc_)) {
        throw std::runtime_error(
            "AMFSystem::setup: family '" + fc_->family
            + "' carries an eikonal residual-mass propagator but reached "
            "eta injection — self-similar family must be ended directly "
            "(EikonalVacuum); amf_system_setup_master interception failed");
    }
    if (opts_.current_depth >= opts_.max_recursion_depth) {
        throw std::runtime_error(
            "AMFSystem::setup: max recursion depth ("
            + std::to_string(opts_.max_recursion_depth)
            + ") reached.  This typically means SingleMass / "
              "FactorizeFamily is needed to terminate the boundary "
              "integrand recursion -- see AUDIT.md.");
    }
    // GRAVITYFLOW_FIXEDPOINT guard: a boundary family whose propagator
    // set reproduces its parent's is a fixed point of the eta->inf
    // recursion — descending further regenerates the same family at
    // every level and grinds to max_recursion_depth without ever
    // producing boundary data.  Detect it here (the recursion
    // continuation point; endings never reach setup) and raise the
    // structured diagnostic instead.  The supported route for such
    // families is the self-consistent boundary closure: probe the fixed
    // point (AMFLOW_FIXEDPOINT_PROBE), then discover/respond/solve/
    // inject with amf_options.explicit_boundary and
    // AMFLOW_FORCE_ENDING_DEPTH (the pmflow pipeline).
    //
    // An identical propagator multiset alone is NOT the fixed-point
    // criterion.  The recursion state is (propagator content, preferred
    // masters), and a boundary child may legitimately re-use the
    // parent's propagator basis while its masters descend to a strictly
    // smaller top sector: the automatic-loop boundary of the massless
    // box (Mma-reference-pinned) produces a child with the parent's 4
    // propagators whose single master is a 2-slot sub-sector bubble —
    // the next level is a 1-slot tadpole ending, and the descent
    // terminates because the top-sector width (active-slot count of
    // qft::get_top_sector over the preferred masters) strictly
    // decreases, a well-founded measure on the naturals.  The true
    // fixed point (b_{k+1} == b_k, the cut-eikonal PM chain) reproduces
    // the parent's masters at the SAME top sector — zero descent.  The
    // guard therefore fires exactly on (signature match) AND (no strict
    // width decrease): the byte-identical fixed point fires at depth 1
    // as before, and any same-propagator chain fires the moment its
    // width stops shrinking (an infinite chain of strictly decreasing
    // naturals is impossible, so nothing self-similar slips through).
    //
    // The guard is scoped to the plain recursion path: it stands down
    // when AMFLOW_FORCE_ENDING_DEPTH or AMFLOW_FIXEDPOINT_PROBE is set
    // (those runs terminate or bypass the descent by construction —
    // exactly the pipeline stages that must traverse self-similar
    // levels), and AMFLOW_ALLOW_FIXEDPOINT_RECURSION=1 restores the
    // silent recursion for debugging.
    if (!opts_.parent_prop_signature.empty()) {
        const char* fe = std::getenv("AMFLOW_FORCE_ENDING_DEPTH");
        const char* fp = std::getenv("AMFLOW_FIXEDPOINT_PROBE");
        const bool pipeline_run = (fe && *fe) || (fp && *fp);
        // Certified descent: the child's preferred masters span a
        // strictly smaller top sector than the parent's.  A missing
        // parent width (< 0) alongside a recorded signature never
        // certifies descent — fire conservatively.
        long child_top_count = 0;
        for (int s : qft::get_top_sector(preferred_)) {
            if (s) ++child_top_count;
        }
        const bool descended =
            opts_.parent_top_sector_count >= 0 &&
            child_top_count < opts_.parent_top_sector_count;
        if (!pipeline_run && !descended &&
            gravityflow_prop_signature(*fc_) == opts_.parent_prop_signature) {
            if (const char* env =
                    std::getenv("AMFLOW_ALLOW_FIXEDPOINT_RECURSION");
                env && *env && std::string(env) != "0") {
                std::cerr << "[setup] WARNING: GRAVITYFLOW_FIXEDPOINT "
                             "guard suppressed by "
                             "AMFLOW_ALLOW_FIXEDPOINT_RECURSION for "
                             "family '" << fc_->family
                          << "' — the recursion will grind to "
                             "max_recursion_depth\n";
            } else {
                std::ostringstream m;
                m << "AMFSystem::setup: GRAVITYFLOW_FIXEDPOINT: boundary "
                     "family '" << fc_->family << "' at depth "
                  << opts_.current_depth
                  << " reproduces its parent's propagator set ("
                  << fc_->propagators_after_conservation.size()
                  << " propagators identical after conservation) — the "
                     "eta->inf boundary recursion is a fixed point and "
                     "would descend to max_recursion_depth without "
                     "producing boundary data. Treat the fixed point as "
                     "the boundary condition instead: probe with "
                     "AMFLOW_FIXEDPOINT_PROBE=1, then close the boundary "
                     "system and inject it via "
                     "amf_options.explicit_boundary with "
                     "AMFLOW_FORCE_ENDING_DEPTH (the pmflow "
                     "detect/discover/respond/solve/inject pipeline). "
                     "Set AMFLOW_ALLOW_FIXEDPOINT_RECURSION=1 to restore "
                     "the silent recursion (debugging only).";
                throw std::runtime_error(m.str());
            }
        }
    }
    fc_with_eta_ = build_injected_fc(*fc_, etac_);

    // Mirrors AMFSystemDirection (AMFlow.m:981-991): pick the
    // continuation direction from the prescriptions of the loops that
    // touch eta-dependent propagators.  Mixed prescriptions abort.
    {
        const auto& fc_eta = *fc_with_eta_;
        const long eta_var = fc_eta.ctx->var_index("eta");
        std::set<long> touched_loops;
        if (eta_var >= 0) {
            const long n_total = fc_eta.ctx->n_vars();
            const long n_loop  = (long)fc_eta.n_loops();
            std::vector<unsigned long> exp((std::size_t)n_total);
            for (const auto& prop : fc_eta.propagators_after_conservation) {
                bool depends_on_eta = false;
                std::vector<bool> seen_loop((std::size_t)n_loop, false);
                const long len = fmpz_mpoly_length(prop.raw(),
                                                     fc_eta.ctx->raw());
                for (long t = 0; t < len; ++t) {
                    fmpz_mpoly_get_term_exp_ui(exp.data(), prop.raw(), t,
                                                fc_eta.ctx->raw());
                    if (exp[(std::size_t)eta_var] > 0) depends_on_eta = true;
                    for (long li = 0; li < n_loop; ++li) {
                        if (exp[(std::size_t)li] > 0) {
                            seen_loop[(std::size_t)li] = true;
                        }
                    }
                }
                if (!depends_on_eta) continue;
                for (long li = 0; li < n_loop; ++li) {
                    if (seen_loop[(std::size_t)li]) {
                        touched_loops.insert(li);
                    }
                }
            }
        }

        std::vector<int> pres_list;
        pres_list.reserve(touched_loops.size());
        for (long li : touched_loops) {
            pres_list.push_back(fc_eta.prescription_of_loop(li));
        }

        // AMFlow.m:534 PrescriptionOf:
        //   all zero -> 0      (returns 0; caller maps to NegIm)
        //   drop zeros, all 1 -> 1   (caller maps to NegIm)
        //   drop zeros, all -1 -> -1 (caller maps to Im)
        //   mixed -> $Failed   (Abort)
        bool all_zero = true;
        for (int p : pres_list) if (p != 0) { all_zero = false; break; }
        if (pres_list.empty() || all_zero) {
            direction_ = numeric::RunningOptions::Direction::NegIm;
        } else {
            int common = 0;
            bool consistent = true;
            for (int p : pres_list) {
                if (p == 0) continue;
                if (common == 0) common = p;
                else if (common != p) { consistent = false; break; }
            }
            if (!consistent) {
                throw std::runtime_error(
                    "AMFSystem::setup: cannot define a self-consistent "
                    "path direction (mixed prescriptions on eta-touching "
                    "loops; mirrors AMFlow.m:987 $Failed branch)");
            }
            direction_ = (common == -1)
                ? numeric::RunningOptions::Direction::Im
                : numeric::RunningOptions::Direction::NegIm;
        }
    }

    build_diffeq();
    build_boundary();
    dump_boundary_json();   // AMFLOW_BOUNDARY_DUMP (env-gated, read-only)
}

// ---------------------------------------------------------------------------
//  build_diffeq
// ---------------------------------------------------------------------------

void AMFSystem::build_diffeq() {
    ibp::ReduceOptions bb = opts_.bb;
    if (!opts_.cache_root.empty()) {
        bb.work_dir = opts_.cache_root + "/system_"
                     + std::to_string(system_id_) + "_diffeq";
    } else if (bb.work_dir.empty()) {
        bb.work_dir = (fs::temp_directory_path() /
            ("amflow_amfsys_" + std::to_string(system_id_) + "_diffeq")).string();
    }
    auto top_pattern = qft::get_top_sector(preferred_);
    diffeq_result_ = ibp::diffeq(*fc_with_eta_, preferred_,
                                        {"eta"}, top_pattern, bb);

    // Map preferred[i] -> sortedmasters index.
    //
    // The throw below enforces `preferred ⊆ sortedmasters` -- the same
    // invariant MMA checks at `AMFSystemCombineSolution`
    // (AMFlow.m:1220 `SubsetQ[Keys[sol], local]`).  MMA stores the
    // Kira preheat output as the system's masters and only catches
    // the violation later, after the system has been solved; we catch
    // it eagerly here so that any divergence between the parent's
    // BlackBoxReduce master list (`preferred_`) and Kira's
    // eta-injected preheat (`sortedmasters`) surfaces with a
    // diagnostic that names the offending integral.  Semantically
    // equivalent to MMA -- just early-failing.  Correct upstream
    // input never trips this; if it fires, the root cause lives in
    // the *top-level* reduce (e.g. symbolic kinematics letting Kira
    // mis-enumerate scaleless sub-sectors -- see
    // `ibp::kira_write_config` for the MMA-faithful numeric
    // substitution that prevents that class of failure).
    std::map<std::string, long> sorted_index;
    for (std::size_t i = 0; i < diffeq_result_.sortedmasters.size(); ++i) {
        sorted_index[jintegral_key(diffeq_result_.sortedmasters[i])] = (long)i;
    }
    pref_to_sorted_.clear();
    pref_to_sorted_.reserve(preferred_.size());
    for (const auto& j : preferred_) {
        auto it = sorted_index.find(jintegral_key(j));
        if (it == sorted_index.end()) {
            throw std::runtime_error(
                "AMFSystem::build_diffeq: preferred master '"
                + jintegral_key(j) + "' not in Kira's master list "
                "(mirrors MMA's AMFSystemCombineSolution SubsetQ check, "
                "early-fired; usually means the parent's top-level Kira "
                "over-enumerated masters)");
        }
        pref_to_sorted_.push_back(it->second);
    }
    AMFLOW_TRACE("AMFLOW_DEBUG_BC") {
        std::cerr << "[build_diffeq] family=" << fc_->family
                  << " preferred_to_sorted=" << pref_to_sorted_.size()
                  << std::endl;
        for (std::size_t i = 0; i < preferred_.size(); ++i) {
            const auto& pref = preferred_[i];
            const long si = pref_to_sorted_[i];
            std::cerr << "  preferred[" << i << "] = J[" << pref.family();
            for (long idx : pref.indices()) std::cerr << "," << idx;
            std::cerr << "] -> sorted[" << si << "]";
            if (si >= 0 && (std::size_t)si < diffeq_result_.sortedmasters.size()) {
                const auto& sm = diffeq_result_.sortedmasters[(std::size_t)si];
                std::cerr << " = J[" << sm.family();
                for (long idx : sm.indices()) std::cerr << "," << idx;
                std::cerr << "]";
            }
            std::cerr << std::endl;
        }
    }
}

// ---------------------------------------------------------------------------
//  build_boundary
// ---------------------------------------------------------------------------
//
//   Mirrors the inner Table loop in AMFSystemBoundaryCondition +
//   ReduceBoundary + AMFSystemsSetup (the recursive part).
//
//   For each non-zero region:
//     1. powers = RegionPower(integrals, region)
//     2. (boundary integrand families, terms) = BoundaryIntegrals(BoundaryIntegrands(...))
//     3. For each family (i.e. each sub-qft::FamilyConfig):
//          a. Build sub_fc by replacing parent's Propagator with the
//             family's prop list, keeping loops/legs/conservation/
//             replacement/numeric/cut/prescription from parent.
//          b. Collect all J integrals appearing in the family's terms.
//          c. Run BlackBoxReduce(those Js, {}) on sub_fc -> (sub_masters,
//             sub_rules).
//          d. Recursively run amf_systems_setup(sub_fc, sub_masters):
//             this either ends (vacuum) or further injects eta.
//          e. Compute the per-(master, term-order) "table" linking
//             boundary integrand terms to sub_masters via sub_rules.

// Helper: evaluate an algebra::Mfrac at ε = 10^-4 (the AMFlow convention) and
// any other vars at 0.  Returns the value as fmpq.  Used by
// build_boundary's adaptive border computation.
namespace {
// The boundary-order ANALYSIS eps: the numeric ε substituted into the
// symbolic-{eta,d,eps} objects (the DE matrix for `deinf`, the region
// `pattern` exponents for `npattern`) when determining boundary orders.
// Default 1/100 — a generic small rational, mirroring AMFlow's 10^-4
// probe (any generic value works because the SAME value is used on
// both sides of the determine_boundary_order comparison).
//
// AMFLOW_DIFFEQ_NUMD EXCEPTION: with the numeric-d substitution active
// the DE matrix arrives with d ALREADY pinned at the run's
// d = 4 - 2*eps_run, so `mfrac_substitute`'s d -> 4-2*(analysis eps)
// is a no-op on it.  If npattern kept the 1/100 frame, deinf (at
// eps_run) and npattern (at 1/100) would sit in DIFFERENT frames and
// the exact-rational matching inside determine_boundary_order breaks
// (observed: a bubble regression case built 1 system instead of 2 and
// returned value 0 — the master=26/order=1 BC-gap class).  The ONLY
// consistent analysis frame for a numeric-d run is eps_run itself:
// eps_run = (4 - P/Q)/2, which the startup assert has already verified
// equals the run's single eps sample.
void analysis_probe_eps(fmpq_t out) {
    fmpq_set_si(out, 1, 100);
    const std::string numd = ibp::diffeq_numd_env();
    if (numd.empty()) return;
    fmpq_t d_q;
    fmpq_init(d_q);
    if (fmpq_set_str(d_q, numd.c_str(), 10) == 0) {
        fmpq_canonicalise(d_q);
        fmpq_t four;
        fmpq_init(four);
        fmpq_set_si(four, 4, 1);
        fmpq_sub(d_q, four, d_q);        // 4 - d
        fmpq_div_2exp(out, d_q, 1);      // (4 - d)/2 = eps_run
        fmpq_clear(four);
    }
    fmpq_clear(d_q);
}

void mfrac_eval_at_small_eps(const algebra::Mfrac& m, fmpq_t out) {
    auto ctx = m.ctx();
    long n_vars = ctx->n_vars();
    long eps_var = ctx->var_index("eps");
    if (eps_var < 0) eps_var = ctx->var_index("__amf_eps");

    auto eval_q = [&](const algebra::Mpoly& p, fmpq_t res) {
        fmpq_t accum, term, eps_pow, base;
        fmpq_init(accum); fmpq_init(term);
        fmpq_init(eps_pow); fmpq_init(base);
        std::vector<unsigned long> exp((std::size_t)n_vars);
        fmpz_t coeff;
        fmpz_init(coeff);
        long len = fmpz_mpoly_length(p.raw(), ctx->raw());
        for (long t = 0; t < len; ++t) {
            fmpz_mpoly_get_term_exp_ui(exp.data(), p.raw(),
                                        t, ctx->raw());
            fmpz_mpoly_get_term_coeff_fmpz(coeff, p.raw(), t, ctx->raw());
            // Skip if any non-eps variable has positive exp.
            bool drop = false;
            for (long v = 0; v < n_vars; ++v) {
                if (v == eps_var) continue;
                if (exp[(std::size_t)v] != 0) { drop = true; break; }
            }
            if (drop) continue;
            fmpq_set_fmpz(term, coeff);
            if (eps_var >= 0 && exp[(std::size_t)eps_var] > 0) {
                fmpq_set_si(eps_pow, 1, 1);
                analysis_probe_eps(base);       // matches num_q (numd-aware)
                for (unsigned long k = 0;
                        k < exp[(std::size_t)eps_var]; ++k) {
                    fmpq_mul(eps_pow, eps_pow, base);
                }
                fmpq_mul(term, term, eps_pow);
            }
            fmpq_add(accum, accum, term);
        }
        fmpz_clear(coeff);
        fmpq_set(res, accum);
        fmpq_clear(accum); fmpq_clear(term);
        fmpq_clear(eps_pow); fmpq_clear(base);
    };

    fmpq_t qn, qd;
    fmpq_init(qn); fmpq_init(qd);
    eval_q(m.numerator(),   qn);
    eval_q(m.denominator(), qd);
    if (fmpq_is_zero(qd)) {
        fmpq_zero(out);
    } else {
        fmpq_div(out, qn, qd);
    }
    fmpq_clear(qn); fmpq_clear(qd);
}
}  // namespace

void AMFSystem::build_boundary() {
    auto rctx = qft::make_region_context(*fc_with_eta_);
    // Numeric substitution map for parent-only variables (kinematic
    // invariants and mass scales) that may appear in the boundary's
    // Laporta coefficient but get dropped from a sub-family's reduced
    // context when `to_complete_explicit` rank-filters linearly
    // dependent propagators.  Mirrors MMA's `/. Numeric` at
    // AMFlow.m:817 inside `ReduceBoundary`.
    auto numeric_q = build_numeric_q(opts_.bb);
    auto top_pattern = qft::get_top_sector(preferred_);
    std::vector<std::size_t> top_posi;
    for (std::size_t i = 0; i < top_pattern.size(); ++i) {
        if (top_pattern[i]) top_posi.push_back(i);
    }
    if (top_posi.empty()) return;   // no boundary needed

    auto regions = qft::find_all_region(*fc_with_eta_, rctx, top_posi);
    std::vector<qft::Region> non_zero;
    for (auto& r : regions) {
        if (!qft::zero_region_q(*fc_with_eta_, rctx, r, top_posi)) {
            non_zero.push_back(std::move(r));
        }
    }
    if (non_zero.empty()) {
        // GRAVITYFLOW_CUTREGION guard.
        //
        // On a family carrying a CUT quadratic propagator (an on-shell
        // graviton/photon line delta+((sum l - q)^2), reverse-unitarity
        // cut, ...) `zero_region_q` is cut-blind: it can classify every
        // eta->inf region as scaleless/zero even though the true
        // boundary masters are the NONZERO soft phase-space-cut
        // integrals (the delta-functions pin a scale the region test
        // never sees).  Returning silently here solves the eta-DE with
        // an empty BC, so every target comes back as a single hard
        // eps^0 = 0 with rc=0 — a FALSE ZERO (measured on a production
        // 3-loop radiation family).
        // For an UNCUT family an all-zero region set is a legitimate
        // scaleless classification, so the silent return is kept there.
        //
        // Fix for cut families: supply the cut-region boundary via
        // amf_options.explicit_boundary (+ AMFLOW_FORCE_ENDING_DEPTH) —
        // the gravityflow discover/respond/solve/inject pipeline; this
        // error string is its MISS surface.
        bool has_cut_quadratic = false;
        std::size_t first_cq = 0;
        {
            const long n_loop = (long)fc_->n_loops();
            const auto& cuts  = fc_->cut;
            const auto& props = fc_->propagators_after_conservation;
            for (std::size_t i = 0;
                 i < props.size() && i < cuts.size(); ++i) {
                if (cuts[i] == 0) continue;
                // Loop-total-degree of prop i: max over monomials of the
                // summed exponents of the loop variables.  Worldline /
                // eikonal 2u.k cut props have degree 1 (kept on their
                // fixed-point path); a cut graviton line (k1+k2-q)^2
                // has degree 2.
                const algebra::Mpoly& d = props[i];
                const long len = fmpz_mpoly_length(d.raw(),
                                                    fc_->ctx->raw());
                std::vector<unsigned long> exp(
                    (std::size_t)fc_->ctx->n_vars());
                long ldeg = 0;
                for (long t = 0; t < len; ++t) {
                    fmpz_mpoly_get_term_exp_ui(exp.data(), d.raw(), t,
                                                fc_->ctx->raw());
                    long s = 0;
                    for (long li = 0; li < n_loop; ++li) {
                        s += (long)exp[(std::size_t)li];
                    }
                    if (s > ldeg) ldeg = s;
                }
                if (ldeg >= 2) {
                    has_cut_quadratic = true;
                    first_cq = i;
                    break;
                }
            }
        }
        if (has_cut_quadratic) {
            if (const char* env =
                    std::getenv("AMFLOW_ALLOW_CUT_FALSE_ZERO");
                env && *env && std::string(env) != "0") {
                std::cerr << "[build_boundary] WARNING: "
                             "GRAVITYFLOW_CUTREGION guard suppressed by "
                             "AMFLOW_ALLOW_CUT_FALSE_ZERO for family '"
                          << fc_->family
                          << "' — targets WILL be silent false zeros\n";
                return;
            }
            std::ostringstream m;
            m << "AMFSystem::build_boundary: GRAVITYFLOW_CUTREGION: "
                 "family '" << fc_->family
              << "' carries a cut quadratic propagator (first at "
                 "position " << first_cq << ", 0-based) and all "
              << regions.size()
              << " eta->inf regions were dropped by zero_region_q — "
                 "the boundary would be silently empty and every "
                 "target would solve to a FALSE ZERO. Inject the "
                 "cut-region boundary masters (soft phase-space-cut "
                 "integrals) via amf_options.explicit_boundary "
                 "(+ AMFLOW_FORCE_ENDING_DEPTH; gravityflow "
                 "discover/respond/solve/inject). eta positions "
                 "(etac != 0):";
            bool any_eta = false;
            for (std::size_t i = 0; i < etac_.size(); ++i) {
                if (etac_[i] != 0) { m << ' ' << i; any_eta = true; }
            }
            if (!any_eta) m << " none";
            m << ". Set AMFLOW_ALLOW_CUT_FALSE_ZERO=1 to bypass "
                 "(debugging only).";
            throw std::runtime_error(m.str());
        }
        return;
    }

    auto pctx = qft::make_powers_context(*fc_with_eta_);

    // We need a BC for *every* master in the diffeq matrix, not just
    // the user-supplied preferred[].  Mathematica builds the BC over
    // sortedmasters.  We use the same set.
    const auto& diff_masters = diffeq_result_.sortedmasters;

    // Step 1 (AMFSystemBoundaryOrder, .m line 922-926):
    //   Compute powers per region, then compute pattern via
    //   BoundaryPattern.  pattern is shape [group_id][integral_id].
    std::vector<std::vector<algebra::Mfrac>> all_powers;
    all_powers.reserve(non_zero.size());
    for (const auto& region : non_zero) {
        all_powers.push_back(qft::region_power(*fc_with_eta_, rctx, pctx,
                                            region, diff_masters));
    }
    auto pattern = qft::boundary_pattern(all_powers);
    // Stash the full pattern list; solve_one_eps will use it to append
    // (mu_pattern_i -> 0) entries to bc_sorted.  Mirrors MMA's final
    //   bc = Join[#1, Thread[(pattern/.epsrule) -> 0]] step.
    bc_pattern_.clear();
    bc_pattern_.reserve(pattern.size());
    for (const auto& pg : pattern) {
        std::vector<algebra::Mfrac> row;
        row.reserve(pg.size());
        for (const auto& m : pg) row.push_back(m.clone());
        bc_pattern_.push_back(std::move(row));
    }

    // Step 2 (AMFSystemBoundaryOrder wlscript template):
    //   Compute deinf = -η^-2 (diffeq /. eta -> 1/eta), evaluated at
    //   ε = 10^-4 (the AMFlow convention).
    //
    //   Then npattern = -pattern /. epsrule (per-pattern, per-integral).
    //   Then orders[g] = DetermineBoundaryOrder[deinf, npattern[g]] for
    //   each pattern group.
    long n_masters = (long)diff_masters.size();
    if ((long)diffeq_result_.diffeq.size() != 1) {
        throw std::runtime_error(
            "AMFSystem::build_boundary: expected single-variable diffeq");
    }
    long eta_in_red = diffeq_result_.red_ctx.ctx->var_index("eta");
    if (eta_in_red < 0) {
        throw std::runtime_error(
            "AMFSystem::build_boundary: red_ctx is missing 'eta'");
    }

    // Build deinf as a numeric::RationalMatrix(η) by:
    //   (a) for each (i, j), evaluate diffeq[i][j] at ε = 10^-4 and
    //       any numeric_values
    //   (b) substitute eta -> 1/eta (as a polynomial transform)
    //   (c) multiply by -1/η²
    //
    // For step (a) we need numeric values for any non-eta vars.
    // Use opts_.bb.numeric_values + ε = 1e-4 fallback.
    // Mathematica uses ε = 10^-4 in the wlscript; we use a smaller-
    // denominator fraction to avoid exact-rational arithmetic blowup
    // (Mathematica falls back to machine-precision floats internally).
    std::map<std::string, FmpqHolder> num_q;
    {
        FmpqHolder h;
        // Analysis eps for the deinf evaluation.  MUST be the same
        // value `mfrac_eval_at_small_eps` uses for npattern — the two
        // sit on opposite sides of determine_boundary_order's exact-
        // rational matching.  With AMFLOW_DIFFEQ_NUMD active this is
        // the run eps (the DE already carries numeric d = 4-2*eps_run);
        // otherwise the stock 1/100 probe.  See analysis_probe_eps.
        analysis_probe_eps(h.raw());
        num_q.emplace("eps", h);
        FmpqHolder h2;
        fmpq_set(h2.raw(), h.raw());
        num_q.emplace("__amf_eps", h2);
    }
    auto parse_str = [](const std::string& s, fmpq_t v) -> bool {
        std::size_t slash = s.find('/');
        try {
            if (slash != std::string::npos) {
                long p = std::stol(s.substr(0, slash));
                long q = std::stol(s.substr(slash + 1));
                fmpq_set_si(v, p, q); return true;
            }
            std::size_t pos;
            long p = std::stol(s, &pos);
            if (pos == s.size()) {
                fmpq_set_si(v, p, 1); return true;
            }
            std::size_t dot = s.find('.');
            if (dot == std::string::npos) return false;
            std::string sint = s.substr(0, dot) + s.substr(dot + 1);
            long denom = 1;
            for (std::size_t i = dot + 1; i < s.size(); ++i) denom *= 10;
            long num = std::stol(sint);
            fmpq_set_si(v, num, denom); return true;
        } catch (...) { return false; }
    };
    for (const auto& [name, val_str] : opts_.bb.numeric_values) {
        if (num_q.count(name)) continue;
        FmpqHolder h;
        if (parse_str(val_str, h.raw())) num_q.emplace(name, std::move(h));
    }

    std::set<std::string> keep_eta = {"eta"};
    numeric::RationalMatrix deinf((std::size_t)n_masters, (std::size_t)n_masters);
    const auto& dmat = diffeq_result_.diffeq[0];
    for (long i = 0; i < n_masters; ++i) {
        for (long j = 0; j < n_masters; ++j) {
            // Evaluate dmat[i][j] at numeric values (keep eta).
            algebra::Mfrac sub = mfrac_substitute(dmat[(std::size_t)i][(std::size_t)j],
                                            num_q, keep_eta);
            // Convert to numeric::RationalFunction in eta.
            numeric::RationalFunction r = mfrac_to_eta_rational(sub, eta_in_red);
            deinf(i, j) = std::move(r);
        }
    }
    // Apply  eta -> 1/eta  and multiply by -1/η²:
    //   deinf'(η) = -1/η^2 * deinf(1/η)
    //
    // numeric::RationalFunction::substitute_inverse_eta implements only r(1/η);
    // we must explicitly multiply by -1/η² afterwards (mirrors AMFlow's
    // `deinf = -eta^-2 * (de /. eta -> 1/eta)`).  Missing the -1/η² factor
    // caused DetermineBoundaryOrder to see a wrong matrix (e.g. a non-
    // proper ODE with an η in the numerator instead of a 1/η pole), which
    // in turn produced wrong `border[i]` values and left `rb.integrand_terms`
    // empty for the non-zero-sector regions.  See Stage B notes.
    fmpq_t neg_one;
    fmpq_init(neg_one);
    fmpq_set_si(neg_one, -1, 1);
    numeric::RationalFunction prefactor = numeric::RationalFunction::monomial(-2, neg_one);
    fmpq_clear(neg_one);
    for (std::size_t i = 0; i < (std::size_t)n_masters; ++i) {
        for (std::size_t j = 0; j < (std::size_t)n_masters; ++j) {
            deinf(i, j) = deinf(i, j).substitute_inverse_eta() * prefactor;
        }
    }


    // Step 3: For each pattern group g, compute orders[g] (length n_masters).
    //
    //   npattern[g][i] = -pattern[g][i] /. epsrule
    //
    //   orders[g] = ode::determine_boundary_order(deinf, npattern[g])
    std::vector<std::vector<long>> orders;
    orders.reserve(pattern.size());
    for (const auto& pat : pattern) {
        if ((long)pat.size() != n_masters) {
            throw std::runtime_error(
                "AMFSystem::build_boundary: pattern size mismatch");
        }
        // Build npat_q: per-integral ini value for DetermineBoundaryOrder.
        // We DO NOT round to integer — MMA accepts rational ini and gives
        // different answers for -1 vs -99/100 (see Stage B notes).
        std::vector<numeric::RationalFunction> npat_q;
        npat_q.reserve(pat.size());
        for (const auto& m : pat) {
            fmpq_t q;
            fmpq_init(q);
            mfrac_eval_at_small_eps(m, q);
            fmpq_neg(q, q);    // npattern = -pattern /. epsrule
            npat_q.push_back(numeric::RationalFunction::from_fmpq(q));
            fmpq_clear(q);
        }

        std::vector<long> ord = ode::determine_boundary_order(deinf, npat_q);
        AMFLOW_TRACE("AMFLOW_DEBUG_BC") {
            std::cerr << "[build_boundary] family=" << fc_->family
                      << " pattern_group " << (orders.size())
                      << " npat_q=[";
            for (std::size_t i = 0; i < npat_q.size(); ++i) {
                if (i) std::cerr << ",";
                std::cerr << npat_q[i].to_string();
            }
            std::cerr << "] orders=[";
            for (std::size_t i = 0; i < ord.size(); ++i) {
                if (i) std::cerr << ",";
                std::cerr << ord[i];
            }
            std::cerr << "]" << std::endl;
        }
        orders.push_back(std::move(ord));
    }

    // Step 4: For each region, compute per-master border:
    //   k = first pattern group containing region's powers (matched by
    //       integer-shift on integral 0)
    //   border[r][i] = orders[k][i] - (pattern[k][i] - powers[r][i])
    //                 (intent: how many extra orders we need to expand
    //                  to capture all the boundary information for
    //                  master i in region r)
    //   border = Map[If[#<0, -1, #]&, border, {2}]
    //
    // We still need to find which group each region belongs to.
    // Local mfrac_is_rational: returns true iff `m` is a rational
    // constant (no algebra::Mpoly variables in num or den).  Sets `out` to the
    // rational value.
    auto mfrac_is_rat = [](const algebra::Mfrac& m, fmpq_t out) -> bool {
        algebra::Mpoly num = m.numerator();
        algebra::Mpoly den = m.denominator();
        auto ctx = m.ctx();
        std::vector<unsigned long> z((std::size_t)ctx->n_vars(), 0);
        // Check that num is purely the const term (no other monomials).
        long len_n = fmpz_mpoly_length(num.raw(), ctx->raw());
        long len_d = fmpz_mpoly_length(den.raw(), ctx->raw());
        if (len_n > 1 || len_d > 1) return false;
        if (len_n == 1) {
            std::vector<unsigned long> e((std::size_t)ctx->n_vars());
            fmpz_mpoly_get_term_exp_ui(e.data(), num.raw(), 0, ctx->raw());
            for (long v = 0; v < ctx->n_vars(); ++v) {
                if (e[(std::size_t)v] != 0) return false;
            }
        }
        if (len_d == 1) {
            std::vector<unsigned long> e((std::size_t)ctx->n_vars());
            fmpz_mpoly_get_term_exp_ui(e.data(), den.raw(), 0, ctx->raw());
            for (long v = 0; v < ctx->n_vars(); ++v) {
                if (e[(std::size_t)v] != 0) return false;
            }
        }
        fmpz_t cn, cd;
        fmpz_init(cn); fmpz_init(cd);
        fmpz_mpoly_get_coeff_fmpz_ui(cn, num.raw(), z.data(), ctx->raw());
        fmpz_mpoly_get_coeff_fmpz_ui(cd, den.raw(), z.data(), ctx->raw());
        if (fmpz_is_zero(cd)) {
            fmpz_clear(cn); fmpz_clear(cd);
            return false;
        }
        fmpq_set_fmpz_frac(out, cn, cd);
        fmpz_clear(cn); fmpz_clear(cd);
        return true;
    };

    auto find_group = [&](std::size_t r) -> long {
        const algebra::Mfrac& this_pow0 = all_powers[r][0];
        for (std::size_t g = 0; g < pattern.size(); ++g) {
            // Test IntegerQ[pattern[g][0] - powers[r][0]]
            algebra::Mfrac diff = pattern[g][0].clone();
            diff -= this_pow0;
            fmpq_t q;
            fmpq_init(q);
            bool ok = mfrac_is_rat(diff, q);
            bool is_int = ok && fmpz_is_one(fmpq_denref(q));
            fmpq_clear(q);
            if (is_int) return (long)g;
        }
        return -1;
    };

    std::vector<std::vector<long>> region_borders;
    region_borders.reserve(non_zero.size());
    for (std::size_t r = 0; r < non_zero.size(); ++r) {
        long k = find_group(r);
        if (k < 0) {
            throw std::runtime_error(
                "AMFSystem::build_boundary: region " + std::to_string(r)
                + " has no matching pattern group");
        }
        std::vector<long> b((std::size_t)n_masters, -1);
        for (long i = 0; i < n_masters; ++i) {
            // pattern[k][i] - powers[r][i] should be a rational number
            // (small integer + fractional from ε).  We round.
            algebra::Mfrac diff = pattern[(std::size_t)k][(std::size_t)i].clone();
            diff -= all_powers[r][(std::size_t)i];
            fmpq_t q;
            fmpq_init(q);
            mfrac_eval_at_small_eps(diff, q);
            // border = orders[k][i] - (pattern[k][i] - powers[r][i])
            fmpq_t order_q;
            fmpq_init(order_q);
            fmpq_set_si(order_q, orders[(std::size_t)k][(std::size_t)i], 1);
            fmpq_sub(order_q, order_q, q);
            // Round to nearest long.
            fmpz_t fz, abs_num, abs_den, two_abs_num, two_abs_den;
            fmpz_init(fz); fmpz_init(abs_num); fmpz_init(abs_den);
            fmpz_init(two_abs_num); fmpz_init(two_abs_den);
            fmpz_abs(abs_num, fmpq_numref(order_q));
            fmpz_abs(abs_den, fmpq_denref(order_q));
            fmpz_mul_si(two_abs_num, abs_num, 2);
            fmpz_add(two_abs_num, two_abs_num, abs_den);
            fmpz_mul_si(two_abs_den, abs_den, 2);
            fmpz_fdiv_q(fz, two_abs_num, two_abs_den);
            long bv = fmpz_get_si(fz);
            int bv_sgn = fmpq_sgn(order_q);
            if (bv_sgn < 0) bv = -bv;
            fmpz_clear(fz); fmpz_clear(abs_num); fmpz_clear(abs_den);
            fmpz_clear(two_abs_num); fmpz_clear(two_abs_den);
            fmpq_clear(order_q); fmpq_clear(q);
            // Map[If[#<0, -1, #]&, border, {2}]  (clamp negative to -1)
            if (bv < 0) bv = -1;
            b[(std::size_t)i] = bv;
        }
        AMFLOW_TRACE("AMFLOW_DEBUG_BC") {
            std::cerr << "[build_boundary] family=" << fc_->family
                      << " region " << r << " border=[";
            for (std::size_t i = 0; i < b.size(); ++i) {
                if (i) std::cerr << ",";
                std::cerr << b[i];
            }
            std::cerr << "]" << std::endl;
        }
        region_borders.push_back(std::move(b));
    }

    // For each region: build boundary integrand decomposition + sub-system.
    for (std::size_t reg_i = 0; reg_i < non_zero.size(); ++reg_i) {
        const auto& region = non_zero[reg_i];

        // Use the per-master border computed above.
        std::vector<long> border = region_borders[reg_i];

        qft::BoundaryIntegrandsResult bi;
        std::vector<qft::BoundaryFamily> fams;
        try {
            bi = qft::boundary_integrands(*fc_with_eta_, rctx, diff_masters,
                                       border, region);
            fams = qft::boundary_integrals(*fc_with_eta_, rctx, bi);
        } catch (const std::exception& e) {
            throw std::runtime_error(
                "AMFSystem::build_boundary: region " + std::to_string(reg_i)
                + ": " + e.what());
        }
        AMFLOW_TRACE("AMFLOW_DEBUG_BC") {
            std::cerr << "[build_boundary] family=" << fc_->family
                      << " reg_i=" << reg_i
                      << " border=[";
            for (std::size_t i = 0; i < border.size(); ++i) {
                if (i) std::cerr << ",";
                std::cerr << border[i];
            }
            std::cerr << "] fams.size=" << fams.size() << std::endl;
            for (std::size_t fi = 0; fi < fams.size(); ++fi) {
                std::cerr << "  fam[" << fi << "] prop.size="
                          << fams[fi].prop.size()
                          << " terms.size=" << fams[fi].terms.size()
                          << std::endl;
                for (std::size_t ti = 0; ti < fams[fi].terms.size(); ++ti) {
                    std::size_t nonempty = 0, total = 0;
                    for (const auto& tk : fams[fi].terms[ti]) {
                        total += tk.size();
                        if (!tk.empty()) nonempty++;
                    }
                    std::cerr << "    terms[" << ti << "] #ks="
                              << fams[fi].terms[ti].size()
                              << " nonempty_ks=" << nonempty
                              << " total_lts=" << total << std::endl;
                }
            }
        }

        auto make_region_boundary = [&]() {
            RegionBoundary rb;
            rb.powers.reserve(diff_masters.size());
            for (auto& m : all_powers[reg_i]) rb.powers.push_back(m.clone());

            rb.integrand_terms.resize(diff_masters.size());
            for (std::size_t i = 0; i < diff_masters.size(); ++i) {
                for (long oi = 0; oi <= border[i]; ++oi) {
                    rb.integrand_terms[i].emplace_back();
                }
            }
            return rb;
        };

        // Mathematica's ReduceBoundary returns one entry per boundary family:
        //   {powers, config, masters, table}
        // AMFSystemBoundaryCondition then flattens all entries across regions.
        // Store each family as its own RegionBoundary to keep that outer
        // multiplicity distinct from SingleMass multi-root systems.
        //
        // Upstream `ReduceBoundary` (AMFlow.m:790-803) projects the parent's
        // `Cut` propagators (after applying the region's bare LoopTransform)
        // onto the boundary sub-family's propagator basis and emits the
        // corresponding `cut_propagators` to Kira:
        //
        //   cutde = Pick[ReducedPropagator /. region[[1]], Cut, 1];
        //   cut   = If[MemberQ[Expand[# - cutde] /. ReducedReplacement, 0],
        //              1, 0]& /@ prop;
        //   If[Count[cut, 1] =!= Count[Cut, 1], Abort[]];
        //
        // The parent `ReducedPropagator` here is the *eta-injected* version
        // (the AMFSystemSetup config), so any cut prop that received an η
        // injection (etac[k] != 0) carries an `eta` term.  After the bare
        // LoopTransform the η term survives, and no eta-free fam.prop entry
        // can match it under reduced_replacement → the count assertion
        // fires, abort.  This is the upstream-intended behaviour ("eta may
        // have been inserted to cut denominators").
        std::vector<algebra::Mfrac> cutde_in_fc_eta;
        long parent_cut_count = 0;
        bool parent_has_cut = false;
        for (int c : fc_with_eta_->cut) if (c == 1) ++parent_cut_count;
        if (parent_cut_count > 0) {
            parent_has_cut = true;
            cutde_in_fc_eta.reserve((std::size_t)parent_cut_count);
            for (std::size_t k = 0; k < fc_with_eta_->cut.size(); ++k) {
                if (fc_with_eta_->cut[k] != 1) continue;
                const auto& parent_prop =
                    fc_with_eta_->propagators_after_conservation[k];
                // Apply the *bare* transform (region.transform.map), NOT
                // the scaled rule from `region_rule(...)` — upstream uses
                // `ReducedPropagator /. region[[1]]` where `region[[1]]` is
                // just the LoopTransform without the half_eta scaling.
                algebra::Mfrac transformed = qft::apply_region_rule(
                    *fc_with_eta_, rctx, parent_prop, region.transform.map);
                // The bare transform only references rctx loops/legs/
                // invariants/eta (not __amf_eta, __amf_half_eta), so the
                // result projects cleanly back to fc_with_eta_->ctx.
                algebra::Mfrac in_eta_ctx =
                    project_mfrac_by_name(transformed, fc_with_eta_->ctx);
                cutde_in_fc_eta.push_back(std::move(in_eta_ctx));
            }
        }

        for (auto& fam : fams) {
            // 1. Build sub-fc.
            std::vector<std::string> prop_strs;
            prop_strs.reserve(fam.prop.size());
            for (auto& p : fam.prop) prop_strs.push_back(p.to_string());

            // 1a. Project parent's `cut` onto fam.prop (mirrors
            //     AMFlow.m:800-801).  An `i` in the new propagator list
            //     is "cut" iff some transformed parent cut prop matches it
            //     modulo `reduced_replacement`.  If the matched count
            //     differs from the parent count, η was injected into a
            //     cut denominator → abort.
            std::vector<int> sub_cut(fam.prop.size(), 0);
            if (parent_has_cut) {
                long sub_cut_count = 0;
                for (std::size_t i = 0; i < fam.prop.size(); ++i) {
                    algebra::Mfrac fp_mfrac =
                        algebra::Mfrac::from_mpoly(fam.prop[i].clone());
                    for (std::size_t j = 0; j < cutde_in_fc_eta.size(); ++j) {
                        algebra::Mfrac diff = fp_mfrac.clone();
                        diff -= cutde_in_fc_eta[j].clone();
                        algebra::Mfrac simplified =
                            fc_with_eta_->apply_replacement(diff);
                        if (simplified.is_zero()) {
                            sub_cut[i] = 1;
                            ++sub_cut_count;
                            break;
                        }
                    }
                }
                if (sub_cut_count != parent_cut_count) {
                    std::ostringstream m;
                    m << "AMFSystem::build_boundary: parent cut propagator "
                         "count mismatch in boundary sub-family (matched "
                      << sub_cut_count << ", expected "
                      << parent_cut_count
                      << ").  Mirrors AMFlow.m:801: 'eta may have been "
                         "inserted to cut denominators'.  This typically "
                         "means the chosen η-injection point landed on a "
                         "cut propagator; pick a different `EndingScheme` "
                         "or rearrange the family so the cut props are "
                         "η-free.";
                    throw std::runtime_error(m.str());
                }
            }

            // The sub-fc inherits everything from parent except prop + cut.
            // We give it a unique family name; users don't observe it.
            std::string sub_fam_name = fc_->family + "_b"
                + std::to_string(system_id_) + "_r"
                + std::to_string(reg_i) + "_f"
                + std::to_string(subsystems_.size());

            qft::FamilyConfig sub_fc = qft::FamilyConfig::build(
                sub_fam_name,
                fc_->loops, fc_->legs,
                extract_conservation_strings(*fc_),
                extract_replacement_strings(*fc_),
                prop_strs,
                std::move(sub_cut),
                fc_->prescription);

            // 2. Collect all J integrals appearing in this family's terms.
            //    Translate from the parent's family name to sub_fam_name.
            std::set<std::string> seen_js;
            std::vector<qft::JIntegral> all_js;
            for (std::size_t i = 0; i < fam.terms.size(); ++i) {
                for (std::size_t k = 0; k < fam.terms[i].size(); ++k) {
                    for (const auto& lt : fam.terms[i][k]) {
                        qft::JIntegral j(sub_fam_name, lt.indices);
                        std::string key = jintegral_key(j);
                        if (seen_js.insert(key).second) all_js.push_back(j);
                    }
                }
            }

            AMFLOW_TRACE("AMFLOW_DEBUG_BC") {
                std::cerr << "[build_boundary] parent=" << fc_->family
                          << " reg_i=" << reg_i
                          << " all_js.size=" << all_js.size() << std::endl;
                for (const auto& j : all_js) {
                    std::cerr << "    j=J[" << j.family();
                    for (auto x : j.indices()) std::cerr << "," << x;
                    std::cerr << "]" << std::endl;
                }
            }
            if (all_js.empty()) continue;

            // 3. Run BlackBoxReduce(all_js, {}) on sub_fc to find
            //    sub-masters + reduction rules.
            //
            //    The Kira call lets it choose its own masters (passing
            //    {} for preferred).
            ibp::ReduceOptions reduce_opts = opts_.bb;
            // Mark this as boundary-recursion staging so ibp::reduce()
            // honors AMFLOW_BOUNDARY_PREHEAT_DOT/RANK caps on the
            // Masters-mode preheat seed (the boundary-staging fix
            // for a production 4.57M-mandatory bad_alloc; see
            // ReduceOptions::boundary_stage).  Env unset => no-op.
            reduce_opts.boundary_stage = true;
            reduce_opts.work_dir = (fs::temp_directory_path() /
                ("amflow_amfsys_" + std::to_string(system_id_)
                  + "_b" + std::to_string(reg_i)
                  + "_f" + std::to_string(subsystems_.size())
                  + "_p" + std::to_string(getpid())
                  + "_u"
                  + std::to_string(
                        g_boundary_reduce_unique.fetch_add(1,
                            std::memory_order_relaxed)))).string();
            fs::remove_all(reduce_opts.work_dir);

            std::vector<int> top_pat = qft::get_top_sector(all_js);
            ibp::ReduceResult red;
            try {
                // stage-2 endings: denominator-corner bias on
                // eikonal-vacuum boundary families (empty otherwise).
                red = ibp::reduce(sub_fc, all_js,
                                  build_eikonal_corner_preferred(sub_fc,
                                                                 all_js),
                                          top_pat, reduce_opts);
            } catch (const std::exception& e) {
                throw std::runtime_error(
                    "AMFSystem::build_boundary: BlackBoxReduce failed: "
                    + std::string(e.what()));
            }

            if (red.masters.empty()) {
                // Boundary integrand reduces to zero -- skip family.
                continue;
            }

            RegionBoundary rb = make_region_boundary();

            // 4. Recursively setup the sub-system.  amf_systems_setup
            //    decides ending vs. inject based on red.masters.
            AMFSystemOptions sub_opts = opts_;
            sub_opts.current_depth = opts_.current_depth + 1;
            // Record this family's propagator signature so the child's
            // setup() can fire the GRAVITYFLOW_FIXEDPOINT guard when the
            // boundary family reproduces it (self-similar recursion),
            // together with this system's top-sector width (top_posi is
            // the active-slot set of get_top_sector(preferred_), computed
            // at the head of build_boundary): the guard stands down only
            // for a child whose own top-sector width is strictly smaller
            // — certified descent — and fires on the width-preserving
            // reproduction that is the true fixed point.
            sub_opts.parent_prop_signature =
                gravityflow_prop_signature(*fc_);
            sub_opts.parent_top_sector_count = (long)top_posi.size();
            AMFLOW_TRACE("AMFLOW_DEBUG_BC") {
                std::cerr << "[build_boundary] parent=" << fc_->family
                          << " region=" << reg_i
                          << " sub_fam=" << sub_fam_name
                          << " red.masters.size=" << red.masters.size()
                          << std::endl;
                for (const auto& m : red.masters) {
                    std::cerr << "  master J[" << m.family();
                    for (auto x : m.indices()) std::cerr << "," << x;
                    std::cerr << "]" << std::endl;
                }
                std::cerr << "  sub_fc loops=" << sub_fc.loops.size()
                          << " props=" << sub_fc.propagators_after_conservation.size()
                          << std::endl;
            }
            // FIXED-POINT PROBE: skip the recursive sub-system setup and
            // record the sub-family master basis instead.  solve_one_eps
            // will source the sub-value vector from fp_injected_.  Only at
            // depth 0 (the probe runs the self-similar family AS the top
            // system).  Multiple regions/factors must share ONE basis —
            // mismatch means the family is not a clean fixed point.
            const char* fp_env = std::getenv("AMFLOW_FIXEDPOINT_PROBE");
            if (fp_env && *fp_env && opts_.current_depth == 0) {
                fp_probe_active_ = true;
                if (fp_sub_masters_.empty()) {
                    for (const auto& m : red.masters) fp_sub_masters_.push_back(m);
                    std::cerr << "[amflow][fp-probe] recursion skipped at "
                              << sub_fam_name << "; sub-basis = "
                              << red.masters.size() << " masters\n";
                } else {
                    if (fp_sub_masters_.size() != red.masters.size()) {
                        throw std::runtime_error(
                            "fixedpoint_probe: multiple boundary families with "
                            "DIFFERENT master bases (" +
                            std::to_string(fp_sub_masters_.size()) + " vs " +
                            std::to_string(red.masters.size()) +
                            ") — not a clean fixed point");
                    }
                    for (std::size_t m = 0; m < red.masters.size(); ++m) {
                        if (jintegral_key(fp_sub_masters_[m]) !=
                            jintegral_key(red.masters[m])) {
                            throw std::runtime_error(
                                "fixedpoint_probe: boundary-family master basis "
                                "mismatch at index " + std::to_string(m));
                        }
                    }
                }
                // rb.sub_system_indices stays EMPTY — probe path keys off it.
            } else {
            auto sub_systems = amf_systems_setup(sub_fc, red.masters, sub_opts);

            // Register every root system belonging to this boundary family.
            rb.sub_system_indices.reserve(sub_systems.size());
            for (auto& sub_sys : sub_systems) {
                long sub_idx = (long)subsystems_.size();
                subsystems_.push_back(std::move(sub_sys));
                rb.sub_system_indices.push_back(sub_idx);
            }
            }

            // 5. For each (master, integrand_term), look up its
            //    sub-master coefficient via red.rules.
            //
            // red.rules[k] = list of ibp::DerivTerm{coef, integ} for all_js[k].
            // We build a key->index map for all_js, then use that to
            // resolve LaportaTerm indices; multiply LaportaTerm.coef
            // (in fc_with_eta_.ctx) by red.rules's per-master coef
            // (in red.red_ctx); accumulate into rb.integrand_terms.

            std::map<std::string, std::size_t> all_js_index;
            for (std::size_t k = 0; k < all_js.size(); ++k) {
                all_js_index[jintegral_key(all_js[k])] = k;
            }
            std::map<std::string, std::size_t> sub_master_idx;
            for (std::size_t m = 0; m < red.masters.size(); ++m) {
                sub_master_idx[jintegral_key(red.masters[m])] = m;
            }

            // After the 2024-12 fix to boundary_integrals, fam.terms is
            // shaped [input_i][order_k][laporta_term], with input_i
            // ranging over diff_masters and order_k ranging over
            // border[input_i]+1.  This means fam.terms.size() ==
            // rb.integrand_terms.size() and fam.terms[i].size() ==
            // rb.integrand_terms[i].size().  We assert to catch any
            // regression.
            if (fam.terms.size() != rb.integrand_terms.size()) {
                throw std::runtime_error(
                    "AMFSystem::build_boundary: fam.terms.size "
                    + std::to_string(fam.terms.size()) + " != "
                    + "rb.integrand_terms.size "
                    + std::to_string(rb.integrand_terms.size())
                    + " — boundary_integrals layout drift; see Stage B notes.");
            }
            for (std::size_t i = 0; i < fam.terms.size(); ++i) {
                if (fam.terms[i].size() != rb.integrand_terms[i].size()) {
                    throw std::runtime_error(
                        "AMFSystem::build_boundary: fam.terms[" +
                        std::to_string(i) + "].size " +
                        std::to_string(fam.terms[i].size()) + " != " +
                        "rb.integrand_terms[" + std::to_string(i) +
                        "].size " + std::to_string(rb.integrand_terms[i].size()));
                }
                for (std::size_t k = 0; k < fam.terms[i].size(); ++k) {
                    for (const auto& lt : fam.terms[i][k]) {
                        qft::JIntegral j(sub_fam_name, lt.indices);
                        std::string jkey = jintegral_key(j);
                        auto idx_it = all_js_index.find(jkey);
                        if (idx_it == all_js_index.end()) continue;

                        const auto& rule = red.rules[idx_it->second];
                        // rule is per-master: each entry is (coef, sub_J).
                        // For each entry, accumulate
                        //   contribution = lt.coef (project to fc_->ctx)
                        //                * rule_term.coef (already on red.red_ctx)
                        //                * sub_solution[sub_idx]
                        // We multiply lt.coef and rule_term.coef
                        // symbolically (after projecting both into a
                        // common ctx = sub-system's red_ctx).

                        // lt.coef is on bi.dctx (= parent fc + DListSymbol
                        // + eta + half_eta).  Project it down to
                        // fc_with_eta_.ctx (drop __amf_*).
                        algebra::Mfrac lt_coef_in_fc;
                        try {
                            lt_coef_in_fc = project_mfrac_by_name(
                                lt.coeff, fc_with_eta_->ctx);
                        } catch (const std::exception& e) {
                            throw std::runtime_error(
                                "AMFSystem::build_boundary: lt.coef projection: "
                                + std::string(e.what()));
                        }

                        // When `to_complete_explicit` rank-filters linearly
                        // dependent propagators of the boundary sub-family
                        // (e.g. two `l1^2 - m_i sq` props with the same
                        // SP signature), the sub-family's reduced context
                        // loses one of the mass scales.  `lt.coeff` may
                        // still carry that mass symbolically.  Substitute
                        // its numeric value here, mirroring MMA's `/.
                        // Numeric` at AMFlow.m:817 inside `ReduceBoundary`,
                        // so the projection to `red.red_ctx.ctx` below
                        // never sees a residual non-droppable variable.
                        std::set<std::string> sub_red_keep_names;
                        for (long v = 0; v < red.red_ctx.ctx->n_vars(); ++v) {
                            sub_red_keep_names.insert(
                                red.red_ctx.ctx->var_name(v));
                        }
                        try {
                            lt_coef_in_fc = mfrac_substitute(
                                lt_coef_in_fc, numeric_q,
                                sub_red_keep_names);
                        } catch (const std::exception& e) {
                            throw std::runtime_error(
                                "AMFSystem::build_boundary: lt.coef numeric "
                                "substitution for sub-red projection: "
                                + std::string(e.what()));
                        }

                        for (const auto& rt : rule) {
                            auto mit = sub_master_idx.find(
                                jintegral_key(rt.integ));
                            if (mit == sub_master_idx.end()) continue;
                            // rt.coef is on red.red_ctx (= sub_fc.ctx + 'd').
                            // Project lt_coef_in_fc to that ctx.
                            algebra::Mfrac lt_in_sub;
                            try {
                                lt_in_sub = project_mfrac_by_name(
                                    lt_coef_in_fc, red.red_ctx.ctx);
                            } catch (const std::exception& e) {
                                throw std::runtime_error(
                                    "AMFSystem::build_boundary: lt-to-sub-red projection: "
                                    + std::string(e.what()));
                            }
                            algebra::Mfrac product = lt_in_sub.clone();
                            product *= rt.coef;
                            BoundaryTerm bt;
                            bt.coef = std::move(product);
                            bt.sub_master_idx = (long)mit->second;
                            bt.sub_master_key = jintegral_key(rt.integ);
                            rb.integrand_terms[i][k].push_back(std::move(bt));
                        }
                    }
                }
            }

            regions_.push_back(std::move(rb));
        }
    }
}

// ---------------------------------------------------------------------------
//  dump_boundary_json (AMFLOW_BOUNDARY_DUMP, env-gated, read-only)
// ---------------------------------------------------------------------------

void AMFSystem::dump_boundary_json() const {
    const std::string& dir = boundary_dump_dir();
    if (dir.empty() || is_ending_) return;
    try {
        {
            std::error_code ec;
            fs::create_directories(dir, ec);
        }
        const std::string path = dir + "/system_"
            + std::to_string(system_id_) + "_boundary.json";
        std::ofstream os(path);
        if (!os) {
            std::cerr << "[amflow][boundary-dump] WARNING: cannot open '"
                      << path << "' -- symbolic dump skipped\n";
            return;
        }
        const auto& sm = diffeq_result_.sortedmasters;

        os << "{\n";
        os << "  \"system_id\": " << system_id_ << ",\n";
        os << "  \"family\": \"" << json_escape(fc_->family) << "\",\n";
        os << "  \"scheme\": \"" << ending_scheme_name(ending_) << "\",\n";
        os << "  \"etac\": [";
        for (std::size_t i = 0; i < etac_.size(); ++i) {
            if (i) os << ", ";
            os << etac_[i];
        }
        os << "],\n";
        os << "  \"grammar\": {\n"
              "    \"expression_printer\": \"FLINT fmpz_mpoly_get_str_pretty:"
              " integer coefficients, '*' for products, '^' for integer"
              " powers, terms joined by ' + '/' - '; an Mfrac prints as NUM"
              " when its denominator is 1, else (NUM)/(DEN)\",\n"
              "    \"mu_convention\": \"mu strings are exact Mfracs on the"
              " powers context whose regulator variable is '__amf_eps' =="
              " eps (d = 4 - 2*eps); when mu is linear, mu_a/mu_b are exact"
              " rationals with mu = mu_a + mu_b*eps\",\n"
              "    \"coef_convention\": \"each term coef is an exact Mfrac"
              " on the boundary sub-family reduction context (variables"
              " listed in the region's coef_ctx_vars, normally ['eta','d']);"
              " parent kinematic invariants were substituted numerically at"
              " build time; at solve time d -> 4 - 2*eps and eps -> the run"
              " eps rational, and the bc value for (master i, order k) is"
              " sum_terms coef(eps) * sub_value[sub_master_key] attached at"
              " mu_i - k\",\n"
              "    \"sm0_prefactor\": \"prefactor = (-1)^n * Gamma(2-eps-m)"
              " * Gamma(-2+eps+m+n) / (Gamma(2-eps)*Gamma(n)) with"
              " m = m0 + m_eps_coef*eps; global_value[i] = prefactor(gamma_"
              "params[i]) * master_values[jpre_idx[i]], multiplied over all"
              " roots of a multi-root child\"\n"
              "  },\n";

        // sortedmasters: the eta-DE row order (diffeq_result_.sortedmasters),
        // plus the Kira masters.final order (preferred_) and the map.
        os << "  \"sortedmasters\": [";
        for (std::size_t i = 0; i < sm.size(); ++i) {
            if (i) os << ", ";
            os << "\"" << json_escape(jintegral_key(sm[i])) << "\"";
        }
        os << "],\n";
        os << "  \"preferred\": [";
        for (std::size_t i = 0; i < preferred_.size(); ++i) {
            if (i) os << ", ";
            os << "\"" << json_escape(jintegral_key(preferred_[i])) << "\"";
        }
        os << "],\n";
        os << "  \"pref_to_sorted\": [";
        for (std::size_t i = 0; i < pref_to_sorted_.size(); ++i) {
            if (i) os << ", ";
            os << pref_to_sorted_[i];
        }
        os << "],\n";

        // Subsystems.
        os << "  \"subsystems\": [\n";
        for (std::size_t si = 0; si < subsystems_.size(); ++si) {
            const AMFSystem& sub = *subsystems_[si];
            os << "    {\"index\": " << si
               << ", \"system_id\": " << sub.system_id_
               << ", \"family\": \"" << json_escape(sub.fc_->family) << "\""
               << ", \"scheme\": \"" << ending_scheme_name(sub.ending_) << "\""
               << ", \"kind\": \""
               << (sub.is_ending_ ? "ending" : "recursion") << "\""
               << ", \"sm0\": " << (sub.sm_cfg_.has_value() ? "true" : "false");
            os << ", \"preferred\": [";
            for (std::size_t i = 0; i < sub.preferred_.size(); ++i) {
                if (i) os << ", ";
                os << "\"" << json_escape(jintegral_key(sub.preferred_[i]))
                   << "\"";
            }
            os << "]";
            if (sub.sm_cfg_.has_value()) {
                const auto& cfg = *sub.sm_cfg_;
                os << ", \"single_mass\": {\"original_preferred\": [";
                for (std::size_t i = 0; i < cfg.original_preferred.size(); ++i) {
                    if (i) os << ", ";
                    os << "\""
                       << json_escape(jintegral_key(cfg.original_preferred[i]))
                       << "\"";
                }
                os << "], \"jpre_idx\": [";
                for (std::size_t i = 0; i < cfg.jpre_idx.size(); ++i) {
                    if (i) os << ", ";
                    os << cfg.jpre_idx[i];
                }
                os << "], \"gamma_params\": [";
                for (std::size_t i = 0; i < cfg.gamma_params.size(); ++i) {
                    const auto& gp = cfg.gamma_params[i];
                    if (i) os << ", ";
                    os << "{\"n\": " << gp.n << ", \"m0\": " << gp.m0
                       << ", \"m_eps_coef\": " << gp.m_eps_coef
                       << ", \"loop_num\": " << gp.loop_num << "}";
                }
                os << "]}";
            }
            if (sub.ck_cfg_.has_value()) {
                os << ", \"cutkosky\": {\"phase_loop_num\": "
                   << sub.ck_cfg_->phase_loop_num << "}";
            }
            if (sub.is_ending_) {
                os << ", \"ending_masters\": {";
                for (std::size_t i = 0; i < sub.preferred_.size(); ++i) {
                    if (i) os << ", ";
                    os << "\"" << json_escape(jintegral_key(sub.preferred_[i]))
                       << "\": "
                       << classify_ending_master_json(*sub.fc_,
                                                      sub.preferred_[i],
                                                      sub.opts_);
                }
                os << "}";
            }
            os << "}" << (si + 1 < subsystems_.size() ? "," : "") << "\n";
        }
        os << "  ],\n";

        // Regions.  One entry per (asy region x boundary family) pair, in
        // the exact order solve_one_eps walks regions_.
        os << "  \"regions\": [\n";
        for (std::size_t r = 0; r < regions_.size(); ++r) {
            const RegionBoundary& rb = regions_[r];
            os << "    {\"region\": " << r << ", \"sub_system_indices\": [";
            for (std::size_t i = 0; i < rb.sub_system_indices.size(); ++i) {
                if (i) os << ", ";
                os << rb.sub_system_indices[i];
            }
            os << "], \"root_multiplicity\": " << rb.sub_system_indices.size();
            // Variables of the ctx the term coefs live on (first term found).
            os << ", \"coef_ctx_vars\": [";
            {
                const algebra::Mfrac* first_coef = nullptr;
                for (const auto& per_master : rb.integrand_terms) {
                    for (const auto& per_k : per_master) {
                        if (!per_k.empty()) { first_coef = &per_k[0].coef; break; }
                    }
                    if (first_coef) break;
                }
                if (first_coef && first_coef->ctx()) {
                    auto cctx = first_coef->ctx();
                    for (long v = 0; v < cctx->n_vars(); ++v) {
                        if (v) os << ", ";
                        os << "\"" << json_escape(cctx->var_name(v)) << "\"";
                    }
                }
            }
            os << "],\n     \"masters\": [\n";
            for (std::size_t i = 0; i < sm.size(); ++i) {
                os << "      {\"i\": " << i << ", \"master\": \""
                   << json_escape(jintegral_key(sm[i])) << "\", ";
                bdump_emit_mu(os, rb.powers[i]);
                const long border =
                    (long)rb.integrand_terms[i].size() - 1;
                os << ", \"border\": " << border << ", \"orders\": [";
                for (std::size_t k = 0; k < rb.integrand_terms[i].size(); ++k) {
                    if (k) os << ", ";
                    os << "[";
                    const auto& terms = rb.integrand_terms[i][k];
                    for (std::size_t t = 0; t < terms.size(); ++t) {
                        if (t) os << ", ";
                        os << "{\"coef\": \""
                           << json_escape(mfrac_str(terms[t].coef))
                           << "\", \"sub_master_key\": \""
                           << json_escape(terms[t].sub_master_key)
                           << "\", \"sub_master_idx\": "
                           << terms[t].sub_master_idx << "}";
                    }
                    os << "]";
                }
                os << "]}" << (i + 1 < sm.size() ? "," : "") << "\n";
            }
            os << "     ]}" << (r + 1 < regions_.size() ? "," : "") << "\n";
        }
        os << "  ],\n";

        // Pattern-zero groups: per group, one mu per sortedmaster; each
        // (mu -> 0) entry is appended to the bc at solve time.
        os << "  \"pattern_zero_groups\": [\n";
        for (std::size_t g = 0; g < bc_pattern_.size(); ++g) {
            os << "    [";
            for (std::size_t i = 0; i < bc_pattern_[g].size(); ++i) {
                if (i) os << ", ";
                os << "{";
                bdump_emit_mu(os, bc_pattern_[g][i]);
                os << "}";
            }
            os << "]" << (g + 1 < bc_pattern_.size() ? "," : "") << "\n";
        }
        os << "  ]\n";
        os << "}\n";
        os.flush();
        std::cerr << "[amflow][boundary-dump] wrote symbolic boundary tables"
                  << " (system_id=" << system_id_
                  << ", family=" << fc_->family
                  << ", sortedmasters=" << sm.size()
                  << ", regions=" << regions_.size()
                  << ", pattern_groups=" << bc_pattern_.size()
                  << ") -> " << path << "\n";
    } catch (const std::exception& e) {
        std::cerr << "[amflow][boundary-dump] WARNING: symbolic dump failed"
                  << " for system_id=" << system_id_ << ": " << e.what()
                  << " -- continuing\n";
    } catch (...) {
        std::cerr << "[amflow][boundary-dump] WARNING: symbolic dump failed"
                  << " for system_id=" << system_id_ << " -- continuing\n";
    }
}

// ===========================================================================
//  Solve
// ===========================================================================

// Compute the SingleMass prefactor numerically:
//
//   prefactor = (-1)^(-n) * Gamma(2 - eps - m) * Gamma(-2 + eps + m + n) /
//               (Gamma(2 - eps) * Gamma(n))
//
// where m = m0 + m_eps_coef * eps.
namespace {
void compute_single_mass_prefactor(const AMFSystem::GammaParams& gp,
                                      const acb_t eps,
                                      acb_t out,
                                      long prec) {
    acb_t m, two_minus_eps, gamma1, gamma2, gamma3, gamma4;
    acb_t arg, sign;
    acb_init(m); acb_init(two_minus_eps);
    acb_init(gamma1); acb_init(gamma2);
    acb_init(gamma3); acb_init(gamma4);
    acb_init(arg); acb_init(sign);

    // m = m0 + m_eps_coef * eps
    acb_set_si(m, gp.m0);
    if (gp.m_eps_coef != 0) {
        acb_t tmp; acb_init(tmp);
        acb_set_si(tmp, gp.m_eps_coef);
        acb_mul(tmp, tmp, eps, prec);
        acb_add(m, m, tmp, prec);
        acb_clear(tmp);
    }

    // Gamma(2 - eps - m)
    acb_set_si(arg, 2);
    acb_sub(arg, arg, eps, prec);
    acb_sub(arg, arg, m, prec);
    acb_gamma(gamma1, arg, prec);

    // Gamma(-2 + eps + m + n)
    acb_set_si(arg, -2 + (long)gp.n);
    acb_add(arg, arg, eps, prec);
    acb_add(arg, arg, m, prec);
    acb_gamma(gamma2, arg, prec);

    // Gamma(2 - eps)
    acb_set_si(arg, 2);
    acb_sub(arg, arg, eps, prec);
    acb_gamma(gamma3, arg, prec);

    // Gamma(n)
    acb_set_si(arg, gp.n);
    acb_gamma(gamma4, arg, prec);

    // sign = (-1)^(-n) = (-1)^n
    acb_one(sign);
    if (gp.n % 2 != 0) acb_neg(sign, sign);

    // out = sign * gamma1 * gamma2 / (gamma3 * gamma4)
    acb_mul(out, sign, gamma1, prec);
    acb_mul(out, out, gamma2, prec);
    acb_div(out, out, gamma3, prec);
    acb_div(out, out, gamma4, prec);

    acb_clear(m); acb_clear(two_minus_eps);
    acb_clear(gamma1); acb_clear(gamma2);
    acb_clear(gamma3); acb_clear(gamma4);
    acb_clear(arg); acb_clear(sign);
}

void compute_cutkosky_prefactor(long phase_loop_num,
                                  const acb_t eps,
                                  acb_t out,
                                  long prec) {
    if (phase_loop_num < 0) {
        throw std::runtime_error(
            "compute_cutkosky_prefactor: phase loop number is negative");
    }

    acb_t pi, two_pi, exp1, exp2, four, pi_pow, two_pi_pow, base;
    acb_init(pi); acb_init(two_pi);
    acb_init(exp1); acb_init(exp2); acb_init(four);
    acb_init(pi_pow); acb_init(two_pi_pow);
    acb_init(base);

    acb_const_pi(pi, prec);
    acb_mul_2exp_si(two_pi, pi, 1);

    // Pi^(2 - eps)
    acb_set_si(exp1, 2);
    acb_sub(exp1, exp1, eps, prec);
    acb_pow(pi_pow, pi, exp1, prec);

    // (2 Pi)^(2 eps - 4)
    acb_mul_2exp_si(exp2, eps, 1);
    acb_set_si(four, 4);
    acb_sub(exp2, exp2, four, prec);
    acb_pow(two_pi_pow, two_pi, exp2, prec);

    acb_mul(base, pi_pow, two_pi_pow, prec);
    acb_pow_ui(out, base, static_cast<unsigned long>(phase_loop_num), prec);
    acb_mul_2exp_si(out, out, 1);
    if ((phase_loop_num + 1) % 2 != 0) acb_neg(out, out);

    acb_clear(pi); acb_clear(two_pi);
    acb_clear(exp1); acb_clear(exp2); acb_clear(four);
    acb_clear(pi_pow); acb_clear(two_pi_pow);
    acb_clear(base);
}
}  // namespace

void AMFSystem::solve(const std::vector<numeric::AcbValue>& epslist) {
    // ------------------------------------------------------------------
    // CHECKPOINT / RESUME (per completed system DE solve).
    //
    // A node's solution fully incorporates its sub-systems, so a HIT here
    // lets us SKIP the recursion into sub-systems entirely.  The key is a
    // content hash of this node's identity + the ε-grid + precision, so a
    // re-run with the same config resumes past every already-completed
    // system; any config change re-keys and re-solves.  The ibp_cache
    // (which persists across crashes) already protects the Kira reduction;
    // this protects the numeric solve.  Disabled unless a root resolves.
    // ------------------------------------------------------------------
    const std::string ckpt_root = system_checkpoint_resolve_root(
        ibp::ibp_cache_resolve_root(opts_.bb.ibp_cache_dir));
    std::string ckpt_key;
    if (!ckpt_root.empty()) {
        std::vector<std::string> pref_keys;
        pref_keys.reserve(preferred_.size());
        for (const auto& j : preferred_) pref_keys.push_back(jintegral_key(j));
        ckpt_key = system_checkpoint_key(
            fc_->family, pref_keys, etac_, epslist,
            numeric::working_prec_bits(), static_cast<int>(direction_),
            is_ending_);
        std::vector<AMFSystemSolution> loaded;
        if (system_checkpoint_load(ckpt_root, ckpt_key, loaded)
            && loaded.size() == epslist.size()) {
            solutions_ = std::move(loaded);
            return;  // sub-systems are subsumed in this node's solution
        }
    }

    solutions_.clear();
    solutions_.resize(epslist.size());
    std::map<std::string, ibp::ReduceResult> ending_reduction_cache;

    if (is_ending_) {
        for (std::size_t k = 0; k < epslist.size(); ++k) {
            std::vector<numeric::AcbValue> values;
            values.reserve(preferred_.size());
            auto numeric_q =
                build_numeric_q_for_eps(opts_.bb, epslist[k], numeric::working_prec_bits());
            std::map<std::string, numeric::AcbValue> ending_value_cache;
            std::set<std::string> ending_active_keys;
            for (const auto& j : preferred_) {
                values.push_back(solve_ending_master_value(
                    *fc_, j, system_id_, opts_, k, epslist[k], numeric_q,
                    ending_reduction_cache, ending_value_cache,
                    ending_active_keys, numeric::working_prec_bits()));
            }
            solutions_[k].master_values = std::move(values);
        }
        // NOTE: do NOT return here -- we still need to apply SingleMass
        // prefactor below if sm_cfg_ is set (an ending SingleMass system
        // computes preferred = jpre, multiplies by Γ-prefactor, stores
        // result in global_values for the parent's boundary integrals).
    } else {
        // Non-ending: solve sub-systems first, then per-ε solve.
        //
        // ------------------------------------------------------------------
        // LEVER 2: node-tree parallelism.
        //
        // The sub-system solves below are mutually independent — each
        // writes only its own `solutions_` (and its own per-node
        // checkpoint key), so sibling sub-trees are dispatched to a
        // thread pool exactly the way the ε-loop is below.  Default OFF
        // (node_parallel = 1): the serial loop is byte-identical to the
        // pre-parallel code path.  Lever 1 already saturates the cores
        // for typical ε-grid sizes ~10-30; this lever is for wide
        // shallow trees / short ε-grids — see `tree_stats()` for the
        // sibling-width diagnostic that quantifies the headroom.
        //
        // Concurrency audit — ending sub-systems may invoke Kira via
        // `get_ending_reduction`; that is thread-safe here because
        //   (a) all reduce/cache coordination is file-based and already
        //       multi-process hardened — the ibp_cache uses a per-key
        //       flock stampede lock plus write-then-rename publication,
        //       and flock excludes per open file description, so two
        //       THREADS of one process (each lock opens its own fd)
        //       exclude each other exactly like two processes;
        //   (b) kira_run forks and execs immediately, with the chdir /
        //       setpgid / rlimit calls in the CHILD — no process-global
        //       mutation in the parent, and the child env is built
        //       pre-fork and passed via execve;
        //   (c) each ending reduce works in a directory keyed by
        //       (system_id, integral key), each node's
        //       `ending_reduction_cache` is local to its own solve()
        //       call, and checkpoint temp files are keyed per node —
        //       distinct siblings never share a path;
        //   (d) `numeric::*_options()` are thread-local; each worker
        //       installs the parent snapshot before recursing, so nested
        //       ε-grid pools inside a sub-tree inherit correctly.
        // Note the ε-grid parallelism's "only one node runs its ε-loop at any instant"
        // invariant is relaxed while this knob is on — the caller sizes
        // node_parallel × eps_parallel (see AMFSystemOptions).
        // ------------------------------------------------------------------
        const long n_node = std::min<long>(resolve_node_parallel(opts_),
                                           (long)subsystems_.size());
        if (n_node <= 1) {
            for (auto& sub : subsystems_) {
                sub->solve(epslist);
            }
        } else {
            const auto snap = numeric::capture_options();
            std::atomic<std::size_t> next_sub{0};
            std::exception_ptr first_err = nullptr;
            std::mutex err_mu;

            auto worker = [&]() {
                numeric::install_options(snap);
                for (;;) {
                    const std::size_t i =
                        next_sub.fetch_add(1, std::memory_order_relaxed);
                    if (i >= subsystems_.size()) break;
                    try {
                        subsystems_[i]->solve(epslist);
                    } catch (...) {
                        std::lock_guard<std::mutex> lk(err_mu);
                        if (!first_err) first_err = std::current_exception();
                        // Drain remaining work so other threads exit quickly.
                        next_sub.store(subsystems_.size(),
                                       std::memory_order_relaxed);
                        return;
                    }
                }
                flint_cleanup();
            };

            std::vector<std::thread> pool;
            pool.reserve((std::size_t)n_node);
            for (long t = 0; t < n_node; ++t) pool.emplace_back(worker);
            for (auto& t : pool) t.join();
            if (first_err) std::rethrow_exception(first_err);
        }

        // ------------------------------------------------------------------
        // LEVER 1: ε-grid parallelism.
        //
        // After setup() every Kira/IBP call has fired (or hit the
        // ibp_cache).  `solve_one_eps` is pure FLINT/Arb numerics:
        // substitute ε into the symbolic DE, read pre-computed sub-system
        // values, run `ode::amflow` (calc_inf/calc_taylor/calc_zero).  No
        // fork(), no filesystem writes, no shared mutable state except
        // `numeric::*_options()`, which is now thread-local with explicit
        // snapshot inheritance.  Distinct ε write to distinct
        // `solutions_[k]` slots that are pre-sized above.
        //
        // The recursion is depth-first: at any instant only ONE node is
        // running its ε-loop, so threading it here uses ≤ eps_parallel
        // workers regardless of tree depth — no oversubscription.
        // ------------------------------------------------------------------
        const long n_workers =
            std::min<long>(resolve_eps_parallel(opts_), (long)epslist.size());

        if (n_workers <= 1) {
            for (std::size_t k = 0; k < epslist.size(); ++k) {
                solutions_[k] = solve_one_eps(k, epslist[k]);
            }
        } else {
            const auto snap = numeric::capture_options();
            std::atomic<std::size_t> next{0};
            std::exception_ptr first_err = nullptr;
            std::mutex err_mu;

            auto worker = [&]() {
                // Inherit the parent thread's working_pre / x_order /
                // run_direction etc. (thread_local singletons start at
                // package defaults otherwise).
                numeric::install_options(snap);
                for (;;) {
                    const std::size_t k =
                        next.fetch_add(1, std::memory_order_relaxed);
                    if (k >= epslist.size()) break;
                    try {
                        solutions_[k] = solve_one_eps(k, epslist[k]);
                    } catch (...) {
                        std::lock_guard<std::mutex> lk(err_mu);
                        if (!first_err) first_err = std::current_exception();
                        // Drain remaining work so other threads exit quickly.
                        next.store(epslist.size(), std::memory_order_relaxed);
                        return;
                    }
                }
                flint_cleanup();
            };

            std::vector<std::thread> pool;
            pool.reserve((std::size_t)n_workers);
            for (long t = 0; t < n_workers; ++t) pool.emplace_back(worker);
            for (auto& t : pool) t.join();
            if (first_err) std::rethrow_exception(first_err);
        }
    }

    // Apply SingleMass prefactor to produce global_values.
    if (sm_cfg_.has_value()) {
        long prec = numeric::working_prec_bits();
        const auto& cfg = *sm_cfg_;
        for (std::size_t k = 0; k < epslist.size(); ++k) {
            auto& sol = solutions_[k];
            sol.global_preferred = cfg.original_preferred;
            sol.global_values.clear();
            sol.global_values.reserve(cfg.original_preferred.size());
            for (std::size_t i = 0; i < cfg.original_preferred.size(); ++i) {
                long jpre_i = cfg.jpre_idx[i];
                if (jpre_i < 0 || (std::size_t)jpre_i >= sol.master_values.size()) {
                    throw std::runtime_error(
                        "AMFSystem::solve: jpre_idx out of range");
                }
                numeric::AcbValue prefactor;
                compute_single_mass_prefactor(cfg.gamma_params[i],
                                                 epslist[k].raw(),
                                                 prefactor.raw(), prec);
                numeric::AcbValue gv;
                acb_mul(gv.raw(), prefactor.raw(),
                          sol.master_values[(std::size_t)jpre_i].raw(), prec);
                sol.global_values.push_back(std::move(gv));
            }
            AMFLOW_TRACE("AMFLOW_DEBUG_BC") {
                std::cerr << "[single_mass_solve] family=" << fc_->family
                          << " eps_index=" << k
                          << " global_values.size=" << sol.global_values.size()
                          << std::endl;
                for (std::size_t i = 0; i < sol.global_values.size(); ++i) {
                    const auto& jp = sol.global_preferred[i];
                    std::cerr << "  J[" << jp.family();
                    for (long idx : jp.indices()) std::cerr << "," << idx;
                    std::cerr << "] = "
                              << sol.global_values[i].to_string(20) << std::endl;
                }
            }
        }
    }

    if (ck_cfg_.has_value()) {
        long prec = numeric::working_prec_bits();
        const auto& cfg = *ck_cfg_;
        for (std::size_t k = 0; k < epslist.size(); ++k) {
            auto& sol = solutions_[k];
            if (sol.master_values.size() < cfg.original_preferred.size()) {
                throw std::runtime_error(
                    "AMFSystem::solve: Cutkosky master value count mismatch");
            }

            numeric::AcbValue prefactor;
            compute_cutkosky_prefactor(cfg.phase_loop_num,
                                         epslist[k].raw(),
                                         prefactor.raw(), prec);

            sol.global_preferred = cfg.original_preferred;
            sol.global_values.clear();
            sol.global_values.reserve(cfg.original_preferred.size());
            for (std::size_t i = 0; i < cfg.original_preferred.size(); ++i) {
                numeric::AcbValue im_value;
                arb_set(acb_realref(im_value.raw()),
                        acb_imagref(sol.master_values[i].raw()));
                arb_zero(acb_imagref(im_value.raw()));

                numeric::AcbValue gv;
                acb_mul(gv.raw(), prefactor.raw(), im_value.raw(), prec);
                sol.global_values.push_back(std::move(gv));
            }
        }
    }

    // Persist this completed system's solution so a crashed/restarted flow
    // resumes past it.  Key was computed at entry; empty ⇒ checkpointing
    // disabled.  Stored only after all of master_values + SingleMass /
    // Cutkosky global_* are populated above, so a HIT restores a complete
    // node.
    if (!ckpt_key.empty()) {
        system_checkpoint_store(ckpt_root, ckpt_key, solutions_);
    }
}

// ---------------------------------------------------------------------------
//  solve_one_eps
// ---------------------------------------------------------------------------

AMFSystemSolution
AMFSystem::solve_one_eps(std::size_t eps_index, const numeric::AcbValue& eps) {
    long prec = numeric::working_prec_bits();
    AMFSystemSolution out;

    // ---- Build numeric_q (ε + family numeric values) ----
    auto numeric_q = build_numeric_q_for_eps(opts_.bb, eps, prec);

    // AMFLOW_BOUNDARY_DUMP (env-gated, read-only): when set, snapshot the
    // per-region sub-master values and the final bc_sorted for this eps.
    const bool bdump_on = !boundary_dump_dir().empty();
    std::vector<std::string> bdump_region_vals;

    // ---- Step 1.  Build numeric::RationalMatrix(η) by evaluating diffeq. ----
    if (diffeq_result_.diffeq.size() != 1) {
        throw std::runtime_error(
            "AMFSystem::solve_one_eps: expected exactly one diffeq variable");
    }
    const auto& dmat = diffeq_result_.diffeq[0];
    std::size_t Nm = diffeq_result_.sortedmasters.size();
    long eta_in_red = diffeq_result_.red_ctx.ctx->var_index("eta");
    if (eta_in_red < 0) {
        throw std::runtime_error(
            "AMFSystem::solve_one_eps: red_ctx is missing 'eta'");
    }
    std::set<std::string> keep_eta = {"eta"};

    numeric::RationalMatrix de((std::size_t)Nm, (std::size_t)Nm);
    for (std::size_t i = 0; i < Nm; ++i) {
        for (std::size_t j = 0; j < Nm; ++j) {
            try {
                algebra::Mfrac sub = mfrac_substitute(dmat[i][j], numeric_q, keep_eta);
                de(i, j) = mfrac_to_eta_rational(sub, eta_in_red);
            } catch (const std::exception& e) {
                throw std::runtime_error(
                    "AMFSystem::solve_one_eps: diffeq("
                    + std::to_string(i) + "," + std::to_string(j) + "): "
                    + e.what());
            }
        }
    }

    // ---- Step 2.  Build ode::BoundarySpec for each sortedmaster. ----
    //
    //   For each region:
    //     For each sortedmaster i:
    //       evaluate mu(ε), then for each integrand-order k:
    //         evaluate Σ_term (coef(ε) * sub_solution_at_eps),
    //         append (μ - k, value) to BC[i].
    std::vector<ode::BoundarySpec> bc_sorted(Nm);
    struct CombinedSubValues {
        std::vector<numeric::AcbValue> values;
        std::map<std::string, std::size_t> key_to_index;
    };
    auto combined_sub_values =
        [&](const std::vector<long>& sub_indices) -> CombinedSubValues {
            CombinedSubValues out_vals;
            if (sub_indices.empty()) return out_vals;

            const auto& first_sys = *subsystems_[(std::size_t)sub_indices.front()];
            if (eps_index >= first_sys.solutions().size()) {
                throw std::runtime_error(
                    "AMFSystem::solve_one_eps: sub-system not solved at this ε");
            }
            const auto& first_sol = first_sys.solutions()[eps_index];

            // SingleMass multi-root children expose the physically relevant
            // boundary-family masters via global_values/global_preferred, and
            // the final answer is the product over roots.  This mirrors
            // Mathematica's AMFSystemsSolution / AMFSystemCombineSolution.
            if (!first_sol.global_values.empty()) {
                std::vector<numeric::AcbValue> combined;
                combined.reserve(first_sol.global_values.size());
                for (std::size_t i = 0; i < first_sol.global_values.size(); ++i) {
                    numeric::AcbValue one;
                    acb_one(one.raw());
                    combined.push_back(std::move(one));
                }
                for (long sub_idx : sub_indices) {
                    const auto& sub_sys = *subsystems_[(std::size_t)sub_idx];
                    if (eps_index >= sub_sys.solutions().size()) {
                        throw std::runtime_error(
                            "AMFSystem::solve_one_eps: sub-system not solved at this ε");
                    }
                    const auto& sub_sol = sub_sys.solutions()[eps_index];
                    if (sub_sol.global_values.size() != combined.size()) {
                        throw std::runtime_error(
                            "AMFSystem::solve_one_eps: multi-root sub-system "
                            "global_values size mismatch");
                    }
                    for (std::size_t i = 0; i < combined.size(); ++i) {
                        numeric::AcbValue tmp;
                        acb_mul(tmp.raw(), combined[i].raw(),
                                sub_sol.global_values[i].raw(), prec);
                        combined[i] = std::move(tmp);
                    }
                }
                for (std::size_t i = 0; i < combined.size(); ++i) {
                    combined[i] = midpoint_only(combined[i]);
                    if (i < first_sol.global_preferred.size()) {
                        out_vals.key_to_index[jintegral_key(first_sol.global_preferred[i])] = i;
                    }
                }
                out_vals.values = std::move(combined);
                return out_vals;
            }

            if (sub_indices.size() != 1) {
                throw std::runtime_error(
                    "AMFSystem::solve_one_eps: multi-root sub-system without "
                    "global_values is unsupported");
            }

            out_vals.values.reserve(first_sol.master_values.size());
            for (const auto& v : first_sol.master_values) {
                out_vals.values.push_back(midpoint_only(v));
            }
            for (std::size_t i = 0; i < first_sys.preferred().size(); ++i) {
                out_vals.key_to_index[jintegral_key(first_sys.preferred()[i])] = i;
            }
            return out_vals;
        };

    for (const auto& rb : regions_) {
        CombinedSubValues sub_vals;
        if (rb.sub_system_indices.empty()) {
            // FIXED-POINT PROBE path: regions carry no sub-systems (the
            // recursion was skipped in build_boundary); the sub-value
            // vector is the injected one on the recorded fp basis.
            if (!(fp_probe_active_ && fp_injected_)) continue;
            if (fp_injected_->size() != fp_sub_masters_.size()) {
                throw std::runtime_error(
                    "fixedpoint_probe: injected vector size mismatch");
            }
            sub_vals.values.reserve(fp_injected_->size());
            for (const auto& v : *fp_injected_) {
                sub_vals.values.push_back(midpoint_only(v));
            }
            for (std::size_t m = 0; m < fp_sub_masters_.size(); ++m) {
                sub_vals.key_to_index[jintegral_key(fp_sub_masters_[m])] = m;
            }
        } else {
            sub_vals = combined_sub_values(rb.sub_system_indices);
        }
        AMFLOW_TRACE("AMFLOW_DEBUG_BC") {
            std::cerr << "[solve_one_eps] family=" << fc_->family
                      << " sub_roots=[";
            for (std::size_t ridx = 0; ridx < rb.sub_system_indices.size(); ++ridx) {
                if (ridx) std::cerr << ",";
                const auto& sub_sys = *subsystems_[(std::size_t)rb.sub_system_indices[ridx]];
                std::cerr << sub_sys.family().family;
            }
            std::cerr << "] combined_sub_vals.size=" << sub_vals.values.size()
                      << std::endl;
            for (std::size_t z = 0; z < sub_vals.values.size(); ++z) {
                std::cerr << "    sub_vals[" << z << "] = "
                          << sub_vals.values[z].to_string(200) << std::endl;
            }
            std::cerr << "  rb.powers.size=" << rb.powers.size()
                      << " rb.integrand_terms.size=" << rb.integrand_terms.size()
                      << std::endl;
            for (std::size_t i = 0; i < Nm; ++i) {
                std::cerr << "  i=" << i << " integrand_terms[i].size="
                          << rb.integrand_terms[i].size() << std::endl;
                for (std::size_t k = 0; k < rb.integrand_terms[i].size(); ++k) {
                    std::cerr << "    k=" << k << " terms.size="
                              << rb.integrand_terms[i][k].size() << std::endl;
                    for (const auto& bt : rb.integrand_terms[i][k]) {
                        std::cerr << "      sub_master_idx=" << bt.sub_master_idx
                                  << " coef=" << bt.coef.to_string() << std::endl;
                    }
                }
            }
        }

        for (std::size_t i = 0; i < Nm; ++i) {
            // mu(ε)
            numeric::AcbValue mu_acb;
            try {
                evaluate_const_mfrac(rb.powers[i], numeric_q, prec, mu_acb);
            } catch (const std::exception& e) {
                throw std::runtime_error(
                    "AMFSystem::solve_one_eps: mu eval: " + std::string(e.what()));
            }
            for (std::size_t k = 0; k < rb.integrand_terms[i].size(); ++k) {
                const auto& terms = rb.integrand_terms[i][k];
                if (terms.empty()) continue;
                numeric::AcbValue value;
                acb_zero(value.raw());
                for (const auto& bt : terms) {
                    numeric::AcbValue coef_acb;
                    try {
                        evaluate_const_mfrac(bt.coef, numeric_q, prec, coef_acb);
                    } catch (const std::exception& e) {
                        throw std::runtime_error(
                            "AMFSystem::solve_one_eps: coef eval: "
                            + std::string(e.what()));
                    }
                    std::size_t sub_idx = (std::size_t)bt.sub_master_idx;
                    auto key_it = sub_vals.key_to_index.find(bt.sub_master_key);
                    if (key_it != sub_vals.key_to_index.end()) {
                        sub_idx = key_it->second;
                    }
                    if (sub_idx >= sub_vals.values.size()) {
                        throw std::runtime_error(
                            "AMFSystem::solve_one_eps: sub-master lookup out of range");
                    }
                    numeric::AcbValue prod;
                    acb_mul(prod.raw(), coef_acb.raw(),
                             sub_vals.values[sub_idx].raw(), prec);
                    acb_add(value.raw(), value.raw(), prod.raw(), prec);
                }
                ode::BoundaryEntry be;
                if (k == 0) {
                    acb_set(be.mu.raw(), mu_acb.raw());
                } else {
                    numeric::AcbValue offset;
                    acb_set_si(offset.raw(), -(long)k);
                    acb_add(be.mu.raw(), mu_acb.raw(), offset.raw(), prec);
                }
                acb_set(be.value.raw(), value.raw());
                bc_sorted[i].push_back(std::move(be));
            }
        }
    }

    // ---- Step 2.5. Append pattern-zero BC entries. ----
    //
    // Mirrors MMA's AMFSystemSolution final line:
    //   bc = MapThread[Join[#1,#2]&,
    //                  {bc, Transpose[Thread /@ Thread[(pattern/.epsrule) -> 0]]}];
    //
    // For every pattern group and every sortedmaster i, add a trivial
    // (mu = pattern[g][i] | eps, value = 0) boundary entry.  Without this,
    // a sub-system whose regions all have border=-1 ends up with an empty
    // bc_sorted and amflow produces 0 -- see Stage B notes.
    for (const auto& pg : bc_pattern_) {
        for (std::size_t i = 0; i < Nm && i < pg.size(); ++i) {
            numeric::AcbValue mu_acb;
            try {
                evaluate_const_mfrac(pg[i], numeric_q, prec, mu_acb);
            } catch (const std::exception& e) {
                throw std::runtime_error(
                    "AMFSystem::solve_one_eps: pattern mu eval: "
                    + std::string(e.what()));
            }
            ode::BoundaryEntry be;
            acb_set(be.mu.raw(), mu_acb.raw());
            acb_zero(be.value.raw());
            bc_sorted[i].push_back(std::move(be));
        }
    }

    // ---- Step 2.6. Finiteness gate on the assembled boundary. ----
    //
    // A non-finite (NaN/Inf) boundary value here means a leaf boundary
    // constant diverged at the requested epsilon -- the canonical case being a
    // SingleMass / vacuum-type leaf whose ending reduction carries a 1/eps
    // pole at eps=0 (its assembled integral is eps^0-finite only after
    // subsector recombination).  Caught here, this produces a clear, actionable
    // error naming the system; left to propagate, the NaN reaches calcx00 and
    // resurfaces as the misleading "rationalized boundary matching ... interval
    // too wide" abort, which no working-precision increase can ever clear.
    for (std::size_t i = 0; i < bc_sorted.size(); ++i) {
        for (const auto& be : bc_sorted[i]) {
            if (!acb_is_finite(be.value.raw())) {
                std::ostringstream msg;
                msg << "AMFSystem::solve_one_eps: non-finite boundary value at "
                       "system_id=" << system_id_
                    << "; family=" << fc_->family
                    << "; eps_index=" << eps_index
                    << "; sorted_master=" << i
                    << ".  A leaf boundary constant is divergent at this "
                       "epsilon (typically a 1/eps pole at eps=0); evaluating "
                       "the target at eps=0 exactly is ill-posed.  Use a small "
                       "nonzero-eps grid and extrapolate to eps=0 (e.g. a "
                       "base-2 grid eps = 2^-4 .. 2^-18) instead of "
                       "eps_samples=[0].";
                throw std::runtime_error(msg.str());
            }
        }
    }

    // ---- Step 2.7 (AMFLOW_BOUNDARY_DUMP, env-gated, read-only). ----
    //
    // Dump exactly the bc_sorted handed to ode::amflow below, plus the
    // sub-master values consumed per region.  Distinct eps indices write
    // distinct files, so this is safe under AMFLOW_EPS_PARALLEL.  Write
    // failures warn and continue; they never abort the solve.
    if (bdump_on) {
        try {
            const std::string& dir = boundary_dump_dir();
            {
                std::error_code ec;
                fs::create_directories(dir, ec);
            }
            const std::string path = dir + "/system_"
                + std::to_string(system_id_) + "_bc_eps"
                + std::to_string(eps_index) + ".json";
            std::ofstream os(path);
            if (!os) {
                std::cerr << "[amflow][boundary-dump] WARNING: cannot open '"
                          << path << "' -- bc dump skipped\n";
            } else {
                os << "{\n";
                os << "  \"system_id\": " << system_id_ << ",\n";
                os << "  \"family\": \"" << json_escape(fc_->family)
                   << "\",\n";
                os << "  \"eps_index\": " << eps_index << ",\n";
                {
                    auto eq = numeric_q.find("eps");
                    if (eq != numeric_q.end()) {
                        // The exact rational actually substituted for eps
                        // (acb_to_fmpq_rational of the sample).
                        os << "  \"eps\": \""
                           << bdump_fmpq_str(eq->second.raw()) << "\",\n";
                    }
                }
                os << "  \"eps_acb\": [\""
                   << bdump_arb_str40(acb_realref(eps.raw())) << "\", \""
                   << bdump_arb_str40(acb_imagref(eps.raw())) << "\"],\n";
                os << "  \"bc\": [\n";
                for (std::size_t i = 0; i < bc_sorted.size(); ++i) {
                    os << "    [";
                    for (std::size_t e = 0; e < bc_sorted[i].size(); ++e) {
                        const auto& be = bc_sorted[i][e];
                        if (e) os << ", ";
                        os << "[\""
                           << bdump_arb_str40(acb_realref(be.mu.raw()))
                           << "\", \""
                           << bdump_arb_str40(acb_realref(be.value.raw()))
                           << "\", \""
                           << bdump_arb_str40(acb_imagref(be.value.raw()))
                           << "\"]";
                    }
                    os << "]" << (i + 1 < bc_sorted.size() ? "," : "")
                       << "\n";
                }
                os << "  ],\n";
                os << "  \"sub_values\": [\n";
                for (std::size_t r = 0; r < bdump_region_vals.size(); ++r) {
                    os << bdump_region_vals[r]
                       << (r + 1 < bdump_region_vals.size() ? "," : "")
                       << "\n";
                }
                os << "  ]\n";
                os << "}\n";
                os.flush();
            }
        } catch (const std::exception& e) {
            std::cerr << "[amflow][boundary-dump] WARNING: bc dump failed"
                      << " (system_id=" << system_id_ << ", eps_index="
                      << eps_index << "): " << e.what()
                      << " -- continuing\n";
        } catch (...) {
            std::cerr << "[amflow][boundary-dump] WARNING: bc dump failed"
                      << " (system_id=" << system_id_ << ", eps_index="
                      << eps_index << ") -- continuing\n";
        }
    }

    // ---- Step 3. amflow ----
    AMFLOW_TRACE("AMFLOW_DEBUG_BC") {
        std::cerr << "[solve_one_eps] family=" << fc_->family
                  << " sortedmasters=" << diffeq_result_.sortedmasters.size()
                  << " boundary_entries=" << regions_.size() << std::endl;
        for (std::size_t i = 0; i < diffeq_result_.sortedmasters.size(); ++i) {
            const auto& m = diffeq_result_.sortedmasters[i];
            std::cerr << "  smaster[" << i << "] = J[" << m.family();
            for (auto x : m.indices()) std::cerr << "," << x;
            std::cerr << "]" << std::endl;
        }
        for (std::size_t i = 0; i < bc_sorted.size(); ++i) {
            std::cerr << "  bc_sorted[" << i << "]: " << bc_sorted[i].size()
                      << " entries" << std::endl;
            for (const auto& be : bc_sorted[i]) {
                std::cerr << "    mu=" << be.mu.to_string(200)
                          << "  value=" << be.value.to_string(200) << std::endl;
            }
        }
    }
    AMFLOW_TRACE("AMFLOW_DEBUG_BC") {
        std::cerr << "[solve_one_eps] family=" << fc_->family
                  << " de matrix " << de.rows() << "x" << de.cols() << ":\n";
        for (std::size_t i = 0; i < de.rows(); ++i) {
            for (std::size_t j = 0; j < de.cols(); ++j) {
                std::cerr << "    de[" << i << "," << j << "] = "
                          << de(i, j).to_string() << std::endl;
            }
        }
    }
    std::vector<numeric::AcbValue> sol_acb;
    try {
        // Per-system path direction (AMFlow.m:981-991).  We override
        // numeric::run_direction() for the duration of this call so
        // ode::amflow → ode::run_eta picks up the system-specific
        // contour orientation derived from the η-touching loops'
        // prescriptions.
        numeric::GlobalScope dir_scope;
        dir_scope.running.run_direction = direction_;
        dir_scope.commit();

        sol_acb = ode::amflow(de, bc_sorted, prec);
    } catch (const std::exception& e) {
        std::ostringstream msg;
        msg << "AMFSystem::solve_one_eps: amflow failed"
            << "; system_id=" << system_id_
            << "; family=" << fc_->family
            << "; eps_index=" << eps_index
            << "; masters=" << Nm
            << ": " << e.what();
        throw std::runtime_error(msg.str());
    }
    AMFLOW_TRACE("AMFLOW_DEBUG_BC") {
        std::cerr << "[solve_one_eps] family=" << fc_->family
                  << " sol_acb (sorted) values:" << std::endl;
        for (std::size_t i = 0; i < sol_acb.size(); ++i) {
            std::cerr << "    sol[" << i << "] = "
                      << sol_acb[i].to_string(200) << std::endl;
        }
    }

    out.master_values.reserve(preferred_.size());
    for (std::size_t i = 0; i < preferred_.size(); ++i) {
        long si = pref_to_sorted_[i];
        out.master_values.push_back(sol_acb[(std::size_t)si].clone());
    }
    AMFLOW_TRACE("AMFLOW_DEBUG_BC") {
        std::cerr << "[solve_one_eps] family=" << fc_->family
                  << " preferred_values=" << out.master_values.size()
                  << std::endl;
        for (std::size_t i = 0; i < preferred_.size(); ++i) {
            const auto& pref = preferred_[i];
            std::cerr << "  preferred[" << i << "] = J[" << pref.family();
            for (long idx : pref.indices()) std::cerr << "," << idx;
            std::cerr << "] value="
                      << out.master_values[i].to_string(200) << std::endl;
        }
    }
    return out;
}

// ===========================================================================
//  fp_probe_solve_eps_batched  (fp-probe column batching)
// ===========================================================================
//
//  AMFLOW_FP_PROBE_BATCH=1: instead of one full eta-transport per injected
//  unit column (64 x 12 = 768 on the PX2 tops), build the numeric eta-DE
//  and the boundary-entry mu/coef structure ONCE per eps and push ALL unit
//  columns through ode::amflow_multi (one fundamental-solution transport:
//  shared contour, shared CalcInf region plans, shared CalcRun
//  re-expansions, ONE NormalizeMat/ToFuchsian per eps).  The per-column
//  value arithmetic below mirrors solve_one_eps' operation sequence exactly
//  (coef * {1,0} products and in-order accumulation), so the batched A(eps)
//  columns are bit-identical to the unbatched ones.

std::vector<std::vector<numeric::AcbValue>>
AMFSystem::fp_probe_solve_eps_batched(std::size_t eps_index,
                                      const numeric::AcbValue& eps) {
    long prec = numeric::working_prec_bits();
    if (!fp_probe_active_) {
        throw std::runtime_error(
            "fp_probe_solve_eps_batched: probe not active");
    }
    const std::size_t ncols = fp_sub_masters_.size();

    // ---- numeric_q + eta-DE: built ONCE per eps (solve_one_eps step 1) ----
    auto numeric_q = build_numeric_q_for_eps(opts_.bb, eps, prec);

    if (diffeq_result_.diffeq.size() != 1) {
        throw std::runtime_error(
            "fp_probe_solve_eps_batched: expected exactly one diffeq variable");
    }
    const auto& dmat = diffeq_result_.diffeq[0];
    std::size_t Nm = diffeq_result_.sortedmasters.size();
    long eta_in_red = diffeq_result_.red_ctx.ctx->var_index("eta");
    if (eta_in_red < 0) {
        throw std::runtime_error(
            "fp_probe_solve_eps_batched: red_ctx is missing 'eta'");
    }
    std::set<std::string> keep_eta = {"eta"};

    numeric::RationalMatrix de((std::size_t)Nm, (std::size_t)Nm);
    for (std::size_t i = 0; i < Nm; ++i) {
        for (std::size_t j = 0; j < Nm; ++j) {
            try {
                algebra::Mfrac sub = mfrac_substitute(dmat[i][j], numeric_q, keep_eta);
                de(i, j) = mfrac_to_eta_rational(sub, eta_in_red);
            } catch (const std::exception& e) {
                throw std::runtime_error(
                    "fp_probe_solve_eps_batched: diffeq("
                    + std::to_string(i) + "," + std::to_string(j) + "): "
                    + e.what());
            }
        }
    }

    // ---- boundary entries: mu/coef evaluated ONCE, values scattered over
    //      the unit columns with the same {coef*1, coef*0} op sequence as
    //      solve_one_eps ----
    std::vector<std::vector<ode::BoundarySpec>> bc_cols(ncols);
    for (auto& col : bc_cols) col.resize(Nm);

    std::map<std::string, std::size_t> key_to_index;
    for (std::size_t m = 0; m < ncols; ++m) {
        key_to_index[jintegral_key(fp_sub_masters_[m])] = m;
    }

    numeric::AcbValue one_acb;
    acb_one(one_acb.raw());

    for (const auto& rb : regions_) {
        if (!rb.sub_system_indices.empty()) {
            throw std::runtime_error(
                "fp_probe_solve_eps_batched: region with recursive "
                "sub-systems on the probe path — unsupported");
        }
        for (std::size_t i = 0; i < Nm; ++i) {
            numeric::AcbValue mu_acb;
            try {
                evaluate_const_mfrac(rb.powers[i], numeric_q, prec, mu_acb);
            } catch (const std::exception& e) {
                throw std::runtime_error(
                    "fp_probe_solve_eps_batched: mu eval: "
                    + std::string(e.what()));
            }
            for (std::size_t k = 0; k < rb.integrand_terms[i].size(); ++k) {
                const auto& terms = rb.integrand_terms[i][k];
                if (terms.empty()) continue;

                std::vector<numeric::AcbValue> coefs;
                std::vector<std::size_t>       idxs;
                coefs.reserve(terms.size());
                idxs.reserve(terms.size());
                for (const auto& bt : terms) {
                    numeric::AcbValue coef_acb;
                    try {
                        evaluate_const_mfrac(bt.coef, numeric_q, prec, coef_acb);
                    } catch (const std::exception& e) {
                        throw std::runtime_error(
                            "fp_probe_solve_eps_batched: coef eval: "
                            + std::string(e.what()));
                    }
                    std::size_t sub_idx = (std::size_t)bt.sub_master_idx;
                    auto key_it = key_to_index.find(bt.sub_master_key);
                    if (key_it != key_to_index.end()) {
                        sub_idx = key_it->second;
                    }
                    if (sub_idx >= ncols) {
                        throw std::runtime_error(
                            "fp_probe_solve_eps_batched: sub-master lookup "
                            "out of range");
                    }
                    coefs.push_back(std::move(coef_acb));
                    idxs.push_back(sub_idx);
                }

                numeric::AcbValue mu_entry;
                if (k == 0) {
                    acb_set(mu_entry.raw(), mu_acb.raw());
                } else {
                    numeric::AcbValue offset;
                    acb_set_si(offset.raw(), -(long)k);
                    acb_add(mu_entry.raw(), mu_acb.raw(), offset.raw(), prec);
                }

                for (std::size_t c = 0; c < ncols; ++c) {
                    numeric::AcbValue value;
                    acb_zero(value.raw());
                    for (std::size_t t = 0; t < coefs.size(); ++t) {
                        numeric::AcbValue prod;
                        if (idxs[t] == c) {
                            // coef * midpoint_only(1)  (exact one)
                            acb_mul(prod.raw(), coefs[t].raw(),
                                    one_acb.raw(), prec);
                        } else {
                            // coef * midpoint_only(0) == exact zero
                            acb_zero(prod.raw());
                        }
                        acb_add(value.raw(), value.raw(), prod.raw(), prec);
                    }
                    ode::BoundaryEntry be;
                    acb_set(be.mu.raw(), mu_entry.raw());
                    acb_set(be.value.raw(), value.raw());
                    bc_cols[c][i].push_back(std::move(be));
                }
            }
        }
    }

    // ---- pattern-zero BC entries (solve_one_eps step 2.5) ----
    for (const auto& pg : bc_pattern_) {
        for (std::size_t i = 0; i < Nm && i < pg.size(); ++i) {
            numeric::AcbValue mu_acb;
            try {
                evaluate_const_mfrac(pg[i], numeric_q, prec, mu_acb);
            } catch (const std::exception& e) {
                throw std::runtime_error(
                    "fp_probe_solve_eps_batched: pattern mu eval: "
                    + std::string(e.what()));
            }
            for (std::size_t c = 0; c < ncols; ++c) {
                ode::BoundaryEntry be;
                acb_set(be.mu.raw(), mu_acb.raw());
                acb_zero(be.value.raw());
                bc_cols[c][i].push_back(std::move(be));
            }
        }
    }

    // ---- one batched fundamental-solution transport ----
    std::vector<std::vector<numeric::AcbValue>> sol_cols;
    try {
        numeric::GlobalScope dir_scope;
        dir_scope.running.run_direction = direction_;
        dir_scope.commit();

        sol_cols = ode::amflow_multi(de, bc_cols, prec);
    } catch (const std::exception& e) {
        std::ostringstream msg;
        msg << "AMFSystem::fp_probe_solve_eps_batched: amflow_multi failed"
            << "; system_id=" << system_id_
            << "; family=" << fc_->family
            << "; eps_index=" << eps_index
            << "; masters=" << Nm
            << "; cols=" << ncols
            << ": " << e.what();
        throw std::runtime_error(msg.str());
    }

    // ---- extract preferred values per column (solve_one_eps tail) ----
    std::vector<std::vector<numeric::AcbValue>> out(ncols);
    for (std::size_t c = 0; c < ncols; ++c) {
        out[c].reserve(preferred_.size());
        for (std::size_t i = 0; i < preferred_.size(); ++i) {
            long si = pref_to_sorted_[i];
            out[c].push_back(sol_cols[c][(std::size_t)si].clone());
        }
    }
    return out;
}

// ===========================================================================
//  Fixed-point probe driver (see header).  Assembles M(0) = A(ε)·v
//  column-by-column via unit-vector injections and writes a JSON table.
// ===========================================================================
void AMFSystem::fixedpoint_probe(const std::vector<numeric::AcbValue>& epslist,
                                 const std::string& out_path) {
    if (!fp_probe_active_) {
        throw std::runtime_error(
            "fixedpoint_probe: probe not active — boundary recursion was not "
            "skipped (no regions, or AMFLOW_FIXEDPOINT_PROBE unset at setup)");
    }
    const std::size_t n = fp_sub_masters_.size();
    const std::size_t np = preferred_.size();
    auto jkey = [](const qft::JIntegral& j) {
        std::ostringstream k;
        k << j.family();
        for (long v : j.indices()) k << "|" << v;
        return k.str();
    };

    // fp-probe column batching: AMFLOW_FP_PROBE_BATCH=1 solves
    // all n unit columns of one eps point through ONE batched
    // fundamental-solution transport (see fp_probe_solve_eps_batched).
    const char* batch_env = std::getenv("AMFLOW_FP_PROBE_BATCH");
    const bool  batch = batch_env && *batch_env
                        && std::string(batch_env) != "0";

    // columns[j][i] = preferred_i value with v = e_j, at eps index k.
    // EXACT operation sequence of the pre-shard code in both modes —
    // eps sharding only SELECTS grid indices, it never changes the
    // per-eps work (bit-identity of shard vs monolithic by construction).
    auto compute_cols = [&](std::size_t k)
        -> std::vector<std::vector<std::string>> {
        std::vector<std::vector<std::string>> cols(n);
        if (batch) {
            auto vals = fp_probe_solve_eps_batched(k, epslist[k]);
            if (vals.size() != n) {
                throw std::runtime_error(
                    "fixedpoint_probe: batched column count mismatch");
            }
            for (std::size_t j = 0; j < n; ++j) {
                cols[j].reserve(np);
                for (std::size_t i = 0; i < np; ++i) {
                    cols[j].push_back(vals[j][i].to_string(200));
                }
            }
            std::cerr << "[amflow][fp-probe] eps " << (k+1) << "/"
                      << epslist.size() << " BATCHED " << n
                      << " cols done (one fundamental transport)\n";
        } else {
        for (std::size_t j = 0; j < n; ++j) {
            std::vector<numeric::AcbValue> v;
            v.reserve(n);
            for (std::size_t m = 0; m < n; ++m) {
                numeric::AcbValue z;
                if (m == j) acb_one(z.raw()); else acb_zero(z.raw());
                v.push_back(std::move(z));
            }
            fp_injected_ = &v;
            AMFSystemSolution s;
            try {
                s = solve_one_eps(k, epslist[k]);
            } catch (...) {
                fp_injected_ = nullptr;
                throw;
            }
            fp_injected_ = nullptr;
            cols[j].reserve(np);
            for (std::size_t i = 0; i < np; ++i) {
                cols[j].push_back(s.master_values[i].to_string(200));
            }
            std::cerr << "[amflow][fp-probe] eps " << (k+1) << "/"
                      << epslist.size() << " col " << (j+1) << "/" << n
                      << " done\n";
        }
        }
        return cols;
    };

    // Shared metadata head: family, sub_masters, preferred, eps grid.
    // IDENTICAL text in the canonical dump and in every eps fragment
    // (fragments append their shard keys after it), so a merge step can
    // validate grid consistency by exact string comparison.
    auto write_head = [&](std::ostream& os) {
        os << "{\n  \"family\": \"" << fc_->family << "\",\n  \"sub_masters\": [";
        for (std::size_t m = 0; m < n; ++m) {
            os << (m ? "," : "") << "\"" << jkey(fp_sub_masters_[m]) << "\"";
        }
        os << "],\n  \"preferred\": [";
        for (std::size_t i = 0; i < np; ++i) {
            os << (i ? "," : "") << "\"" << jkey(preferred_[i]) << "\"";
        }
        os << "],\n  \"eps\": [";
        for (std::size_t k = 0; k < epslist.size(); ++k) {
            os << (k ? "," : "") << "\"" << epslist[k].to_string(200) << "\"";
        }
        os << "]";
    };

    // One A-row (all np x n values of one eps) in the canonical format:
    // [[cols of pref_0],[cols of pref_1],...], indented like the dump.
    auto write_row = [&](std::ostream& os,
                         const std::vector<std::vector<std::string>>& cols) {
        os << "    [";
        for (std::size_t i = 0; i < np; ++i) {
            os << (i ? "," : "") << "[";
            for (std::size_t j = 0; j < n; ++j) {
                os << (j ? "," : "") << "\"" << cols[j][i] << "\"";
            }
            os << "]";
        }
        os << "]";
    };

    // ---- eps sharding: AMFLOW_FP_PROBE_EPS_INDICES ----
    //
    // Comma list of 0-based eps-grid indices (e.g. "0,1,2"): run ONLY those
    // eps points and write one self-contained fragment per index,
    // <out_path minus .json>.eps<k>.json, carrying the same row schema as
    // the canonical dump plus the grid metadata needed for merging.
    // The eps grid is ALWAYS the full grid handed in by the
    // pipeline (identical regardless of sharding) — indices only select.
    // Each fragment is computed fully BEFORE its file is opened, so a
    // killed run leaves complete fragments only.  Unset => canonical
    // single-file behavior, byte-untouched.
    const char* shard_env = std::getenv("AMFLOW_FP_PROBE_EPS_INDICES");
    if (shard_env && *shard_env) {
        std::vector<std::size_t> idx;
        {
            const std::string s(shard_env);
            std::size_t pos = 0;
            while (pos <= s.size()) {
                const std::size_t c = s.find(',', pos);
                std::string tok = (c == std::string::npos)
                                      ? s.substr(pos)
                                      : s.substr(pos, c - pos);
                const auto a = tok.find_first_not_of(" \t");
                const auto b = tok.find_last_not_of(" \t");
                tok = (a == std::string::npos)
                          ? std::string()
                          : tok.substr(a, b - a + 1);
                if (!tok.empty()) {
                    if (tok.find_first_not_of("0123456789") !=
                        std::string::npos) {
                        throw std::runtime_error(
                            "fixedpoint_probe: bad AMFLOW_FP_PROBE_EPS_INDICES "
                            "token '" + tok + "' (0-based decimal indices)");
                    }
                    const unsigned long k = std::stoul(tok);
                    if (k >= epslist.size()) {
                        throw std::runtime_error(
                            "fixedpoint_probe: AMFLOW_FP_PROBE_EPS_INDICES "
                            "index " + tok + " out of range (n_eps = "
                            + std::to_string(epslist.size()) + ")");
                    }
                    idx.push_back((std::size_t)k);
                }
                if (c == std::string::npos) break;
                pos = c + 1;
            }
        }
        if (idx.empty()) {
            throw std::runtime_error(
                "fixedpoint_probe: AMFLOW_FP_PROBE_EPS_INDICES set but empty");
        }
        std::sort(idx.begin(), idx.end());
        if (std::adjacent_find(idx.begin(), idx.end()) != idx.end()) {
            throw std::runtime_error(
                "fixedpoint_probe: duplicate index in "
                "AMFLOW_FP_PROBE_EPS_INDICES");
        }
        std::string stem = out_path;
        if (stem.size() >= 5 &&
            stem.compare(stem.size() - 5, 5, ".json") == 0) {
            stem.resize(stem.size() - 5);
        }
        std::cerr << "[amflow][fp-probe] eps-shard mode: " << idx.size()
                  << " of " << epslist.size() << " eps indices selected\n";
        for (std::size_t k : idx) {
            auto cols = compute_cols(k);
            const std::string frag_path =
                stem + ".eps" + std::to_string(k) + ".json";
            std::ofstream fo(frag_path);
            if (!fo) {
                throw std::runtime_error(
                    "fixedpoint_probe: cannot open " + frag_path);
            }
            write_head(fo);
            fo << ",\n  \"n_eps\": " << epslist.size()
               << ",\n  \"eps_index\": " << k
               << ",\n  \"eps_value\": \"" << epslist[k].to_string(200)
               << "\",\n  \"A\": [\n";
            write_row(fo, cols);
            fo << "\n  ]\n}\n";
            fo.close();
            std::cerr << "[amflow][fp-probe] eps " << (k+1) << "/"
                      << epslist.size() << " fragment written: "
                      << frag_path << "\n";
        }
        std::cerr << "[amflow][fp-probe] eps-shard complete (" << idx.size()
                  << " fragments; merge on the shared metadata head)\n";
        return;
    }

    // ---- canonical single-file path (byte-identical to pre-shard code) ----
    std::ofstream out(out_path);
    if (!out) {
        throw std::runtime_error("fixedpoint_probe: cannot open " + out_path);
    }
    write_head(out);
    out << ",\n  \"A\": [\n";
    for (std::size_t k = 0; k < epslist.size(); ++k) {
        auto cols = compute_cols(k);
        out << (k ? ",\n" : "");
        write_row(out, cols);
        out.flush();
    }
    out << "\n  ]\n}\n";
    out.close();
    std::cerr << "[amflow][fp-probe] wrote " << out_path << "\n";
}

// ===========================================================================
//  Top-level orchestrators
// ===========================================================================

// ===========================================================================
//  SingleMass setup
// ===========================================================================
//
//   Mirrors AMFSystemSetupMaster[..., "SingleMass"] (.m line 1039-1061).
//
//   The single_mass_setup_master function:
//     1. Apply legs->0 to ReducedPropagator[topposi].
//     2. Call FactorizeFamily(props, [FromJ[i][topposi] for i in preferred]).
//     3. For each component (loop, prop, mas):
//          a. Find `drop` = position where mass = -1 (numerically).
//          b. droploop = the loop momentum involved in prop[drop].
//          c. Build sub-qft::FamilyConfig with droploop as leg, replacement
//             {droploop² = -1}, propagators = ToCompleteExplicit[Drop[prop, drop]].
//          d. jpre[i] = ToJ[PadRight[Drop[mas[i], {drop}], n_orig_props]]
//             on the new sub-fc family.
//          e. prefactor[i] = analytic Gamma function expression.
//          f. Recursively setup the sub-system on jpre.
//     4. Return list of AMFSystems (one per component).

namespace {
std::vector<std::unique_ptr<AMFSystem>>
single_mass_setup_master(const qft::FamilyConfig& fc,
                           const std::vector<qft::JIntegral>& preferred,
                           const AMFSystemOptions& opts);

// Find Position[mass_list, -1] (under numeric sub).  Returns -1 if not found.
//
// Same C++ enhancement pattern as `single_mass_q_numeric` above
// (locked by SingleMassEnhancement_* tests): upstream's
// `Position[ToSquareAll[prop][[2]], -1]` tests the mass list
// literally, while this routine first applies `numeric_q` so a
// symbolic mass like `msq` with `Numeric = {msq -> -1}` (or
// effectively -1 after the SingleMass loop-promotion's sign flip on
// a `Numeric = {msq -> 1}` original family) resolves to the expected
// literal `-1`.  Locked by `tests/test_amflow_amfsystem.cpp`.
// `FactorizeFamilyMassMinusOne_NumericResolvesSymbolic`.
long find_mass_minus_one(const std::vector<algebra::Mfrac>& masses,
                            const std::map<std::string, FmpqHolder>& numeric_q) {
    for (std::size_t k = 0; k < masses.size(); ++k) {
        algebra::Mfrac sub = apply_numeric(masses[k], numeric_q);
        // mass=-1 means num = -1 * den, i.e. (num + den).is_zero().
        algebra::Mpoly sum = sub.numerator() + sub.denominator();
        if (sum.is_zero()) return (long)k;
    }
    return -1;
}

// Convert ToSquareResult.masses (algebra::Mfrac on fc.ctx) to a vector and
// return.  We need this for FactorizeFamily output's prop list.
std::vector<algebra::Mfrac>
to_square_masses(const qft::FamilyConfig& fc, const std::vector<algebra::Mpoly>& props) {
    auto sq = to_square_all(fc, props);
    return std::move(sq.masses);
}

// Find the loop momentum involved in `prop` (the chosen squared-loop
// variable index).  Returns -1 if no loop appears.
long find_loop_in_prop(const qft::FamilyConfig& fc, const algebra::Mpoly& prop) {
    long n_loop = (long)fc.n_loops();
    for (long j = 0; j < n_loop; ++j) {
        if (!prop.coeff_of(j, 2).is_zero()) return j;
    }
    return -1;
}
}  // namespace

// SingleMass scheme implementation.
//
// For each FactorizeFamily component, build a new sub-qft::FamilyConfig
// where one loop momentum is "promoted to a leg" (with on-shell
// constraint loop² = -1) and one propagator is dropped.  The new
// system is on (n_loops - 1, n_props - 1) form.  An analytic
// Gamma-function prefactor links the original master to the new
// (sub-)master.
namespace {

std::vector<std::unique_ptr<AMFSystem>>
single_mass_setup_master(const qft::FamilyConfig& fc,
                           const std::vector<qft::JIntegral>& preferred,
                           const AMFSystemOptions& opts) {
    auto top = qft::get_top_position(preferred);
    if (top.empty()) return {};

    // Step 1: extract top-sector propagators with legs->0 substitution.
    std::vector<algebra::Mpoly> top_props;
    top_props.reserve(top.size());
    long n_loop_orig = (long)fc.n_loops();
    long n_red_orig  = (long)fc.n_red_legs();
    fmpz_t zero;
    fmpz_init(zero);
    for (std::size_t pos : top) {
        algebra::Mpoly p = fc.propagators_after_conservation[pos].clone();
        for (long k = 0; k < n_red_orig; ++k) {
            algebra::Mpoly out(fc.ctx);
            fmpz_mpoly_evaluate_one_fmpz(out.raw(), p.raw(),
                                          n_loop_orig + k, zero,
                                          fc.ctx->raw());
            p = std::move(out);
        }
        top_props.push_back(std::move(p));
    }
    fmpz_clear(zero);

    // Step 2: build patts (per preferred, indices at top positions).
    std::vector<std::vector<long>> patts;
    patts.reserve(preferred.size());
    for (const auto& j : preferred) {
        std::vector<long> p;
        p.reserve(top.size());
        for (std::size_t pos : top) {
            if (pos < (std::size_t)j.indices().size()) {
                p.push_back(j.indices()[pos]);
            } else {
                p.push_back(0);
            }
        }
        patts.push_back(std::move(p));
    }

    // Step 3: factorize.
    auto faminfo = factorize_family(fc, top_props, patts);
    if (faminfo.empty()) {
        throw std::runtime_error(
            "single_mass_setup_master: FactorizeFamily returned no components");
    }

    auto numeric_q = build_numeric_q(opts.bb);

    std::vector<std::unique_ptr<AMFSystem>> out;

    for (auto& fcomp : faminfo) {
        // FactorizeFamily output: fcomp has loops (subset of fc.loops),
        // propagators (legs-zeroed), patterns (per master, projected
        // to component positions), original_positions (where in
        // top_props each prop comes from).
        //
        // For SingleMass, we operate on the COMPONENT's propagators and
        // loops.  Build a temporary qft::FamilyConfig for this component.
        // Note: legs are empty (legs were zeroed in FactorizeFamily input).
        std::vector<std::string> comp_prop_strs;
        comp_prop_strs.reserve(fcomp.propagators.size());
        for (const auto& pp : fcomp.propagators) {
            comp_prop_strs.push_back(pp.to_string());
        }
        qft::FamilyConfig comp_fc = qft::FamilyConfig::build(
            fc.family + "_smcomp", fcomp.loops, /*legs=*/{},
            /*conservation=*/{},
            /*replacement=*/{},
            comp_prop_strs);

        // Re-parse propagators on comp_fc.ctx (string round-trip).
        std::vector<algebra::Mpoly> comp_props;
        comp_props.reserve(fcomp.propagators.size());
        for (const auto& s : comp_prop_strs) {
            comp_props.push_back(algebra::Mpoly::from_string(comp_fc.ctx, s));
        }

        // Find drop = position where mass = -1 numerically.
        auto comp_masses = to_square_masses(comp_fc, comp_props);
        long drop = find_mass_minus_one(comp_masses, numeric_q);
        if (drop < 0) {
            throw std::runtime_error(
                "single_mass_setup_master: no mass=-1 propagator after "
                "numeric substitution; check that AMFlowInfo[\"Numeric\"] "
                "sets exactly one propagator's mass to 1");
        }

        // droploop = loop variable in dropped propagator (use comp_props on comp_fc.ctx).
        long droploop = find_loop_in_prop(comp_fc, comp_props[(std::size_t)drop]);
        if (droploop < 0) {
            throw std::runtime_error(
                "single_mass_setup_master: cannot find loop in dropped propagator");
        }
        std::string droploop_name = fcomp.loops[(std::size_t)droploop];

        // Build sub-fc: loops = component loops minus droploop;
        //              legs = {droploop};
        //              replacement = {droploop² -> -1};
        //              propagators = Drop[fcomp.propagators, drop],
        //              then ToCompleteExplicit (added implicitly when fc is
        //              built? No -- we do it manually).
        std::vector<std::string> sub_loops;
        for (const auto& l : fcomp.loops) {
            if (l != droploop_name) sub_loops.push_back(l);
        }
        std::vector<std::string> sub_legs = {droploop_name};
        std::vector<std::pair<std::string, std::string>> sub_replacement = {
            {droploop_name + "^2", "-1"}
        };

        // Build sub propagator list: drop the `drop`-th propagator,
        // keep the others.  Use comp_props (already on comp_fc.ctx).
        std::vector<algebra::Mpoly> dropped_props;
        for (std::size_t i = 0; i < comp_props.size(); ++i) {
            if ((long)i == drop) continue;
            dropped_props.push_back(comp_props[i].clone());
        }

        // Apply ToCompleteExplicit to dropped_props.  We need a
        // qft::FamilyConfig with the new loop/leg structure to do this.
        qft::FamilyConfig pre_sub_fc = qft::FamilyConfig::build(
            fc.family + "_sm" + std::to_string(out.size()),
            sub_loops, sub_legs,
            /*conservation=*/{}, sub_replacement,
            /*propagators=*/{},
            /*cut=*/{}, /*prescription=*/{});
        // ToCompleteExplicit expects propagators on this fc's ctx.
        // dropped_props live on comp_fc.ctx; translate by variable name
        // into pre_sub_fc.ctx after promoting droploop to an external leg.
        //
        // Build a name->index map for pre_sub_fc.ctx.
        std::map<std::string, long> name_to_dst;
        for (long i = 0; i < pre_sub_fc.ctx->n_vars(); ++i) {
            name_to_dst[pre_sub_fc.ctx->var_name(i)] = i;
        }
        std::vector<algebra::Mpoly> dropped_in_sub;
        dropped_in_sub.reserve(dropped_props.size());
        for (const auto& dp : dropped_props) {
            auto src_ctx = dp.ctx();
            // Walk monomials and re-encode on pre_sub_fc.ctx.
            algebra::Mpoly out_p(pre_sub_fc.ctx);
            long len = fmpz_mpoly_length(dp.raw(), src_ctx->raw());
            std::vector<unsigned long> sexp((std::size_t)src_ctx->n_vars());
            std::vector<unsigned long> dexp(
                (std::size_t)pre_sub_fc.ctx->n_vars(), 0);
            fmpz_t coeff;
            fmpz_init(coeff);
            for (long t = 0; t < len; ++t) {
                fmpz_mpoly_get_term_exp_ui(sexp.data(), dp.raw(),
                                            t, src_ctx->raw());
                fmpz_mpoly_get_term_coeff_fmpz(coeff, dp.raw(),
                                                t, src_ctx->raw());
                std::fill(dexp.begin(), dexp.end(), 0);
                for (long sv = 0; sv < src_ctx->n_vars(); ++sv) {
                    if (sexp[(std::size_t)sv] == 0) continue;
                    const std::string& name = src_ctx->var_name(sv);
                    auto it = name_to_dst.find(name);
                    if (it == name_to_dst.end()) {
                        fmpz_clear(coeff);
                        throw std::runtime_error(
                            "single_mass_setup_master: variable '" + name
                            + "' not in SingleMass sub-family context");
                    }
                    dexp[(std::size_t)it->second] = sexp[(std::size_t)sv];
                }
                fmpz_mpoly_set_coeff_fmpz_ui(out_p.raw(), coeff,
                                              dexp.data(), pre_sub_fc.ctx->raw());
            }
            fmpz_clear(coeff);
            dropped_in_sub.push_back(std::move(out_p));
        }

        // ToCompleteExplicit on dropped_in_sub.
        auto completede = qft::to_complete_explicit(pre_sub_fc, dropped_in_sub);
        std::vector<std::string> sub_prop_strs;
        sub_prop_strs.reserve(completede.size());
        for (const auto& pp : completede) {
            sub_prop_strs.push_back(pp.to_string());
        }
        // Rebuild final sub-fc with the completed propagator list.
        qft::FamilyConfig sub_fc = qft::FamilyConfig::build(
            pre_sub_fc.family,
            sub_loops, sub_legs,
            /*conservation=*/{}, sub_replacement,
            sub_prop_strs);

        // jpre[i] = ToJ[PadRight[Drop[mas[i], {drop}], n_orig_props]]
        // mas[i] = fcomp.patterns[i] (length = top.size())
        // We drop position `drop` from mas[i], then PadRight to
        // n_orig_props (i.e. pad with zeros on the right).
        //
        // Wait -- Mathematica's PadRight is to Length[ReducedPropagator]
        // = original family size.  But in our sub-system, the J integral
        // lives on sub_fc which has Length(sub_prop_strs) propagators.
        // We use sub_prop_strs.size() as the index length.
        std::vector<qft::JIntegral> jpre;
        jpre.reserve(fcomp.patterns.size());
        for (const auto& patt : fcomp.patterns) {
            std::vector<long> idx;
            idx.reserve(sub_prop_strs.size());
            for (std::size_t k = 0; k < patt.size(); ++k) {
                if ((long)k == drop) continue;
                idx.push_back(patt[k]);
            }
            // Pad with zeros up to sub_prop_strs.size().
            while (idx.size() < sub_prop_strs.size()) idx.push_back(0);
            jpre.emplace_back(sub_fc.family, idx);
        }

        // Compute prefactor[i] symbolically as algebra::Mfrac on a context that
        // includes 'd' / 'eps' (for runtime substitution).
        //
        // prefactor = (-1)^(-n) * Gamma[2-eps-m] * Gamma[-2+eps+m+n] /
        //              (Gamma[2-eps] * Gamma[n])
        //
        // where n = patt[drop], m = Total[patt] - patt[drop] - (loop_num-1)*(2-eps)
        //
        // We can't represent Gamma symbolically with algebra::Mfrac (which is
        // polynomial-rational only).  Instead we evaluate prefactor at
        // runtime in solve_one_eps using acb_gamma.  We store the
        // *integer* parameters (n, total, loop_num) per master and
        // evaluate symbolically when ε is given.
        //
        // For now, we represent the prefactor as a symbolic placeholder
        // using a small data structure attached to AMFSystem via
        // SingleMassConfig.  The fields n_param, m_int, total_int,
        // loop_num let us compute the Gamma expression at runtime.
        //
        // algebra::Mfrac is only used for *coefficient* multiplication; for
        // Gamma we extend the representation.

        AMFSystem::SingleMassConfig sm_cfg;
        sm_cfg.original_preferred = preferred;
        sm_cfg.jpre_idx.reserve(preferred.size());
        sm_cfg.gamma_params.reserve(preferred.size());

        for (std::size_t i = 0; i < preferred.size(); ++i) {
            const auto& patt = fcomp.patterns[i];
            long n_param = patt[(std::size_t)drop];
            long total = 0;
            for (long v : patt) total += v;
            long loop_num = (long)fcomp.loops.size();
            // m = total - n - (loop_num-1) * (2 - eps)
            //   = (total - n - 2*(loop_num-1)) + (loop_num-1) * eps
            //   = m0 + m_eps_coef * eps
            AMFSystem::GammaParams gp;
            gp.n            = n_param;
            gp.m0           = total - n_param - 2 * (loop_num - 1);
            gp.m_eps_coef   = loop_num - 1;
            gp.loop_num     = loop_num;
            sm_cfg.gamma_params.push_back(gp);

            // Find which jpre index this preferred maps to.  Since
            // jpre is built per-preferred from fcomp.patterns,
            // jpre[i] (before dedup) corresponds to preferred[i].
            sm_cfg.jpre_idx.push_back((long)i);     // placeholder; fixed below
        }

        // DeleteDuplicates[jpre] and remap jpre_idx accordingly.
        std::map<std::string, long> jpre_dedup;
        std::vector<qft::JIntegral> jpre_unique;
        std::vector<long> orig_to_unique(jpre.size(), -1);
        for (std::size_t i = 0; i < jpre.size(); ++i) {
            std::string key = jintegral_key(jpre[i]);
            auto it = jpre_dedup.find(key);
            if (it == jpre_dedup.end()) {
                jpre_dedup[key] = (long)jpre_unique.size();
                jpre_unique.push_back(jpre[i]);
                orig_to_unique[i] = (long)jpre_unique.size() - 1;
            } else {
                orig_to_unique[i] = it->second;
            }
        }
        // Remap sm_cfg.jpre_idx to point into jpre_unique.
        for (std::size_t i = 0; i < sm_cfg.jpre_idx.size(); ++i) {
            sm_cfg.jpre_idx[i] = orig_to_unique[(std::size_t)sm_cfg.jpre_idx[i]];
        }

        // Compute etac on jpre_unique and DIRECTLY construct the
        // sub-system (mirrors Mathematica's AMFSystemSetup call at .m
        // line 1057 -- NOT AMFSystemSetupMaster, so the sub-system
        // does not re-enter scheme selection).
        AMFSystemOptions sub_opts = opts;
        sub_opts.current_depth = opts.current_depth + 1;
        // Sub-fc has empty Numeric (Mathematica .m line 1052 clears it).
        // Drop all numeric values; ε is supplied at solve-time.
        sub_opts.bb.numeric_values.clear();
        // A factorized (SingleMass) sub-family is a different structure
        // from the eta->inf boundary chain the GRAVITYFLOW_FIXEDPOINT
        // guard tracks — drop the inherited signature (and its top-sector
        // width) so the guard never compares a factor against an ancestor
        // it did not descend from.
        sub_opts.parent_prop_signature.clear();
        sub_opts.parent_top_sector_count = -1;

        // Compute etac for the sub-system directly:
        //   etac = AMFEtaC[AMFPosition[GetTopPosition[jpre_unique], $qft::AMFMode]]
        // If GetTopPosition is empty (vacuum), etac is all zeros and
        // the sub-system is ending.
        auto sub_top = qft::get_top_position(jpre_unique);
        std::vector<int> sub_etac(
            sub_fc.propagators_after_conservation.size(), 0);
        if (!sub_top.empty()) {
            auto sub_pos = qft::amf_position(sub_fc, sub_top, opts.amf_modes);
            if (!sub_pos.empty()) {
                sub_etac = qft::amf_eta_c(sub_fc, sub_pos);
            }
        }

        // Build sub_fc copy (qft::FamilyConfig is move-only).
        qft::FamilyConfig sub_fc_copy = qft::FamilyConfig::build(
            sub_fc.family, sub_fc.loops, sub_fc.legs,
            extract_conservation_strings(sub_fc),
            extract_replacement_strings(sub_fc),
            extract_propagator_strings(sub_fc),
            sub_fc.cut, sub_fc.prescription);

        AMFLOW_TRACE("AMFLOW_DEBUG_BC") {
            std::cerr << "[single_mass_setup] parent fc family=" << fc.family
                      << " parent_n_props="
                      << fc.propagators_after_conservation.size()
                      << " fcomp.loops.size=" << fcomp.loops.size()
                      << " fcomp.props.size=" << fcomp.propagators.size()
                      << " fcomp.patterns.size=" << fcomp.patterns.size()
                      << std::endl;
            for (std::size_t pi = 0; pi < fcomp.patterns.size(); ++pi) {
                std::cerr << "  patt[" << pi << "] = [";
                for (auto v : fcomp.patterns[pi]) std::cerr << v << ",";
                std::cerr << "]" << std::endl;
            }
            std::cerr << "  drop=" << drop
                      << " droploop_name=" << droploop_name
                      << " comp_props.size=" << comp_props.size()
                      << " dropped_props.size=" << dropped_props.size()
                      << " completede.size=" << completede.size()
                      << std::endl;
            std::cerr << "[single_mass_setup] component fam="
                      << sub_fc.family
                      << " sub_loops=" << sub_loops.size()
                      << " sub_legs=" << sub_legs.size()
                      << " props=" << sub_fc.propagators_after_conservation.size()
                      << " sub_top.size=" << sub_top.size()
                      << " sub_etac=[";
            for (auto v : sub_etac) std::cerr << v << ",";
            std::cerr << "]" << std::endl;
            std::cerr << "[single_mass_setup] jpre.size=" << jpre.size()
                      << " jpre_unique.size=" << jpre_unique.size() << std::endl;
            for (const auto& jp : jpre_unique) {
                std::cerr << "  J[" << jp.family();
                for (auto x : jp.indices()) std::cerr << "," << x;
                std::cerr << "]" << std::endl;
            }
            std::cerr << "[single_mass_setup] gamma_params:" << std::endl;
            for (std::size_t i = 0; i < sm_cfg.gamma_params.size(); ++i) {
                const auto& gp = sm_cfg.gamma_params[i];
                std::cerr << "  gp[" << i << "] n=" << gp.n
                          << " m0=" << gp.m0
                          << " m_eps=" << gp.m_eps_coef
                          << " loop_num=" << gp.loop_num << std::endl;
            }
        }
        auto sub_sys = std::make_unique<AMFSystem>(
            std::move(sub_fc_copy), jpre_unique, std::move(sub_etac),
            EndingScheme::SingleMass, sub_opts);
        sub_sys->set_single_mass_config(std::move(sm_cfg));
        out.push_back(std::move(sub_sys));
    }

    return out;
}

}  // namespace

std::vector<std::unique_ptr<AMFSystem>>
amf_system_setup_master(const qft::FamilyConfig& fc,
                          const std::vector<qft::JIntegral>& preferred,
                          const AMFSystemOptions& opts) {
    if (preferred.empty()) return {};

    // EikonalVacuum interception (region engine): a family carrying
    // an eikonal propagator with residual mass (2 l.u - 1) is SELF-SIMILAR
    // under eta re-injection — the l~eta region of (2l.u - 1 - eta)
    // regenerates the family, so Tradition/SingleMass would recurse
    // forever (or misclassify).  End it directly: force the Trivial
    // scheme; solve() then values every master via the EikonalVacuum
    // Gamma forms in try_builtin_ending_value.  MUST run before any
    // ending_q / analyze_top_sector call (those are not defined on
    // eikonal propagators).
    if (eikonal_regions_on() && has_eikonal_mass_prop(fc)) {
        AMFLOW_TRACE("AMFLOW_DEBUG_SCHEME") {
            std::cerr << "[scheme] EikonalVacuum (forced Trivial ending): "
                      << "family=" << fc.family << std::endl;
        }
        std::vector<int> zero_etac(
            fc.propagators_after_conservation.size(), 0);
        qft::FamilyConfig fc_copy = qft::FamilyConfig::build(
            fc.family, fc.loops, fc.legs,
            extract_conservation_strings(fc),
            extract_replacement_strings(fc),
            extract_propagator_strings(fc),
            fc.cut, fc.prescription);
        std::vector<std::unique_ptr<AMFSystem>> out_eik;
        out_eik.push_back(std::make_unique<AMFSystem>(
            std::move(fc_copy), preferred, std::move(zero_etac),
            EndingScheme::Trivial, opts));
        return out_eik;
    }

    auto numeric_q = build_family_numeric_q(opts.bb);
    // Amf-scheme decisions (AMFSystemEndingQ, AMFPosition,
    // AnalyzeTopSector) run on the symbolic `fc` — MMA does the same
    // (`AnalyzeTopology` at `AMFlow.m:503-510` calls `EvaluateUF`
    // without `/.Numeric`).  Substituting `numeric_q` into the
    // replacement map can collapse external-momentum invariants
    // (`p1^2 -> 0` under `ssq=0`), spuriously flipping `info.vacQ`
    // and forcing `amf_candidate_component` into Branch-mode fallback.
    //
    // Numeric substitution still flows downstream via
    // `opts.bb.numeric_values`: Kira's yaml writer (`NumericSubs` in
    // `kira_yaml.cpp`) applies it when emitting yaml; the Cutkosky
    // mass-positivity check below applies it via
    // `apply_numeric(..., numeric_q)`; SingleMass's enhanced
    // predicates apply it per mass entry.

    // Mirrors AMFSystemSetupMaster (AMFlow.m:1010 + 1034): auto-append
    // the Trivial scheme as a final fallback, then pick the first
    // scheme whose `ending_q` returns false.  Since
    // `ending_q(Trivial)` is unconditionally false, the for-loop
    // always finds something — there is no "all ending, fall back"
    // edge case to handle separately.
    std::vector<EndingScheme> schemes_with_trivial = opts.ending_schemes;
    bool has_trivial = false;
    for (auto s : schemes_with_trivial) {
        if (s == EndingScheme::Trivial) { has_trivial = true; break; }
    }
    if (!has_trivial) schemes_with_trivial.push_back(EndingScheme::Trivial);

    EndingScheme scheme = EndingScheme::Trivial;  // safe default
    for (auto s : schemes_with_trivial) {
        if (!ending_q(fc, preferred, s, opts)) { scheme = s; break; }
    }

    // FORCE-ENDING (fixed-point boundary injection): when
    // AMFLOW_FORCE_ENDING_DEPTH=<n> is set and this node sits at
    // depth >= n, take the Trivial branch unconditionally — the node's
    // masters resolve via amf_options.explicit_boundary (e.g. the
    // fixed-point eigenvector values from AMFLOW_FIXEDPOINT_PROBE) or
    // the Vacuum table.  Terminates self-similar boundary towers.
    if (const char* fe = std::getenv("AMFLOW_FORCE_ENDING_DEPTH");
        fe && *fe) {
        long lim = std::atol(fe);
        if (opts.current_depth >= lim && scheme != EndingScheme::Trivial) {
            std::cerr << "[amflow][force-ending] depth "
                      << opts.current_depth << " >= " << lim
                      << ": forcing Trivial for family=" << fc.family << "\n";
            scheme = EndingScheme::Trivial;
        }
    }

    std::vector<std::unique_ptr<AMFSystem>> out;

    if (scheme == EndingScheme::Trivial) {
        // Mirrors upstream `AMFSystemSetupMaster[preferred_, "Trivial"]`
        // (AMFlow.m:1086-1095): set up a system with NO eta injection
        // (etac all-zero).  AMFSystem::setup() will see `is_ending_ =
        // true` and the solver will look up the Vacuum table or use
        // explicit boundary conditions provided by the caller.
        AMFLOW_TRACE("AMFLOW_DEBUG_SCHEME") {
            std::cerr << "[scheme] Trivial fired: family="
                      << fc.family << std::endl;
        }
        std::vector<int> zero_etac(fc.propagators_after_conservation.size(), 0);
        qft::FamilyConfig fc_copy = qft::FamilyConfig::build(
            fc.family, fc.loops, fc.legs,
            extract_conservation_strings(fc),
            extract_replacement_strings(fc),
            extract_propagator_strings(fc),
            fc.cut, fc.prescription);
        auto sys = std::make_unique<AMFSystem>(
            std::move(fc_copy), preferred, std::move(zero_etac),
            scheme, opts);
        out.push_back(std::move(sys));
        return out;
    }

    if (scheme == EndingScheme::SingleMass) {
        AMFLOW_TRACE("AMFLOW_DEBUG_SCHEME") {
            auto top = qft::get_top_position(preferred);
            std::cerr << "[scheme] SingleMass fired: family="
                      << fc.family
                      << " top_position={";
            for (std::size_t i = 0; i < top.size(); ++i) {
                if (i) std::cerr << ',';
                std::cerr << top[i];
            }
            std::cerr << "}" << std::endl;
        }
        return single_mass_setup_master(fc, preferred, opts);
    }

    if (scheme == EndingScheme::Cutkosky) {
        auto top = qft::get_top_position(preferred);
        auto info = qft::analyze_top_sector(fc, top);
        long phase_loop_num = -1;
        const qft::TopSectorComponentInfo* phase_comp = nullptr;
        for (const auto& comp : info) {
            if (qft::phase_volume_q(comp)) {
                if (phase_loop_num >= 0) {
                    throw std::runtime_error(
                        "amf_system_setup_master: Cutkosky found multiple "
                        "phase-volume components");
                }
                phase_loop_num = comp.loopnum;
                phase_comp = &comp;
            }
        }
        if (phase_loop_num < 0) {
            throw std::runtime_error(
                "amf_system_setup_master: Cutkosky found no "
                "phase-volume component");
        }

        // Mirrors AMFlow.m:1050: abort if any cutcom mass is negative
        // after `Numeric` substitution.  Negative squared masses in the
        // cut component are unphysical and cannot yield a real cut
        // rate; allowing them produces meaningless numbers.
        for (std::size_t mi = 0; mi < phase_comp->mass.size(); ++mi) {
            algebra::Mfrac evaluated =
                apply_numeric(phase_comp->mass[mi], numeric_q);
            const auto* mctx = evaluated.ctx().get();
            algebra::Mpoly num = evaluated.numerator();
            algebra::Mpoly den = evaluated.denominator();
            if (!fmpz_mpoly_is_fmpz(num.raw(), mctx->raw()) ||
                !fmpz_mpoly_is_fmpz(den.raw(), mctx->raw())) {
                throw std::runtime_error(
                    "amf_system_setup_master: Cutkosky cannot validate "
                    "physical-mass condition because a phase-volume "
                    "component mass remains symbolic after Numeric "
                    "substitution; supply a numeric value for every "
                    "mass and kinematic invariant via "
                    "amf_options.blackbox.numeric_values.  "
                    "(Mirrors AMFlow.m:1050.)");
            }
            fmpz_t num_z, den_z;
            fmpz_init(num_z); fmpz_init(den_z);
            fmpz_mpoly_get_fmpz(num_z, num.raw(), mctx->raw());
            fmpz_mpoly_get_fmpz(den_z, den.raw(), mctx->raw());
            const int sgn_num = fmpz_sgn(num_z);
            const int sgn_den = fmpz_sgn(den_z);
            fmpz_clear(num_z); fmpz_clear(den_z);
            if (sgn_den == 0) {
                throw std::runtime_error(
                    "amf_system_setup_master: Cutkosky mass denominator "
                    "is zero after Numeric substitution");
            }
            if (sgn_num * sgn_den < 0) {
                std::ostringstream m;
                m << "amf_system_setup_master: Cutkosky cannot proceed: "
                  << "phase-volume component mass[" << mi
                  << "] is negative after Numeric substitution.  "
                     "Negative squared masses are unphysical for the "
                     "Cutkosky scheme.  (Mirrors AMFlow.m:1050.)";
                throw std::runtime_error(m.str());
            }
        }

        AMFLOW_TRACE("AMFLOW_DEBUG_SCHEME") {
            std::cerr << "[scheme] Cutkosky fired: family="
                      << fc.family
                      << " phase_loop_num=" << phase_loop_num
                      << " top_position={";
            for (std::size_t i = 0; i < top.size(); ++i) {
                if (i) std::cerr << ',';
                std::cerr << top[i];
            }
            std::cerr << "}" << std::endl;
        }

        qft::FamilyConfig fc_cutless = qft::FamilyConfig::build(
            fc.family, fc.loops, fc.legs,
            extract_conservation_strings(fc),
            extract_replacement_strings(fc),
            extract_propagator_strings(fc),
            /*cut=*/{}, /*prescription=*/{});

        auto pos = qft::amf_position(fc_cutless, top, opts.amf_modes);
        if (pos.empty() && !top.empty()) {
            throw std::runtime_error(
                "amf_system_setup_master: Cutkosky could not find an "
                "eta-injection position after clearing cut/prescription");
        }
        auto etac = qft::amf_eta_c(fc_cutless, pos);

        AMFSystem::CutkoskyConfig ck_cfg;
        ck_cfg.original_preferred = preferred;
        ck_cfg.phase_loop_num = phase_loop_num;

        auto sys = std::make_unique<AMFSystem>(
            std::move(fc_cutless), preferred, std::move(etac),
            scheme, opts);
        sys->set_cutkosky_config(std::move(ck_cfg));
        out.push_back(std::move(sys));
        return out;
    }

    if (scheme != EndingScheme::Tradition) {
        throw std::runtime_error(
            "amf_system_setup_master: unexpected scheme " +
            std::string(ending_scheme_name(scheme)));
    }

    auto top = qft::get_top_position(preferred);
    auto pos = qft::amf_position(fc, top, opts.amf_modes);
    if (pos.empty()) {
        // This shouldn't happen since ending_q(Tradition) returned false.
        throw std::runtime_error(
            "amf_system_setup_master: no eta-injection position found");
    }
    auto etac = qft::amf_eta_c(fc, pos);

    // ETAC OVERRIDE (multi-flow fixed-point closure): force a custom
    // η-injection vector at the TOP system (depth 0) only — used with
    // AMFLOW_FIXEDPOINT_PROBE to generate INDEPENDENT flows A_β on the
    // same master basis (different prop subsets), whose joint
    // (A_i − I)-constraints close the self-similar boundary system.
    // Format: comma-separated ints, length = #propagators, entries in
    // {0,-1}; cut positions must be 0.
    if (const char* eo = std::getenv("AMFLOW_ETAC_OVERRIDE");
        eo && *eo && opts.current_depth == 0) {
        std::vector<int> ov;
        std::string s(eo);
        std::size_t p = 0;
        while (p <= s.size()) {
            std::size_t c = s.find(',', p);
            if (c == std::string::npos) c = s.size();
            ov.push_back(std::atoi(s.substr(p, c - p).c_str()));
            p = c + 1;
        }
        if (ov.size() != etac.size()) {
            throw std::runtime_error(
                "AMFLOW_ETAC_OVERRIDE: length " + std::to_string(ov.size())
                + " != #propagators " + std::to_string(etac.size()));
        }
        for (std::size_t i = 0; i < ov.size(); ++i) {
            if (ov[i] != 0 && (i < fc.cut.size() && fc.cut[i])) {
                throw std::runtime_error(
                    "AMFLOW_ETAC_OVERRIDE: eta on cut propagator "
                    + std::to_string(i));
            }
        }
        std::cerr << "[amflow][etac-override] applying custom etac at depth 0\n";
        etac = std::move(ov);
    }
    AMFLOW_TRACE("AMFLOW_DEBUG_SCHEME") {
        std::cerr << "[scheme] Tradition fired: family=" << fc.family
                  << " top={";
        for (std::size_t i = 0; i < top.size(); ++i) {
            if (i) std::cerr << ',';
            std::cerr << top[i];
        }
        std::cerr << "} pos={";
        for (std::size_t i = 0; i < pos.size(); ++i) {
            if (i) std::cerr << ',';
            std::cerr << pos[i];
        }
        std::cerr << "} etac=[";
        for (std::size_t i = 0; i < etac.size(); ++i) {
            if (i) std::cerr << ',';
            std::cerr << etac[i];
        }
        std::cerr << "]" << std::endl;

        // Detailed dump for scheme-dispatch diagnosis: top-sector
        // propagators + their U-polynomial component decomposition.
        std::cerr << "    top-sector props (after conservation):" << std::endl;
        for (std::size_t k = 0; k < top.size(); ++k) {
            std::cerr << "      [" << top[k] << "] "
                      << fc.propagators_after_conservation[top[k]].to_string()
                      << std::endl;
        }
        auto comps = qft::analyze_top_sector(fc, top);
        std::cerr << "    components (" << comps.size() << "):" << std::endl;
        for (std::size_t ci = 0; ci < comps.size(); ++ci) {
            const auto& comp = comps[ci];
            std::cerr << "      comp[" << ci << "] loopnum=" << comp.loopnum
                      << " prop_index={";
            for (std::size_t k = 0; k < comp.prop_index.size(); ++k) {
                if (k) std::cerr << ',';
                std::cerr << comp.prop_index[k];
            }
            std::cerr << "} mass=[";
            for (std::size_t k = 0; k < comp.mass.size(); ++k) {
                if (k) std::cerr << ',';
                std::cerr << comp.mass[k].to_string();
            }
            std::cerr << "]" << std::endl;
        }
    }
    qft::FamilyConfig fc_copy = qft::FamilyConfig::build(
        fc.family, fc.loops, fc.legs,
        extract_conservation_strings(fc),
        extract_replacement_strings(fc),
        extract_propagator_strings(fc),
        fc.cut, fc.prescription);
    auto sys = std::make_unique<AMFSystem>(
        std::move(fc_copy), preferred, std::move(etac), scheme, opts);
    out.push_back(std::move(sys));
    return out;
}

// ---------------------------------------------------------------------------
//  tree_stats — node-tree parallelism diagnostic.
// ---------------------------------------------------------------------------
void AMFSystem::tree_stats(long& total, long& leaves,
                           long& max_siblings) const {
    ++total;
    if (is_ending_ || subsystems_.empty()) {
        ++leaves;
        return;
    }
    if ((long)subsystems_.size() > max_siblings) {
        max_siblings = (long)subsystems_.size();
    }
    for (const auto& sub : subsystems_) {
        sub->tree_stats(total, leaves, max_siblings);
    }
}

std::vector<std::unique_ptr<AMFSystem>>
amf_systems_setup(const qft::FamilyConfig& fc,
                  const std::vector<qft::JIntegral>& preferred,
                  const AMFSystemOptions& opts) {
    auto roots = amf_system_setup_master(fc, preferred, opts);
    for (auto& r : roots) r->setup();
    return roots;
}

std::vector<AMFSystemSolution>
amf_systems_solution(std::vector<std::unique_ptr<AMFSystem>>& systems,
                     const std::vector<numeric::AcbValue>& epslist) {
    if (systems.empty()) return {};

    // Node-tree diagnostic: report the AMFSystem dependency-DAG shape so
    // the caller can see how much node-level parallelism is on the
    // table (independent siblings within one ε), alongside the active
    // ε-grid parallelism and node-tree parallelism.
    {
        long total = 0, leaves = 0, max_sib = 0;
        for (const auto& sys : systems) sys->tree_stats(total, leaves, max_sib);
        const long ep = resolve_eps_parallel(systems.front()->options());
        const long np = resolve_node_parallel(systems.front()->options());
        std::cerr << "[amflow] node-tree: " << total << " systems, "
                  << leaves << " leaves, max " << max_sib
                  << " independent siblings; eps-grid=" << epslist.size()
                  << " samples, eps_parallel=" << ep
                  << (ep > 1 ? " (threaded)" : " (serial)")
                  << ", node_parallel=" << np
                  << (np > 1 ? " (threaded)" : " (serial)") << "\n";
    }

    for (auto& sys : systems) {
        sys->solve(epslist);
    }

    // Mirror Mathematica's AMFSystemCombineSolution (.m line 1188-1197):
    //   For each sysid, read GlobalPreferred = (global, coe, local, op).
    //   Compute  values[sysid][i] = coe[i] * op(local[i] /. sol[sysid]).
    //   All sysids share the same `global` list (Keys); the final
    //   answer per global[i] is the PRODUCT over sysids:
    //       result[i] = Product_sysid values[sysid][i]
    //
    // Tradition systems have an empty SingleMassConfig, so their
    // global_values is empty and they directly contribute their
    // master_values for the global preferred.
    std::vector<AMFSystemSolution> out;
    const auto& src0 = systems.front()->solutions();
    out.reserve(src0.size());
    long prec = numeric::working_prec_bits();

    for (std::size_t k = 0; k < src0.size(); ++k) {
        AMFSystemSolution cs;
        // master_values: root's solution (kept for inspection).
        cs.master_values.reserve(src0[k].master_values.size());
        for (const auto& v : src0[k].master_values) {
            cs.master_values.push_back(v.clone());
        }
        // global_preferred: take from the FIRST system's
        // global_preferred (all sysids share the same global list).
        // If the first system has no global_preferred (i.e. Tradition
        // single root), use master_values directly as the answer.
        if (!systems.front()->solutions()[k].global_preferred.empty()) {
            cs.global_preferred = systems.front()->solutions()[k].global_preferred;
            // Initialize global_values to all 1s (multiplicative identity).
            cs.global_values.reserve(cs.global_preferred.size());
            for (std::size_t i = 0; i < cs.global_preferred.size(); ++i) {
                numeric::AcbValue one;
                acb_one(one.raw());
                cs.global_values.push_back(std::move(one));
            }
            // Multiply per-system contributions.
            for (auto& sys : systems) {
                const auto& sol_k = sys->solutions()[k];
                if (sol_k.global_values.size() != cs.global_values.size()) {
                    throw std::runtime_error(
                        "amf_systems_solution: global_values size mismatch "
                        "across roots (FactorizeFamily components must "
                        "agree on global count)");
                }
                for (std::size_t i = 0; i < cs.global_values.size(); ++i) {
                    numeric::AcbValue tmp;
                    acb_mul(tmp.raw(), cs.global_values[i].raw(),
                              sol_k.global_values[i].raw(), prec);
                    cs.global_values[i] = std::move(tmp);
                }
            }
        }
        out.push_back(std::move(cs));
    }
    return out;
}

// stage-2 endings: subsector corners + one-dot variants of the given
// targets, as a Kira preferred-masters list, ONLY for eikonal-vacuum
// families with the region engine on (else empty: default paths untouched).
// The stage-2 ending values are wired for denominator indices only, so the
// reduction must not park masters on numerator (ISP) integrals.
std::vector<qft::JIntegral>
build_eikonal_corner_preferred(const qft::FamilyConfig& fc,
                               const std::vector<qft::JIntegral>& targets) {
    std::vector<qft::JIntegral> preferred;
    if (!eikonal_regions_on() || !has_eikonal_mass_prop(fc)) {
        return preferred;
    }
    std::set<std::string> seen;
    for (const auto& target : targets) {
        const auto& tidx = target.indices();
        std::vector<std::size_t> pos;
        for (std::size_t i = 0; i < tidx.size(); ++i) {
            if (tidx[i] > 0) pos.push_back(i);
        }
        if (pos.size() > 10) continue;
        const std::size_t n = pos.size();
        auto push = [&](const std::vector<long>& idx) {
            std::string k;
            for (long v : idx) k += std::to_string(v) + "|";
            if (seen.insert(k).second) {
                preferred.emplace_back(target.family(), idx);
            }
        };
        for (std::size_t mask = 1; mask < (1ull << n); ++mask) {
            std::vector<long> idx(tidx.size(), 0);
            for (std::size_t b = 0; b < n; ++b) {
                if (mask & (1ull << b)) idx[pos[b]] = 1;
            }
            push(idx);
            for (std::size_t b = 0; b < n; ++b) {
                if (!(mask & (1ull << b))) continue;
                std::vector<long> idx2 = idx;
                idx2[pos[b]] = 2;
                push(idx2);
            }
        }
    }
    return preferred;
}

}  // namespace amflow::pipeline
