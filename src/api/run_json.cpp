// SPDX-License-Identifier: MIT
// amflow::api::run_json — implementation.
//

#include "amflow/api/run_json.hpp"

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <flint/acb.h>
#include <flint/arb.h>
#include <flint/fmpq.h>
#include <flint/fmpq_poly.h>

#include "amflow/pipeline/amfsystem.hpp"
#include "amflow/pipeline/solve_integrals.hpp"
#include "amflow/ibp/reduce.hpp"
#include "amflow/numeric/acb_wrap.hpp"
#include "amflow/numeric/options.hpp"
#include "amflow/numeric/rational.hpp"
#include "amflow/ode/amflow.hpp"
#include "amflow/ode/inf.hpp"
#include "amflow/qft/amfmode.hpp"
#include "amflow/qft/family_config.hpp"
#include "amflow/qft/jintegral.hpp"
#include "amflow/qft/region.hpp"

namespace amflow::api {

using nlohmann::json;
using algebra::Mfrac;
using algebra::Mpoly;
using numeric::AcbValue;
using numeric::FmpqPoly;
using numeric::RationalComplex;
using numeric::RationalFunction;
using numeric::RationalMatrix;

namespace {

// Forward declaration: defined below, used by parse_amf_options for
// explicit_boundary parsing.
AcbValue parse_complex(const json& obj, long prec);

std::string json_scalar_to_string(const json& value, const std::string& name) {
    if (value.is_string()) return value.get<std::string>();
    if (value.is_number_integer()) {
        return std::to_string(value.get<long long>());
    }
    if (value.is_number_unsigned()) {
        return std::to_string(value.get<unsigned long long>());
    }
    if (value.is_number_float()) {
        std::ostringstream os;
        os << value.get<double>();
        return os.str();
    }
    throw std::runtime_error(name + " must be a string or number");
}

std::vector<std::string> parse_string_array(const json& arr,
                                              const std::string& name) {
    if (!arr.is_array()) {
        throw std::runtime_error(name + " must be an array");
    }
    std::vector<std::string> out;
    out.reserve(arr.size());
    for (std::size_t i = 0; i < arr.size(); ++i) {
        if (!arr[i].is_string()) {
            throw std::runtime_error(name + "[" + std::to_string(i)
                                       + "] must be a string");
        }
        out.push_back(arr[i].get<std::string>());
    }
    return out;
}

std::vector<int> parse_int_array(const json& arr, const std::string& name) {
    if (!arr.is_array()) {
        throw std::runtime_error(name + " must be an array");
    }
    std::vector<int> out;
    out.reserve(arr.size());
    for (std::size_t i = 0; i < arr.size(); ++i) {
        if (!arr[i].is_number_integer()) {
            throw std::runtime_error(name + "[" + std::to_string(i)
                                       + "] must be an integer");
        }
        out.push_back(arr[i].get<int>());
    }
    return out;
}

std::vector<std::pair<std::string, std::string>>
parse_rule_list(const json& node, const std::string& name) {
    std::vector<std::pair<std::string, std::string>> out;
    if (node.is_null()) return out;

    if (node.is_object()) {
        for (auto it = node.begin(); it != node.end(); ++it) {
            out.emplace_back(it.key(),
                              json_scalar_to_string(it.value(),
                                                     name + "." + it.key()));
        }
        return out;
    }

    if (!node.is_array()) {
        throw std::runtime_error(name + " must be an object or array");
    }

    out.reserve(node.size());
    for (std::size_t i = 0; i < node.size(); ++i) {
        const auto& entry = node[i];
        if (entry.is_array()) {
            if (entry.size() != 2 || !entry[0].is_string()) {
                throw std::runtime_error(name + "[" + std::to_string(i)
                                           + "] must be [lhs, rhs]");
            }
            out.emplace_back(entry[0].get<std::string>(),
                              json_scalar_to_string(entry[1],
                                                     name + "[" + std::to_string(i) + "][1]"));
            continue;
        }
        if (!entry.is_object() || !entry.contains("lhs") || !entry.contains("rhs")
                || !entry.at("lhs").is_string()) {
            throw std::runtime_error(name + "[" + std::to_string(i)
                                       + "] must contain string lhs/rhs");
        }
        out.emplace_back(entry.at("lhs").get<std::string>(),
                          json_scalar_to_string(entry.at("rhs"),
                                                 name + "[" + std::to_string(i) + "].rhs"));
    }

    return out;
}

// ---------------------------------------------------------------------------
// Config hygiene (this fork): unrecognized keys in
// any config object WARN LOUDLY on stderr instead of being silently ignored.
// Failure mode this prevents: a knob placed in the wrong JSON object (the canonical
// case: `eps_parallel` under `options` instead of `amf_options`) parses
// cleanly, does nothing, and costs a serial night.
//
// Behaviour:
//   - default: one stderr warning per (object, key) per process; parsing
//     continues unchanged (legacy configs must keep running).
//   - env AMFLOW_STRICT_CONFIG=1 (any non-empty value != "0"): upgrade every
//     such warning to an abort (std::runtime_error).
//   - keys beginning with '_' are reserved for comments and never flagged.
//
// If the stray key IS a known knob of another object, the warning names the
// location(s) where it WOULD be read.  The accepted-key list at each call
// site and `config_key_locations()` below are enumerated from the parsing
// code in THIS file — when adding a knob, update both in the same edit.

bool strict_config_enabled() {
    const char* v = std::getenv("AMFLOW_STRICT_CONFIG");
    return v != nullptr && *v != '\0' && std::string(v) != "0";
}

// key -> where that key IS read (the "known-elsewhere" hint table).
const std::map<std::string, std::string>& config_key_locations() {
    static const std::map<std::string, std::string> m = {
        // options.* (apply_options)
        {"working_pre",       "options.working_pre"},
        {"chop_pre",          "options.chop_pre"},
        {"rationalize_pre",   "options.rationalize_pre"},
        {"silent_mode",       "options.silent_mode"},
        {"d0",                "options.d0"},
        {"x_order",           "options.x_order"},
        {"extra_x_order",     "options.extra_x_order"},
        {"learn_x_order",     "options.learn_x_order"},
        {"test_x_order",      "options.test_x_order"},
        {"run_radius",        "options.run_radius"},
        {"run_length",        "options.run_length"},
        {"run_candidate",     "options.run_candidate"},
        {"run_direction",     "options.run_direction"},
        // amf_options.* (parse_amf_options)
        {"amf_modes",           "amf_options.amf_modes"},
        {"ending_schemes",      "amf_options.ending_schemes"},
        {"cache_root",          "amf_options.cache_root"},
        {"direction",           "amf_options.direction"},
        {"max_recursion_depth", "amf_options.max_recursion_depth"},
        {"eps_parallel",
         "amf_options.eps_parallel or AMFLOW_EPS_PARALLEL"},
        {"node_parallel",
         "amf_options.node_parallel or AMFLOW_NODE_PARALLEL"},
        {"explicit_boundary",   "amf_options.explicit_boundary"},
        {"blackbox",            "amf_options.blackbox"},
        {"bb",                  "amf_options.bb (alias of blackbox)"},
        // amf_options.blackbox.* (apply_blackbox_options)
        {"ibp_rank",          "amf_options.blackbox.ibp_rank"},
        {"ibp_dot",           "amf_options.blackbox.ibp_dot"},
        {"n_thread",          "amf_options.blackbox.n_thread"},
        {"integral_order",    "amf_options.blackbox.integral_order"},
        {"kira_executable",   "amf_options.blackbox.kira_executable"},
        {"fermat_executable", "amf_options.blackbox.fermat_executable"},
        {"run_firefly",       "amf_options.blackbox.run_firefly"},
        {"firefly_bunch",     "amf_options.blackbox.firefly_bunch"},
        {"auto_firefly",
         "amf_options.blackbox.auto_firefly "
         "(env AMFLOW_NO_AUTO_FIREFLY=1 forces the Fermat path)"},
        {"kira3_truncate_sp", "amf_options.blackbox.kira3_truncate_sp"},
        {"kira3_sectorwise",  "amf_options.blackbox.kira3_sectorwise"},
        {"symbolic_ibp",
         "amf_options.blackbox.symbolic_ibp or AMFLOW_SYMBOLIC_IBP"},
        {"ibp_cache_dir",
         "amf_options.blackbox.ibp_cache_dir or AMFLOW_IBP_CACHE"},
        {"work_dir",
         "top-level work_dir or amf_options.blackbox.work_dir"},
        {"log_file",          "amf_options.blackbox.log_file"},
        {"numeric_values",    "amf_options.blackbox.numeric_values"},
        // family.* (parse_family_config)
        {"name",          "family.name"},
        {"loops",         "family.loops"},
        {"legs",          "family.legs"},
        {"conservation",  "family.conservation"},
        {"replacement",   "family.replacement"},
        {"propagators",   "family.propagators"},
        {"cut",           "family.cut"},
        {"prescription",  "family.prescription"},
        {"family",
         "top-level family (the family object), family.name alias, "
         "or integrals[i].family"},
        // per-integral entries (parse_jintegral)
        {"indices",       "integrals[i].indices"},
        // top-level (run_json)
        {"mode",             "top-level mode"},
        {"options",          "top-level options"},
        {"amf_options",      "top-level amf_options"},
        {"matrix",           "top-level matrix (mode:amflow)"},
        {"boundaries",       "top-level boundaries (mode:amflow)"},
        {"integrals",        "top-level integrals"},
        {"jints",            "top-level jints (alias of integrals)"},
        {"targets",          "top-level targets (alias of integrals)"},
        {"goal_digits",      "top-level goal_digits (mode:solve_integrals)"},
        {"goal",             "top-level goal (alias of goal_digits)"},
        {"eps_order",        "top-level eps_order (mode:solve_integrals)"},
        {"order",            "top-level order (alias of eps_order)"},
        {"target_eps_power",
         "top-level target_eps_power (mode:solve_integrals)"},
        {"eps_samples",
         "top-level eps_samples (mode:black_box_amflow)"},
        {"epslist",
         "top-level epslist (alias of eps_samples)"},
        {"variables",        "top-level variables (mode:diffeq)"},
        {"top_pattern",      "top-level top_pattern (mode:diffeq)"},
    };
    return m;
}

// Warn (or abort under AMFLOW_STRICT_CONFIG=1) about keys of `obj` that are
// not in `accepted`.  `obj_path` is the human-readable location used in the
// message; `dedup_path` (defaults to obj_path) is the key under which the
// once-per-process suppression is recorded — pass a normalized path for
// per-element objects (e.g. "integrals[i]") so a repeated mistake across
// hundreds of entries warns once, not hundreds of times.
void warn_unknown_keys(const json& obj,
                       const std::string& obj_path,
                       const std::vector<std::string>& accepted,
                       const std::string& dedup_path = std::string{}) {
    if (!obj.is_object()) return;
    static std::set<std::string> warned_once;
    const std::string& dedup =
        dedup_path.empty() ? obj_path : dedup_path;
    for (auto it = obj.begin(); it != obj.end(); ++it) {
        const std::string& key = it.key();
        if (!key.empty() && key[0] == '_') continue;  // comment convention
        if (std::find(accepted.begin(), accepted.end(), key)
                != accepted.end()) {
            continue;
        }
        std::string msg = "'" + key + "' found in " + obj_path
                          + " — not read; ";
        const auto& locs = config_key_locations();
        const auto hit = locs.find(key);
        if (hit != locs.end()) {
            msg += "use " + hit->second;
        } else {
            std::string keys;
            for (const auto& a : accepted) {
                if (!keys.empty()) keys += ", ";
                keys += a;
            }
            msg += "not a known key anywhere; accepted keys in " + obj_path
                   + ": " + keys;
        }
        if (strict_config_enabled()) {
            throw std::runtime_error(
                "unrecognized config key (fatal because "
                "AMFLOW_STRICT_CONFIG=1): " + msg);
        }
        if (warned_once.insert(dedup + '\x1f' + key).second) {
            std::fprintf(stderr, "[amflow][config] WARNING: %s\n",
                         msg.c_str());
            std::fflush(stderr);
        }
    }
}

long parse_required_long(const json& obj,
                            const std::string& primary,
                            const std::string& alias = "") {
    if (obj.contains(primary)) {
        if (!obj.at(primary).is_number_integer()) {
            throw std::runtime_error(primary + " must be an integer");
        }
        return obj.at(primary).get<long>();
    }
    if (!alias.empty() && obj.contains(alias)) {
        if (!obj.at(alias).is_number_integer()) {
            throw std::runtime_error(alias + " must be an integer");
        }
        return obj.at(alias).get<long>();
    }
    throw std::runtime_error("missing required integer field '" + primary + "'");
}

qft::FamilyConfig parse_family_config(const json& obj) {
    if (!obj.is_object()) {
        throw std::runtime_error("family must be an object");
    }
    warn_unknown_keys(obj, "family",
                      {"name", "family", "loops", "legs", "conservation",
                       "replacement", "propagators", "cut", "prescription"});

    const std::string name = obj.contains("name")
        ? obj.at("name").get<std::string>()
        : obj.at("family").get<std::string>();
    const auto loops = parse_string_array(obj.at("loops"), "family.loops");
    const auto legs = parse_string_array(obj.at("legs"), "family.legs");
    const auto conservation = obj.contains("conservation")
        ? parse_rule_list(obj.at("conservation"), "family.conservation")
        : std::vector<std::pair<std::string, std::string>>{};
    const auto replacement = obj.contains("replacement")
        ? parse_rule_list(obj.at("replacement"), "family.replacement")
        : std::vector<std::pair<std::string, std::string>>{};
    const auto propagators =
        parse_string_array(obj.at("propagators"), "family.propagators");
    const auto cut = obj.contains("cut")
        ? parse_int_array(obj.at("cut"), "family.cut")
        : std::vector<int>{};
    const auto prescription = obj.contains("prescription")
        ? parse_int_array(obj.at("prescription"), "family.prescription")
        : std::vector<int>{};

    qft::FamilyConfig fc =
        qft::FamilyConfig::build(name, loops, legs, conservation, replacement,
                                   propagators, cut, prescription);

    // Auto-complete the propagator basis with irreducible scalar products
    // (ISPs) when the user supplies an INCOMPLETE family.  Kira (and the IBP
    // machinery generally) requires a basis that spans every loop-dependent
    // scalar product: n_sp = L(L+1)/2 + L*E_indep denominators.  A graph
    // like the massless planar double box has only its 7 physical lines but
    // needs 9 denominators (2 ISPs), and without completion Kira aborts with
    // "N propagators are needed and M are provided".  The single-mass /
    // factorize sub-family code paths already complete internally
    // (amfsystem.cpp), and 1-loop topologies happen to be complete as given,
    // which is why this only bit fully-massless >=2-loop graphs.
    //
    // `to_complete_explicit` reuses the physical denominators FIRST (they are
    // the leading rows of its RREF and the returned picks are sorted), so the
    // original propagators keep their positions and the appended ISP slots
    // come strictly after.  Target integral indices are therefore made valid
    // by zero-padding on the right (see parse_jintegral), which is what an
    // ISP power of 0 means.  Completing an already-complete family is a no-op
    // (picks == 0..M-1 in order), so this is backward-compatible.
    if (fc.n_propagators() < fc.n_sp()) {
        auto completed =
            qft::to_complete_explicit(fc, fc.propagators_after_conservation);
        std::vector<std::string> completed_strs;
        completed_strs.reserve(completed.size());
        for (const auto& p : completed) completed_strs.push_back(p.to_string());
        // Rebuild on the completed denominator list.  The completed
        // denominators are already expressed post-conservation (they never
        // contain the eliminated leg), but `legs` and `replacement` are still
        // stated in terms of the full external set, so we MUST re-supply the
        // same `conservation` here: that is what eliminates the dependent leg
        // (e.g. p4 = -p1-p2-p3) and reduces the kinematics to the independent
        // legs.  Re-applying conservation to denominators that do not mention
        // the dependent leg is a no-op, so this is safe.  A user-supplied
        // `cut` only makes sense against the physical propagators, so extend
        // it with zeros for the appended ISP slots.
        std::vector<int> cut_ext = cut;
        if (!cut_ext.empty()) {
            cut_ext.resize(completed_strs.size(), 0);
        }
        fc = qft::FamilyConfig::build(name, loops, legs,
                                        conservation, replacement,
                                        completed_strs, cut_ext, prescription);
    }

    return fc;
}

qft::JIntegral parse_jintegral(const json& entry,
                                  const qft::FamilyConfig& fc,
                                  const std::string& name) {
    std::string family = fc.family;
    const json* idx_node = &entry;

    if (entry.is_object()) {
        warn_unknown_keys(entry, name, {"family", "indices"},
                          /*dedup_path=*/"integrals[i]");
        if (entry.contains("family")) {
            if (!entry.at("family").is_string()) {
                throw std::runtime_error(name + ".family must be a string");
            }
            family = entry.at("family").get<std::string>();
        }
        if (!entry.contains("indices")) {
            throw std::runtime_error(name + " must contain 'indices'");
        }
        idx_node = &entry.at("indices");
    }

    if (!idx_node->is_array()) {
        throw std::runtime_error(name + ".indices must be an array");
    }

    std::vector<long> indices;
    indices.reserve(idx_node->size());
    for (std::size_t i = 0; i < idx_node->size(); ++i) {
        if (!(*idx_node)[i].is_number_integer()) {
            throw std::runtime_error(name + ".indices[" + std::to_string(i)
                                       + "] must be an integer");
        }
        indices.push_back((*idx_node)[i].get<long>());
    }

    // When the family was auto-completed with ISP denominators (see
    // parse_family_config), a user-supplied index vector that matches the
    // original physical propagator count is shorter than n_propagators().
    // The appended slots are ISPs, which a physical target carries at power
    // 0, so zero-pad on the right.  An over-long vector is still an error.
    if (indices.size() > fc.n_propagators()) {
        throw std::runtime_error(name + " has " + std::to_string(indices.size())
                                   + " indices, expected at most "
                                   + std::to_string(fc.n_propagators()));
    }
    indices.resize(fc.n_propagators(), 0);
    return qft::JIntegral(family, std::move(indices));
}

std::vector<qft::JIntegral> parse_jintegral_list(const json& obj,
                                                    const qft::FamilyConfig& fc) {
    const json* arr = nullptr;
    if (obj.contains("integrals")) arr = &obj.at("integrals");
    else if (obj.contains("jints")) arr = &obj.at("jints");
    else if (obj.contains("targets")) arr = &obj.at("targets");
    else throw std::runtime_error("missing required array 'integrals'");

    if (!arr->is_array()) {
        throw std::runtime_error("integrals must be an array");
    }

    std::vector<qft::JIntegral> out;
    out.reserve(arr->size());
    for (std::size_t i = 0; i < arr->size(); ++i) {
        out.push_back(parse_jintegral((*arr)[i], fc,
                                         "integrals[" + std::to_string(i) + "]"));
    }
    return out;
}

pipeline::EndingScheme parse_ending_scheme_string(const std::string& s) {
    if (s == "Tradition") return pipeline::EndingScheme::Tradition;
    if (s == "Cutkosky") return pipeline::EndingScheme::Cutkosky;
    if (s == "SingleMass") return pipeline::EndingScheme::SingleMass;
    if (s == "Trivial") return pipeline::EndingScheme::Trivial;
    throw std::runtime_error("unknown ending scheme '" + s + "'");
}

void apply_blackbox_options(const json& obj, ibp::ReduceOptions& bb,
                            const std::string& obj_path) {
    warn_unknown_keys(obj, obj_path,
                      {"ibp_rank", "ibp_dot", "n_thread", "integral_order",
                       "kira_executable", "fermat_executable", "run_firefly",
                       "firefly_bunch", "auto_firefly", "kira3_truncate_sp",
                       "kira3_sectorwise", "symbolic_ibp", "ibp_cache_dir",
                       "work_dir", "log_file", "numeric_values"});
    if (obj.contains("ibp_rank")) bb.ibp_rank = obj.at("ibp_rank").get<long>();
    if (obj.contains("ibp_dot")) bb.ibp_dot = obj.at("ibp_dot").get<long>();
    if (obj.contains("n_thread")) bb.n_thread = obj.at("n_thread").get<long>();
    if (obj.contains("integral_order")) {
        bb.integral_order = obj.at("integral_order").get<long>();
    }
    if (obj.contains("kira_executable")) {
        bb.kira_executable = obj.at("kira_executable").get<std::string>();
    }
    if (obj.contains("fermat_executable")) {
        bb.fermat_executable = obj.at("fermat_executable").get<std::string>();
    }
    if (obj.contains("run_firefly")) {
        bb.run_firefly = obj.at("run_firefly").get<bool>();
    }
    if (obj.contains("firefly_bunch")) {
        bb.firefly_bunch = obj.at("firefly_bunch").get<long>();
    }
    // Auto-route the aux-mass (η) Reduce to FireFly (default ON).  See
    // ReduceOptions::auto_firefly: the Fermat backend SIG6-asserts on the
    // η-deformed master tower; FireFly is the memory-light fix, enabled
    // automatically for any reduction whose family carries `eta`.  Set
    // false here (or AMFLOW_NO_AUTO_FIREFLY=1) to force the Fermat path.
    if (obj.contains("auto_firefly")) {
        bb.auto_firefly = obj.at("auto_firefly").get<bool>();
    }
    // Kira 3 seeding levers.  BOTH default OFF (`truncate_sp` regressed on
    // the box1 integral family, `sectorwise` hangs on the aux-mass DE
    // reduce -- see
    // ReduceOptions).  But on the 3-loop crossed-box / two-mass 5-point top, `truncate_sp:true` is
    // the load-bearing fix.  This is a per-family decision, so surface it
    // as an explicit JSON knob (was env-only, which is what bit the 5-pt
    // probe with a Fermat SIG11).
    if (obj.contains("kira3_truncate_sp")) {
        bb.kira3_truncate_sp = obj.at("kira3_truncate_sp").get<bool>();
    }
    if (obj.contains("kira3_sectorwise")) {
        bb.kira3_sectorwise = obj.at("kira3_sectorwise").get<bool>();
    }
    // Layer-B symbolic-kinematics IBP (multi-point cache accelerator).
    // See `ibp_cache.hpp` / `KiraConfig::symbolic_ibp`.  Also overridable
    // via env `AMFLOW_SYMBOLIC_IBP`.
    if (obj.contains("symbolic_ibp")) {
        bb.symbolic_ibp = obj.at("symbolic_ibp").get<bool>();
    }
    // Layer-A content-addressed Kira-results cache root.  Also
    // overridable via env `AMFLOW_IBP_CACHE`; empty ⇒ caching off.
    if (obj.contains("ibp_cache_dir")) {
        bb.ibp_cache_dir = obj.at("ibp_cache_dir").get<std::string>();
    }
    if (obj.contains("work_dir")) {
        bb.work_dir = obj.at("work_dir").get<std::string>();
    }
    if (obj.contains("log_file")) {
        bb.log_file = obj.at("log_file").get<std::string>();
    }
    if (obj.contains("numeric_values")) {
        const auto& nv = obj.at("numeric_values");
        if (!nv.is_object()) {
            throw std::runtime_error(
                "amf_options.blackbox.numeric_values must be an object");
        }
        for (auto it = nv.begin(); it != nv.end(); ++it) {
            const auto& v = it.value();
            const std::string ctx_path =
                "amf_options.blackbox.numeric_values." + it.key();
            // Complex-numeric kinematics is intentionally not
            // implemented; the upstream `ComplexMode` /
            // `CompensateRule` machinery (`Kira/interface.m:50-57`,
            // `488`, `519`) is out of scope.  Reject the object form
            // loudly rather than silently truncating to the real part
            // or returning a wrong answer.
            if (v.is_object()) {
                throw std::runtime_error(
                    ctx_path + ": complex-numeric kinematics "
                    "(object form `{\"re\":..,\"im\":..}`) is not "
                    "implemented; see `docs/AUDIT_MMA_PARITY.md`. "
                    "All numeric_values must be real scalars "
                    "(integer, rational like \"1/100\", or decimal).");
            }
            bb.numeric_values[it.key()] = json_scalar_to_string(v, ctx_path);
        }
    }
}

pipeline::AMFSystemOptions parse_amf_options(const json& obj) {
    if (!obj.is_object()) {
        throw std::runtime_error("amf_options must be an object");
    }
    warn_unknown_keys(obj, "amf_options",
                      {"amf_modes", "ending_schemes", "cache_root",
                       "direction", "max_recursion_depth", "eps_parallel",
                       "node_parallel",
                       "explicit_boundary", "blackbox", "bb"});

    pipeline::AMFSystemOptions out;
    out.bb.kira_executable = "/usr/local/bin/kira";
    out.bb.fermat_executable = "/usr/share/Ferl7/fer64";

    if (obj.contains("amf_modes")) {
        out.amf_modes.clear();
        for (const auto& mode : obj.at("amf_modes")) {
            if (!mode.is_string()) {
                throw std::runtime_error("amf_options.amf_modes entries must be strings");
            }
            out.amf_modes.push_back(qft::parse_amf_mode(mode.get<std::string>()));
        }
    }

    if (obj.contains("ending_schemes")) {
        out.ending_schemes.clear();
        for (const auto& scheme : obj.at("ending_schemes")) {
            if (!scheme.is_string()) {
                throw std::runtime_error(
                    "amf_options.ending_schemes entries must be strings");
            }
            out.ending_schemes.push_back(
                parse_ending_scheme_string(scheme.get<std::string>()));
        }
    }

    if (obj.contains("cache_root")) {
        out.cache_root = obj.at("cache_root").get<std::string>();
    }
    if (obj.contains("direction")) {
        out.direction = obj.at("direction").get<std::string>();
    }
    if (obj.contains("max_recursion_depth")) {
        out.max_recursion_depth = obj.at("max_recursion_depth").get<long>();
    }
    // ε-grid parallelism (Lever 1).  Also overridable via env
    // `AMFLOW_EPS_PARALLEL`.  See AMFSystemOptions::eps_parallel.
    if (obj.contains("eps_parallel")) {
        out.eps_parallel = obj.at("eps_parallel").get<long>();
    }
    // Node-tree parallelism (Lever 2).  Also overridable via env
    // `AMFLOW_NODE_PARALLEL`.  See AMFSystemOptions::node_parallel.
    if (obj.contains("node_parallel")) {
        out.node_parallel = obj.at("node_parallel").get<long>();
    }

    // Optional explicit boundary values for ending-system masters that are
    // not in the Vacuum table (e.g. cut/maximal-cut vacua).  Shape:
    //   "explicit_boundary": {
    //       "<family>|<i1>|...|<iN>": <scalar-or-complex>            // broadcast
    //       "<family>|<i1>|...|<iN>": [ <val>, <val>, ... ]          // per-eps
    //   }
    // A single value (scalar, {"re":..,"im":..}, or a length-1 array) is
    // broadcast to every internal eps sample; a longer array is indexed by
    // eps sample order.  See amfsystem.cpp::try_builtin_ending_value.
    if (obj.contains("explicit_boundary")) {
        const long prec = numeric::working_prec_bits();
        const auto& eb = obj.at("explicit_boundary");
        if (!eb.is_object()) {
            throw std::runtime_error(
                "amf_options.explicit_boundary must be an object "
                "(key -> value or key -> [values])");
        }
        for (auto it = eb.begin(); it != eb.end(); ++it) {
            auto vec = std::make_shared<std::vector<numeric::AcbValue>>();
            const auto& v = it.value();
            if (v.is_array()) {
                if (v.empty()) {
                    throw std::runtime_error(
                        "amf_options.explicit_boundary['" + it.key()
                        + "'] array must be non-empty");
                }
                for (const auto& e : v) vec->push_back(parse_complex(e, prec));
            } else {
                vec->push_back(parse_complex(v, prec));
            }
            out.explicit_boundary[it.key()] = std::move(vec);
        }
    }

    if (obj.contains("blackbox")) {
        apply_blackbox_options(obj.at("blackbox"), out.bb,
                               "amf_options.blackbox");
    }
    if (obj.contains("bb")) {
        apply_blackbox_options(obj.at("bb"), out.bb, "amf_options.bb");
    }
    return out;
}

pipeline::AMFSystemOptions default_solve_integrals_options() {
    pipeline::AMFSystemOptions out;
    out.bb.kira_executable = "/usr/local/bin/kira";
    out.bb.fermat_executable = "/usr/share/Ferl7/fer64";
    return out;
}

void ensure_cli_blackbox_log_file(pipeline::AMFSystemOptions& opts) {
    // A shared default path (`/tmp/amflow_cli_kira.log`) interleaves/
    // truncates across concurrent amflow runs, so monitors read the
    // wrong job's events.  Leave empty: `ibp::kira_run` then defaults
    // to the per-run `<working_dir>/kira_stdout.log` (the wd is
    // job-unique).  A user
    // can still force a path via `amf_options.blackbox.log_file`.
    (void)opts;
}

json jintegral_to_json(const qft::JIntegral& j) {
    return json{
        {"family", j.family()},
        {"indices", j.indices()},
    };
}

bool parse_fmpq(fmpq_t out, const std::string& s) {
    int rc = fmpq_set_str(out, s.c_str(), 10);
    return rc == 0;
}

void load_real_component_from_json(arb_t out,
                                     const json& value,
                                     const std::string& name,
                                     long prec) {
    const std::string text = json_scalar_to_string(value, name);

    fmpq_t q;
    fmpq_init(q);
    if (parse_fmpq(q, text)) {
        arb_set_fmpq(out, q, prec);
        fmpq_clear(q);
        return;
    }
    fmpq_clear(q);

    if (arb_set_str(out, text.c_str(), prec) == 0) {
        return;
    }
    throw std::runtime_error("failed to parse real value '" + text
                               + "' for " + name);
}

FmpqPoly parse_coeff_list(const json& arr) {
    FmpqPoly p;
    long n = 0;
    fmpq_t q;
    fmpq_init(q);
    for (const auto& c : arr) {
        std::string s;
        if (c.is_string()) s = c.get<std::string>();
        else if (c.is_number_integer()) s = std::to_string(c.get<long>());
        else if (c.is_number_float()) s = std::to_string(c.get<double>());
        else throw std::runtime_error("coefficient must be a string or integer");
        if (!parse_fmpq(q, s)) {
            fmpq_clear(q);
            throw std::runtime_error("failed to parse fmpq from '" + s + "'");
        }
        p.set_coeff_fmpq(n++, q);
    }
    fmpq_clear(q);
    return p;
}

RationalFunction parse_rational(const json& obj) {
    if (obj.contains("coeffs")) {
        return RationalFunction::from_polynomial(parse_coeff_list(obj.at("coeffs")));
    }
    if (!obj.contains("num") || !obj.contains("den")) {
        throw std::runtime_error(
            "rational entry must have 'num' and 'den' (or just 'coeffs')");
    }
    FmpqPoly num = parse_coeff_list(obj.at("num"));
    FmpqPoly den = parse_coeff_list(obj.at("den"));
    return RationalFunction(std::move(num), std::move(den));
}

RationalMatrix parse_matrix(const json& arr) {
    if (!arr.is_array() || arr.empty()) {
        throw std::runtime_error("matrix must be a non-empty array of rows");
    }
    std::size_t rows = arr.size();
    std::size_t cols = arr[0].size();
    if (cols != rows) {
        throw std::runtime_error("matrix must be square");
    }
    RationalMatrix out(rows, cols);
    for (std::size_t i = 0; i < rows; ++i) {
        const auto& row = arr[i];
        if (row.size() != cols) {
            throw std::runtime_error("matrix is not rectangular");
        }
        for (std::size_t j = 0; j < cols; ++j) {
            out(i, j) = parse_rational(row[j]);
        }
    }
    return out;
}

AcbValue parse_complex(const json& obj, long prec) {
    AcbValue out;
    if (obj.is_object()) {
        load_real_component_from_json(acb_realref(out.raw()), obj.at("re"),
                                         "complex.re", prec);
        if (obj.contains("im")) {
            load_real_component_from_json(acb_imagref(out.raw()), obj.at("im"),
                                             "complex.im", prec);
        } else {
            arb_zero(acb_imagref(out.raw()));
        }
        return out;
    }

    load_real_component_from_json(acb_realref(out.raw()), obj, "scalar", prec);
    arb_zero(acb_imagref(out.raw()));
    return out;
}

std::vector<AcbValue> parse_eps_samples(const json& obj, long prec) {
    const json* arr = nullptr;
    if (obj.contains("eps_samples")) arr = &obj.at("eps_samples");
    else if (obj.contains("epslist")) arr = &obj.at("epslist");
    else throw std::runtime_error("missing required array 'eps_samples'");

    if (!arr->is_array()) {
        throw std::runtime_error("eps_samples must be an array");
    }

    std::vector<AcbValue> out;
    out.reserve(arr->size());
    for (std::size_t i = 0; i < arr->size(); ++i) {
        out.push_back(parse_complex((*arr)[i], prec));
    }
    return out;
}

std::vector<ode::BoundarySpec> parse_boundaries(const json& arr, long prec) {
    if (!arr.is_array()) throw std::runtime_error("boundaries must be an array");
    std::vector<ode::BoundarySpec> out(arr.size());
    for (std::size_t i = 0; i < arr.size(); ++i) {
        const auto& list = arr[i];
        if (!list.is_array()) throw std::runtime_error("each boundary must be an array");
        for (const auto& entry : list) {
            ode::BoundaryEntry e;
            e.mu    = parse_complex(entry.at("mu"),    prec);
            e.value = parse_complex(entry.at("value"), prec);
            out[i].push_back(std::move(e));
        }
    }
    return out;
}

void apply_options(const json& obj) {
    warn_unknown_keys(obj, "options",
                      {"working_pre", "chop_pre", "rationalize_pre",
                       "silent_mode", "d0", "x_order", "extra_x_order",
                       "learn_x_order", "test_x_order", "run_radius",
                       "run_length", "run_candidate", "run_direction"});
    if (obj.contains("working_pre"))    { auto g = numeric::global_options(); g.working_pre = obj.at("working_pre").get<int>(); numeric::set_global_options(g); }
    if (obj.contains("chop_pre"))       { auto g = numeric::global_options(); g.chop_pre    = obj.at("chop_pre").get<int>();    numeric::set_global_options(g); }
    if (obj.contains("rationalize_pre")){ auto g = numeric::global_options(); g.rationalize_pre = obj.at("rationalize_pre").get<int>(); numeric::set_global_options(g); }
    if (obj.contains("silent_mode"))    { auto g = numeric::global_options(); g.silent_mode = obj.at("silent_mode").get<bool>(); numeric::set_global_options(g); }
    if (obj.contains("d0")) {
        auto g = numeric::global_options();
        const auto& v = obj.at("d0");
        if (v.is_string())             g.d0 = v.get<std::string>();
        else if (v.is_number_integer()) g.d0 = std::to_string(v.get<long>());
        else throw std::runtime_error("d0 must be a rational string ('p/q') or integer");
        fmpq_t q; fmpq_init(q);
        const bool ok = parse_fmpq(q, g.d0);
        fmpq_clear(q);
        if (!ok) throw std::runtime_error("d0: failed to parse rational '" + g.d0 + "'");
        numeric::set_global_options(g);
    }

    if (obj.contains("x_order"))        { auto e = numeric::expansion_options(); e.x_order       = obj.at("x_order").get<int>();       numeric::set_expansion_options(e); }
    if (obj.contains("extra_x_order"))  { auto e = numeric::expansion_options(); e.extra_x_order = obj.at("extra_x_order").get<int>(); numeric::set_expansion_options(e); }
    if (obj.contains("learn_x_order"))  { auto e = numeric::expansion_options(); e.learn_x_order = obj.at("learn_x_order").get<int>(); numeric::set_expansion_options(e); }
    if (obj.contains("test_x_order"))   { auto e = numeric::expansion_options(); e.test_x_order  = obj.at("test_x_order").get<int>();  numeric::set_expansion_options(e); }

    if (obj.contains("run_radius"))     { auto r = numeric::running_options(); r.run_radius    = obj.at("run_radius").get<int>();    numeric::set_running_options(r); }
    if (obj.contains("run_length"))     { auto r = numeric::running_options(); r.run_length    = obj.at("run_length").get<int>();    numeric::set_running_options(r); }
    if (obj.contains("run_candidate"))  { auto r = numeric::running_options(); r.run_candidate = obj.at("run_candidate").get<int>(); numeric::set_running_options(r); }
    if (obj.contains("run_direction"))  {
        auto r = numeric::running_options();
        std::string s = obj.at("run_direction").get<std::string>();
        if      (s == "Re")    r.run_direction = numeric::RunningOptions::Direction::Re;
        else if (s == "Im")    r.run_direction = numeric::RunningOptions::Direction::Im;
        else if (s == "NegRe") r.run_direction = numeric::RunningOptions::Direction::NegRe;
        else if (s == "NegIm") r.run_direction = numeric::RunningOptions::Direction::NegIm;
        else throw std::runtime_error("unknown run_direction '" + s + "'");
        numeric::set_running_options(r);
    }
}

// Precision-output fix: display digits must not silently
// truncate high-precision results.  arb_get_str with flags=0 is radius-aware
// (never prints digits the ball does not certify), so over-asking is safe.
// Default was a hardcoded 30, which capped every CLI result at ~30 digits
// regardless of goal_digits / working_pre.
long output_display_digits = 30;

json acb_to_json(const AcbValue& v) {
    const long digits = std::max<long>(output_display_digits, 30);
    json o;
    char* re = arb_get_str(acb_realref(v.raw()), digits, 0);
    char* im = arb_get_str(acb_imagref(v.raw()), digits, 0);
    o["re"] = re ? re : "";
    o["im"] = im ? im : "";
    if (re) flint_free(re);
    if (im) flint_free(im);
    return o;
}

json options_to_json() {
    json o;
    o["working_pre"]     = numeric::global_options().working_pre;
    o["chop_pre"]        = numeric::global_options().chop_pre;
    o["rationalize_pre"] = numeric::global_options().rationalize_pre;
    o["silent_mode"]     = numeric::global_options().silent_mode;
    o["d0"]              = numeric::global_options().d0;
    o["x_order"]         = numeric::expansion_options().x_order;
    o["extra_x_order"]   = numeric::expansion_options().extra_x_order;
    o["learn_x_order"]   = numeric::expansion_options().learn_x_order;
    o["test_x_order"]    = numeric::expansion_options().test_x_order;
    o["run_radius"]      = numeric::running_options().run_radius;
    o["run_length"]      = numeric::running_options().run_length;
    o["run_candidate"]   = numeric::running_options().run_candidate;
    switch (numeric::running_options().run_direction) {
        case numeric::RunningOptions::Direction::Re:    o["run_direction"] = "Re"; break;
        case numeric::RunningOptions::Direction::Im:    o["run_direction"] = "Im"; break;
        case numeric::RunningOptions::Direction::NegRe: o["run_direction"] = "NegRe"; break;
        case numeric::RunningOptions::Direction::NegIm: o["run_direction"] = "NegIm"; break;
        case numeric::RunningOptions::Direction::Custom: o["run_direction"] = "Custom"; break;
    }
    return o;
}

json laurent_solutions_to_json(
        const std::vector<pipeline::LaurentIntegralSolution>& sol) {
    json out = json::array();
    for (const auto& row : sol) {
        json item;
        item["integral"] = jintegral_to_json(row.integral);
        item["leading_order"] = row.leading_order;
        item["coefficients"] = json::array();
        for (std::size_t i = 0; i < row.coefficients.size(); ++i) {
            item["coefficients"].push_back(json{
                {"order", row.leading_order + static_cast<long>(i)},
                {"value", acb_to_json(row.coefficients[i])},
            });
        }
        out.push_back(std::move(item));
    }
    return out;
}

json sampled_solutions_to_json(
        const std::vector<pipeline::SampledIntegralSolution>& sol,
        const std::vector<AcbValue>& eps_samples) {
    json out = json::array();
    for (const auto& row : sol) {
        if (row.values.size() != eps_samples.size()) {
            throw std::runtime_error(
                "sampled_solutions_to_json: epsilon/value grid size mismatch");
        }

        json item;
        item["integral"] = jintegral_to_json(row.integral);
        item["samples"] = json::array();
        for (std::size_t i = 0; i < row.values.size(); ++i) {
            item["samples"].push_back(json{
                {"eps", acb_to_json(eps_samples[i])},
                {"value", acb_to_json(row.values[i])},
            });
        }
        out.push_back(std::move(item));
    }
    return out;
}

}  // namespace

json run_json(const json& input) {
    if (input.contains("options")) apply_options(input.at("options"));
    std::string mode = input.value("mode", std::string("amflow"));

    // Config hygiene: top-level unrecognized-key check, mode-aware (a key
    // that only another mode reads gets the known-elsewhere hint).  Unknown
    // modes fall through to the loud throw at the end of the dispatch.
    {
        std::vector<std::string> accepted = {"mode", "options"};
        bool known_mode = true;
        if (mode == "amflow") {
            accepted.insert(accepted.end(), {"matrix", "boundaries"});
        } else if (mode == "solve_integrals") {
            accepted.insert(accepted.end(),
                            {"family", "integrals", "jints", "targets",
                             "goal_digits", "goal", "eps_order", "order",
                             "target_eps_power", "amf_options", "work_dir"});
        } else if (mode == "black_box_amflow") {
            accepted.insert(accepted.end(),
                            {"family", "integrals", "jints", "targets",
                             "eps_samples", "epslist", "amf_options",
                             "work_dir"});
        } else if (mode == "diffeq") {
            accepted.insert(accepted.end(),
                            {"family", "integrals", "jints", "targets",
                             "variables", "top_pattern", "amf_options",
                             "work_dir"});
        } else {
            known_mode = false;
        }
        if (known_mode) {
            warn_unknown_keys(input, "top-level (mode:" + mode + ")",
                              accepted, /*dedup_path=*/"top-level");
        }
    }

    json out;
    out["mode"] = mode;

    if (mode == "amflow") {
        long prec = numeric::working_prec_bits();
        RationalMatrix de = parse_matrix(input.at("matrix"));
        std::vector<ode::BoundarySpec> bcs =
            parse_boundaries(input.at("boundaries"), prec);
        if (bcs.size() != de.rows()) {
            throw std::runtime_error(
                "boundaries length must match matrix size; got "
                + std::to_string(bcs.size()) + " vs "
                + std::to_string(de.rows()));
        }

        auto sol = ode::amflow(de, bcs, prec);
        out["result"] = json::array();
        for (const auto& v : sol) out["result"].push_back(acb_to_json(v));
    } else if (mode == "solve_integrals") {
        auto fc = parse_family_config(input.at("family"));
        auto targets = parse_jintegral_list(input, fc);
        const long goal_digits =
            parse_required_long(input, "goal_digits", "goal");
        // Either the AMFlow-convention `eps_order` (= K + 2L, highest fitted
        // power eps^{eps_order-2L}) OR the physical `target_eps_power` K.
        // The latter is the unambiguous one; if both are given they must
        // agree.  This removes the eps_order = K+2L pitfall entirely.
        const long L = static_cast<long>(fc.n_loops());
        long eps_order;
        if (input.contains("target_eps_power")) {
            const long K = input.at("target_eps_power").get<long>();
            eps_order = K + 2 * L;
            if (input.contains("eps_order") || input.contains("order")) {
                const long given =
                    parse_required_long(input, "eps_order", "order");
                if (given != eps_order) {
                    throw std::runtime_error(
                        "solve_integrals: target_eps_power=" + std::to_string(K)
                        + " (=> eps_order=" + std::to_string(eps_order)
                        + " for " + std::to_string(L) + "-loop) conflicts with "
                        "eps_order=" + std::to_string(given));
                }
            }
        } else {
            eps_order = parse_required_long(input, "eps_order", "order");
        }
        // Precision-output fix: print what was asked for (radius-aware, see acb_to_json)
        output_display_digits =
            std::max<long>(goal_digits + 10,
                           numeric::global_options().working_pre + 10);
        pipeline::AMFSystemOptions opts = input.contains("amf_options")
            ? parse_amf_options(input.at("amf_options"))
            : default_solve_integrals_options();
        ensure_cli_blackbox_log_file(opts);
        const std::string work_dir =
            input.value("work_dir", std::string{});

        auto sol = pipeline::solve_integrals(fc, targets, goal_digits, eps_order,
                                              opts, work_dir);
        out["result"] = laurent_solutions_to_json(sol);
    } else if (mode == "black_box_amflow") {
        const long prec = numeric::working_prec_bits();
        // Precision-output fix: see acb_to_json
        output_display_digits = numeric::global_options().working_pre + 10;
        auto fc = parse_family_config(input.at("family"));
        auto targets = parse_jintegral_list(input, fc);
        auto eps_samples = parse_eps_samples(input, prec);
        pipeline::AMFSystemOptions opts = input.contains("amf_options")
            ? parse_amf_options(input.at("amf_options"))
            : default_solve_integrals_options();
        ensure_cli_blackbox_log_file(opts);
        const std::string work_dir =
            input.value("work_dir", std::string{});

        auto sol = pipeline::black_box_amflow(fc, targets, eps_samples, opts, work_dir);
        out["result"] = sampled_solutions_to_json(sol, eps_samples);
    } else if (mode == "diffeq") {
        // Kinematic-variable differential-equation matrix at a numeric
        // point.  Exposes `ibp::diffeq` (src/ibp/reduce.cpp) — which now,
        // via the `LIBPDerivivative` momentum-derivative chain rule
        // (src/ibp/libp_deriv.cpp), supports replacement-defined
        // invariants `s, t, sᵢⱼ, …` and not just `eta`.  Consumed by
        // `tools/wayfinder/sampler.py` (`backend="amflow_diffeq"`) so the
        // connection-once route can sample ∂J/∂x_i exactly (one Kira IBP
        // per point) instead of central-differencing 2·|kinvars|+1
        // `solve_integrals` calls per point at ~⅔·goal_digits.
        //
        //   input:  family, integrals (= preferred masters), variables,
        //           amf_options.blackbox.{numeric_values, work_dir, …}
        //   output: masters (sorted Kira basis), vars, diffeq[var][i][j]
        //           = (∂Iᵢ/∂var)|_masters as Mfrac.to_string() in `d`.
        auto fc = parse_family_config(input.at("family"));
        auto preferred = parse_jintegral_list(input, fc);
        if (!input.contains("variables")) {
            throw std::runtime_error("mode:diffeq requires 'variables'");
        }
        auto vars = parse_string_array(input.at("variables"), "variables");
        pipeline::AMFSystemOptions amf_opts = input.contains("amf_options")
            ? parse_amf_options(input.at("amf_options"))
            : default_solve_integrals_options();
        ensure_cli_blackbox_log_file(amf_opts);
        ibp::ReduceOptions bb = amf_opts.bb;
        if (bb.work_dir.empty()) {
            bb.work_dir = input.value("work_dir", std::string{"/tmp/amflow_diffeq"});
        }
        std::vector<int> top = input.contains("top_pattern")
            ? parse_int_array(input.at("top_pattern"), "top_pattern")
            : qft::get_top_sector(preferred);

        auto de = ibp::diffeq(fc, preferred, vars, top, bb);

        out["masters"] = json::array();
        for (const auto& m : de.sortedmasters) {
            out["masters"].push_back(jintegral_to_json(m));
        }
        out["vars"] = de.vars;
        json de_json = json::object();
        for (std::size_t v = 0; v < de.vars.size(); ++v) {
            json mat = json::array();
            for (std::size_t r = 0; r < de.diffeq[v].size(); ++r) {
                json row = json::array();
                for (std::size_t c = 0; c < de.diffeq[v][r].size(); ++c) {
                    row.push_back(de.diffeq[v][r][c].to_string());
                }
                mat.push_back(std::move(row));
            }
            de_json[de.vars[v]] = std::move(mat);
        }
        out["diffeq"] = std::move(de_json);
    } else {
        throw std::runtime_error("unknown mode '" + mode + "'");
    }
    out["options"] = options_to_json();

    return out;
}

}  // namespace amflow::api
