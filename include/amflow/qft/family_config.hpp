// SPDX-License-Identifier: MIT
// qft::family_config — FamilyConfig and propagator algebra.
//
//
// Mirrors AMFlow.m:
//   * AMFlowInfo[...] global state                      (lines 188-212)
//   * SPList / DListSymbol / CheckCompleteness          (lines 398-422)
//   * ToCompleteExplicit / ToSquare / ToSquareAll       (lines 425-443)
//
// Design
// ------
//   The Mathematica package keeps everything in process-global state via
//   AMFlowInfo[].  We bundle the per-family data into a single immutable
//   FamilyConfig object the user constructs once and passes around.
//
//   The single shared MpolyContext exposes, in this order:
//     loop names         (n_loop)
//     reduced-leg names  (n_red_leg)            (legs minus conservation keys)
//     invariant names    (n_inv)                (e.g. s, t, msq)

#ifndef AMFLOW_QFT_FAMILY_CONFIG_HPP
#define AMFLOW_QFT_FAMILY_CONFIG_HPP

#include <cstddef>
#include <map>
#include <memory>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "amflow/algebra/mpoly.hpp"

namespace amflow::qft {

class FamilyConfig {
public:
    std::shared_ptr<algebra::MpolyContext> ctx;

    std::string                       family;
    std::vector<std::string>          loops;
    std::vector<std::string>          legs;
    std::vector<std::string>          reduced_legs;
    std::vector<std::string>          invariants;

    std::map<std::string, algebra::Mpoly> conservation;
    std::map<std::string, algebra::Mfrac> reduced_replacement;

    std::vector<algebra::Mpoly>       propagators_after_conservation;
    std::vector<std::string>          propagator_denominators;

    std::vector<algebra::Mpoly>       sp_list;

    std::vector<int>                  cut;            // length = n_propagators
    std::vector<int>                  prescription;   // length = n_loops

    static FamilyConfig build(
        std::string family,
        std::vector<std::string> loops,
        std::vector<std::string> legs,
        std::vector<std::pair<std::string, std::string>> conservation,
        std::vector<std::pair<std::string, std::string>> replacement,
        std::vector<std::string> propagator_strings,
        std::vector<int> cut          = {},
        std::vector<int> prescription = {});

    FamilyConfig(const FamilyConfig&) = delete;
    FamilyConfig& operator=(const FamilyConfig&) = delete;
    FamilyConfig(FamilyConfig&&) noexcept = default;
    FamilyConfig& operator=(FamilyConfig&&) noexcept = default;

    std::size_t n_loops()       const noexcept { return loops.size(); }
    std::size_t n_red_legs()    const noexcept { return reduced_legs.size(); }
    std::size_t n_invariants()  const noexcept { return invariants.size(); }
    std::size_t n_propagators() const noexcept { return propagators_after_conservation.size(); }
    std::size_t n_sp()          const noexcept { return sp_list.size(); }

    long var_index(const std::string& name) const { return ctx->var_index(name); }

    int prescription_of_loop(long loop_idx) const noexcept;
    std::optional<int> prescription_of_prop(const algebra::Mpoly& prop) const;

    algebra::Mpoly apply_conservation(const algebra::Mpoly& p) const;
    algebra::Mfrac apply_replacement(const algebra::Mfrac& p) const;
    algebra::Mfrac apply_replacement(const algebra::Mpoly& p) const;

private:
    FamilyConfig() = default;

    void build_context();
    void build_sp_list();
};

}  // namespace amflow::qft

#endif  // AMFLOW_QFT_FAMILY_CONFIG_HPP
