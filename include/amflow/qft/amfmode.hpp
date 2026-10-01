// SPDX-License-Identifier: MIT
// qft::amfmode — eta-injection strategies (AMFMode).
//
//
// Mirrors AMFlow.m lines 528-622.

#ifndef AMFLOW_QFT_AMFMODE_HPP
#define AMFLOW_QFT_AMFMODE_HPP

#include <cstddef>
#include <string>
#include <vector>

#include "amflow/algebra/mpoly.hpp"
#include "amflow/qft/family_config.hpp"
#include "amflow/qft/topology.hpp"

namespace amflow::qft {

enum class AMFMode {
    Prescription,
    Mass,
    Propagator,
    Branch,
    Loop,
    All,
    // C++-port extension (not in upstream AMFlow.m): Prescription-style
    // union, but EXCLUDING propagators that form a 1-loop bubble
    // sub-integral (a loop momentum carried by exactly two top-sector
    // propagators).  Keeps the bubble at fixed mass so the η-DE Kira
    // reduction stays at (L-1)-loop size — the same mechanism that
    // already makes a 1-prop tadpole sub-loop cheap (there the bubble
    // analogue is a U-polynomial factor and the vacuum-component skip
    // in `amf_candidate` drops it automatically; the 2-prop bubble's U
    // does NOT factor so it needs this explicit exclusion).
    Bubble,
};

AMFMode parse_amf_mode(const std::string& name);
const char* amf_mode_name(AMFMode m) noexcept;

// ComponentInfo augmented with per-Feynman-parameter cut & prescription.
struct TopSectorComponentInfo {
    algebra::Mfrac           u0;
    std::vector<long>        var;        // uf_ctx indices (Feynman params)
    std::vector<std::size_t> prop_index; // global propagator index per var
    long                     loopnum;
    std::vector<algebra::Mfrac> mass;
    bool                     vacQ;
    std::vector<int>         cut;
    std::vector<int>         pres;
};

std::vector<TopSectorComponentInfo>
analyze_top_sector(const FamilyConfig& fc,
                   const std::vector<std::size_t>& topposi);

bool vacuum_q     (const TopSectorComponentInfo& info);
bool single_mass_q(const TopSectorComponentInfo& info);
bool phase_volume_q(const TopSectorComponentInfo& info);
bool ending_q     (const TopSectorComponentInfo& info);

std::vector<std::vector<std::size_t>>
all_possible_position(const FamilyConfig& fc,
                      const TopSectorComponentInfo& info,
                      AMFMode mode);

std::vector<std::size_t>
amf_candidate_component(const FamilyConfig& fc,
                        const TopSectorComponentInfo& info,
                        AMFMode mode);

std::vector<std::size_t>
amf_candidate(const FamilyConfig& fc,
              const std::vector<TopSectorComponentInfo>& info_list,
              AMFMode mode);

std::vector<std::size_t>
amf_position(const FamilyConfig& fc,
             const std::vector<std::size_t>& topposi,
             AMFMode mode);

std::vector<std::size_t>
amf_position(const FamilyConfig& fc,
             const std::vector<std::size_t>& topposi,
             const std::vector<AMFMode>& modes);

std::vector<int>
amf_eta_c(const FamilyConfig& fc,
          const std::vector<std::size_t>& posi);

}  // namespace amflow::qft

#endif  // AMFLOW_QFT_AMFMODE_HPP
