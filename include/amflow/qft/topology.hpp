// SPDX-License-Identifier: MIT
// qft::topology — AnalyzeComponent / AnalyzeTopology / ZeroSectorQ.
//
//
// Mirrors AMFlow.m: ZeroSectorQ, AnalyzeComponent, AnalyzeTopology.

#ifndef AMFLOW_QFT_TOPOLOGY_HPP
#define AMFLOW_QFT_TOPOLOGY_HPP

#include <memory>
#include <vector>

#include "amflow/algebra/mpoly.hpp"
#include "amflow/qft/family_config.hpp"
#include "amflow/qft/family_uf.hpp"

namespace amflow::qft {

// Scaleless sector test (Euler-scaling consistency).  Returns true iff the
// sector is scaleless (i.e. all integrals vanish).
bool zero_sector_q(const FamilyConfig& fc,
                   const std::vector<algebra::Mpoly>& denominators);

// One factor of U after factorization.
struct ComponentInfo {
    algebra::Mfrac    u0;       // factor itself, lifted into Mfrac
    std::vector<long> var;      // Feynman-parameter indices (uf_ctx-relative)
    long              loopnum;  // # distinct vars in leading monomial
    std::vector<algebra::Mfrac> mass;   // Coefficient[f0, x_var[k]]
    bool              vacQ;     // f restricted to var has no high-rank monomial
};

std::vector<ComponentInfo>
analyze_topology(const FamilyConfig& fc,
                 const std::vector<algebra::Mpoly>& denominators);

// Distinct Feynman-parameter indices in `p` (sorted ascending).
std::vector<long>
feynman_vars_in(const algebra::Mfrac& p, long first_x_var, long n_x);
std::vector<long>
feynman_vars_in(const algebra::Mpoly& p, long first_x_var, long n_x);

// Factorize Mpoly. Constant factors come first (only emitted when != 1).
// Each non-constant factor with multiplicity k appears k times.
std::vector<algebra::Mpoly>
factor_mpoly(const algebra::Mpoly& p);

}  // namespace amflow::qft

#endif  // AMFLOW_QFT_TOPOLOGY_HPP
