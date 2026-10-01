// SPDX-License-Identifier: MIT
// amflow::factorize — FactorizeFamily decomposition.
//
//
// Mirrors AMFlow.m lines 488-525.
//
// The Symanzik U-polynomial of a multi-loop family factors over the
// ring of Feynman parameters into independent connected components.
// Each component corresponds to a vacuum-like sub-family that shares
// no propagator with the others.

#ifndef AMFLOW_PIPELINE_FACTORIZE_HPP
#define AMFLOW_PIPELINE_FACTORIZE_HPP

#include <memory>
#include <string>
#include <vector>

#include "amflow/algebra/mpoly.hpp"
#include "amflow/qft/family_config.hpp"

namespace amflow::pipeline {

struct FactorizedComponent {
    std::vector<std::string>           loops;
    std::vector<algebra::Mpoly>        propagators;
    std::vector<std::vector<long>>     patterns;
    std::vector<long>                   original_positions;
};

std::vector<FactorizedComponent>
factorize_family(const qft::FamilyConfig& fc,
                 const std::vector<algebra::Mpoly>& denominators,
                 const std::vector<std::vector<long>>& patts);

}  // namespace amflow::pipeline

#endif  // AMFLOW_PIPELINE_FACTORIZE_HPP
