// SPDX-License-Identifier: MIT
// qft::boundary — BoundaryPattern + ApartRationals + BoundaryIntegrands +
//                  LaportaIntegrals + BoundaryIntegrals.
//
//
// Mirrors AMFlow.m lines 680-812.

#ifndef AMFLOW_QFT_BOUNDARY_HPP
#define AMFLOW_QFT_BOUNDARY_HPP

#include <memory>
#include <string>
#include <vector>

#include "amflow/algebra/mpoly.hpp"
#include "amflow/qft/family_config.hpp"
#include "amflow/qft/jintegral.hpp"
#include "amflow/qft/region.hpp"

namespace amflow::qft {

std::vector<std::vector<algebra::Mfrac>>
boundary_pattern(const std::vector<std::vector<algebra::Mfrac>>& powers);

std::vector<algebra::Mfrac>
apart_one_var(const algebra::Mpoly& poly, long var);

std::vector<std::vector<algebra::Mfrac>>
apart_rationals(const std::vector<algebra::Mfrac>& rationals,
                const DListContext& dctx);

struct BoundaryIntegrandsResult {
    std::vector<algebra::Mpoly>              completede;
    std::vector<std::vector<algebra::Mfrac>> integrands;
    std::shared_ptr<algebra::MpolyContext>   dctx;
    long                                      first_d_var;
};

BoundaryIntegrandsResult
boundary_integrands(const FamilyConfig&             fc,
                    const RegionContext&            rctx,
                    const std::vector<JIntegral>&   integrals,
                    const std::vector<long>&        border,
                    const Region&                   region);

struct LaportaTerm {
    std::vector<long>  indices;
    algebra::Mfrac     coeff;
};

std::vector<LaportaTerm>
laporta_integrals(const algebra::Mfrac& term, const DListContext& dctx);

struct BoundaryFamily {
    std::vector<algebra::Mpoly>                          prop;
    std::vector<std::vector<std::vector<LaportaTerm>>>   terms;
};

std::vector<BoundaryFamily>
boundary_integrals(const FamilyConfig& fc,
                   const RegionContext& rctx,
                   const BoundaryIntegrandsResult& bi);

}  // namespace amflow::qft

#endif  // AMFLOW_QFT_BOUNDARY_HPP
