// SPDX-License-Identifier: MIT
// qft::region — regions, region rules, FindAllRegion, RegionPower,
//                ToCompleteExplicit, ToSquareAll.
//

#ifndef AMFLOW_QFT_REGION_HPP
#define AMFLOW_QFT_REGION_HPP

#include <memory>
#include <string>
#include <vector>

#include "amflow/algebra/mpoly.hpp"
#include "amflow/algebra/mpoly_matrix.hpp"
#include "amflow/qft/family_config.hpp"
#include "amflow/qft/jintegral.hpp"

namespace amflow::qft {

// ---------------------------------------------------------------------------
//  RegionContext
// ---------------------------------------------------------------------------

struct RegionContext {
    std::shared_ptr<algebra::MpolyContext> ctx;
    long eta_var;
    long half_eta_var;
    std::vector<algebra::Mpoly> lift_gens;
};

RegionContext make_region_context(const FamilyConfig& fc);

algebra::Mpoly lift_to_region(const algebra::Mpoly& src, const RegionContext& rctx);

// ---------------------------------------------------------------------------
//  branch_momenta / BranchToLoop / RegionRule / BranchScale
// ---------------------------------------------------------------------------

std::vector<algebra::Mpoly>
branch_momenta(const FamilyConfig& fc,
               const std::vector<algebra::Mpoly>& denominators);

struct LoopTransform {
    bool ok;
    std::vector<algebra::Mfrac> map;
};

LoopTransform branch_to_loop(const FamilyConfig& fc,
                              const RegionContext& rctx,
                              const std::vector<algebra::Mpoly>& candidates);

std::vector<algebra::Mfrac>
region_rule(const RegionContext& rctx,
            const LoopTransform& tr,
            const std::vector<int>& scale);

algebra::Mfrac apply_region_rule(const FamilyConfig& fc,
                                  const RegionContext& rctx,
                                  const algebra::Mfrac& expr,
                                  const std::vector<algebra::Mfrac>& rule);

algebra::Mfrac apply_region_rule(const FamilyConfig& fc,
                                  const RegionContext& rctx,
                                  const algebra::Mpoly& expr,
                                  const std::vector<algebra::Mfrac>& rule);

std::vector<int>
branch_scale(const FamilyConfig& fc,
             const std::vector<algebra::Mfrac>& transformed_branches,
             const std::vector<int>& scale);

// ---------------------------------------------------------------------------
//  Region + FindAllRegion + ZeroRegionQ + RegionPower
// ---------------------------------------------------------------------------

struct Region {
    LoopTransform     transform;
    std::vector<int>  scale;
};

std::vector<Region>
find_all_region(const FamilyConfig& fc,
                const RegionContext& rctx,
                const std::vector<std::size_t>& topposi);

bool zero_region_q(const FamilyConfig& fc,
                   const RegionContext& rctx,
                   const Region& region,
                   const std::vector<std::size_t>& topposi);

struct PowersContext {
    std::shared_ptr<algebra::MpolyContext> ctx;
    long eps_var;
    std::vector<algebra::Mpoly> lift_gens;
};

PowersContext make_powers_context(const FamilyConfig& fc);

std::vector<algebra::Mfrac>
region_power(const FamilyConfig& fc,
             const RegionContext& rctx,
             const PowersContext& pctx,
             const Region& region,
             const std::vector<JIntegral>& integrals);

// ---------------------------------------------------------------------------
//  Layer 14d: ToCompleteExplicit + SPListToDListSymbol + SquaredDenominators
// ---------------------------------------------------------------------------

struct DListContext {
    std::shared_ptr<algebra::MpolyContext> ctx;
    long                                    first_d_var;
    long                                    n_d;
    std::vector<algebra::Mpoly>             lift_gens;
};

DListContext make_dlist_context(const FamilyConfig& fc,
                                  std::size_t n_denominators);

std::vector<algebra::Mfrac>
sp_list_to_dlist_symbol(const FamilyConfig& fc,
                          const DListContext& dctx,
                          const std::vector<algebra::Mpoly>& denominators);

std::vector<algebra::Mpoly>
to_complete_explicit(const FamilyConfig& fc,
                       const std::vector<algebra::Mpoly>& denominators);

std::vector<algebra::Mpoly>
squared_denominators(const FamilyConfig& fc,
                       const std::vector<algebra::Mpoly>& denominators);

struct ToSquareResult {
    std::vector<algebra::Mfrac> momenta;
    std::vector<algebra::Mfrac> masses;
};

ToSquareResult
to_square_all(const FamilyConfig& fc,
                const std::vector<algebra::Mpoly>& denominators);

// Internal helper used between region.cpp and complete.cpp / findregion.cpp.
void coeff_over_splist(const algebra::Mpoly& den,
                        const std::vector<algebra::Mpoly>& sp_list,
                        std::vector<fmpq*>& out);

}  // namespace amflow::qft

#endif  // AMFLOW_QFT_REGION_HPP
