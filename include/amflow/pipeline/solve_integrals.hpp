// SPDX-License-Identifier: MIT
#ifndef AMFLOW_PIPELINE_SOLVE_INTEGRALS_HPP
#define AMFLOW_PIPELINE_SOLVE_INTEGRALS_HPP

#include <string>
#include <vector>

#include "amflow/numeric/acb_wrap.hpp"
#include "amflow/pipeline/amfsystem.hpp"
#include "amflow/qft/family_config.hpp"
#include "amflow/qft/jintegral.hpp"

namespace amflow::pipeline {

// Layer 17 (partial): top-level AMFlow helpers that sit above AMFSystem.
//
// We start with the two pieces that are independent of Kira orchestration:
//   * GenerateNumericalConfig[goal, order]
//   * FitEps[epslist, values, leading]
//
// The full SolveIntegrals / BlackBoxAMFlow surface will build on these.

struct NumericalConfig {
    std::vector<numeric::AcbValue> eps_samples;
    long working_pre = 0;   // decimal digits
    long x_order     = 0;
};

struct SampledIntegralSolution {
    qft::JIntegral integral;
    std::vector<numeric::AcbValue> values;  // one entry per requested eps sample
};

struct LaurentIntegralSolution {
    qft::JIntegral integral;
    long leading_order = 0;

    // coefficients[k] multiplies eps^(leading_order + k).
    std::vector<numeric::AcbValue> coefficients;
};

// Generate Mathematica-style epsilon sampling for a family with
// `loop_count` loops, targeting `goal_digits` digits up to Laurent order
// `eps_order`.
//
// Mirrors:
//   GenerateNumericalConfig[goal_, order_]
//
// in AMFlow.m, except that the loop count is explicit instead of being
// pulled from Mathematica's global Loop symbol.
NumericalConfig
generate_numerical_config(long loop_count,
                          long goal_digits,
                          long eps_order);

// Fit a Laurent expansion
//   sum_{k=0}^{n-1} c_k * eps^(leading_order + k)
// through the supplied sample values.
//
// `eps_samples.size()` must match `values.size()`.  The number of returned
// coefficients equals the number of samples, mirroring Mathematica's
// Fit[..., eps^(leading + Range[0, n-1]), eps].
std::vector<numeric::AcbValue>
fit_eps(const std::vector<numeric::AcbValue>& eps_samples,
        const std::vector<numeric::AcbValue>& values,
        long leading_order,
        long prec = numeric::working_prec_bits());

// Layer-17 orchestration mirroring AMFlow.m:
//
//   * BlackBoxAMFlowSingle[jints, epslist, dir]
//   * BlackBoxAMFlow[jints, epslist, dir]
//   * SolveIntegrals[jints, goal, order]
//
// `work_dir` plays the role of Mathematica's cache directory.  When left
// empty, the implementation derives a temporary root automatically.
std::vector<SampledIntegralSolution>
black_box_amflow_single(const qft::FamilyConfig& fc,
                        const std::vector<qft::JIntegral>& jints,
                        const std::vector<numeric::AcbValue>& eps_samples,
                        const AMFSystemOptions& opts,
                        const std::string& work_dir = "");

std::vector<SampledIntegralSolution>
black_box_amflow(const qft::FamilyConfig& fc,
                 const std::vector<qft::JIntegral>& jints,
                 const std::vector<numeric::AcbValue>& eps_samples,
                 const AMFSystemOptions& opts,
                 const std::string& work_dir = "");

std::vector<LaurentIntegralSolution>
solve_integrals(const qft::FamilyConfig& fc,
                const std::vector<qft::JIntegral>& jints,
                long goal_digits,
                long eps_order,
                const AMFSystemOptions& opts,
                const std::string& work_dir = "");

}  // namespace amflow::pipeline

#endif  // AMFLOW_PIPELINE_SOLVE_INTEGRALS_HPP
