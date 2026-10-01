// SPDX-License-Identifier: MIT
// ode::path — integration contour generation.
//
//
// Mirrors AMFlow.m / DESolver.m lines 519-586.

#ifndef AMFLOW_ODE_PATH_HPP
#define AMFLOW_ODE_PATH_HPP

#include <cstddef>
#include <vector>

#include "amflow/numeric/options.hpp"
#include "amflow/numeric/rational.hpp"

namespace amflow::ode {

numeric::RationalComplex
first_step(const std::vector<numeric::RationalComplex>& poles,
           long prec = numeric::working_prec_bits());

numeric::RationalComplex
last_step(const std::vector<numeric::RationalComplex>& poles,
          long prec = numeric::working_prec_bits());

std::vector<numeric::RationalComplex>
run_unit(const std::vector<numeric::RationalComplex>& poles,
         long prec = numeric::working_prec_bits());

std::vector<numeric::RationalComplex>
run_segment(const std::vector<numeric::RationalComplex>& poles,
            const numeric::RationalComplex& ini,
            const numeric::RationalComplex& fin,
            long prec = numeric::working_prec_bits());

std::vector<numeric::RationalComplex>
run_eta_direction_positive(const std::vector<numeric::RationalComplex>& poles,
                           long prec = numeric::working_prec_bits());

std::vector<numeric::RationalComplex>
run_eta_direction(const std::vector<numeric::RationalComplex>& poles,
                  const numeric::RationalComplex& direction,
                  long prec = numeric::working_prec_bits());

std::vector<numeric::RationalComplex>
run_eta_direction(const std::vector<numeric::RationalComplex>& poles,
                  numeric::RunningOptions::Direction mode,
                  long prec = numeric::working_prec_bits());

std::vector<numeric::RationalComplex>
run_eta(const std::vector<numeric::RationalComplex>& poles,
        long prec = numeric::working_prec_bits());

}  // namespace amflow::ode

#endif  // AMFLOW_ODE_PATH_HPP
