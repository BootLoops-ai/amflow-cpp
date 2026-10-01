// SPDX-License-Identifier: MIT
// ode::amflow — System handle and the top-level AMFlow pipeline.
//
//
// Mirrors AMFlow.m / DESolver.m lines 985-1138.
//
// We model the per-system state (DE, BC, P, AsyExp) as a single object
// instead of Mathematica's `LoadSystem[sysid, ...]` global lookup-table
// pattern.  The methods inf_to_regular / regular_run / regular_interpolation
// / solve_asy_exp / pick_zero_solution mirror their Mathematica counterparts.
//
// The top-level pipeline `amflow(de, bcs)` wires them together: build path
// with run_eta(get_poles(de)), solve at infinity, propagate along the path,
// solve at zero, and pick the leading numerical value per integral.

#ifndef AMFLOW_ODE_AMFLOW_HPP
#define AMFLOW_ODE_AMFLOW_HPP

#include <cstddef>
#include <stdexcept>
#include <vector>

#include "amflow/numeric/acb_wrap.hpp"
#include "amflow/numeric/options.hpp"
#include "amflow/numeric/rational.hpp"
#include "amflow/ode/asy.hpp"
#include "amflow/ode/inf.hpp"

namespace amflow::ode {

// ---------------------------------------------------------------------------
//  System
// ---------------------------------------------------------------------------
//
//  The boundary condition starts as an asymptotic specification (one
//  BoundarySpec per integral, as accepted by CalcInf) when point == infinity.
//  After inf_to_regular() it switches to a numerical vector (one AcbValue
//  per integral) and `point` becomes a regular complex point.  After
//  solve_asy_exp() the asymptotic expansion at zero is recorded.

class System {
public:
    enum class BCKind {
        Singular,
        Regular
    };

    System(numeric::RationalMatrix de,
           std::vector<BoundarySpec> bc,
           bool point_at_infinity = true);

    System(numeric::RationalMatrix de,
           std::vector<numeric::AcbValue> bc_regular,
           numeric::RationalComplex point);

    void inf_to_regular(const numeric::RationalComplex& x0,
                        long prec = numeric::working_prec_bits());

    void regular_run(const std::vector<numeric::RationalComplex>& run,
                     long prec = numeric::working_prec_bits());

    std::vector<std::vector<numeric::AcbValue>>
    regular_interpolation(const std::vector<numeric::RationalComplex>& samples,
                          long prec = numeric::working_prec_bits());

    void solve_asy_exp(long prec = numeric::working_prec_bits());

    std::vector<numeric::AcbValue>
    pick_zero_solution(long prec = numeric::working_prec_bits()) const;

    BCKind                                  bc_kind()  const { return bc_kind_; }
    const std::vector<BoundarySpec>&        bc_singular() const { return bc_singular_; }
    const std::vector<numeric::AcbValue>&   bc_regular()  const { return bc_regular_; }
    const numeric::RationalMatrix&          de() const { return de_; }
    bool                                     point_at_infinity() const { return point_at_infinity_; }
    const numeric::RationalComplex&         point() const { return point_; }
    const std::vector<AsyExpansion>&        asy_exp() const { return asy_exp_; }

private:
    numeric::RationalMatrix              de_;
    std::vector<BoundarySpec>            bc_singular_;
    std::vector<numeric::AcbValue>       bc_regular_;
    BCKind                                bc_kind_;
    bool                                  point_at_infinity_;
    numeric::RationalComplex             point_;
    std::vector<AsyExpansion>            asy_exp_;
};

// ---------------------------------------------------------------------------
//  amflow
// ---------------------------------------------------------------------------
//
//  Top-level entry point that mirrors AMFlow[de, bcs].  Returns one acb per
//  integral: the leading numerical value at eta = 0.
std::vector<numeric::AcbValue>
amflow(const numeric::RationalMatrix& de,
       const std::vector<BoundarySpec>& bcs,
       long prec = numeric::working_prec_bits());

// ---------------------------------------------------------------------------
//  amflow_multi  (fp-probe column batching)
// ---------------------------------------------------------------------------
//
//  Batched top-level pipeline for the AMFLOW_FP_PROBE_BATCH=1 fixed-point
//  probe: ONE fundamental-solution transport per (de, eps) instead of one
//  full transport per injected unit column.  All columns must share the
//  same de and the same BoundarySpec mu structure (only the boundary
//  VALUES differ — exactly what the probe's unit-vector injections
//  produce).  Column-independent work (contour, CalcInf region plans, the
//  CalcRun re-expansions, CalcZero's NormalizeMat/ToFuchsian chain) is
//  done once; per-column arithmetic is the exact operation sequence of
//  amflow(), so results are bit-identical to per-column amflow() calls.
//
//  Returns result[col][integral].
std::vector<std::vector<numeric::AcbValue>>
amflow_multi(const numeric::RationalMatrix& de,
             const std::vector<std::vector<BoundarySpec>>& bcs_cols,
             long prec = numeric::working_prec_bits());

}  // namespace amflow::ode

#endif  // AMFLOW_ODE_AMFLOW_HPP
