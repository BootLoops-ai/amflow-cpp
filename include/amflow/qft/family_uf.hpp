// SPDX-License-Identifier: MIT
// qft::family_uf — Symanzik U/F polynomials via EvaluateABC / EvaluateUF.
//
//
// Mirrors AMFlow.m lines 446-473.

#ifndef AMFLOW_QFT_FAMILY_UF_HPP
#define AMFLOW_QFT_FAMILY_UF_HPP

#include <memory>
#include <vector>

#include "amflow/algebra/mpoly.hpp"
#include "amflow/algebra/mpoly_matrix.hpp"
#include "amflow/qft/family_config.hpp"

namespace amflow::qft {

// ---------------------------------------------------------------------------
//  ABC decomposition
// ---------------------------------------------------------------------------

struct ABCResult {
    std::shared_ptr<algebra::MpolyContext> uf_ctx;
    long                                    first_x_var;
    algebra::MpolyMatrix                    Ax;
    std::vector<algebra::Mpoly>             Bx;
    algebra::Mpoly                          Cx;
};

ABCResult evaluate_abc(const FamilyConfig& fc,
                       const std::vector<algebra::Mpoly>& denominators);

// ---------------------------------------------------------------------------
//  Symanzik U / F polynomials
// ---------------------------------------------------------------------------

struct UFResult {
    std::shared_ptr<algebra::MpolyContext> uf_ctx;
    long                                    first_x_var;

    algebra::Mfrac u;
    algebra::Mfrac f;
    algebra::Mfrac f0;

    bool degenerate;
};

UFResult evaluate_uf(const FamilyConfig& fc,
                     const std::vector<algebra::Mpoly>& denominators);

// ---------------------------------------------------------------------------
//  apply_uf_replacement — fc.reduced_replacement on extended ctx
// ---------------------------------------------------------------------------

algebra::Mfrac
apply_uf_replacement(const FamilyConfig& fc,
                     const std::shared_ptr<algebra::MpolyContext>& uf_ctx,
                     const algebra::Mfrac& p);

}  // namespace amflow::qft

#endif  // AMFLOW_QFT_FAMILY_UF_HPP
