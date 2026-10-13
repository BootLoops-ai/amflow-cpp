// SPDX-License-Identifier: MIT
// ibp::libp_deriv — symbolic derivative of J integrals via chain rule.
//
//
// Mirrors AMFlow.m / Kira interface lines 282-405.
//
// Provenance
// --------------------------------
//   The block of AMFlow's `ibp_interface/Kira/interface.m` mirrored here is
//   from LiteIBP.m by Tiziano Peraro (FiniteFlow MathTools,
//   https://github.com/peraro/finiteflow-mathtools, packages/LiteIBP.m;
//   T. Peraro, arXiv:1905.08019).  AMFlow marks it "compute derivatives
//   (code from LiteIBP.m provided by T.Peraro)" (interface.m line 280 at
//   AMFlow commit efda1db, line 400 at tag 2.0).  This file and
//   `src/ibp/libp_deriv.cpp` are upstream AMFlow.cpp's C++ port of that
//   block.  The momentum-derivative chain rule of scope item 2 below was
//   added by the BootLoops fork of AMFlow.cpp; it takes the place of
//   LiteIBP's `LIBPDerivivative` (LiteIBP.m line 565 at commit f77e2b5;
//   interface.m line 383 at efda1db, line 503 at tag 2.0).
//   AMFlow and AMFlow.cpp are MIT.  The finiteflow-mathtools repository
//   carries no license file (read at commit f77e2b5).  The Anthropic
//   copyright and the MIT grant on the BootLoops modifications do not cover
//   the portion of this file and of `libp_deriv.cpp` that derives from
//   LiteIBP.  See NOTICE and PATCHES.md.
//
//   For J(a) = Integral[ Π_k 1 / D_k^{a_k} ]
//
//      ∂J(a)/∂s = Σ_k (∂D_k/∂s) · (-a_k · J(a + e_k))
//
//   where (∂D_k/∂s) is decomposed over the propagator basis as
//
//      ∂D_k/∂s = Σ_j coef_kj · D_j  +  const_k
//
//   so that  D_j · J(a + e_k) = J(a + e_k - e_j).
//
// Scope
// --------------------------------
//   Variables `s_name` for which differentiation is supported:
//
//     1. variables that appear directly as a polynomial-coefficient
//        in `fc.propagators_after_conservation[k]` — i.e. free
//        symbols (typical examples: `eta` after η-insertion, mass
//        parameters like `msq`, `m1sq`, `m2sq`).
//
//     2. Replacement-defined kinematic invariants (e.g. `s` defined
//        via `(p1+p2)^2 -> s` in `fc.replacement`).  Propagators are
//        stored pre-replacement, so the direct ∂D_k/∂s vanishes; the
//        derivative is instead computed via the momentum-derivative
//        chain rule of T. Peraro's LiteIBP (`LIBPDerivivative`, as
//        carried in AMFlow's `Kira/interface.m`), which this routine
//        follows in what it computes; the coefficients are obtained in
//        closed form: with G_{ij} = p_i·p_j = R_{ij}(s,t,…)
//        the external Gram matrix and J_{ij} = ∂R_{ij}/∂s, the
//        per-propagator effective derivative is
//
//           ∂D_k/∂s|_eff  =  ∂D_k/∂s|_direct
//                          + Σ_{a,i} (½ G⁻¹ J)_{ai} · (p_a · ∂D_k/∂p_i)
//
//        decomposed onto the propagator basis as in (1).  The chain
//        term is skipped (`J ≡ 0`) for any variable not appearing in
//        `fc.reduced_replacement`, so case-1 callers (the production
//        η-DE path, `src/pipeline/amfsystem.cpp:1012`) take the
//        identical pre-patch code path.  Requires the external Gram
//        matrix to be non-singular (generic kinematics) and every
//        leg·leg scalar product to have a replacement rule.
//
//   Tests exercising both scopes: `test_ibp_libp_deriv.cpp`
//   `*TwoMassBubble*` cases (mass + kinematic invariants) and
//   `*Box1_Wrt_s*` cases (multi-leg kinematic chain).

#ifndef AMFLOW_IBP_LIBP_DERIV_HPP
#define AMFLOW_IBP_LIBP_DERIV_HPP

#include <map>
#include <memory>
#include <string>
#include <vector>

#include "amflow/algebra/mpoly.hpp"
#include "amflow/qft/family_config.hpp"
#include "amflow/qft/jintegral.hpp"
#include "amflow/qft/region.hpp"   // DListContext

namespace amflow::ibp {

struct DerivTerm {
    algebra::Mfrac    coef;
    qft::JIntegral    integ;
};

struct LibpDenomsDerivResult {
    std::vector<algebra::Mpoly>               completede;
    std::vector<std::vector<algebra::Mfrac>>  coef;
    std::vector<algebra::Mfrac>               constant;
};

LibpDenomsDerivResult
libp_denoms_deriv(const qft::FamilyConfig& fc, const std::string& s_name);

std::vector<DerivTerm>
libp_deriv(const qft::FamilyConfig& fc,
           const qft::JIntegral& j,
           const std::string& s_name);

// Narrow-context overload: same as `libp_denoms_deriv(fc, s_name)`
// but applies `numeric_values` substitution to every intermediate
// Mfrac AND reprojects the result onto `target_ctx`.  Use this when
// you have a narrow target context (e.g. `{eta, d}` from
// `make_reduction_context`) and want to avoid FLINT's
// multivariate-GCD overhead on the wide `fc.ctx` polynomial ring.
//
// Preconditions: `target_ctx` must contain `s_name` as one of its
// variables, and `numeric_values` must supply rationals for every
// `fc.ctx` variable that is NOT in `target_ctx`.
//
// Returned `coef[k][j]` and `constant[k]` Mfracs live on `target_ctx`.
LibpDenomsDerivResult
libp_denoms_deriv(const qft::FamilyConfig& fc, const std::string& s_name,
                  const std::map<std::string, std::string>& numeric_values,
                  const std::shared_ptr<algebra::MpolyContext>& target_ctx);

// Narrow-context overload of libp_deriv mirroring the above:
// returns DerivTerms whose `.coef` lives on `target_ctx`.
std::vector<DerivTerm>
libp_deriv(const qft::FamilyConfig& fc,
           const qft::JIntegral& j,
           const std::string& s_name,
           const std::map<std::string, std::string>& numeric_values,
           const std::shared_ptr<algebra::MpolyContext>& target_ctx);

std::vector<DerivTerm>
compute_derivative(const qft::FamilyConfig& fc,
                   const std::vector<DerivTerm>& expr,
                   const std::string& x_name);

std::vector<DerivTerm>
simplify_terms(std::vector<DerivTerm> terms);

}  // namespace amflow::ibp

#endif  // AMFLOW_IBP_LIBP_DERIV_HPP
