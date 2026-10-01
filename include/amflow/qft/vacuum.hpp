// SPDX-License-Identifier: MIT
// qft::vacuum — single-mass vacuum integrals (Layer 14f).
//
//
// Mirrors AMFlow.m lines 819-824, but the (loop, prop) coverage is extended
// beyond the original fixed 5-entry table so that auxiliary-mass flow has an
// analytic eta->infinity boundary for FULLY-MASSLESS multi-loop graphs.
// Available (loop, prop) pairs:
//   * the general single-mass "banana" family  n = L+1  for every L >= 1
//     (closed product of Gammas; subsumes (1,1),(2,3),(3,4),(4,5)), and
//   * the non-banana three-loop five-line vacuum (3,5).
// See src/qft/vacuum.cpp for the closed form and conventions.

#ifndef AMFLOW_QFT_VACUUM_HPP
#define AMFLOW_QFT_VACUUM_HPP

#include <flint/acb.h>

#include "amflow/numeric/acb_wrap.hpp"

namespace amflow::qft {

void vacuum(long loop_count, long prop_count,
            const acb_t eps, acb_t out, long prec);

inline numeric::AcbValue vacuum(long loop_count, long prop_count,
                                 const acb_t eps, long prec) {
    numeric::AcbValue v;
    vacuum(loop_count, prop_count, eps, v.raw(), prec);
    return v;
}

bool vacuum_known(long loop_count, long prop_count) noexcept;

}  // namespace amflow::qft

#endif  // AMFLOW_QFT_VACUUM_HPP
