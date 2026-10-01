// SPDX-License-Identifier: MIT
// ode::qqbar — exact algebraic numbers via FLINT qqbar (Layer 7 helper).
//
//
// Why this exists
// ---------------
//
//   The Layer-7 normalisation pipeline must carry exact algebraic eigen-data
//   far enough that NormalizeMat stays algorithmically aligned with
//   Mathematica before the final algebraic Jordan step.  FLINT exposes
//   qqbar_t (exact algebraic number); this wrapper gives Layer 7 a minimal
//   C++ lifetime/utility surface.
//
// Scope
// -----
//
//   * exact lifetime / copy helpers for qqbar_t
//   * exact conversion to acb for the existing numeric plumbing
//   * exact root/eigenvalue extraction for rational polynomials/matrices
//   * exact floor(Re(lambda)) used by algebraic shearing decisions

#ifndef AMFLOW_ODE_QQBAR_HPP
#define AMFLOW_ODE_QQBAR_HPP

#include <string>
#include <vector>

#include <flint/fmpq_mat.h>
#include <flint/fmpq_poly.h>
#include <flint/qqbar.h>

#include "amflow/numeric/acb_wrap.hpp"
#include "amflow/numeric/options.hpp"

namespace amflow::ode {

class QqbarValue {
public:
    QqbarValue() noexcept                 { qqbar_init(handle_); }
    ~QqbarValue() noexcept                { qqbar_clear(handle_); }

    QqbarValue(const QqbarValue&)            = delete;
    QqbarValue& operator=(const QqbarValue&) = delete;

    QqbarValue(QqbarValue&& other) noexcept {
        qqbar_init(handle_);
        qqbar_swap(handle_, other.handle_);
    }
    QqbarValue& operator=(QqbarValue&& other) noexcept {
        if (this != &other) qqbar_swap(handle_, other.handle_);
        return *this;
    }

    QqbarValue clone() const {
        QqbarValue out;
        qqbar_set(out.handle_, handle_);
        return out;
    }

    qqbar_ptr    raw()       { return handle_; }
    qqbar_srcptr raw() const { return handle_; }

    void set_zero()               { qqbar_zero(handle_); }
    void set_si(long x)           { qqbar_set_si(handle_, x); }
    void set_fmpq(const fmpq_t q) { qqbar_set_fmpq(handle_, q); }

    bool is_zero() const     { return qqbar_is_zero(handle_); }
    bool is_rational() const { return qqbar_is_rational(handle_); }
    bool is_real() const     { return qqbar_is_real(handle_); }

    bool operator==(const QqbarValue& other) const {
        return qqbar_equal(handle_, other.handle_) != 0;
    }
    bool operator!=(const QqbarValue& other) const { return !(*this == other); }

    long              floor_real_si() const;
    numeric::AcbValue to_acb(long prec = numeric::working_prec_bits()) const;
    std::string       to_string(long display_digits = 30) const;

private:
    qqbar_t handle_;
};

struct QqbarEigenvalue {
    QqbarValue value;
    long       multiplicity = 0;
};

std::vector<QqbarValue>
algebraic_roots_from_fmpq_poly(const fmpq_poly_t poly);

std::vector<QqbarEigenvalue>
algebraic_eigenvalues_from_fmpq_mat(const fmpq_mat_t A);

}  // namespace amflow::ode

#endif  // AMFLOW_ODE_QQBAR_HPP
