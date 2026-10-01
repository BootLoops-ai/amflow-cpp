// SPDX-License-Identifier: MIT
// ode::acb_rational — rational functions over Acb (complex arbitrary-prec)
// coefficients.  Layer 7 internal helper.
//
//
// Why this exists
// ---------------
//
//   The DESolver normalisation pipeline produces matrices whose coefficients
//   are *complex arbitrary-precision numbers* (eigenvalues of Jordan blocks,
//   shifts induced by ShearingTransformation, etc.).  These cannot be
//   represented in fmpq.  AcbRationalFunction / AcbRationalMatrix are
//   deliberately minimal: they hold (numerator, denominator) pairs of
//   acb_poly_t with NO automatic reduction (no GCD, since acb has no exact
//   GCD).  Call `value_at_zero` (which strips the leading eta-prefix) before
//   reading off coefficients.
//
//   These types are only needed inside Layer 7.  If you can stay in fmpq
//   land, use numeric::RationalFunction instead.

#ifndef AMFLOW_ODE_ACB_RATIONAL_HPP
#define AMFLOW_ODE_ACB_RATIONAL_HPP

#include <cstddef>
#include <utility>
#include <vector>

#include <flint/acb.h>
#include <flint/acb_mat.h>
#include <flint/acb_poly.h>
#include <flint/fmpq.h>

#include "amflow/numeric/acb_wrap.hpp"
#include "amflow/numeric/options.hpp"
#include "amflow/numeric/rational.hpp"

namespace amflow::ode {

// ---------------------------------------------------------------------------
//  AcbRationalFunction  -- (acb_poly num) / (acb_poly den).  No reduction.
// ---------------------------------------------------------------------------

class AcbRationalFunction {
public:
    AcbRationalFunction();   // default = 0 / 1
    ~AcbRationalFunction() = default;

    AcbRationalFunction(const AcbRationalFunction& other);
    AcbRationalFunction& operator=(const AcbRationalFunction& other);
    AcbRationalFunction(AcbRationalFunction&& other) noexcept = default;
    AcbRationalFunction& operator=(AcbRationalFunction&& other) noexcept = default;

    static AcbRationalFunction
    from_rational(const numeric::RationalFunction& r,
                  long prec = numeric::working_prec_bits());

    static AcbRationalFunction from_acb(acb_srcptr c);
    static AcbRationalFunction from_si(long n);

    // c * eta^k for any (positive or negative) k.
    static AcbRationalFunction monomial(long k);
    static AcbRationalFunction monomial(long k, acb_srcptr c);

    acb_poly_struct*       num()       { return num_; }
    const acb_poly_struct* num() const { return num_; }
    acb_poly_struct*       den()       { return den_; }
    const acb_poly_struct* den() const { return den_; }

    bool is_zero() const;
    bool is_one()  const;

    AcbRationalFunction& operator+=(const AcbRationalFunction& other);
    AcbRationalFunction& operator-=(const AcbRationalFunction& other);
    AcbRationalFunction& operator*=(const AcbRationalFunction& other);
    AcbRationalFunction& operator/=(const AcbRationalFunction& other);
    AcbRationalFunction operator-() const;

    // Multiply by eta^k (k can be negative) in place.
    AcbRationalFunction& multiply_by_eta_power(long k);

    // Multiply by an acb scalar.
    AcbRationalFunction& multiply_by_acb(acb_srcptr c,
                                         long prec = numeric::working_prec_bits());

    AcbRationalFunction derivative(long prec = numeric::working_prec_bits()) const;

    // Evaluate at eta = 0.  Throws std::domain_error if denom vanishes there.
    numeric::AcbValue
    value_at_zero(long prec = numeric::working_prec_bits()) const;

    // Strip a leading common factor of eta from numerator/denominator.
    // Trivial detectable case only.
    void strip_eta_prefix();

    // Round num / den to `prec` bits.
    void round(long prec);

    long prec_used() const { return prec_; }

private:
    acb_poly_t num_;
    acb_poly_t den_;
    long       prec_;
};

// ---------------------------------------------------------------------------
//  AcbRationalMatrix
// ---------------------------------------------------------------------------

class AcbRationalMatrix {
public:
    AcbRationalMatrix() = default;
    AcbRationalMatrix(std::size_t rows, std::size_t cols);

    static AcbRationalMatrix identity(std::size_t n);
    static AcbRationalMatrix
    from_rational(const numeric::RationalMatrix& m,
                  long prec = numeric::working_prec_bits());

    std::size_t rows() const { return rows_; }
    std::size_t cols() const { return cols_; }
    bool        empty() const { return rows_ == 0 || cols_ == 0; }

    AcbRationalFunction&       operator()(std::size_t i, std::size_t j)       { return data_[i * cols_ + j]; }
    const AcbRationalFunction& operator()(std::size_t i, std::size_t j) const { return data_[i * cols_ + j]; }

    AcbRationalMatrix submatrix(const std::vector<std::size_t>& rows,
                                const std::vector<std::size_t>& cols) const;

    AcbRationalMatrix matmul(const AcbRationalMatrix& other,
                             long prec = numeric::working_prec_bits()) const;

    AcbRationalMatrix operator-(const AcbRationalMatrix& other) const;
    AcbRationalMatrix operator+(const AcbRationalMatrix& other) const;

private:
    std::size_t rows_ = 0;
    std::size_t cols_ = 0;
    std::vector<AcbRationalFunction> data_;
};

// ---------------------------------------------------------------------------
//  Numeric matrix helpers
// ---------------------------------------------------------------------------

// Evaluate `m` at eta = 0 (after stripping eta prefixes) into a complex
// matrix `out`.  Throws if any entry has a pole at eta=0.
void acb_rational_matrix_value_at_zero(acb_mat_t out,
                                       const AcbRationalMatrix& m,
                                       long prec = numeric::working_prec_bits());

// Evaluate `eta * m` at eta = 0  (the matrix L0 in NormalizeEigenQ /
// AsymptoticBehavior).
void acb_rational_matrix_eta_at_zero(acb_mat_t out,
                                     const AcbRationalMatrix& m,
                                     long prec = numeric::working_prec_bits());

}  // namespace amflow::ode

#endif  // AMFLOW_ODE_ACB_RATIONAL_HPP
