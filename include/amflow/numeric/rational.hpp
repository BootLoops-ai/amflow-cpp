// SPDX-License-Identifier: MIT
// numeric::Rational — single-variable rational types over ℚ in η, plus
// rational complex (poles) and dense rational matrices.
//
//
// Coefficients are exact rationals (fmpz / fmpq).  This matches the
// Mathematica side, where the ODE matrix lives in ℚ(η) and only the
// evaluation step touches arbitrary-precision floats.
//
// Types:
//   * FmpqPoly         — RAII over fmpq_poly_t, implicit "eta" variable.
//                          Move-only with explicit clone().
//   * RationalComplex  — (fmpq, fmpq) pair for poles / contour points.
//                          Has both copy and move (lightweight).
//   * RationalFunction — num/den always reduced (gcd-cancelled, denom monic).
//                          Has both copy and move.
//   * RationalMatrix   — dense vector<RationalFunction>, row-major.
//
// `get_poles(matrix)` is the bridge from rational matrix to rational complex
// poles: LCM of denominators → factor → complex roots → rationalise.

#ifndef AMFLOW_NUMERIC_RATIONAL_HPP
#define AMFLOW_NUMERIC_RATIONAL_HPP

#include <cstddef>
#include <iosfwd>
#include <string>
#include <utility>
#include <vector>

#include <flint/fmpq.h>
#include <flint/fmpq_poly.h>
#include <flint/fmpz.h>
#include <flint/fmpz_poly.h>

#include "amflow/numeric/acb_wrap.hpp"
#include "amflow/numeric/options.hpp"

namespace amflow::numeric {

// ---------------------------------------------------------------------------
//  FmpqPoly  -- RAII wrapper around fmpq_poly_t.
// ---------------------------------------------------------------------------
//
// Move-only; explicit clone().  Implicit variable name is "eta" (only matters
// for to_string).  FLINT's fmpq_poly_t is "integer coefficients + one shared
// denominator" internally even though we treat it as "fmpq per coefficient"
// semantically.

class FmpqPoly {
public:
    FmpqPoly() noexcept                 { fmpq_poly_init(handle_); }
    ~FmpqPoly() noexcept                { fmpq_poly_clear(handle_); }

    FmpqPoly(const FmpqPoly&)            = delete;
    FmpqPoly& operator=(const FmpqPoly&) = delete;

    FmpqPoly(FmpqPoly&& other) noexcept {
        fmpq_poly_init(handle_);
        fmpq_poly_swap(handle_, other.handle_);
    }
    FmpqPoly& operator=(FmpqPoly&& other) noexcept {
        if (this != &other) fmpq_poly_swap(handle_, other.handle_);
        return *this;
    }

    FmpqPoly clone() const {
        FmpqPoly out;
        fmpq_poly_set(out.handle_, handle_);
        return out;
    }

    // ---- accessors ----
    fmpq_poly_struct*       raw()       { return handle_; }
    const fmpq_poly_struct* raw() const { return handle_; }

    long length() const  { return fmpq_poly_length(handle_); }
    long degree() const  { return fmpq_poly_degree(handle_); }   // -1 for zero
    bool is_zero() const { return fmpq_poly_is_zero(handle_); }
    bool is_one()  const { return fmpq_poly_is_one(handle_); }

    // Smallest k such that coefficient of η^k is nonzero.  Returns -1 for
    // the zero polynomial.  Equivalent to MMA's Exponent[poly, eta, Min].
    long valuation() const;

    // ---- assignment helpers ----
    void set_zero()             { fmpq_poly_zero(handle_); }
    void set_one()              { fmpq_poly_one(handle_); }
    void set_si(long c)         { fmpq_poly_set_si(handle_, c); }
    void set_fmpq(const fmpq_t q)               { fmpq_poly_set_fmpq(handle_, q); }
    void set_coeff_si(long n, long c)           { fmpq_poly_set_coeff_si(handle_, n, c); }
    void set_coeff_fmpq(long n, const fmpq_t c) { fmpq_poly_set_coeff_fmpq(handle_, n, c); }

    // Copy a coefficient out into a freshly-initialised fmpq_t (caller clears).
    void coeff(long n, fmpq_t out) const { fmpq_poly_get_coeff_fmpq(out, handle_, n); }

    // Pretty-print as "1 + 2*x - 3*x^2".
    std::string to_string(const char* var = "eta") const;

private:
    fmpq_poly_t handle_;
};

std::ostream& operator<<(std::ostream& os, const FmpqPoly& p);

// Free helpers operating on raw fmpq_poly handles.
void fmpq_poly_make_denominator_monic_pair(fmpq_poly_t num, fmpq_poly_t den);
long fmpq_poly_valuation(const fmpq_poly_t p);

// ---------------------------------------------------------------------------
//  RationalComplex — (re, im) of fmpq.
// ---------------------------------------------------------------------------
//
// Used as the return type of get_poles, the contour-point list, and any
// place a "rationalised complex number" is needed.  Lightweight — has both
// copy and move.

class RationalComplex {
public:
    RationalComplex();
    ~RationalComplex();
    RationalComplex(const RationalComplex& other);
    RationalComplex(RationalComplex&& other) noexcept;
    RationalComplex& operator=(const RationalComplex& other);
    RationalComplex& operator=(RationalComplex&& other) noexcept;

    static RationalComplex from_fractions(long re_num, long re_den,
                                          long im_num, long im_den);

    // The members are mutable on purpose — get_poles writes into them.
    fmpq*       re()       { return re_; }
    const fmpq* re() const { return re_; }
    fmpq*       im()       { return im_; }
    const fmpq* im() const { return im_; }

    bool is_real() const;
    bool is_zero() const;

    AcbValue to_acb(long prec = working_prec_bits()) const;
    std::string to_string() const;

    bool operator==(const RationalComplex& other) const;
    bool operator!=(const RationalComplex& other) const { return !(*this == other); }

private:
    fmpq_t re_;
    fmpq_t im_;
};

std::ostream& operator<<(std::ostream& os, const RationalComplex& z);

// ---------------------------------------------------------------------------
//  RationalFunction — num(η) / den(η), always reduced.
// ---------------------------------------------------------------------------
//
// THE working type for the differential-equation matrix.  All arithmetic
// produces a reduced result (i.e. Together is implicit).  Has both copy and
// move (FmpqPoly clone is cheap).

class RationalFunction {
public:
    RationalFunction();                                  // = 0/1

    // num / den, reduced.  Throws std::invalid_argument if den is zero.
    RationalFunction(FmpqPoly num, FmpqPoly den);

    static RationalFunction from_polynomial(FmpqPoly num);
    static RationalFunction from_fmpq(const fmpq_t q);
    static RationalFunction from_si(long n);
    static RationalFunction from_si_si(long p, long q);
    static RationalFunction monomial(long k);             // η^k
    static RationalFunction monomial(long k, const fmpq_t c);

    RationalFunction(const RationalFunction& other);
    RationalFunction& operator=(const RationalFunction& other);
    RationalFunction(RationalFunction&& other) noexcept = default;
    RationalFunction& operator=(RationalFunction&& other) noexcept = default;

    const FmpqPoly& numerator()   const { return num_; }
    const FmpqPoly& denominator() const { return den_; }

    bool is_zero()       const { return num_.is_zero(); }
    bool is_one()        const { return num_.is_one() && den_.is_one(); }
    bool is_polynomial() const;
    bool is_constant()   const;

    long degree_num() const { return num_.degree(); }
    long degree_den() const { return den_.degree(); }
    long valuation_num() const { return num_.valuation(); }
    long valuation_den() const { return den_.valuation(); }
    long net_valuation() const;

    // ---- arithmetic ----
    RationalFunction operator-() const;
    RationalFunction& operator+=(const RationalFunction& other);
    RationalFunction& operator-=(const RationalFunction& other);
    RationalFunction& operator*=(const RationalFunction& other);
    RationalFunction& operator/=(const RationalFunction& other);

    RationalFunction& multiply_by_polynomial(const FmpqPoly& p);
    RationalFunction& multiply_by_eta_power(long k);
    RationalFunction& multiply_by_fmpq(const fmpq_t q);

    RationalFunction reciprocal() const;
    RationalFunction derivative() const;

    // η -> 1/η, then reduce.  Used by calc_inf:
    //   deinf = -η^-2 * (de /. η -> 1/η)
    RationalFunction substitute_inverse_eta() const;

    AcbValue evaluate(acb_srcptr x, long prec = working_prec_bits()) const;
    AcbValue evaluate(const AcbValue& x, long prec = working_prec_bits()) const {
        return evaluate(x.raw(), prec);
    }
    AcbValue evaluate(const RationalComplex& x, long prec = working_prec_bits()) const;

    std::string to_string(const char* var = "eta") const;

    bool operator==(const RationalFunction& other) const;
    bool operator!=(const RationalFunction& other) const { return !(*this == other); }

private:
    void reduce();
    FmpqPoly num_;
    FmpqPoly den_;
};

RationalFunction operator+(RationalFunction a, const RationalFunction& b);
RationalFunction operator-(RationalFunction a, const RationalFunction& b);
RationalFunction operator*(RationalFunction a, const RationalFunction& b);
RationalFunction operator/(RationalFunction a, const RationalFunction& b);

std::ostream& operator<<(std::ostream& os, const RationalFunction& r);

// ---------------------------------------------------------------------------
//  RationalMatrix — dense, row-major.
// ---------------------------------------------------------------------------

class RationalMatrix {
public:
    RationalMatrix() = default;
    RationalMatrix(std::size_t rows, std::size_t cols);

    static RationalMatrix identity(std::size_t n);
    static RationalMatrix zeros(std::size_t rows, std::size_t cols);

    std::size_t rows() const noexcept { return rows_; }
    std::size_t cols() const noexcept { return cols_; }
    bool        empty() const noexcept { return rows_ == 0 || cols_ == 0; }

    RationalFunction&       operator()(std::size_t i, std::size_t j)       { return data_[i * cols_ + j]; }
    const RationalFunction& operator()(std::size_t i, std::size_t j) const { return data_[i * cols_ + j]; }

    RationalFunction&       at(std::size_t i, std::size_t j);
    const RationalFunction& at(std::size_t i, std::size_t j) const;

    double density() const;
    long   max_numerator_degree() const;
    long   max_denominator_degree() const;

    RationalMatrix submatrix(const std::vector<std::size_t>& rows,
                             const std::vector<std::size_t>& cols) const;
    RationalMatrix substituted_inverse_eta() const;

    FmpqPoly denominator_lcm(const std::vector<std::size_t>& row_indices,
                             const std::vector<std::size_t>& col_indices) const;
    FmpqPoly denominator_lcm() const;

    std::string to_string(const char* var = "eta") const;
    bool operator==(const RationalMatrix& other) const;
    bool operator!=(const RationalMatrix& other) const { return !(*this == other); }

private:
    std::size_t                   rows_ = 0;
    std::size_t                   cols_ = 0;
    std::vector<RationalFunction> data_;
};

std::ostream& operator<<(std::ostream& os, const RationalMatrix& m);

// ---------------------------------------------------------------------------
//  get_poles
// ---------------------------------------------------------------------------
//
//  Steps:
//    1. LCM of all entry denominators across the matrix.
//    2. Convert to primitive integer polynomial.
//    3. Square-free factorise (fmpz_poly_factor); discard repeats.
//    4. Find complex roots of each remaining factor at the current
//       working precision.
//    5. Rationalise each root via best-rational-approximation at tolerance
//       10^-RationalizePre.
//    6. Deduplicate.
//
//  Returns a deduplicated list of rational complex poles.  Order not
//  guaranteed.

std::vector<RationalComplex> get_poles(const RationalMatrix& matrix);

}  // namespace amflow::numeric

#endif  // AMFLOW_NUMERIC_RATIONAL_HPP
