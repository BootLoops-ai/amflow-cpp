// SPDX-License-Identifier: MIT
// numeric::acb_wrap — thin RAII wrappers around FLINT's acb_t / acb_mat_t /
// acb_poly_t.
//
//
// Design goals
// ------------
//   * Zero overhead beyond what FLINT already costs.
//   * Move-only.  Copies are *explicit* via clone().  This avoids accidental
//     O(N^2) deep copies hidden inside container manipulations.
//   * Direct access to the underlying FLINT handles via raw() / data(), so
//     we can drop into FLINT C API at any time.
//   * NO global precision is implied: every arithmetic helper that takes
//     a precision either accepts it explicitly or pulls it from
//     amflow::numeric::working_prec_bits().
//
// We deliberately do NOT overload operator+ / operator* / etc.  FLINT
// routines all want a destination handle plus operands plus precision; hiding
// that behind operators leads to lots of temporaries and makes the precision
// argument disappear.  Inner-loop code in higher domains calls FLINT
// directly via raw().

#ifndef AMFLOW_NUMERIC_ACB_WRAP_HPP
#define AMFLOW_NUMERIC_ACB_WRAP_HPP

#include <cstddef>
#include <iosfwd>
#include <string>
#include <utility>

#include <flint/acb.h>
#include <flint/acb_mat.h>
#include <flint/acb_poly.h>
#include <flint/arb.h>
#include <flint/arf.h>
#include <flint/fmpq.h>
#include <flint/fmpz.h>

#include "amflow/numeric/options.hpp"

namespace amflow::numeric {

// ---------------------------------------------------------------------------
//  AcbValue  -- single complex ball
// ---------------------------------------------------------------------------

class AcbValue {
public:
    AcbValue() noexcept           { acb_init(handle_); }
    ~AcbValue() noexcept          { acb_clear(handle_); }

    AcbValue(const AcbValue&)            = delete;
    AcbValue& operator=(const AcbValue&) = delete;

    AcbValue(AcbValue&& other) noexcept {
        acb_init(handle_);
        acb_swap(handle_, other.handle_);
    }
    AcbValue& operator=(AcbValue&& other) noexcept {
        if (this != &other) acb_swap(handle_, other.handle_);
        return *this;
    }

    // Explicit deep copy.
    AcbValue clone() const {
        AcbValue out;
        acb_set(out.handle_, handle_);
        return out;
    }

    // ---- accessors ----
    acb_ptr        raw()       { return handle_; }
    acb_srcptr     raw() const { return handle_; }
    acb_t&         data()       { return handle_; }
    const acb_t&   data() const { return handle_; }

    // ---- assignment helpers ----
    void set_zero()                                                 { acb_zero(handle_); }
    void set_one()                                                  { acb_one(handle_); }
    void set_si(long re)                                            { acb_set_si(handle_, re); }
    void set_si_si(long re, long im)                                { acb_set_si_si(handle_, re, im); }
    void set_d_d(double re, double im, long prec = working_prec_bits()) {
        arb_set_d(acb_realref(handle_), re);
        arb_set_d(acb_imagref(handle_), im);
        // d->arb is exact if double is representable; round to prec to be safe.
        acb_set_round(handle_, handle_, prec);
    }
    void set_fmpq(const fmpq_t q, long prec = working_prec_bits()) {
        acb_set_fmpq(handle_, q, prec);
    }
    void set_fmpq_fmpq(const fmpq_t re, const fmpq_t im, long prec = working_prec_bits()) {
        arb_set_fmpq(acb_realref(handle_), re, prec);
        arb_set_fmpq(acb_imagref(handle_), im, prec);
    }

    // ---- predicates ----
    bool is_zero()  const { return acb_is_zero(handle_); }
    bool is_one()   const { return acb_is_one(handle_); }
    bool is_real()  const { return acb_is_real(handle_); }
    bool is_finite() const { return acb_is_finite(handle_); }

    // Mathematica's `Chop[x, 10^-chop_pre]` predicate: |x| < 10^-chop_digits
    // using the *midpoint* of the ball (we do not consult the radius —
    // matches MMA semantics, see notes/refactor_design/numeric.md §8c).
    bool is_chop_zero(int chop_digits = chop_pre()) const;

    // Pretty-print as "(re + im*I)".  display_digits controls decimal
    // significand.  When AMFLOW_DEBUG_PRECISE env var is set, prints
    // midpoint digits beyond the ball radius (diagnostic only).
    std::string to_string(long display_digits = 30) const;

private:
    acb_t handle_;
};

std::ostream& operator<<(std::ostream& os, const AcbValue& v);

// ---------------------------------------------------------------------------
//  AcbVector  -- contiguous vector of acb_t (uses _acb_vec_init / _acb_vec_clear)
// ---------------------------------------------------------------------------

class AcbVector {
public:
    AcbVector() noexcept : data_(nullptr), n_(0) {}
    explicit AcbVector(std::size_t n) {
        n_    = n;
        data_ = (n == 0) ? nullptr : _acb_vec_init(static_cast<long>(n));
    }
    ~AcbVector() noexcept { release(); }

    AcbVector(const AcbVector&)            = delete;
    AcbVector& operator=(const AcbVector&) = delete;

    AcbVector(AcbVector&& other) noexcept
        : data_(other.data_), n_(other.n_) {
        other.data_ = nullptr;
        other.n_    = 0;
    }
    AcbVector& operator=(AcbVector&& other) noexcept {
        if (this != &other) {
            release();
            data_       = other.data_;
            n_          = other.n_;
            other.data_ = nullptr;
            other.n_    = 0;
        }
        return *this;
    }

    AcbVector clone() const {
        AcbVector out(n_);
        for (std::size_t i = 0; i < n_; ++i) {
            acb_set(out.data_ + i, data_ + i);
        }
        return out;
    }

    std::size_t size() const noexcept { return n_; }
    bool        empty() const noexcept { return n_ == 0; }

    acb_ptr       raw()       { return data_; }
    acb_srcptr    raw() const { return data_; }
    acb_ptr       at(std::size_t i)       { return data_ + i; }
    acb_srcptr    at(std::size_t i) const { return data_ + i; }

private:
    void release() noexcept {
        if (data_) {
            _acb_vec_clear(data_, static_cast<long>(n_));
            data_ = nullptr;
            n_    = 0;
        }
    }
    acb_ptr     data_;
    std::size_t n_;
};

// ---------------------------------------------------------------------------
//  AcbMatrix  -- wrapper around acb_mat_t
// ---------------------------------------------------------------------------

class AcbMatrix {
public:
    AcbMatrix() noexcept : initialized_(false), rows_(0), cols_(0) {}
    AcbMatrix(long rows, long cols) {
        acb_mat_init(handle_, rows, cols);
        initialized_ = true;
        rows_        = rows;
        cols_        = cols;
    }
    ~AcbMatrix() noexcept { release(); }

    AcbMatrix(const AcbMatrix&)            = delete;
    AcbMatrix& operator=(const AcbMatrix&) = delete;

    AcbMatrix(AcbMatrix&& other) noexcept
        : initialized_(other.initialized_), rows_(other.rows_), cols_(other.cols_) {
        if (initialized_) {
            // acb_mat_swap requires both initialised; instead copy the struct.
            handle_[0]                = other.handle_[0];
            other.initialized_        = false;
            other.rows_ = other.cols_ = 0;
            // other.handle_ left in moved-from state; release() guarded by
            // initialized_ flag so dtor is a no-op.
        }
    }
    AcbMatrix& operator=(AcbMatrix&& other) noexcept {
        if (this != &other) {
            release();
            initialized_ = other.initialized_;
            rows_        = other.rows_;
            cols_        = other.cols_;
            if (initialized_) {
                handle_[0]                = other.handle_[0];
                other.initialized_        = false;
                other.rows_ = other.cols_ = 0;
            }
        }
        return *this;
    }

    AcbMatrix clone() const {
        if (!initialized_) return AcbMatrix();
        AcbMatrix out(rows_, cols_);
        acb_mat_set(out.handle_, handle_);
        return out;
    }

    bool empty() const noexcept { return !initialized_ || rows_ == 0 || cols_ == 0; }
    long rows()  const noexcept { return rows_; }
    long cols()  const noexcept { return cols_; }

    acb_mat_struct*       raw()       { return handle_; }
    const acb_mat_struct* raw() const { return handle_; }

    acb_ptr    entry(long i, long j)       { return acb_mat_entry(handle_, i, j); }
    acb_srcptr entry(long i, long j) const { return acb_mat_entry(handle_, i, j); }

    void set_zero() { if (initialized_) acb_mat_zero(handle_); }
    void set_identity() {
        if (initialized_) {
            acb_mat_zero(handle_);
            long n = rows_ < cols_ ? rows_ : cols_;
            for (long i = 0; i < n; ++i) acb_one(acb_mat_entry(handle_, i, i));
        }
    }

private:
    void release() noexcept {
        if (initialized_) {
            acb_mat_clear(handle_);
            initialized_ = false;
            rows_ = cols_ = 0;
        }
    }
    acb_mat_t handle_;
    bool      initialized_;
    long      rows_;
    long      cols_;
};

// ---------------------------------------------------------------------------
//  AcbPoly  -- wrapper around acb_poly_t  (used for Newton interpolation in
//              api::fit_eps and various power-series ops).
// ---------------------------------------------------------------------------

class AcbPoly {
public:
    AcbPoly() noexcept           { acb_poly_init(handle_); }
    ~AcbPoly() noexcept          { acb_poly_clear(handle_); }

    AcbPoly(const AcbPoly&)            = delete;
    AcbPoly& operator=(const AcbPoly&) = delete;

    AcbPoly(AcbPoly&& other) noexcept {
        acb_poly_init(handle_);
        acb_poly_swap(handle_, other.handle_);
    }
    AcbPoly& operator=(AcbPoly&& other) noexcept {
        if (this != &other) acb_poly_swap(handle_, other.handle_);
        return *this;
    }

    AcbPoly clone() const {
        AcbPoly out;
        acb_poly_set(out.handle_, handle_);
        return out;
    }

    acb_poly_struct*       raw()       { return handle_; }
    const acb_poly_struct* raw() const { return handle_; }

    long length() const  { return acb_poly_length(handle_); }
    long degree() const  { return acb_poly_degree(handle_); }
    bool is_zero() const { return acb_poly_is_zero(handle_); }

    void set_zero() { acb_poly_zero(handle_); }

private:
    acb_poly_t handle_;
};

// ---------------------------------------------------------------------------
//  Free helpers
// ---------------------------------------------------------------------------

// |x| < 10^-digits using the midpoint of x.  Standalone form mirrors the
// member predicate; useful when you have a raw acb_srcptr.
bool acb_is_chop_zero(acb_srcptr x, int chop_digits);

}  // namespace amflow::numeric

#endif  // AMFLOW_NUMERIC_ACB_WRAP_HPP
