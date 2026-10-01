// SPDX-License-Identifier: MIT
// ode::acb_rational — implementation.
//

#include "amflow/ode/acb_rational.hpp"

#include <algorithm>
#include <cstdlib>
#include <stdexcept>
#include <utility>

#include <flint/acb.h>
#include <flint/acb_mat.h>
#include <flint/acb_poly.h>
#include <flint/fmpq.h>

#include "amflow/numeric/options.hpp"

namespace amflow::ode {

using numeric::AcbValue;
using numeric::RationalFunction;
using numeric::RationalMatrix;
using numeric::working_prec_bits;

// ===========================================================================
//  AcbRationalFunction
// ===========================================================================

AcbRationalFunction::AcbRationalFunction()
    : prec_(working_prec_bits()) {
    acb_poly_init(num_);
    acb_poly_init(den_);
    acb_poly_one(den_);
}

AcbRationalFunction::AcbRationalFunction(const AcbRationalFunction& other)
    : prec_(other.prec_) {
    acb_poly_init(num_);
    acb_poly_init(den_);
    acb_poly_set(num_, other.num_);
    acb_poly_set(den_, other.den_);
}

AcbRationalFunction&
AcbRationalFunction::operator=(const AcbRationalFunction& other) {
    if (this != &other) {
        acb_poly_set(num_, other.num_);
        acb_poly_set(den_, other.den_);
        prec_ = other.prec_;
    }
    return *this;
}

AcbRationalFunction
AcbRationalFunction::from_rational(const RationalFunction& r, long prec) {
    AcbRationalFunction out;
    out.prec_ = prec;
    acb_poly_set_fmpq_poly(out.num_, r.numerator().raw(), prec);
    acb_poly_set_fmpq_poly(out.den_, r.denominator().raw(), prec);
    return out;
}

AcbRationalFunction AcbRationalFunction::from_acb(acb_srcptr c) {
    AcbRationalFunction out;
    acb_poly_set_acb(out.num_, c);
    return out;
}

AcbRationalFunction AcbRationalFunction::from_si(long n) {
    AcbRationalFunction out;
    acb_poly_set_si(out.num_, n);
    return out;
}

AcbRationalFunction AcbRationalFunction::monomial(long k) {
    AcbValue one; one.set_one();
    return monomial(k, one.raw());
}

AcbRationalFunction AcbRationalFunction::monomial(long k, acb_srcptr c) {
    AcbRationalFunction out;
    if (k >= 0) {
        acb_poly_set_coeff_acb(out.num_, k, c);
    } else {
        acb_poly_set_acb(out.num_, c);
        // Default ctor leaves den_ = 1.  Without this zero, set_coeff_acb
        // below would leave a residual coeff[0] = 1 -> c / (1 + eta^|k|)
        // instead of c * eta^k.
        acb_poly_zero(out.den_);
        AcbValue one; one.set_one();
        acb_poly_set_coeff_acb(out.den_, -k, one.raw());
    }
    return out;
}

bool AcbRationalFunction::is_zero() const { return acb_poly_is_zero(num_); }
bool AcbRationalFunction::is_one()  const {
    return acb_poly_equal(num_, den_) && !acb_poly_is_zero(den_);
}

namespace {

// (a/b) +/- (c/d) = (a*d +/- b*c) / (b*d).  No reduction.
inline void rat_add_sub(acb_poly_t out_num,
                        acb_poly_t out_den,
                        const acb_poly_struct* a_num,
                        const acb_poly_struct* a_den,
                        const acb_poly_struct* b_num,
                        const acb_poly_struct* b_den,
                        bool subtract,
                        long prec) {
    acb_poly_t ad, bc, bd;
    acb_poly_init(ad);
    acb_poly_init(bc);
    acb_poly_init(bd);
    acb_poly_mul(ad, a_num, b_den, prec);
    acb_poly_mul(bc, a_den, b_num, prec);
    acb_poly_mul(bd, a_den, b_den, prec);
    if (subtract) {
        acb_poly_sub(out_num, ad, bc, prec);
    } else {
        acb_poly_add(out_num, ad, bc, prec);
    }
    acb_poly_set(out_den, bd);
    acb_poly_clear(ad);
    acb_poly_clear(bc);
    acb_poly_clear(bd);
}

}  // namespace

AcbRationalFunction&
AcbRationalFunction::operator+=(const AcbRationalFunction& other) {
    rat_add_sub(num_, den_, num_, den_, other.num_, other.den_, false, prec_);
    return *this;
}
AcbRationalFunction&
AcbRationalFunction::operator-=(const AcbRationalFunction& other) {
    rat_add_sub(num_, den_, num_, den_, other.num_, other.den_, true, prec_);
    return *this;
}
AcbRationalFunction&
AcbRationalFunction::operator*=(const AcbRationalFunction& other) {
    acb_poly_mul(num_, num_, other.num_, prec_);
    acb_poly_mul(den_, den_, other.den_, prec_);
    return *this;
}
AcbRationalFunction&
AcbRationalFunction::operator/=(const AcbRationalFunction& other) {
    if (acb_poly_is_zero(other.num_)) {
        throw std::domain_error("AcbRationalFunction divide by zero");
    }
    acb_poly_mul(num_, num_, other.den_, prec_);
    acb_poly_mul(den_, den_, other.num_, prec_);
    return *this;
}

AcbRationalFunction AcbRationalFunction::operator-() const {
    AcbRationalFunction out = *this;
    acb_poly_neg(out.num_, out.num_);
    return out;
}

AcbRationalFunction& AcbRationalFunction::multiply_by_eta_power(long k) {
    if (k == 0) return *this;
    acb_poly_t shift;
    acb_poly_init(shift);
    AcbValue one; one.set_one();
    acb_poly_set_coeff_acb(shift, std::abs(k), one.raw());
    if (k > 0) acb_poly_mul(num_, num_, shift, prec_);
    else        acb_poly_mul(den_, den_, shift, prec_);
    acb_poly_clear(shift);
    return *this;
}

AcbRationalFunction&
AcbRationalFunction::multiply_by_acb(acb_srcptr c, long prec) {
    acb_poly_scalar_mul(num_, num_, c, prec);
    return *this;
}

AcbRationalFunction AcbRationalFunction::derivative(long prec) const {
    // d/deta (n/d) = (n' d - n d') / d^2
    acb_poly_t nprime, dprime, lhs, rhs, top, bottom;
    acb_poly_init(nprime); acb_poly_init(dprime);
    acb_poly_init(lhs);    acb_poly_init(rhs);
    acb_poly_init(top);    acb_poly_init(bottom);

    acb_poly_derivative(nprime, num_, prec);
    acb_poly_derivative(dprime, den_, prec);
    acb_poly_mul(lhs, nprime, den_, prec);
    acb_poly_mul(rhs, num_,   dprime, prec);
    acb_poly_sub(top, lhs, rhs, prec);
    acb_poly_mul(bottom, den_, den_, prec);

    AcbRationalFunction out;
    out.prec_ = prec;
    acb_poly_set(out.num_, top);
    acb_poly_set(out.den_, bottom);

    acb_poly_clear(nprime); acb_poly_clear(dprime);
    acb_poly_clear(lhs); acb_poly_clear(rhs);
    acb_poly_clear(top); acb_poly_clear(bottom);
    return out;
}

void AcbRationalFunction::strip_eta_prefix() {
    long Ln = acb_poly_length(num_);
    long Ld = acb_poly_length(den_);
    long k_num = 0;
    while (k_num < Ln) {
        AcbValue c;
        acb_poly_get_coeff_acb(c.raw(), num_, k_num);
        if (!c.is_zero()) break;
        ++k_num;
    }
    long k_den = 0;
    while (k_den < Ld) {
        AcbValue c;
        acb_poly_get_coeff_acb(c.raw(), den_, k_den);
        if (!c.is_zero()) break;
        ++k_den;
    }
    long k = std::min(k_num, k_den);
    if (k == 0) return;

    acb_poly_shift_right(num_, num_, k);
    acb_poly_shift_right(den_, den_, k);
}

AcbValue AcbRationalFunction::value_at_zero(long prec) const {
    AcbRationalFunction tmp = *this;
    tmp.strip_eta_prefix();

    AcbValue d0;
    acb_poly_get_coeff_acb(d0.raw(), tmp.den(), 0);
    if (d0.is_zero()) {
        throw std::domain_error("AcbRationalFunction::value_at_zero: denominator vanishes");
    }
    AcbValue n0;
    acb_poly_get_coeff_acb(n0.raw(), tmp.num(), 0);
    AcbValue out;
    acb_div(out.raw(), n0.raw(), d0.raw(), prec);
    return out;
}

void AcbRationalFunction::round(long prec) {
    prec_ = prec;
    acb_poly_set_round(num_, num_, prec);
    acb_poly_set_round(den_, den_, prec);
}

// ===========================================================================
//  AcbRationalMatrix
// ===========================================================================

AcbRationalMatrix::AcbRationalMatrix(std::size_t rows, std::size_t cols)
    : rows_(rows), cols_(cols), data_(rows * cols) {}

AcbRationalMatrix AcbRationalMatrix::identity(std::size_t n) {
    AcbRationalMatrix m(n, n);
    for (std::size_t i = 0; i < n; ++i) {
        m(i, i) = AcbRationalFunction::from_si(1);
    }
    return m;
}

AcbRationalMatrix
AcbRationalMatrix::from_rational(const RationalMatrix& m, long prec) {
    AcbRationalMatrix out(m.rows(), m.cols());
    for (std::size_t i = 0; i < m.rows(); ++i) {
        for (std::size_t j = 0; j < m.cols(); ++j) {
            out(i, j) = AcbRationalFunction::from_rational(m(i, j), prec);
        }
    }
    return out;
}

AcbRationalMatrix
AcbRationalMatrix::submatrix(const std::vector<std::size_t>& rows,
                             const std::vector<std::size_t>& cols) const {
    AcbRationalMatrix out(rows.size(), cols.size());
    for (std::size_t i = 0; i < rows.size(); ++i) {
        for (std::size_t j = 0; j < cols.size(); ++j) {
            out(i, j) = (*this)(rows[i], cols[j]);
        }
    }
    return out;
}

AcbRationalMatrix
AcbRationalMatrix::matmul(const AcbRationalMatrix& other, long prec) const {
    if (cols_ != other.rows_) {
        throw std::invalid_argument("AcbRationalMatrix::matmul: dimension mismatch");
    }
    AcbRationalMatrix out(rows_, other.cols_);
    for (std::size_t i = 0; i < rows_; ++i) {
        for (std::size_t j = 0; j < other.cols_; ++j) {
            AcbRationalFunction acc;
            for (std::size_t k = 0; k < cols_; ++k) {
                AcbRationalFunction prod = (*this)(i, k);
                prod *= other(k, j);
                acc  += prod;
            }
            out(i, j) = std::move(acc);
        }
    }
    return out;
}

AcbRationalMatrix
AcbRationalMatrix::operator-(const AcbRationalMatrix& other) const {
    if (rows_ != other.rows_ || cols_ != other.cols_) {
        throw std::invalid_argument("AcbRationalMatrix::operator-: shape mismatch");
    }
    AcbRationalMatrix out(rows_, cols_);
    for (std::size_t i = 0; i < rows_; ++i) {
        for (std::size_t j = 0; j < cols_; ++j) {
            out(i, j) = (*this)(i, j);
            out(i, j) -= other(i, j);
        }
    }
    return out;
}

AcbRationalMatrix
AcbRationalMatrix::operator+(const AcbRationalMatrix& other) const {
    if (rows_ != other.rows_ || cols_ != other.cols_) {
        throw std::invalid_argument("AcbRationalMatrix::operator+: shape mismatch");
    }
    AcbRationalMatrix out(rows_, cols_);
    for (std::size_t i = 0; i < rows_; ++i) {
        for (std::size_t j = 0; j < cols_; ++j) {
            out(i, j) = (*this)(i, j);
            out(i, j) += other(i, j);
        }
    }
    return out;
}

// ===========================================================================
//  Numeric matrix helpers
// ===========================================================================

void acb_rational_matrix_value_at_zero(acb_mat_t out,
                                       const AcbRationalMatrix& m,
                                       long prec) {
    if (acb_mat_nrows(out) != static_cast<long>(m.rows())
        || acb_mat_ncols(out) != static_cast<long>(m.cols())) {
        throw std::invalid_argument("acb_rational_matrix_value_at_zero: out shape mismatch");
    }
    for (std::size_t i = 0; i < m.rows(); ++i) {
        for (std::size_t j = 0; j < m.cols(); ++j) {
            AcbValue v;
            try {
                v = m(i, j).value_at_zero(prec);
            } catch (const std::domain_error&) {
                throw std::domain_error(
                    "acb_rational_matrix_value_at_zero: entry has a pole at eta=0");
            }
            acb_set(acb_mat_entry(out, i, j), v.raw());
        }
    }
}

void acb_rational_matrix_eta_at_zero(acb_mat_t out,
                                     const AcbRationalMatrix& m,
                                     long prec) {
    if (acb_mat_nrows(out) != static_cast<long>(m.rows())
        || acb_mat_ncols(out) != static_cast<long>(m.cols())) {
        throw std::invalid_argument("acb_rational_matrix_eta_at_zero: out shape mismatch");
    }
    AcbRationalFunction eta_fn = AcbRationalFunction::monomial(1);
    for (std::size_t i = 0; i < m.rows(); ++i) {
        for (std::size_t j = 0; j < m.cols(); ++j) {
            AcbRationalFunction tmp = m(i, j);
            tmp *= eta_fn;
            AcbValue v = tmp.value_at_zero(prec);
            acb_set(acb_mat_entry(out, i, j), v.raw());
        }
    }
}

}  // namespace amflow::ode
