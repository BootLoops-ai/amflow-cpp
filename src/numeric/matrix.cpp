// SPDX-License-Identifier: MIT
// numeric::Rational — RationalMatrix and get_poles.
//
// amflow::numeric.

#include "amflow/numeric/rational.hpp"

#include <ostream>
#include <sstream>
#include <stdexcept>
#include <utility>
#include <vector>

#include <flint/acb.h>
#include <flint/arb.h>
#include <flint/arb_fmpz_poly.h>
#include <flint/arf.h>
#include <flint/fmpq.h>
#include <flint/fmpq_poly.h>
#include <flint/fmpz.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz_poly_factor.h>

namespace amflow::numeric {

// ===========================================================================
//  RationalMatrix
// ===========================================================================

RationalMatrix::RationalMatrix(std::size_t rows, std::size_t cols)
    : rows_(rows), cols_(cols), data_(rows * cols) {}

RationalMatrix RationalMatrix::identity(std::size_t n) {
    RationalMatrix m(n, n);
    for (std::size_t i = 0; i < n; ++i) m(i, i) = RationalFunction::from_si(1);
    return m;
}

RationalMatrix RationalMatrix::zeros(std::size_t rows, std::size_t cols) {
    return RationalMatrix(rows, cols);
}

RationalFunction& RationalMatrix::at(std::size_t i, std::size_t j) {
    if (i >= rows_ || j >= cols_) throw std::out_of_range("RationalMatrix::at");
    return data_[i * cols_ + j];
}
const RationalFunction& RationalMatrix::at(std::size_t i, std::size_t j) const {
    if (i >= rows_ || j >= cols_) throw std::out_of_range("RationalMatrix::at");
    return data_[i * cols_ + j];
}

double RationalMatrix::density() const {
    if (rows_ == 0 || cols_ == 0) return 0.0;
    std::size_t nz = 0;
    for (const auto& e : data_) if (!e.is_zero()) ++nz;
    return static_cast<double>(nz) / static_cast<double>(rows_ * cols_);
}

long RationalMatrix::max_numerator_degree() const {
    long m = -1;
    for (const auto& e : data_) {
        long d = e.degree_num();
        if (d > m) m = d;
    }
    return m;
}

long RationalMatrix::max_denominator_degree() const {
    long m = 0;
    for (const auto& e : data_) {
        long d = e.degree_den();
        if (d > m) m = d;
    }
    return m;
}

RationalMatrix RationalMatrix::submatrix(const std::vector<std::size_t>& rows,
                                         const std::vector<std::size_t>& cols) const {
    RationalMatrix out(rows.size(), cols.size());
    for (std::size_t i = 0; i < rows.size(); ++i) {
        for (std::size_t j = 0; j < cols.size(); ++j) {
            out(i, j) = at(rows[i], cols[j]);
        }
    }
    return out;
}

RationalMatrix RationalMatrix::substituted_inverse_eta() const {
    RationalMatrix out(rows_, cols_);
    for (std::size_t i = 0; i < rows_; ++i) {
        for (std::size_t j = 0; j < cols_; ++j) {
            out(i, j) = at(i, j).substitute_inverse_eta();
        }
    }
    return out;
}

FmpqPoly RationalMatrix::denominator_lcm(const std::vector<std::size_t>& row_indices,
                                         const std::vector<std::size_t>& col_indices) const {
    FmpqPoly lcm;
    lcm.set_one();
    bool initialised = false;
    for (auto i : row_indices) {
        for (auto j : col_indices) {
            const auto& d = at(i, j).denominator();
            if (!initialised) {
                fmpq_poly_set(lcm.raw(), d.raw());
                initialised = true;
            } else {
                fmpq_poly_lcm(lcm.raw(), lcm.raw(), d.raw());
            }
        }
    }
    if (!initialised) lcm.set_one();
    if (!lcm.is_zero() && lcm.degree() >= 0) {
        long d = lcm.degree();
        fmpq_t lead;
        fmpq_init(lead);
        fmpq_poly_get_coeff_fmpq(lead, lcm.raw(), d);
        if (!fmpq_is_zero(lead) && !fmpq_is_one(lead)) {
            fmpq_t inv;
            fmpq_init(inv);
            fmpq_inv(inv, lead);
            fmpq_poly_scalar_mul_fmpq(lcm.raw(), lcm.raw(), inv);
            fmpq_clear(inv);
        }
        fmpq_clear(lead);
    }
    return lcm;
}

FmpqPoly RationalMatrix::denominator_lcm() const {
    std::vector<std::size_t> rows(rows_);
    std::vector<std::size_t> cols(cols_);
    for (std::size_t i = 0; i < rows_; ++i) rows[i] = i;
    for (std::size_t j = 0; j < cols_; ++j) cols[j] = j;
    return denominator_lcm(rows, cols);
}

std::string RationalMatrix::to_string(const char* var) const {
    std::ostringstream oss;
    oss << "{\n";
    for (std::size_t i = 0; i < rows_; ++i) {
        oss << "  {";
        for (std::size_t j = 0; j < cols_; ++j) {
            if (j != 0) oss << ", ";
            oss << at(i, j).to_string(var);
        }
        oss << "}";
        if (i + 1 < rows_) oss << ",";
        oss << "\n";
    }
    oss << "}";
    return oss.str();
}

bool RationalMatrix::operator==(const RationalMatrix& other) const {
    if (rows_ != other.rows_ || cols_ != other.cols_) return false;
    for (std::size_t k = 0; k < data_.size(); ++k) {
        if (!(data_[k] == other.data_[k])) return false;
    }
    return true;
}

std::ostream& operator<<(std::ostream& os, const RationalMatrix& m) {
    return os << m.to_string();
}

// ===========================================================================
//  get_poles
// ===========================================================================

namespace {

// Round x to the nearest fraction p/q such that |x - p/q| < 10^-digits.
// Multiply midpoint by 10^digits, round to integer, divide.  Mirrors
// Mathematica's Rationalize.
void arb_to_rational(fmpq_t out, arb_srcptr x, int digits) {
    arf_t mid;
    arf_init(mid);
    arf_set(mid, arb_midref(x));

    if (arf_is_zero(mid)) {
        fmpq_zero(out);
        arf_clear(mid);
        return;
    }

    long bits = decimal_digits_to_bits(digits) + 32;
    arb_t scale, scaled;
    arb_init(scale);
    arb_init(scaled);
    arb_set_si(scale, 10);
    arb_pow_ui(scale, scale, static_cast<unsigned long>(digits), bits);

    arb_t xb;
    arb_init(xb);
    arb_set_arf(xb, mid);
    arb_mul(scaled, xb, scale, bits);

    fmpz_t numer, pow;
    fmpz_init(numer);
    fmpz_init(pow);

    arf_get_fmpz(numer, arb_midref(scaled), ARF_RND_NEAR);
    fmpz_set_si(pow, 10);
    fmpz_pow_ui(pow, pow, static_cast<unsigned long>(digits));

    fmpq_set_fmpz_frac(out, numer, pow);

    fmpz_clear(numer);
    fmpz_clear(pow);
    arb_clear(scale);
    arb_clear(scaled);
    arb_clear(xb);
    arf_clear(mid);
}

// Convert fmpq_poly to primitive fmpz_poly (clear common denominator + content).
void fmpq_poly_to_primitive_fmpz_poly(fmpz_poly_t out, const fmpq_poly_t in) {
    fmpz_poly_t numerator;
    fmpz_t denominator;
    fmpz_poly_init(numerator);
    fmpz_init(denominator);
    fmpq_poly_get_numerator(numerator, in);
    fmpq_poly_get_denominator(denominator, in);
    fmpz_poly_set(out, numerator);
    if (!fmpz_poly_is_zero(out)) {
        fmpz_t content;
        fmpz_init(content);
        fmpz_poly_content(content, out);
        if (!fmpz_is_one(content)) fmpz_poly_scalar_divexact_fmpz(out, out, content);
        fmpz_clear(content);
    }
    fmpz_poly_clear(numerator);
    fmpz_clear(denominator);
}

}  // namespace

std::vector<RationalComplex> get_poles(const RationalMatrix& matrix) {
    FmpqPoly lcm = matrix.denominator_lcm();
    if (lcm.is_zero()) return {};
    if (lcm.degree() <= 0) return {};

    fmpz_poly_t int_poly;
    fmpz_poly_init(int_poly);
    fmpq_poly_to_primitive_fmpz_poly(int_poly, lcm.raw());

    if (fmpz_poly_is_zero(int_poly) || fmpz_poly_degree(int_poly) <= 0) {
        fmpz_poly_clear(int_poly);
        return {};
    }

    fmpz_poly_factor_t factors;
    fmpz_poly_factor_init(factors);
    fmpz_poly_factor(factors, int_poly);

    std::vector<RationalComplex> result;

    long prec = working_prec_bits();
    int rdigits = rationalize_pre();
    if (rdigits <= 0) rdigits = 20;

    // PRECISION-MISMATCH FIX (see src/ode/inf.cpp:acb_real_to_fmpq_local):
    // cap rdigits at what working_prec can resolve so the rationalization of
    // pole locations doesn't capture binary representation noise.
    {
        long working_digits = static_cast<long>(prec * 0.301029995663981195L);
        if (working_digits >= 5) working_digits -= 5;
        if (working_digits < 1) working_digits = 1;
        if (rdigits > working_digits) rdigits = static_cast<int>(working_digits);
    }

    for (long fi = 0; fi < factors->num; ++fi) {
        fmpz_poly_struct* fac = factors->p + fi;
        long fac_deg = fmpz_poly_degree(fac);
        if (fac_deg <= 0) continue;

        acb_ptr roots = _acb_vec_init(fac_deg);
        arb_fmpz_poly_complex_roots(roots, fac, 0, prec);

        for (long ri = 0; ri < fac_deg; ++ri) {
            RationalComplex z;
            arb_to_rational(z.re(), acb_realref(roots + ri), rdigits);
            arb_to_rational(z.im(), acb_imagref(roots + ri), rdigits);

            bool seen = false;
            for (const auto& existing : result) {
                if (existing == z) { seen = true; break; }
            }
            if (!seen) result.push_back(std::move(z));
        }

        _acb_vec_clear(roots, fac_deg);
    }

    fmpz_poly_factor_clear(factors);
    fmpz_poly_clear(int_poly);

    return result;
}

}  // namespace amflow::numeric
