// SPDX-License-Identifier: MIT
// ode::jordan — implementation.
//

#include "amflow/ode/jordan.hpp"

#include <cstdio>
#include <cstdlib>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <flint/fmpq.h>
#include <flint/fmpq_mat.h>
#include <flint/fmpq_poly.h>
#include <flint/fmpz.h>
#include <flint/fmpz_mat.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz_poly_factor.h>

namespace amflow::ode {

using numeric::FmpqPoly;
using numeric::RationalFunction;
using numeric::RationalMatrix;

namespace {

struct FmpqMat {
    fmpq_mat_t data;
    FmpqMat(slong rows, slong cols) { fmpq_mat_init(data, rows, cols); }
    ~FmpqMat() { fmpq_mat_clear(data); }
    FmpqMat(const FmpqMat&) = delete;
    FmpqMat& operator=(const FmpqMat&) = delete;
};

struct FmpzMat {
    fmpz_mat_t data;
    FmpzMat(slong rows, slong cols) { fmpz_mat_init(data, rows, cols); }
    ~FmpzMat() { fmpz_mat_clear(data); }
    FmpzMat(const FmpzMat&) = delete;
    FmpzMat& operator=(const FmpzMat&) = delete;
};

struct FmpzVal {
    fmpz_t data;
    FmpzVal() { fmpz_init(data); }
    ~FmpzVal() { fmpz_clear(data); }
    FmpzVal(const FmpzVal&) = delete;
    FmpzVal& operator=(const FmpzVal&) = delete;
};

struct FmpqVal {
    fmpq_t data;
    FmpqVal() { fmpq_init(data); }
    ~FmpqVal() { fmpq_clear(data); }
    FmpqVal(const FmpqVal& o) { fmpq_init(data); fmpq_set(data, o.data); }
    FmpqVal(FmpqVal&& o) noexcept { fmpq_init(data); fmpq_swap(data, o.data); }
    FmpqVal& operator=(const FmpqVal& o) {
        if (this != &o) fmpq_set(data, o.data);
        return *this;
    }
    FmpqVal& operator=(FmpqVal&& o) noexcept {
        if (this != &o) fmpq_swap(data, o.data);
        return *this;
    }
};

void fmpq_mat_to_fmpz_mat_scaled(fmpz_mat_t out, const fmpq_mat_t in) {
    slong r = fmpq_mat_nrows(in);
    slong c = fmpq_mat_ncols(in);
    FmpzVal den;
    fmpz_one(den.data);
    for (slong i = 0; i < r; ++i) {
        for (slong j = 0; j < c; ++j) {
            fmpz_lcm(den.data, den.data, fmpq_mat_entry_den(in, i, j));
        }
    }
    FmpzVal scale;
    for (slong i = 0; i < r; ++i) {
        for (slong j = 0; j < c; ++j) {
            fmpz_divexact(scale.data, den.data, fmpq_mat_entry_den(in, i, j));
            fmpz_mul(fmpz_mat_entry(out, i, j),
                     fmpq_mat_entry_num(in, i, j), scale.data);
        }
    }
}

slong fmpq_mat_nullspace_exact(fmpq_mat_t basis, const fmpq_mat_t A) {
    slong nr = fmpq_mat_nrows(A);
    slong nc = fmpq_mat_ncols(A);
    if (nc == 0) return 0;
    if (nr == 0) {
        fmpq_mat_zero(basis);
        for (slong i = 0; i < nc; ++i) fmpq_one(fmpq_mat_entry(basis, i, i));
        return nc;
    }
    FmpzMat Az(nr, nc);
    fmpq_mat_to_fmpz_mat_scaled(Az.data, A);
    FmpzMat N(nc, nc);
    slong nullity = fmpz_mat_nullspace(N.data, Az.data);

    fmpq_mat_zero(basis);
    for (slong j = 0; j < nullity; ++j) {
        for (slong i = 0; i < nc; ++i) {
            fmpz_set(fmpq_numref(fmpq_mat_entry(basis, i, j)),
                     fmpz_mat_entry(N.data, i, j));
            fmpz_one(fmpq_denref(fmpq_mat_entry(basis, i, j)));
        }
        // Normalize each null-space basis vector so its last non-zero entry
        // is 1, matching Mathematica's Eigenvectors / JordanDecomposition
        // convention.  FLINT's fmpz_mat_nullspace returns vectors with
        // denominator-cleared integer entries -- e.g. an eigenvector that
        // Mathematica would emit as (6993/998, 1) is returned as (6993, 998).
        // Without this rescaling, downstream shearing / leading-Jordan T
        // blocks accumulate huge integer scale factors (1/eps-sized) that
        // cascade through the off-diagonal Sylvester step in
        // to_fuchsian_global, producing T entries up to 10^300+ in hard
        // cases (pentabox 2L 5-leg 76-master sub-system) that overwhelm
        // any practical working precision in PSMapRuleS.
        slong last_nz = -1;
        for (slong i = nc - 1; i >= 0; --i) {
            if (!fmpq_is_zero(fmpq_mat_entry(basis, i, j))) {
                last_nz = i;
                break;
            }
        }
        if (last_nz >= 0 && !fmpq_is_one(fmpq_mat_entry(basis, last_nz, j))) {
            fmpq_t pivot;
            fmpq_init(pivot);
            fmpq_set(pivot, fmpq_mat_entry(basis, last_nz, j));
            for (slong i = 0; i < nc; ++i) {
                fmpq_div(fmpq_mat_entry(basis, i, j),
                         fmpq_mat_entry(basis, i, j), pivot);
            }
            fmpq_clear(pivot);
        }
    }
    return nullity;
}

slong fmpq_mat_rank_exact(const fmpq_mat_t A) {
    slong nr = fmpq_mat_nrows(A);
    slong nc = fmpq_mat_ncols(A);
    if (nr == 0 || nc == 0) return 0;
    FmpqMat tmp(nr, nc);
    return fmpq_mat_rref(tmp.data, A);
}

bool column_independent_from_span(const fmpq_mat_t span,
                                  slong             n,
                                  slong             k_span,
                                  const fmpq_mat_t  src,
                                  slong             c) {
    FmpqMat aug(n, k_span + 1);
    for (slong j = 0; j < k_span; ++j) {
        for (slong i = 0; i < n; ++i) {
            fmpq_set(fmpq_mat_entry(aug.data, i, j),
                     fmpq_mat_entry(span, i, j));
        }
    }
    for (slong i = 0; i < n; ++i) {
        fmpq_set(fmpq_mat_entry(aug.data, i, k_span),
                 fmpq_mat_entry(src, i, c));
    }
    return fmpq_mat_rank_exact(aug.data) > k_span;
}

void fmpq_mat_append_col(fmpq_mat_t       out,
                         const fmpq_mat_t in,
                         slong            k_in,
                         const fmpq_mat_t src,
                         slong            c) {
    slong n = fmpq_mat_nrows(out);
    for (slong j = 0; j < k_in; ++j) {
        for (slong i = 0; i < n; ++i) {
            fmpq_set(fmpq_mat_entry(out, i, j),
                     fmpq_mat_entry(in, i, j));
        }
    }
    for (slong i = 0; i < n; ++i) {
        fmpq_set(fmpq_mat_entry(out, i, k_in),
                 fmpq_mat_entry(src, i, c));
    }
}

void fmpq_mat_a_minus_lambda_i(fmpq_mat_t       M,
                               const fmpq_mat_t A,
                               const fmpq_t     lambda) {
    slong n = fmpq_mat_nrows(A);
    fmpq_mat_set(M, A);
    for (slong i = 0; i < n; ++i) {
        fmpq_sub(fmpq_mat_entry(M, i, i),
                 fmpq_mat_entry(M, i, i), lambda);
    }
}

std::vector<std::pair<FmpqVal, long>>
rational_eigenvalues(const fmpq_mat_t A) {
    fmpq_poly_t p;
    fmpq_poly_init(p);
    fmpq_mat_charpoly(p, A);

    fmpz_poly_t pz;
    fmpz_poly_init(pz);
    fmpq_poly_get_numerator(pz, p);

    fmpz_poly_factor_t fac;
    fmpz_poly_factor_init(fac);
    fmpz_poly_factor(fac, pz);

    std::vector<std::pair<FmpqVal, long>> out;
    for (slong i = 0; i < fac->num; ++i) {
        const fmpz_poly_struct* fi = fac->p + i;
        slong deg = fmpz_poly_degree(fi);
        slong e = fac->exp[i];
        if (deg <= 0) continue;
        if (deg > 1) {
            fmpz_poly_factor_clear(fac);
            fmpz_poly_clear(pz);
            fmpq_poly_clear(p);
            throw std::runtime_error(
                "jordan_decomposition_exact: characteristic polynomial has "
                "an irreducible factor of degree > 1 (eigenvalue not in Q)");
        }

        FmpzVal a, b;
        fmpz_poly_get_coeff_fmpz(a.data, fi, 1);
        fmpz_poly_get_coeff_fmpz(b.data, fi, 0);
        fmpz_neg(b.data, b.data);
        FmpqVal root;
        fmpq_set_fmpz_frac(root.data, b.data, a.data);
        out.emplace_back(std::move(root), static_cast<long>(e));
    }

    fmpz_poly_factor_clear(fac);
    fmpz_poly_clear(pz);
    fmpq_poly_clear(p);
    return out;
}

}  // namespace

void jordan_decomposition_exact(fmpq_mat_t         out_S,
                                fmpq_mat_t         out_J,
                                fmpq_mat_t         out_Sinv,
                                std::vector<long>& out_block_sizes,
                                const fmpq_mat_t   A) {
    slong n = fmpq_mat_nrows(A);
    if (fmpq_mat_ncols(A) != n) {
        throw std::invalid_argument(
            "jordan_decomposition_exact: matrix must be square");
    }

    if (const char* dump_path = std::getenv("AMFLOW_JORDAN_DUMP")) {
        std::FILE* fp = std::fopen(dump_path, "a");
        if (fp) {
            std::fprintf(fp, "=== jordan_input n=%ld ===\n", (long)n);
            for (slong i = 0; i < n; ++i) {
                for (slong j = 0; j < n; ++j) {
                    const fmpq* q = fmpq_mat_entry(A, i, j);
                    char* num_s = fmpz_get_str(nullptr, 10, fmpq_numref(q));
                    char* den_s = fmpz_get_str(nullptr, 10, fmpq_denref(q));
                    std::fprintf(fp, "%s %s%s", num_s, den_s,
                                 (j + 1 == n) ? "\n" : "  ");
                    std::free(num_s);
                    std::free(den_s);
                }
            }
            std::fflush(fp);
            std::fclose(fp);
        }
    }

    fmpq_mat_zero(out_S);
    fmpq_mat_zero(out_J);
    fmpq_mat_zero(out_Sinv);
    out_block_sizes.clear();
    if (n == 0) return;

    auto eigs = rational_eigenvalues(A);
    long total_mult = 0;
    for (auto& e : eigs) total_mult += e.second;
    if (total_mult != static_cast<long>(n)) {
        throw std::runtime_error(
            "jordan_decomposition_exact: algebraic multiplicities do not "
            "sum to matrix dimension");
    }

    slong col = 0;

    for (auto& [lambda, mult] : eigs) {
        FmpqMat N(n, n);
        fmpq_mat_a_minus_lambda_i(N.data, A, lambda.data);

        std::vector<FmpqMat*> levels;
        std::vector<slong> dims;

        FmpqMat Npow(n, n);
        fmpq_mat_set(Npow.data, N.data);

        for (slong k = 1; k <= mult; ++k) {
            auto* lvl = new FmpqMat(n, n);
            slong d = fmpq_mat_nullspace_exact(lvl->data, Npow.data);
            levels.push_back(lvl);
            dims.push_back(d);

            if (dims.size() >= 2 &&
                dims[dims.size() - 1] == dims[dims.size() - 2]) {
                delete levels.back();
                levels.pop_back();
                dims.pop_back();
                break;
            }
            if (d >= mult) break;
            if (k == mult) break;
            FmpqMat next(n, n);
            fmpq_mat_mul(next.data, Npow.data, N.data);
            fmpq_mat_swap(Npow.data, next.data);
        }

        slong depth = static_cast<slong>(levels.size());
        if (depth == 0 || dims.back() != mult) {
            for (auto* lv : levels) delete lv;
            throw std::runtime_error(
                "jordan_decomposition_exact: kernel tower top dimension "
                "does not match algebraic multiplicity");
        }

        std::vector<slong> r(depth + 2, 0);
        for (slong s = 1; s <= depth; ++s) {
            slong prev = (s == 1) ? 0 : dims[s - 2];
            r[s] = dims[s - 1] - prev;
        }
        std::vector<slong> blocks_of_size(depth + 2, 0);
        for (slong s = 1; s <= depth; ++s) {
            blocks_of_size[s] = r[s] - r[s + 1];
        }

        FmpqMat chosen_span(n, n);
        slong placed_in_span = 0;

        for (slong s = depth; s >= 1; --s) {
            slong want = blocks_of_size[s];
            if (want <= 0) continue;

            slong base_cols = (s >= 2) ? dims[s - 2] : 0;
            slong excl_cols = base_cols + placed_in_span;
            FmpqMat excl(n, excl_cols);
            for (slong j = 0; j < base_cols; ++j) {
                for (slong i = 0; i < n; ++i) {
                    fmpq_set(fmpq_mat_entry(excl.data, i, j),
                             fmpq_mat_entry(levels[s - 2]->data, i, j));
                }
            }
            for (slong j = 0; j < placed_in_span; ++j) {
                for (slong i = 0; i < n; ++i) {
                    fmpq_set(fmpq_mat_entry(excl.data, i, base_cols + j),
                             fmpq_mat_entry(chosen_span.data, i, j));
                }
            }

            slong ks = dims[s - 1];
            slong got = 0;

            for (slong j = 0; j < ks && got < want; ++j) {
                if (!column_independent_from_span(
                        excl.data, n, excl_cols, levels[s - 1]->data, j)) {
                    continue;
                }

                FmpqMat excl_next(n, excl_cols + 1);
                fmpq_mat_append_col(excl_next.data, excl.data, excl_cols,
                                    levels[s - 1]->data, j);
                fmpq_mat_swap(excl.data, excl_next.data);
                ++excl_cols;

                std::vector<FmpqMat*> chain;
                chain.reserve(s);
                for (slong t = 0; t < s; ++t) chain.push_back(new FmpqMat(n, 1));
                for (slong i = 0; i < n; ++i) {
                    fmpq_set(fmpq_mat_entry(chain[s - 1]->data, i, 0),
                             fmpq_mat_entry(levels[s - 1]->data, i, j));
                }
                for (slong t = s - 1; t >= 1; --t) {
                    fmpq_mat_mul(chain[t - 1]->data, N.data, chain[t]->data);
                }

                for (slong t = 0; t < s; ++t) {
                    slong target = col + t;
                    for (slong i = 0; i < n; ++i) {
                        fmpq_set(fmpq_mat_entry(out_S, i, target),
                                 fmpq_mat_entry(chain[t]->data, i, 0));
                    }
                    fmpq_set(fmpq_mat_entry(out_J, target, target), lambda.data);
                    if (t > 0) {
                        fmpq_one(fmpq_mat_entry(out_J, target - 1, target));
                    }
                    for (slong i = 0; i < n; ++i) {
                        fmpq_set(fmpq_mat_entry(chosen_span.data, i, placed_in_span),
                                 fmpq_mat_entry(chain[t]->data, i, 0));
                    }
                    ++placed_in_span;
                }
                out_block_sizes.push_back(static_cast<long>(s));
                col += s;

                for (auto* c : chain) delete c;
                ++got;
            }

            if (got < want) {
                for (auto* lv : levels) delete lv;
                throw std::runtime_error(
                    "jordan_decomposition_exact: insufficient independent "
                    "Jordan-chain tops at depth " + std::to_string(s));
            }
        }

        for (auto* lv : levels) delete lv;
    }

    if (col != n) {
        throw std::runtime_error(
            "jordan_decomposition_exact: placed column count does not match "
            "matrix dimension");
    }

    int inv_ok = fmpq_mat_inv(out_Sinv, out_S);
    if (!inv_ok) {
        throw std::runtime_error(
            "jordan_decomposition_exact: assembled Jordan basis is singular");
    }
}

void fmpq_mat_from_rational_residue(fmpq_mat_t            out,
                                    const RationalMatrix& mat,
                                    long                   shift) {
    slong r = static_cast<slong>(mat.rows());
    slong c = static_cast<slong>(mat.cols());
    if (fmpq_mat_nrows(out) != r || fmpq_mat_ncols(out) != c) {
        throw std::invalid_argument(
            "fmpq_mat_from_rational_residue: output matrix has wrong shape");
    }
    fmpq_mat_zero(out);

    for (slong i = 0; i < r; ++i) {
        for (slong j = 0; j < c; ++j) {
            RationalFunction tmp = mat(i, j);
            if (shift != 0) tmp.multiply_by_eta_power(shift);

            if (tmp.is_zero()) {
                fmpq_zero(fmpq_mat_entry(out, i, j));
                continue;
            }

            const FmpqPoly& num = tmp.numerator();
            const FmpqPoly& den = tmp.denominator();

            FmpqVal den0;
            den.coeff(0, den0.data);
            if (fmpq_is_zero(den0.data)) {
                throw std::runtime_error(
                    "fmpq_mat_from_rational_residue: entry still has a pole "
                    "at eta = 0 after the shift");
            }
            FmpqVal num0;
            num.coeff(0, num0.data);
            fmpq_div(fmpq_mat_entry(out, i, j), num0.data, den0.data);
        }
    }
}

RationalMatrix rational_matrix_from_fmpq_mat(const fmpq_mat_t in) {
    std::size_t r = static_cast<std::size_t>(fmpq_mat_nrows(in));
    std::size_t c = static_cast<std::size_t>(fmpq_mat_ncols(in));
    RationalMatrix out(r, c);
    for (std::size_t i = 0; i < r; ++i) {
        for (std::size_t j = 0; j < c; ++j) {
            if (fmpq_is_zero(fmpq_mat_entry(in, i, j))) {
                out(i, j) = RationalFunction();
            } else {
                out(i, j) = RationalFunction::from_fmpq(fmpq_mat_entry(in, i, j));
            }
        }
    }
    return out;
}

}  // namespace amflow::ode
