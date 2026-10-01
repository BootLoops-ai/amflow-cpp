// SPDX-License-Identifier: MIT
// ode::normalize — implementation.
//

#include "amflow/ode/normalize.hpp"

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <utility>
#include <vector>

#include <flint/acb.h>
#include <flint/acb_mat.h>
#include <flint/arb.h>
#include <flint/ca.h>
#include <flint/ca_mat.h>
#include <flint/fmpq.h>
#include <flint/fmpq_mat.h>
#include <flint/fmpq_poly.h>
#include <flint/fmpz.h>
#include <flint/fmpz_mat.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz_poly_factor.h>

#include "amflow/numeric/options.hpp"
#include "amflow/ode/acb_rational.hpp"
#include "amflow/ode/blocks.hpp"
#include "amflow/ode/jordan.hpp"
#include "amflow/ode/qqbar.hpp"

namespace amflow::ode {

using numeric::AcbValue;
using numeric::FmpqPoly;
using numeric::RationalFunction;
using numeric::RationalMatrix;
using numeric::acb_is_chop_zero;
using numeric::chop_pre;
using numeric::log_line;
using numeric::resonance_chop_pre;
using numeric::silent_mode;
using numeric::working_prec_bits;

// ===========================================================================
//  Utilities
// ===========================================================================

namespace {

std::string block_to_string(const std::vector<std::size_t>& block) {
    std::ostringstream oss;
    oss << "[";
    for (std::size_t i = 0; i < block.size(); ++i) {
        if (i) oss << ",";
        oss << block[i];
    }
    oss << "]";
    return oss.str();
}

BlockBehavior exact_or_numeric_leading_behavior(const RationalMatrix& sub,
                                                long prec);

// UnionBehavior:  group by mu, take the maximum log_power per group.
BlockBehavior union_behavior(const BlockBehavior& in, long prec) {
    BlockBehavior out;
    std::vector<bool> done(in.size(), false);
    for (std::size_t i = 0; i < in.size(); ++i) {
        if (done[i]) continue;
        BehaviorEntry combined;
        combined.mu = in[i].mu.clone();
        combined.log_power = in[i].log_power;
        done[i] = true;
        for (std::size_t j = i + 1; j < in.size(); ++j) {
            if (done[j]) continue;
            AcbValue diff;
            acb_sub(diff.raw(), in[i].mu.raw(), in[j].mu.raw(), prec);
            if (acb_is_chop_zero(diff.raw(), resonance_chop_pre())) {
                if (in[j].log_power > combined.log_power) combined.log_power = in[j].log_power;
                done[j] = true;
            }
        }
        out.push_back(std::move(combined));
    }
    return out;
}

// UnionBehavior2: group by mu; if a group has a single member keep it,
// otherwise sum log_powers + 1 (subbeh + essbeh contribution).
BlockBehavior union_behavior2(const BlockBehavior& in, long prec) {
    BlockBehavior out;
    std::vector<bool> done(in.size(), false);
    for (std::size_t i = 0; i < in.size(); ++i) {
        if (done[i]) continue;
        std::vector<long> log_powers;
        log_powers.push_back(in[i].log_power);
        AcbValue mu_keep = in[i].mu.clone();
        done[i] = true;
        for (std::size_t j = i + 1; j < in.size(); ++j) {
            if (done[j]) continue;
            AcbValue diff;
            acb_sub(diff.raw(), in[i].mu.raw(), in[j].mu.raw(), prec);
            if (acb_is_chop_zero(diff.raw(), resonance_chop_pre())) {
                log_powers.push_back(in[j].log_power);
                done[j] = true;
            }
        }
        BehaviorEntry e;
        e.mu = std::move(mu_keep);
        if (log_powers.size() == 1) {
            e.log_power = log_powers[0];
        } else {
            long s = 0;
            for (long lp : log_powers) s += lp;
            e.log_power = s + 1;
        }
        out.push_back(std::move(e));
    }
    return out;
}

}  // namespace

AsymptoticBehaviorList asymptotic_behavior(const RationalMatrix& mat, long prec) {
    auto partition = build_partition(mat);
    AsymptoticBehaviorList behavior(partition.size());

    for (std::size_t i = 0; i < partition.size(); ++i) {
        const auto& blk = partition.block(i);
        RationalMatrix sub = mat.submatrix(blk, blk);
        BlockBehavior essbeh = exact_or_numeric_leading_behavior(sub, prec);
        essbeh = union_behavior(essbeh, prec);

        BlockBehavior subbeh;
        for (auto j : partition.lower_dependencies(i)) {
            for (const auto& bh : behavior[j]) {
                BehaviorEntry e;
                e.mu = bh.mu.clone();
                e.log_power = bh.log_power;
                subbeh.push_back(std::move(e));
            }
        }
        subbeh = union_behavior(subbeh, prec);

        BlockBehavior joined;
        for (auto& e : essbeh) joined.push_back(BehaviorEntry{e.mu.clone(), e.log_power});
        for (auto& e : subbeh) joined.push_back(BehaviorEntry{e.mu.clone(), e.log_power});
        behavior[i] = union_behavior2(joined, prec);
    }
    return behavior;
}

// ===========================================================================
//  Sylvester-style solve used by SolveOffDiagonal
// ===========================================================================

namespace {

void vectorise(acb_mat_t v, const acb_mat_t M) {
    long n = acb_mat_nrows(M);
    long m = acb_mat_ncols(M);
    if (acb_mat_nrows(v) != n * m || acb_mat_ncols(v) != 1) {
        throw std::invalid_argument("vectorise: shape mismatch");
    }
    for (long j = 0; j < m; ++j) {
        for (long i = 0; i < n; ++i) {
            acb_set(acb_mat_entry(v, j * n + i, 0),
                    acb_mat_entry(M, i, j));
        }
    }
}

void unvectorise(acb_mat_t M, const acb_mat_t v) {
    long n = acb_mat_nrows(M);
    long m = acb_mat_ncols(M);
    if (acb_mat_nrows(v) != n * m || acb_mat_ncols(v) != 1) {
        throw std::invalid_argument("unvectorise: shape mismatch");
    }
    for (long j = 0; j < m; ++j) {
        for (long i = 0; i < n; ++i) {
            acb_set(acb_mat_entry(M, i, j),
                    acb_mat_entry(v, j * n + i, 0));
        }
    }
}

}  // namespace

namespace internal {

bool solve_off_diagonal_acb(acb_mat_t X,
                            const acb_mat_t a0,
                            const acb_mat_t b0,
                            const acb_mat_t c0,
                            long p,
                            long prec) {
    long m = acb_mat_nrows(a0);
    long n = acb_mat_nrows(c0);
    if (acb_mat_ncols(a0) != m || acb_mat_ncols(c0) != n
        || acb_mat_nrows(b0) != n || acb_mat_ncols(b0) != m) {
        throw std::invalid_argument("solve_off_diagonal_acb: shape mismatch");
    }
    long nm = n * m;

    acb_mat_t M;
    acb_mat_init(M, nm, nm);
    acb_mat_zero(M);

    // I_m ⊗ c0
    for (long jb = 0; jb < m; ++jb) {
        for (long i = 0; i < n; ++i) {
            for (long j = 0; j < n; ++j) {
                AcbValue tmp;
                acb_set(tmp.raw(), acb_mat_entry(M, jb * n + i, jb * n + j));
                acb_add(tmp.raw(), tmp.raw(), acb_mat_entry(c0, i, j), prec);
                acb_set(acb_mat_entry(M, jb * n + i, jb * n + j), tmp.raw());
            }
        }
    }
    // - a0^T ⊗ I_n
    for (long jb1 = 0; jb1 < m; ++jb1) {
        for (long jb2 = 0; jb2 < m; ++jb2) {
            AcbValue scalar;
            acb_neg(scalar.raw(), acb_mat_entry(a0, jb2, jb1));
            for (long i = 0; i < n; ++i) {
                AcbValue tmp;
                acb_set(tmp.raw(), acb_mat_entry(M, jb1 * n + i, jb2 * n + i));
                acb_add(tmp.raw(), tmp.raw(), scalar.raw(), prec);
                acb_set(acb_mat_entry(M, jb1 * n + i, jb2 * n + i), tmp.raw());
            }
        }
    }
    // + p I_{nm}
    for (long i = 0; i < nm; ++i) {
        AcbValue tmp;
        acb_set(tmp.raw(), acb_mat_entry(M, i, i));
        AcbValue add; add.set_si(p);
        acb_add(tmp.raw(), tmp.raw(), add.raw(), prec);
        acb_set(acb_mat_entry(M, i, i), tmp.raw());
    }

    acb_mat_t rhs;
    acb_mat_init(rhs, nm, 1);
    {
        acb_mat_t neg_b;
        acb_mat_init(neg_b, n, m);
        acb_mat_neg(neg_b, b0);
        vectorise(rhs, neg_b);
        acb_mat_clear(neg_b);
    }

    acb_mat_t sol;
    acb_mat_init(sol, nm, 1);
    int ok = acb_mat_solve(sol, M, rhs, prec);

    if (ok) unvectorise(X, sol);

    acb_mat_clear(sol);
    acb_mat_clear(rhs);
    acb_mat_clear(M);
    return ok != 0;
}

bool solve_off_diagonal_fmpq(fmpq_mat_t X,
                             const fmpq_mat_t a0,
                             const fmpq_mat_t b0,
                             const fmpq_mat_t c0,
                             long p) {
    slong m = fmpq_mat_nrows(a0);
    slong n = fmpq_mat_nrows(c0);
    if (fmpq_mat_ncols(a0) != m || fmpq_mat_ncols(c0) != n
        || fmpq_mat_nrows(b0) != n || fmpq_mat_ncols(b0) != m
        || fmpq_mat_nrows(X) != n || fmpq_mat_ncols(X) != m) {
        throw std::invalid_argument("solve_off_diagonal_fmpq: shape mismatch");
    }

    slong nm = n * m;
    fmpq_mat_t M;
    fmpq_mat_init(M, nm, nm);
    fmpq_mat_zero(M);

    for (slong jb = 0; jb < m; ++jb) {
        for (slong i = 0; i < n; ++i) {
            for (slong j = 0; j < n; ++j) {
                fmpq_add(fmpq_mat_entry(M, jb * n + i, jb * n + j),
                         fmpq_mat_entry(M, jb * n + i, jb * n + j),
                         fmpq_mat_entry(c0, i, j));
            }
        }
    }
    for (slong jb1 = 0; jb1 < m; ++jb1) {
        for (slong jb2 = 0; jb2 < m; ++jb2) {
            for (slong i = 0; i < n; ++i) {
                fmpq_sub(fmpq_mat_entry(M, jb1 * n + i, jb2 * n + i),
                         fmpq_mat_entry(M, jb1 * n + i, jb2 * n + i),
                         fmpq_mat_entry(a0, jb2, jb1));
            }
        }
    }
    if (p != 0) {
        fmpq_t pq;
        fmpq_init(pq);
        fmpq_set_si(pq, p, 1);
        for (slong i = 0; i < nm; ++i) {
            fmpq_add(fmpq_mat_entry(M, i, i),
                     fmpq_mat_entry(M, i, i), pq);
        }
        fmpq_clear(pq);
    }

    fmpq_mat_t rhs;
    fmpq_mat_init(rhs, nm, 1);
    for (slong j = 0; j < m; ++j) {
        for (slong i = 0; i < n; ++i) {
            fmpq_neg(fmpq_mat_entry(rhs, j * n + i, 0),
                     fmpq_mat_entry(b0, i, j));
        }
    }

    fmpq_mat_t sol;
    fmpq_mat_init(sol, nm, 1);
    int ok = fmpq_mat_solve_dixon(sol, M, rhs);
    if (ok) {
        for (slong j = 0; j < m; ++j) {
            for (slong i = 0; i < n; ++i) {
                fmpq_set(fmpq_mat_entry(X, i, j),
                         fmpq_mat_entry(sol, j * n + i, 0));
            }
        }
    }

    fmpq_mat_clear(sol);
    fmpq_mat_clear(rhs);
    fmpq_mat_clear(M);
    return ok != 0;
}

}  // namespace internal

// ===========================================================================
//  PBar / Balance / InvBalance + ToFuchsian
// ===========================================================================

namespace {

struct FmpqMatLocal {
    fmpq_mat_t data;
    FmpqMatLocal(slong rows, slong cols) { fmpq_mat_init(data, rows, cols); }
    ~FmpqMatLocal() { fmpq_mat_clear(data); }
    FmpqMatLocal(const FmpqMatLocal&) = delete;
    FmpqMatLocal& operator=(const FmpqMatLocal&) = delete;
    FmpqMatLocal(FmpqMatLocal&& other) noexcept {
        fmpq_mat_init(data, 0, 0);
        fmpq_mat_swap(data, other.data);
    }
    FmpqMatLocal& operator=(FmpqMatLocal&& other) noexcept {
        if (this != &other) fmpq_mat_swap(data, other.data);
        return *this;
    }
};

struct FmpzMatLocal {
    fmpz_mat_t data;
    FmpzMatLocal(slong rows, slong cols) { fmpz_mat_init(data, rows, cols); }
    ~FmpzMatLocal() { fmpz_mat_clear(data); }
    FmpzMatLocal(const FmpzMatLocal&) = delete;
    FmpzMatLocal& operator=(const FmpzMatLocal&) = delete;
};

struct FmpzValLocal {
    fmpz_t data;
    FmpzValLocal() { fmpz_init(data); }
    ~FmpzValLocal() { fmpz_clear(data); }
    FmpzValLocal(const FmpzValLocal&) = delete;
    FmpzValLocal& operator=(const FmpzValLocal&) = delete;
};

void fmpq_mat_to_fmpz_mat_scaled_local(fmpz_mat_t out,
                                       const fmpq_mat_t in) {
    slong r = fmpq_mat_nrows(in);
    slong c = fmpq_mat_ncols(in);
    FmpzValLocal den;
    fmpz_one(den.data);
    for (slong i = 0; i < r; ++i) {
        for (slong j = 0; j < c; ++j) {
            fmpz_lcm(den.data, den.data, fmpq_mat_entry_den(in, i, j));
        }
    }
    FmpzValLocal scale;
    for (slong i = 0; i < r; ++i) {
        for (slong j = 0; j < c; ++j) {
            fmpz_divexact(scale.data, den.data, fmpq_mat_entry_den(in, i, j));
            fmpz_mul(fmpz_mat_entry(out, i, j),
                     fmpq_mat_entry_num(in, i, j), scale.data);
        }
    }
}

slong fmpq_mat_rank_local(const fmpq_mat_t A) {
    slong nr = fmpq_mat_nrows(A);
    slong nc = fmpq_mat_ncols(A);
    if (nr == 0 || nc == 0) return 0;
    FmpqMatLocal tmp(nr, nc);
    return fmpq_mat_rref(tmp.data, A);
}

slong fmpq_mat_nullspace_local(fmpq_mat_t basis, const fmpq_mat_t A) {
    slong nr = fmpq_mat_nrows(A);
    slong nc = fmpq_mat_ncols(A);
    if (nc == 0) return 0;
    if (nr == 0) {
        fmpq_mat_zero(basis);
        for (slong i = 0; i < nc; ++i) fmpq_one(fmpq_mat_entry(basis, i, i));
        return nc;
    }

    FmpzMatLocal Az(nr, nc);
    fmpq_mat_to_fmpz_mat_scaled_local(Az.data, A);
    FmpzMatLocal N(nc, nc);
    slong nullity = fmpz_mat_nullspace(N.data, Az.data);

    fmpq_mat_zero(basis);
    for (slong j = 0; j < nullity; ++j) {
        for (slong i = 0; i < nc; ++i) {
            fmpz_set(fmpq_numref(fmpq_mat_entry(basis, i, j)),
                     fmpz_mat_entry(N.data, i, j));
            fmpz_one(fmpq_denref(fmpq_mat_entry(basis, i, j)));
        }
    }
    return nullity;
}

FmpqMatLocal fmpq_rows_prefix(const fmpq_mat_t A,
                              const std::vector<slong>& rows_1based,
                              slong prefix_cols) {
    FmpqMatLocal out(static_cast<slong>(rows_1based.size()), prefix_cols);
    for (slong i = 0; i < static_cast<slong>(rows_1based.size()); ++i) {
        slong src_i = rows_1based[static_cast<std::size_t>(i)] - 1;
        for (slong j = 0; j < prefix_cols; ++j) {
            fmpq_set(fmpq_mat_entry(out.data, i, j),
                     fmpq_mat_entry(A, src_i, j));
        }
    }
    return out;
}

slong fmpq_rank_rows_prefix(const fmpq_mat_t A,
                            const std::vector<slong>& rows_1based,
                            slong prefix_cols) {
    if (prefix_cols <= 0 || rows_1based.empty()) return 0;
    auto sub = fmpq_rows_prefix(A, rows_1based, prefix_cols);
    return fmpq_mat_rank_local(sub.data);
}

bool contains_slong(const std::vector<slong>& values, slong x) {
    return std::find(values.begin(), values.end(), x) != values.end();
}

void reduce_l0_exact(fmpq_mat_t Delta,
                     slong& k0,
                     std::vector<slong>& selected_blocks,
                     const fmpq_mat_t L0,
                     slong r,
                     const std::vector<long>& lblock) {
    slong nblocks = fmpq_mat_nrows(L0);
    if (fmpq_mat_ncols(L0) != nblocks) {
        throw std::invalid_argument("reduce_l0_exact: L0 must be square");
    }
    if (static_cast<slong>(lblock.size()) != nblocks) {
        throw std::invalid_argument("reduce_l0_exact: block-size mismatch");
    }

    fmpq_mat_zero(Delta);
    FmpqMatLocal L0tmp(nblocks, nblocks);
    fmpq_mat_set(L0tmp.data, L0);

    std::vector<slong> S_1based;
    slong current_n_1based = 1;

    while (current_n_1based > r) {
        std::vector<slong> comS;
        for (slong i = 1; i <= nblocks; ++i) {
            if (!contains_slong(S_1based, i)) comS.push_back(i);
        }
        if (comS.empty()) {
            throw std::runtime_error("reduce_l0_exact: empty complement");
        }

        bool found = false;
        for (slong candidate : comS) {
            slong rank_prev = fmpq_rank_rows_prefix(L0tmp.data, comS,
                                                    candidate - 1);
            slong rank_curr = fmpq_rank_rows_prefix(L0tmp.data, comS,
                                                    candidate);
            if (rank_prev == rank_curr) {
                current_n_1based = candidate;
                found = true;
                break;
            }
        }
        if (!found) {
            throw std::runtime_error(
                "reduce_l0_exact: no dependent prefix column found");
        }

        auto L0t_prefix = fmpq_rows_prefix(L0tmp.data, comS,
                                           current_n_1based);
        FmpqMatLocal null_basis(current_n_1based, current_n_1based);
        slong nullity = fmpq_mat_nullspace_local(null_basis.data,
                                                 L0t_prefix.data);
        if (nullity != 1) {
            throw std::runtime_error(
                "reduce_l0_exact: L0 prefix nullspace is not one-dimensional");
        }

        fmpq_t last;
        fmpq_init(last);
        fmpq_set(last, fmpq_mat_entry(null_basis.data,
                                      current_n_1based - 1, 0));
        if (fmpq_is_zero(last)) {
            fmpq_clear(last);
            throw std::runtime_error(
                "reduce_l0_exact: null vector has zero pivot component");
        }

        FmpqMatLocal Delta0(nblocks, nblocks);
        FmpqMatLocal Delta0t(nblocks, nblocks);
        fmpq_mat_zero(Delta0.data);
        fmpq_mat_zero(Delta0t.data);
        for (slong j = 0; j < current_n_1based - 1; ++j) {
            fmpq_div(fmpq_mat_entry(Delta0.data, j, current_n_1based - 1),
                     fmpq_mat_entry(null_basis.data, j, 0), last);
            if (lblock[static_cast<std::size_t>(j)] ==
                lblock[static_cast<std::size_t>(current_n_1based - 1)]) {
                fmpq_set(fmpq_mat_entry(Delta0t.data, j, current_n_1based - 1),
                         fmpq_mat_entry(Delta0.data, j, current_n_1based - 1));
            }
        }
        fmpq_clear(last);

        FmpqMatLocal I(nblocks, nblocks);
        FmpqMatLocal left(nblocks, nblocks);
        FmpqMatLocal right(nblocks, nblocks);
        FmpqMatLocal tmp(nblocks, nblocks);
        fmpq_mat_one(I.data);
        fmpq_mat_sub(left.data, I.data, Delta0t.data);
        fmpq_mat_add(right.data, I.data, Delta0.data);
        fmpq_mat_mul(tmp.data, left.data, L0tmp.data);
        fmpq_mat_mul(L0tmp.data, tmp.data, right.data);

        FmpqMatLocal delta_prod(nblocks, nblocks);
        fmpq_mat_mul(delta_prod.data, Delta, Delta0.data);
        fmpq_mat_add(Delta, Delta, Delta0.data);
        fmpq_mat_add(Delta, Delta, delta_prod.data);

        S_1based.push_back(current_n_1based);
    }

    k0 = current_n_1based - 1;
    selected_blocks.clear();
    for (slong s : S_1based) {
        if (s != current_n_1based) selected_blocks.push_back(s - 1);
    }
}

std::vector<slong> jordan_block_starts(const std::vector<long>& block_sizes) {
    std::vector<slong> starts;
    starts.reserve(block_sizes.size());
    slong pos = 0;
    for (long size : block_sizes) {
        starts.push_back(pos);
        pos += size;
    }
    return starts;
}

void fmpq_bilinear_entry(fmpq_t out,
                         const fmpq_mat_t row_source,
                         slong row_index,
                         const fmpq_mat_t M,
                         const fmpq_mat_t col_source,
                         slong col_index) {
    slong n = fmpq_mat_ncols(row_source);
    fmpq_zero(out);
    fmpq_t tmp1, tmp2;
    fmpq_init(tmp1);
    fmpq_init(tmp2);
    for (slong a = 0; a < n; ++a) {
        for (slong b = 0; b < n; ++b) {
            fmpq_mul(tmp1, fmpq_mat_entry(row_source, row_index, a),
                     fmpq_mat_entry(M, a, b));
            fmpq_mul(tmp2, tmp1, fmpq_mat_entry(col_source, b, col_index));
            fmpq_add(out, out, tmp2);
        }
    }
    fmpq_clear(tmp2);
    fmpq_clear(tmp1);
}

void fmpq_row_col_entry(fmpq_t out,
                        const fmpq_mat_t row_source,
                        slong row_index,
                        const fmpq_mat_t col_source,
                        slong col_index) {
    slong n = fmpq_mat_ncols(row_source);
    fmpq_zero(out);
    fmpq_t prod;
    fmpq_init(prod);
    for (slong a = 0; a < n; ++a) {
        fmpq_mul(prod, fmpq_mat_entry(row_source, row_index, a),
                 fmpq_mat_entry(col_source, a, col_index));
        fmpq_add(out, out, prod);
    }
    fmpq_clear(prod);
}

RationalMatrix find_projector_exact(const fmpq_mat_t cp,
                                    const fmpq_mat_t cp1) {
    slong n = fmpq_mat_nrows(cp);
    if (fmpq_mat_ncols(cp) != n || fmpq_mat_nrows(cp1) != n ||
        fmpq_mat_ncols(cp1) != n) {
        throw std::invalid_argument(
            "find_projector_exact: matrices must be square with same size");
    }

    FmpqMatLocal S(n, n);
    FmpqMatLocal J(n, n);
    FmpqMatLocal Sinv(n, n);
    std::vector<long> block_sizes;
    // Mirror DESolver.m:347 -- abort on any non-zero eigenvalue.
    try {
        jordan_decomposition_exact(S.data, J.data, Sinv.data, block_sizes, cp);
    } catch (const std::exception&) {
        throw std::runtime_error(
            "find_projector_exact: irreducible Poincare rank "
            "(exact Jordan decomposition failed)");
    }
    for (slong i = 0; i < n; ++i) {
        if (!fmpq_is_zero(fmpq_mat_entry(J.data, i, i))) {
            throw std::runtime_error(
                "find_projector_exact: irreducible Poincare rank "
                "(non-nilpotent leading matrix)");
        }
    }
    if (block_sizes.empty()) {
        return RationalMatrix::zeros(static_cast<std::size_t>(n),
                                     static_cast<std::size_t>(n));
    }

    std::vector<slong> starts = jordan_block_starts(block_sizes);
    slong nb = static_cast<slong>(block_sizes.size());

    FmpqMatLocal L0(nb, nb);
    FmpqMatLocal L1(nb, nb);
    fmpq_t val;
    fmpq_init(val);
    for (slong i = 0; i < nb; ++i) {
        slong left_first_row = starts[static_cast<std::size_t>(i)]
                             + block_sizes[static_cast<std::size_t>(i)] - 1;
        for (slong j = 0; j < nb; ++j) {
            slong right_first_col = starts[static_cast<std::size_t>(j)];
            fmpq_bilinear_entry(val, Sinv.data, left_first_row, cp1,
                                S.data, right_first_col);
            fmpq_set(fmpq_mat_entry(L0.data, i, j), val);
            fmpq_row_col_entry(val, Sinv.data, left_first_row,
                               S.data, right_first_col);
            fmpq_set(fmpq_mat_entry(L1.data, i, j), val);
        }
    }
    fmpq_clear(val);

    slong r = 0;
    for (slong i = 0; i < nb; ++i) {
        if (fmpq_is_zero(fmpq_mat_entry(L1.data, i, i))) ++r;
    }

    FmpqMatLocal Delta(nb, nb);
    slong k0 = 0;
    std::vector<slong> selected;
    reduce_l0_exact(Delta.data, k0, selected, L0.data, r, block_sizes);
    selected.push_back(k0);

    FmpqMatLocal Et(n, n);
    fmpq_mat_zero(Et.data);
    for (slong i = 0; i < nb; ++i) {
        for (slong j = i + 1; j < nb; ++j) {
            const fmpq* d = fmpq_mat_entry(Delta.data, i, j);
            if (fmpq_is_zero(d)) continue;
            long len_i = block_sizes[static_cast<std::size_t>(i)];
            long len_j = block_sizes[static_cast<std::size_t>(j)];
            long common = std::min(len_i, len_j);
            for (long pos = 0; pos < common; ++pos) {
                fmpq_set(fmpq_mat_entry(Et.data,
                                        starts[static_cast<std::size_t>(i)] + pos,
                                        starts[static_cast<std::size_t>(j)] + pos),
                         d);
            }
        }
    }

    FmpqMatLocal I(n, n);
    FmpqMatLocal IplusEt(n, n);
    FmpqMatLocal U(n, n);
    FmpqMatLocal invIplusEt(n, n);
    FmpqMatLocal invU(n, n);
    fmpq_mat_one(I.data);
    fmpq_mat_add(IplusEt.data, I.data, Et.data);
    fmpq_mat_mul(U.data, S.data, IplusEt.data);
    if (!fmpq_mat_inv(invIplusEt.data, IplusEt.data)) {
        throw std::runtime_error(
            "find_projector_exact: I + Et is singular");
    }
    fmpq_mat_mul(invU.data, invIplusEt.data, Sinv.data);

    FmpqMatLocal Q(n, n);
    fmpq_mat_zero(Q.data);
    fmpq_t prod;
    fmpq_init(prod);
    for (slong block : selected) {
        slong idx = starts[static_cast<std::size_t>(block)];
        for (slong i = 0; i < n; ++i) {
            for (slong j = 0; j < n; ++j) {
                fmpq_mul(prod, fmpq_mat_entry(U.data, i, idx),
                         fmpq_mat_entry(invU.data, idx, j));
                fmpq_add(fmpq_mat_entry(Q.data, i, j),
                         fmpq_mat_entry(Q.data, i, j), prod);
            }
        }
    }
    fmpq_clear(prod);

    return rational_matrix_from_fmpq_mat(Q.data);
}

RationalMatrix balance(const RationalMatrix& P) {
    std::size_t n = P.rows();
    RationalMatrix out(n, n);
    fmpq_t qval;
    fmpq_init(qval);
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            const auto& pij = P(i, j);
            FmpqPoly num_poly;
            FmpqPoly den_poly;

            fmpq_t delta_minus_p;
            fmpq_init(delta_minus_p);
            if (i == j) fmpq_one(delta_minus_p);
            else        fmpq_zero(delta_minus_p);
            if (!pij.is_zero() && pij.is_constant()) {
                pij.numerator().coeff(0, qval);
                fmpq_t qden;
                fmpq_init(qden);
                pij.denominator().coeff(0, qden);
                fmpq_div(qval, qval, qden);
                fmpq_clear(qden);
                fmpq_sub(delta_minus_p, delta_minus_p, qval);
            } else if (!pij.is_zero()) {
                throw std::invalid_argument("balance: P has non-constant entry");
            }

            if (pij.is_zero()) {
                num_poly.set_fmpq(delta_minus_p);
                den_poly.set_one();
            } else {
                pij.numerator().coeff(0, qval);
                fmpq_t qden;
                fmpq_init(qden);
                pij.denominator().coeff(0, qden);
                fmpq_div(qval, qval, qden);
                fmpq_clear(qden);
                num_poly.set_coeff_fmpq(0, qval);
                num_poly.set_coeff_fmpq(1, delta_minus_p);
                den_poly.set_coeff_si(1, 1);   // eta
            }
            out(i, j) = RationalFunction(std::move(num_poly), std::move(den_poly));
            fmpq_clear(delta_minus_p);
        }
    }
    fmpq_clear(qval);
    return out;
}

RationalMatrix inv_balance(const RationalMatrix& P) {
    std::size_t n = P.rows();
    RationalMatrix out(n, n);
    fmpq_t qval;
    fmpq_init(qval);
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            const auto& pij = P(i, j);
            FmpqPoly num_poly;
            fmpq_t delta_minus_p;
            fmpq_init(delta_minus_p);
            if (i == j) fmpq_one(delta_minus_p);
            else        fmpq_zero(delta_minus_p);
            if (!pij.is_zero() && pij.is_constant()) {
                pij.numerator().coeff(0, qval);
                fmpq_t qden;
                fmpq_init(qden);
                pij.denominator().coeff(0, qden);
                fmpq_div(qval, qval, qden);
                fmpq_clear(qden);
                fmpq_sub(delta_minus_p, delta_minus_p, qval);
            }
            num_poly.set_coeff_fmpq(0, delta_minus_p);
            if (!pij.is_zero()) {
                pij.numerator().coeff(0, qval);
                fmpq_t qden;
                fmpq_init(qden);
                pij.denominator().coeff(0, qden);
                fmpq_div(qval, qval, qden);
                fmpq_clear(qden);
                num_poly.set_coeff_fmpq(1, qval);
            }
            FmpqPoly den_poly;
            den_poly.set_one();
            out(i, j) = RationalFunction(std::move(num_poly), std::move(den_poly));
            fmpq_clear(delta_minus_p);
        }
    }
    fmpq_clear(qval);
    return out;
}

}  // namespace

// ===========================================================================
//  Rational matrix similarity transforms
// ===========================================================================

namespace {

RationalMatrix matmul_rat(const RationalMatrix& A, const RationalMatrix& B) {
    std::size_t n = A.rows();
    std::size_t m = B.cols();
    std::size_t k = A.cols();
    if (k != B.rows()) throw std::invalid_argument("matmul_rat: shape mismatch");
    RationalMatrix out(n, m);
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < m; ++j) {
            RationalFunction acc;
            for (std::size_t l = 0; l < k; ++l) {
                acc += A(i, l) * B(l, j);
            }
            out(i, j) = std::move(acc);
        }
    }
    return out;
}

RationalMatrix matadd_rat(const RationalMatrix& A, const RationalMatrix& B) {
    if (A.rows() != B.rows() || A.cols() != B.cols()) {
        throw std::invalid_argument("matadd_rat: shape mismatch");
    }
    RationalMatrix out(A.rows(), A.cols());
    for (std::size_t i = 0; i < A.rows(); ++i) {
        for (std::size_t j = 0; j < A.cols(); ++j) {
            out(i, j) = A(i, j) + B(i, j);
        }
    }
    return out;
}

RationalMatrix eta_power_times(const RationalMatrix& A, long power) {
    RationalMatrix out(A.rows(), A.cols());
    for (std::size_t i = 0; i < A.rows(); ++i) {
        for (std::size_t j = 0; j < A.cols(); ++j) {
            RationalFunction entry = A(i, j);
            entry.multiply_by_eta_power(power);
            out(i, j) = std::move(entry);
        }
    }
    return out;
}

void fmpq_mat_from_a1_residue(fmpq_mat_t out,
                              const RationalMatrix& B,
                              long p) {
    RationalMatrix a1(B.rows(), B.cols());
    fmpq_t p_plus_one;
    fmpq_init(p_plus_one);
    fmpq_set_si(p_plus_one, p + 1, 1);
    for (std::size_t i = 0; i < B.rows(); ++i) {
        for (std::size_t j = 0; j < B.cols(); ++j) {
            RationalFunction term = B(i, j).derivative();
            term.multiply_by_eta_power(p + 1);
            RationalFunction leading = B(i, j);
            leading.multiply_by_eta_power(p);
            leading.multiply_by_fmpq(p_plus_one);
            term += leading;
            a1(i, j) = std::move(term);
        }
    }
    fmpq_clear(p_plus_one);
    fmpq_mat_from_rational_residue(out, a1, /*shift=*/0);
}

// ---------------------------------------------------------------------------
//  Env-gated kernel-shear rescue for to_fuchsian_local
//  (AMFLOW_TOFUCHSIAN_KERNEL_RESCUE=1).
//
//  Lee/Barkatou right-kernel shear, moser_prepass-verified on a reference
//  12-block operator pair (regular-singular at eta=0 by the cyclic-vector Fuchs
//  test, yet outside the projector loop's 0<->inf balance candidate set).
//  Behavior is unchanged unless the env knob is set AND find_projector_exact
//  throws "irreducible Poincare rank".
// ---------------------------------------------------------------------------

bool kernel_rescue_enabled() {
    const char* e = std::getenv("AMFLOW_TOFUCHSIAN_KERNEL_RESCUE");
    return e != nullptr && *e == '1';
}

slong fmpq_mat_rank_exact(const fmpq_mat_t A) {
    FmpqMatLocal R(fmpq_mat_nrows(A), fmpq_mat_ncols(A));
    return fmpq_mat_rref(R.data, A);
}

// Exact projector onto a complement of ker(A) along ker(A): columns of
// K = nullspace basis of A (from rref), completed to a basis C = [K | W]
// with W = unit vectors on the pivot columns; Q = C * diag(0_k, I_{n-k}) * C^-1.
void kernel_complement_projector(fmpq_mat_t Q, const fmpq_mat_t A) {
    slong n = fmpq_mat_nrows(A);
    FmpqMatLocal R(n, n);
    slong rank = fmpq_mat_rref(R.data, A);

    std::vector<slong> pivots;
    pivots.reserve(static_cast<std::size_t>(rank));
    {
        slong row = 0;
        for (slong col = 0; col < n && row < rank; ++col) {
            if (!fmpq_is_zero(fmpq_mat_entry(R.data, row, col))) {
                pivots.push_back(col);
                ++row;
            }
        }
    }
    std::vector<char> is_pivot(static_cast<std::size_t>(n), 0);
    for (slong pc : pivots) is_pivot[static_cast<std::size_t>(pc)] = 1;

    FmpqMatLocal C(n, n);
    fmpq_mat_zero(C.data);
    slong kdim = 0;  // nullspace columns go first
    for (slong f = 0; f < n; ++f) {
        if (is_pivot[static_cast<std::size_t>(f)]) continue;
        fmpq_one(fmpq_mat_entry(C.data, f, kdim));
        for (slong r = 0; r < rank; ++r) {
            fmpq_neg(fmpq_mat_entry(C.data,
                                    pivots[static_cast<std::size_t>(r)], kdim),
                     fmpq_mat_entry(R.data, r, f));
        }
        ++kdim;
    }
    for (slong r = 0; r < rank; ++r) {
        fmpq_one(fmpq_mat_entry(C.data,
                                pivots[static_cast<std::size_t>(r)], kdim + r));
    }

    FmpqMatLocal Cinv(n, n);
    if (!fmpq_mat_inv(Cinv.data, C.data)) {
        throw std::runtime_error(
            "kernel_complement_projector: completion matrix is singular");
    }
    FmpqMatLocal CW(n, n);
    fmpq_mat_set(CW.data, C.data);
    for (slong i = 0; i < n; ++i) {
        for (slong j = 0; j < kdim; ++j) {
            fmpq_zero(fmpq_mat_entry(CW.data, i, j));
        }
    }
    fmpq_mat_mul(Q, CW.data, Cinv.data);
}

// One trial application of the balance update (identical formula to the
// main loop below) — computes the trial B without committing anything.
RationalMatrix apply_balance_trial(const RationalMatrix& B,
                                   const RationalMatrix& Q) {
    RationalMatrix Bal    = balance(Q);
    RationalMatrix InvBal = inv_balance(Q);
    RationalMatrix inner  = matmul_rat(B, Bal);
    inner = matadd_rat(inner, eta_power_times(Q, -2));
    return matmul_rat(InvBal, inner);
}

// Attempt the kernel-shear rescue: right-kernel variant on A0 first, then
// the left (cokernel) variant built on A0^T (transposed back so it acts as
// a column-space projector).  A candidate is accepted iff the trial update
// strictly lowers the Poincare rank, or keeps it equal while strictly
// lowering rank(A0).  Returns true with the accepted projector in Q_out;
// false if neither variant is accepted (caller rethrows the original error).
bool try_kernel_shear_rescue(RationalMatrix&       Q_out,
                             const RationalMatrix& B,
                             const fmpq_mat_t      A0,
                             long                  p) {
    slong n     = fmpq_mat_nrows(A0);
    slong rank0 = fmpq_mat_rank_exact(A0);

    for (int variant = 0; variant < 2; ++variant) {
        FmpqMatLocal Qf(n, n);
        if (variant == 0) {
            kernel_complement_projector(Qf.data, A0);
        } else {
            FmpqMatLocal At(n, n);
            FmpqMatLocal Qt(n, n);
            fmpq_mat_transpose(At.data, A0);
            kernel_complement_projector(Qt.data, At.data);
            fmpq_mat_transpose(Qf.data, Qt.data);
        }
        RationalMatrix Q      = rational_matrix_from_fmpq_mat(Qf.data);
        RationalMatrix Btrial = apply_balance_trial(B, Q);
        long  pt     = poincare_rank(Btrial);
        bool  accept = false;
        slong rank_t = -1;
        if (pt < p) {
            accept = true;
        } else if (pt == p) {
            FmpqMatLocal A0t(n, n);
            fmpq_mat_from_rational_residue(A0t.data, Btrial, /*shift=*/pt + 1);
            rank_t = fmpq_mat_rank_exact(A0t.data);
            if (rank_t < rank0) accept = true;
        }
        if (accept) {
            const std::string msg =
                std::string("ToFuchsian: TOFUCHSIAN_KERNEL_RESCUE applied ")
                + (variant == 0 ? "right-kernel" : "left-cokernel")
                + " shear (p " + std::to_string(p) + " -> "
                + std::to_string(pt) + ", rank A0 " + std::to_string(rank0)
                + " -> "
                + (rank_t >= 0 ? std::to_string(rank_t) : std::string("n/a"))
                + ")";
            if (!silent_mode()) log_line(msg);
            else                std::fputs((msg + "\n").c_str(), stderr);
            Q_out = std::move(Q);
            return true;
        }
    }
    return false;
}

NormalizationResult to_fuchsian_local(const RationalMatrix& mat) {
    std::size_t n = mat.rows();
    if (mat.cols() != n) {
        throw std::invalid_argument("to_fuchsian_local: matrix must be square");
    }

    NormalizationResult out;
    out.T = RationalMatrix::identity(n);
    out.invT = RationalMatrix::identity(n);
    out.B = mat;

    const long max_iter = std::max<long>(16, 8 * static_cast<long>(n) *
                                                 static_cast<long>(n) + 8);
    for (long iter = 0; iter < max_iter; ++iter) {
        long p = poincare_rank(out.B);
        if (p <= 0) return out;

        fmpq_mat_t A0, A1;
        fmpq_mat_init(A0, static_cast<slong>(n), static_cast<slong>(n));
        fmpq_mat_init(A1, static_cast<slong>(n), static_cast<slong>(n));
        try {
            fmpq_mat_from_rational_residue(A0, out.B, /*shift=*/p + 1);
            fmpq_mat_from_a1_residue(A1, out.B, p);
            RationalMatrix Q;
            try {
                Q = find_projector_exact(A0, A1);
            } catch (const std::runtime_error& e) {
                if (!kernel_rescue_enabled() ||
                    std::string(e.what()).find("irreducible Poincare rank") ==
                        std::string::npos) {
                    throw;
                }
                // RESCUE (env-gated): right-kernel shear where the projector
                // move dead-ends; accepted only on a measured Poincare-rank /
                // rank(A0) drop, else the original error is rethrown.
                if (!try_kernel_shear_rescue(Q, out.B, A0, p)) throw;
            }
            RationalMatrix Bal = balance(Q);
            RationalMatrix InvBal = inv_balance(Q);

            out.T = matmul_rat(out.T, Bal);
            out.invT = matmul_rat(InvBal, out.invT);

            RationalMatrix inner = matmul_rat(out.B, Bal);
            inner = matadd_rat(inner, eta_power_times(Q, -2));
            out.B = matmul_rat(InvBal, inner);
        } catch (...) {
            fmpq_mat_clear(A1);
            fmpq_mat_clear(A0);
            throw;
        }
        fmpq_mat_clear(A1);
        fmpq_mat_clear(A0);
    }

    throw std::runtime_error("to_fuchsian_local: iteration limit reached");
}

}  // namespace

// ===========================================================================
//  Shearing transformation, NormalizeEigen, leading Jordan
// ===========================================================================

namespace {

struct LocalTransform {
    RationalMatrix T, invT, B;
};

long fmpq_floor_si(const fmpq_t q) {
    fmpz_t fq;
    fmpz_init(fq);
    fmpz_fdiv_q(fq, fmpq_numref(q), fmpq_denref(q));
    long out = fmpz_get_si(fq);
    fmpz_clear(fq);
    return out;
}

struct AlgebraicFactorShift {
    FmpqPoly factor;
    long exponent = 1;
    bool selected = false;
    std::vector<long> floors;
};

struct AlgebraicShearingPlan {
    long ele = 0;
    std::vector<long> floors;
    std::vector<AlgebraicFactorShift> factors;
};

void fmpq_poly_make_monic_local(fmpq_poly_t p) {
    if (fmpq_poly_is_zero(p)) return;
    const slong d = fmpq_poly_degree(p);
    if (d < 0) return;
    fmpq_t lead, inv;
    fmpq_init(lead);
    fmpq_init(inv);
    fmpq_poly_get_coeff_fmpq(lead, p, d);
    if (!fmpq_is_zero(lead) && !fmpq_is_one(lead)) {
        fmpq_inv(inv, lead);
        fmpq_poly_scalar_mul_fmpq(p, p, inv);
    }
    fmpq_clear(inv);
    fmpq_clear(lead);
}

std::vector<long> floors_for_factor(const fmpz_poly_t factor) {
    const slong deg = fmpz_poly_degree(factor);
    std::vector<long> floors;
    if (deg <= 0) return floors;

    fmpq_poly_t factor_q;
    fmpq_poly_init(factor_q);
    fmpq_poly_set_fmpz_poly(factor_q, factor);

    auto roots = algebraic_roots_from_fmpq_poly(factor_q);
    fmpq_poly_clear(factor_q);

    floors.reserve(roots.size());
    for (const auto& root : roots) {
        floors.push_back(root.floor_real_si());
    }
    return floors;
}

AlgebraicShearingPlan algebraic_shearing_plan_or_throw(
        const fmpq_mat_t L0,
        long /*prec*/,
        const std::string& reason) {
    fmpq_poly_t p;
    fmpq_poly_init(p);
    fmpq_mat_charpoly(p, L0);

    fmpz_poly_t pz;
    fmpz_poly_init(pz);
    fmpq_poly_get_numerator(pz, p);

    fmpz_poly_factor_t fac;
    fmpz_poly_factor_init(fac);
    fmpz_poly_factor(fac, pz);

    AlgebraicShearingPlan plan;
    long max_floor = std::numeric_limits<long>::min();
    bool any_nonzero = false;

    for (slong i = 0; i < fac->num; ++i) {
        fmpz_poly_struct* fi = fac->p + i;
        const slong deg = fmpz_poly_degree(fi);
        if (deg <= 0) continue;

        AlgebraicFactorShift info;
        fmpq_poly_set_fmpz_poly(info.factor.raw(), fi);
        fmpq_poly_make_monic_local(info.factor.raw());
        info.exponent = static_cast<long>(fac->exp[i]);
        info.floors = floors_for_factor(fi);
        for (long f : info.floors) {
            plan.floors.push_back(f);
            max_floor = std::max(max_floor, f);
            if (f != 0) any_nonzero = true;
        }
        plan.factors.push_back(std::move(info));
    }

    fmpz_poly_factor_clear(fac);
    fmpz_poly_clear(pz);
    fmpq_poly_clear(p);

    if (!any_nonzero) return plan;
    plan.ele = (max_floor > 0) ? 1 : -1;

    for (auto& f : plan.factors) {
        bool have = false;
        bool selected = false;
        for (long floor : f.floors) {
            bool s = (plan.ele > 0) ? (floor > 0) : (floor < 0);
            if (!have) {
                selected = s;
                have = true;
            } else if (selected != s) {
                std::ostringstream oss;
                oss << "algebraic eigenvalue block has conjugate roots with "
                    << "different integer-floor shifts [";
                for (std::size_t i = 0; i < f.floors.size(); ++i) {
                    if (i) oss << ",";
                    oss << f.floors[i];
                }
                oss << "], requiring algebraic-number shearing: " << reason;
                throw std::runtime_error(oss.str());
            }
        }
        f.selected = selected;
    }

    return plan;
}

void fmpq_mat_poly_eval(fmpq_mat_t out,
                        const fmpq_poly_t poly,
                        const fmpq_mat_t A) {
    const slong n = fmpq_mat_nrows(A);
    if (fmpq_mat_ncols(A) != n) {
        throw std::invalid_argument("fmpq_mat_poly_eval: matrix must be square");
    }

    fmpq_mat_zero(out);
    const slong deg = fmpq_poly_degree(poly);
    if (deg < 0) return;

    fmpq_mat_t tmp;
    fmpq_mat_init(tmp, n, n);
    fmpq_t coeff;
    fmpq_init(coeff);

    for (slong k = deg; k >= 0; --k) {
        if (k != deg) {
            fmpq_mat_mul(tmp, out, A);
            fmpq_mat_swap(out, tmp);
        }
        fmpq_poly_get_coeff_fmpq(coeff, poly, k);
        if (!fmpq_is_zero(coeff)) {
            for (slong i = 0; i < n; ++i) {
                fmpq_add(fmpq_mat_entry(out, i, i),
                         fmpq_mat_entry(out, i, i), coeff);
            }
        }
    }

    fmpq_clear(coeff);
    fmpq_mat_clear(tmp);
}

RationalMatrix algebraic_projector_from_plan(
        const fmpq_mat_t L0,
        const AlgebraicShearingPlan& plan) {
    const slong n = fmpq_mat_nrows(L0);

    FmpqPoly selected;
    FmpqPoly unselected;
    selected.set_one();
    unselected.set_one();

    for (const auto& f : plan.factors) {
        FmpqPoly pow;
        fmpq_poly_pow(pow.raw(), f.factor.raw(),
                      static_cast<ulong>(std::max<long>(1, f.exponent)));
        if (f.selected) {
            fmpq_poly_mul(selected.raw(), selected.raw(), pow.raw());
        } else {
            fmpq_poly_mul(unselected.raw(), unselected.raw(), pow.raw());
        }
    }

    fmpq_poly_t gcd, s, t, proj_poly;
    fmpq_poly_init(gcd);
    fmpq_poly_init(s);
    fmpq_poly_init(t);
    fmpq_poly_init(proj_poly);

    fmpq_poly_xgcd(gcd, s, t, selected.raw(), unselected.raw());
    if (fmpq_poly_degree(gcd) != 0) {
        fmpq_poly_clear(proj_poly);
        fmpq_poly_clear(t);
        fmpq_poly_clear(s);
        fmpq_poly_clear(gcd);
        throw std::runtime_error(
            "algebraic_projector_from_plan: selected factors are not coprime");
    }

    fmpq_t g0;
    fmpq_init(g0);
    fmpq_poly_get_coeff_fmpq(g0, gcd, 0);
    if (!fmpq_is_one(g0)) {
        fmpq_poly_scalar_div_fmpq(t, t, g0);
    }
    fmpq_clear(g0);

    fmpq_poly_mul(proj_poly, t, unselected.raw());

    fmpq_mat_t P;
    fmpq_mat_init(P, n, n);
    fmpq_mat_poly_eval(P, proj_poly, L0);
    RationalMatrix out = rational_matrix_from_fmpq_mat(P);

    fmpq_mat_clear(P);
    fmpq_poly_clear(proj_poly);
    fmpq_poly_clear(t);
    fmpq_poly_clear(s);
    fmpq_poly_clear(gcd);
    return out;
}

LocalTransform algebraic_projector_shearing_transformation(
        const RationalMatrix& mat,
        const fmpq_mat_t L0,
        long prec,
        const std::string& reason) {
    AlgebraicShearingPlan plan =
        algebraic_shearing_plan_or_throw(L0, prec, reason);
    if (plan.ele == 0) {
        LocalTransform out;
        out.T = RationalMatrix::identity(mat.rows());
        out.invT = RationalMatrix::identity(mat.rows());
        out.B = mat;
        return out;
    }

    RationalMatrix P = algebraic_projector_from_plan(L0, plan);
    RationalMatrix IminusP = RationalMatrix::identity(mat.rows());
    for (std::size_t i = 0; i < mat.rows(); ++i) {
        for (std::size_t j = 0; j < mat.cols(); ++j) {
            IminusP(i, j) -= P(i, j);
        }
    }

    LocalTransform out;
    out.T = IminusP;
    out.invT = IminusP;
    for (std::size_t i = 0; i < mat.rows(); ++i) {
        for (std::size_t j = 0; j < mat.cols(); ++j) {
            RationalFunction t = P(i, j);
            t.multiply_by_eta_power(plan.ele);
            out.T(i, j) += t;

            RationalFunction inv = P(i, j);
            inv.multiply_by_eta_power(-plan.ele);
            out.invT(i, j) += inv;
        }
    }

    out.B = matmul_rat(matmul_rat(out.invT, mat), out.T);
    fmpq_t ele_q;
    fmpq_init(ele_q);
    fmpq_set_si(ele_q, plan.ele, 1);
    for (std::size_t i = 0; i < mat.rows(); ++i) {
        for (std::size_t j = 0; j < mat.cols(); ++j) {
            RationalFunction corr = P(i, j);
            corr.multiply_by_eta_power(-1);
            corr.multiply_by_fmpq(ele_q);
            out.B(i, j) -= corr;
        }
    }
    fmpq_clear(ele_q);

    return out;
}

BlockBehavior exact_or_numeric_leading_behavior(const RationalMatrix& sub,
                                                long prec) {
    const std::size_t n = sub.rows();
    fmpq_mat_t L0;
    fmpq_mat_init(L0, static_cast<slong>(n), static_cast<slong>(n));
    fmpq_mat_from_rational_residue(L0, sub, /*shift=*/1);

    fmpq_mat_t S, J, Sinv;
    fmpq_mat_init(S, static_cast<slong>(n), static_cast<slong>(n));
    fmpq_mat_init(J, static_cast<slong>(n), static_cast<slong>(n));
    fmpq_mat_init(Sinv, static_cast<slong>(n), static_cast<slong>(n));

    BlockBehavior out;
    std::vector<long> block_sizes;
    try {
        jordan_decomposition_exact(S, J, Sinv, block_sizes, L0);
        slong offset = 0;
        for (long bs : block_sizes) {
            BehaviorEntry e;
            e.mu.set_fmpq(fmpq_mat_entry(J, offset, offset), prec);
            e.log_power = bs - 1;
            out.push_back(std::move(e));
            offset += bs;
        }
    } catch (...) {
        auto eigs = algebraic_eigenvalues_from_fmpq_mat(L0);
        for (auto& ev : eigs) {
            BehaviorEntry e;
            e.mu = ev.value.to_acb(prec);
            e.log_power = std::max<long>(0, ev.multiplicity - 1);
            out.push_back(std::move(e));
        }
    }

    fmpq_mat_clear(Sinv);
    fmpq_mat_clear(J);
    fmpq_mat_clear(S);
    fmpq_mat_clear(L0);
    return out;
}

LocalTransform shearing_transformation(const RationalMatrix& mat,
                                       long /*prec*/) {
    std::size_t n = mat.rows();

    fmpq_mat_t L0;
    fmpq_mat_init(L0, static_cast<slong>(n), static_cast<slong>(n));
    fmpq_mat_from_rational_residue(L0, mat, /*shift=*/1);

    fmpq_mat_t S, J, Sinv;
    fmpq_mat_init(S, static_cast<slong>(n), static_cast<slong>(n));
    fmpq_mat_init(J, static_cast<slong>(n), static_cast<slong>(n));
    fmpq_mat_init(Sinv, static_cast<slong>(n), static_cast<slong>(n));
    std::vector<long> block_sizes;
    try {
        jordan_decomposition_exact(S, J, Sinv, block_sizes, L0);
    } catch (const std::exception& e) {
        try {
            LocalTransform out = algebraic_projector_shearing_transformation(
                mat, L0, working_prec_bits(), e.what());
            fmpq_mat_clear(Sinv);
            fmpq_mat_clear(J);
            fmpq_mat_clear(S);
            fmpq_mat_clear(L0);
            return out;
        } catch (...) {
            fmpq_mat_clear(Sinv);
            fmpq_mat_clear(J);
            fmpq_mat_clear(S);
            fmpq_mat_clear(L0);
            throw;
        }
    }
    fmpq_mat_clear(L0);

    std::vector<long> int_part(n, 0);
    long max_int = std::numeric_limits<long>::min();
    for (std::size_t i = 0; i < n; ++i) {
        long ip = fmpq_floor_si(fmpq_mat_entry(J, i, i));
        int_part[i] = ip;
        if (ip > max_int) max_int = ip;
    }

    std::vector<long> ele(n, 0);
    if (max_int > 0) {
        for (std::size_t k = 0; k < n; ++k) ele[k] = (int_part[k] > 0 ? 1 : 0);
    } else {
        for (std::size_t k = 0; k < n; ++k) ele[k] = (int_part[k] < 0 ? -1 : 0);
    }

    RationalMatrix u_rat    = rational_matrix_from_fmpq_mat(S);
    RationalMatrix invu_rat = rational_matrix_from_fmpq_mat(Sinv);
    fmpq_mat_clear(Sinv);
    fmpq_mat_clear(J);
    fmpq_mat_clear(S);

    LocalTransform out;
    out.T    = RationalMatrix(n, n);
    out.invT = RationalMatrix(n, n);
    out.B    = RationalMatrix(n, n);
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t k = 0; k < n; ++k) {
            RationalFunction tmp = u_rat(i, k);
            tmp.multiply_by_eta_power(ele[k]);
            out.T(i, k) = std::move(tmp);
        }
        for (std::size_t k = 0; k < n; ++k) {
            RationalFunction tmp = invu_rat(i, k);
            tmp.multiply_by_eta_power(-ele[i]);
            out.invT(i, k) = std::move(tmp);
        }
    }

    RationalMatrix prod = matmul_rat(matmul_rat(out.invT, mat), out.T);
    out.B = prod;
    for (std::size_t k = 0; k < n; ++k) {
        if (ele[k] != 0) {
            RationalFunction corr = RationalFunction::monomial(-1);
            fmpq_t e;
            fmpq_init(e);
            fmpq_set_si(e, ele[k], 1);
            corr.multiply_by_fmpq(e);
            fmpq_clear(e);
            out.B(k, k) -= corr;
        }
    }
    return out;
}

bool normalize_eigen_q(const RationalMatrix& mat, long prec) {
    std::size_t n = mat.rows();
    fmpq_mat_t L0;
    fmpq_mat_init(L0, static_cast<slong>(n), static_cast<slong>(n));
    fmpq_mat_from_rational_residue(L0, mat, /*shift=*/1);

    fmpq_mat_t S, J, Sinv;
    fmpq_mat_init(S, static_cast<slong>(n), static_cast<slong>(n));
    fmpq_mat_init(J, static_cast<slong>(n), static_cast<slong>(n));
    fmpq_mat_init(Sinv, static_cast<slong>(n), static_cast<slong>(n));
    std::vector<long> block_sizes;
    bool any = false;
    try {
        jordan_decomposition_exact(S, J, Sinv, block_sizes, L0);
        for (std::size_t i = 0; i < n; ++i) {
            long ip = fmpq_floor_si(fmpq_mat_entry(J, i, i));
            if (ip != 0) {
                any = true;
                break;
            }
        }
    } catch (const std::exception& e) {
        try {
            AlgebraicShearingPlan plan =
                algebraic_shearing_plan_or_throw(L0, prec, e.what());
            fmpq_mat_clear(Sinv);
            fmpq_mat_clear(J);
            fmpq_mat_clear(S);
            fmpq_mat_clear(L0);
            return plan.ele != 0;
        } catch (...) {
            fmpq_mat_clear(Sinv);
            fmpq_mat_clear(J);
            fmpq_mat_clear(S);
            fmpq_mat_clear(L0);
            throw;
        }
    }
    fmpq_mat_clear(Sinv);
    fmpq_mat_clear(J);
    fmpq_mat_clear(S);
    fmpq_mat_clear(L0);
    return any;
}

LocalTransform normalize_eigen(const RationalMatrix& mat, long prec) {
    std::size_t n = mat.rows();
    LocalTransform st;
    st.T    = RationalMatrix::identity(n);
    st.invT = RationalMatrix::identity(n);
    st.B    = mat;
    while (normalize_eigen_q(st.B, prec)) {
        LocalTransform step = shearing_transformation(st.B, prec);
        st.T    = matmul_rat(st.T, step.T);
        st.invT = matmul_rat(step.invT, st.invT);
        st.B    = step.B;
    }
    return st;
}

struct LeadingJordan {
    RationalMatrix u;
    RationalMatrix invu;
    RationalMatrix jor;
};

std::vector<AcbValue> acb_flat_from_fmpq_mat(const fmpq_mat_t mat, long prec) {
    const slong rows = fmpq_mat_nrows(mat);
    const slong cols = fmpq_mat_ncols(mat);
    std::vector<AcbValue> out(static_cast<std::size_t>(rows * cols));
    for (slong i = 0; i < rows; ++i) {
        for (slong j = 0; j < cols; ++j) {
            out[static_cast<std::size_t>(i * cols + j)].set_fmpq(
                fmpq_mat_entry(mat, i, j), prec);
        }
    }
    return out;
}

std::vector<AcbValue> acb_flat_from_ca_mat(const ca_mat_t mat,
                                           ca_ctx_t       ctx,
                                           long           prec) {
    const slong rows = ca_mat_nrows(mat);
    const slong cols = ca_mat_ncols(mat);
    std::vector<AcbValue> out(static_cast<std::size_t>(rows * cols));
    for (slong i = 0; i < rows; ++i) {
        for (slong j = 0; j < cols; ++j) {
            ca_get_acb(out[static_cast<std::size_t>(i * cols + j)].raw(),
                       ca_mat_entry(mat, i, j), prec, ctx);
        }
    }
    return out;
}

std::vector<AcbValue> acb_flat_matmul(const std::vector<AcbValue>& A,
                                      const std::vector<AcbValue>& B,
                                      std::size_t                  n,
                                      long                         prec) {
    std::vector<AcbValue> out(n * n);
    AcbValue tmp;
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            AcbValue acc;
            for (std::size_t k = 0; k < n; ++k) {
                acb_mul(tmp.raw(),
                        A[i * n + k].raw(),
                        B[k * n + j].raw(),
                        prec);
                acb_add(acc.raw(), acc.raw(), tmp.raw(), prec);
            }
            out[i * n + j] = std::move(acc);
        }
    }
    return out;
}

bool acb_flat_matrix_close(const std::vector<AcbValue>& A,
                           const std::vector<AcbValue>& B) {
    if (A.size() != B.size()) return false;
    for (std::size_t i = 0; i < A.size(); ++i) {
        AcbValue diff;
        acb_sub(diff.raw(), A[i].raw(), B[i].raw(), working_prec_bits());
        if (!acb_is_chop_zero(diff.raw(), chop_pre())) return false;
    }
    return true;
}

BlockJordanRotation block_jordan_rotation_for_calc_zero(
        const std::vector<std::size_t>& block,
        const RationalMatrix&           mat,
        long                            prec) {
    const std::size_t n = mat.rows();
    BlockJordanRotation out;
    out.block = block;

    fmpq_mat_t L0;
    fmpq_mat_init(L0, static_cast<slong>(n), static_cast<slong>(n));
    fmpq_mat_from_rational_residue(L0, mat, /*shift=*/1);

    fmpq_mat_t S, J, Sinv;
    fmpq_mat_init(S, static_cast<slong>(n), static_cast<slong>(n));
    fmpq_mat_init(J, static_cast<slong>(n), static_cast<slong>(n));
    fmpq_mat_init(Sinv, static_cast<slong>(n), static_cast<slong>(n));
    std::vector<long> block_sizes;

    try {
        jordan_decomposition_exact(S, J, Sinv, block_sizes, L0);
        out.u = acb_flat_from_fmpq_mat(S, prec);
        out.invu = acb_flat_from_fmpq_mat(Sinv, prec);
        out.jor = acb_flat_from_fmpq_mat(J, prec);
        fmpq_mat_clear(Sinv);
        fmpq_mat_clear(J);
        fmpq_mat_clear(S);
        fmpq_mat_clear(L0);
        return out;
    } catch (const std::exception& rational_error) {
        ca_ctx_t ctx;
        ca_ctx_init(ctx);

        ca_mat_t A, J_ca, P_ca, Pinv_ca;
        ca_mat_init(A, static_cast<slong>(n), static_cast<slong>(n), ctx);
        ca_mat_init(J_ca, static_cast<slong>(n), static_cast<slong>(n), ctx);
        ca_mat_init(P_ca, static_cast<slong>(n), static_cast<slong>(n), ctx);
        ca_mat_init(Pinv_ca, static_cast<slong>(n), static_cast<slong>(n), ctx);
        ca_mat_zero(A, ctx);
        ca_mat_zero(J_ca, ctx);
        ca_mat_zero(P_ca, ctx);
        ca_mat_zero(Pinv_ca, ctx);

        ca_mat_set_fmpq_mat(A, L0, ctx);
        const int jordan_status = ca_mat_jordan_form(J_ca, P_ca, A, ctx);
        const truth_t inv_status = ca_mat_inv(Pinv_ca, P_ca, ctx);

        bool ok = (inv_status == T_TRUE);
        if (ok) {
            out.u = acb_flat_from_ca_mat(P_ca, ctx, prec);
            out.invu = acb_flat_from_ca_mat(Pinv_ca, ctx, prec);
            out.jor = acb_flat_from_ca_mat(J_ca, ctx, prec);

            const auto A_flat = acb_flat_from_fmpq_mat(L0, prec);
            const auto lhs = acb_flat_matmul(
                acb_flat_matmul(out.invu, A_flat, n, prec), out.u, n, prec);
            if (!acb_flat_matrix_close(lhs, out.jor)) {
                const auto alt = acb_flat_matmul(
                    acb_flat_matmul(out.u, A_flat, n, prec), out.invu, n, prec);
                if (acb_flat_matrix_close(alt, out.jor)) {
                    std::swap(out.u, out.invu);
                } else {
                    ok = false;
                }
            }
        }

        ca_mat_clear(Pinv_ca, ctx);
        ca_mat_clear(P_ca, ctx);
        ca_mat_clear(J_ca, ctx);
        ca_mat_clear(A, ctx);
        ca_ctx_clear(ctx);

        fmpq_mat_clear(Sinv);
        fmpq_mat_clear(J);
        fmpq_mat_clear(S);
        fmpq_mat_clear(L0);

        if (!ok) {
            std::ostringstream msg;
            msg << "block_jordan_rotation_for_calc_zero: algebraic Jordan "
                   "rotation failed for block "
                << block_to_string(block)
                << " after exact-Q failure: " << rational_error.what()
                << "; ca_mat_jordan_form status=" << jordan_status
                << "; ca_mat_inv status=" << static_cast<int>(inv_status);
            throw std::runtime_error(msg.str());
        }

        return out;
    }
}

LeadingJordan leading_jordan(const RationalMatrix& mat, long /*prec*/) {
    std::size_t n = mat.rows();
    fmpq_mat_t L0;
    fmpq_mat_init(L0, static_cast<slong>(n), static_cast<slong>(n));
    fmpq_mat_from_rational_residue(L0, mat, /*shift=*/1);

    fmpq_mat_t S, J, Sinv;
    fmpq_mat_init(S, static_cast<slong>(n), static_cast<slong>(n));
    fmpq_mat_init(J, static_cast<slong>(n), static_cast<slong>(n));
    fmpq_mat_init(Sinv, static_cast<slong>(n), static_cast<slong>(n));
    std::vector<long> block_sizes;

    LeadingJordan out;
    try {
        jordan_decomposition_exact(S, J, Sinv, block_sizes, L0);
        out.u = rational_matrix_from_fmpq_mat(S);
        out.invu = rational_matrix_from_fmpq_mat(Sinv);
        out.jor = rational_matrix_from_fmpq_mat(J);
    } catch (const std::exception& e) {
        fmpq_mat_clear(Sinv);
        fmpq_mat_clear(J);
        fmpq_mat_clear(S);
        fmpq_mat_clear(L0);
        throw std::runtime_error(
            "leading_jordan: algebraic leading residue requires "
            "algebraic Jordan rotation: " + std::string(e.what()));
    }
    fmpq_mat_clear(Sinv);
    fmpq_mat_clear(J);
    fmpq_mat_clear(S);
    fmpq_mat_clear(L0);
    return out;
}

NormalizationResult normalize_diagonal(const RationalMatrix& mat, long prec) {
    auto blocks = analyze_block(mat);
    NormalizationResult out;
    out.T    = RationalMatrix::zeros(mat.rows(), mat.cols());
    out.invT = RationalMatrix::zeros(mat.rows(), mat.cols());
    out.B    = RationalMatrix::zeros(mat.rows(), mat.cols());

    for (const auto& blk : blocks) {
        auto sub = mat.submatrix(blk, blk);
        LocalTransform tf;
        try {
            auto fuchs = to_fuchsian_local(sub);
            tf.T = std::move(fuchs.T);
            tf.invT = std::move(fuchs.invT);
            tf.B = std::move(fuchs.B);
        } catch (const std::exception& e) {
            throw std::runtime_error(
                "normalize_diagonal: ToFuchsian failed for block "
                + block_to_string(blk) + ": " + e.what());
        }

        LocalTransform ne;
        try {
            ne = normalize_eigen(tf.B, prec);
        } catch (const std::exception& e) {
            throw std::runtime_error(
                "normalize_diagonal: NormalizeEigen failed for block "
                + block_to_string(blk) + ": " + e.what());
        }

        LeadingJordan lj;
        try {
            lj = leading_jordan(ne.B, prec);
        } catch (const std::exception& e) {
            throw std::runtime_error(
                "normalize_diagonal: leading Jordan rotation failed for block "
                + block_to_string(blk) + ": " + e.what());
        }
        RationalMatrix T_block    = matmul_rat(matmul_rat(tf.T, ne.T), lj.u);
        RationalMatrix invT_block = matmul_rat(lj.invu,
                                               matmul_rat(ne.invT, tf.invT));

        RationalMatrix B_block = matmul_rat(matmul_rat(lj.invu, ne.B), lj.u);

        for (std::size_t a = 0; a < blk.size(); ++a) {
            for (std::size_t b = 0; b < blk.size(); ++b) {
                out.T(blk[a], blk[b])    = T_block(a, b);
                out.invT(blk[a], blk[b]) = invT_block(a, b);
                out.B(blk[a], blk[b])    = B_block(a, b);
            }
        }
    }

    for (std::size_t i = 0; i < blocks.size(); ++i) {
        for (std::size_t j = 0; j < i; ++j) {
            const auto& rows = blocks[i];
            const auto& cols = blocks[j];
            RationalMatrix invT_diag = out.invT.submatrix(rows, rows);
            RationalMatrix T_diag    = out.T.submatrix(cols, cols);
            RationalMatrix mat_off   = mat.submatrix(rows, cols);
            RationalMatrix block_out = matmul_rat(matmul_rat(invT_diag, mat_off), T_diag);
            for (std::size_t a = 0; a < rows.size(); ++a) {
                for (std::size_t b = 0; b < cols.size(); ++b) {
                    out.B(rows[a], cols[b]) = block_out(a, b);
                }
            }
        }
    }

    return out;
}

struct CalcZeroDiagonalNormalization {
    NormalizationResult              base;
    std::vector<BlockJordanRotation> rotations;
};

CalcZeroDiagonalNormalization
normalize_diagonal_for_calc_zero(const RationalMatrix& mat, long prec) {
    auto blocks = analyze_block(mat);
    CalcZeroDiagonalNormalization out;
    out.base.T    = RationalMatrix::zeros(mat.rows(), mat.cols());
    out.base.invT = RationalMatrix::zeros(mat.rows(), mat.cols());
    out.base.B    = RationalMatrix::zeros(mat.rows(), mat.cols());
    out.rotations.reserve(blocks.size());

    for (const auto& blk : blocks) {
        auto sub = mat.submatrix(blk, blk);
        LocalTransform tf;
        try {
            auto fuchs = to_fuchsian_local(sub);
            tf.T = std::move(fuchs.T);
            tf.invT = std::move(fuchs.invT);
            tf.B = std::move(fuchs.B);
        } catch (const std::exception& e) {
            throw std::runtime_error(
                "normalize_diagonal_for_calc_zero: ToFuchsian failed for block "
                + block_to_string(blk) + ": " + e.what());
        }

        LocalTransform ne;
        try {
            ne = normalize_eigen(tf.B, prec);
        } catch (const std::exception& e) {
            const std::string msg = e.what();
            if (msg.find("algebraic-number shearing") != std::string::npos) {
                ne.T = RationalMatrix::identity(tf.B.rows());
                ne.invT = RationalMatrix::identity(tf.B.rows());
                ne.B = tf.B;
            } else {
            throw std::runtime_error(
                "normalize_diagonal_for_calc_zero: NormalizeEigen failed for block "
                + block_to_string(blk) + ": " + e.what());
            }
        }

        BlockJordanRotation rot;
        try {
            rot = block_jordan_rotation_for_calc_zero(blk, ne.B, prec);
        } catch (const std::exception& e) {
            throw std::runtime_error(
                "normalize_diagonal_for_calc_zero: leading Jordan rotation failed for block "
                + block_to_string(blk) + ": " + e.what());
        }

        RationalMatrix T_block    = matmul_rat(tf.T, ne.T);
        RationalMatrix invT_block = matmul_rat(ne.invT, tf.invT);
        RationalMatrix B_block    = ne.B;

        for (std::size_t a = 0; a < blk.size(); ++a) {
            for (std::size_t b = 0; b < blk.size(); ++b) {
                out.base.T(blk[a], blk[b])    = T_block(a, b);
                out.base.invT(blk[a], blk[b]) = invT_block(a, b);
                out.base.B(blk[a], blk[b])    = B_block(a, b);
            }
        }
        out.rotations.push_back(std::move(rot));
    }

    for (std::size_t i = 0; i < blocks.size(); ++i) {
        for (std::size_t j = 0; j < i; ++j) {
            const auto& rows = blocks[i];
            const auto& cols = blocks[j];
            RationalMatrix invT_diag = out.base.invT.submatrix(rows, rows);
            RationalMatrix T_diag    = out.base.T.submatrix(cols, cols);
            RationalMatrix mat_off   = mat.submatrix(rows, cols);
            RationalMatrix block_out =
                matmul_rat(matmul_rat(invT_diag, mat_off), T_diag);
            for (std::size_t a = 0; a < rows.size(); ++a) {
                for (std::size_t b = 0; b < cols.size(); ++b) {
                    out.base.B(rows[a], cols[b]) = block_out(a, b);
                }
            }
        }
    }

    return out;
}

NormalizationResult to_fuchsian_global(const RationalMatrix& B_in, long /*prec*/) {
    auto blocks = analyze_block(B_in);
    std::size_t Nblk = blocks.size();
    std::size_t total = B_in.rows();

    NormalizationResult out;
    out.T    = RationalMatrix::identity(total);
    out.invT = RationalMatrix::identity(total);
    out.B    = B_in;

    std::vector<fmpq_mat_struct> diag(Nblk);
    for (std::size_t i = 0; i < Nblk; ++i) {
        std::size_t Nb = blocks[i].size();
        fmpq_mat_init(&diag[i], static_cast<slong>(Nb),
                      static_cast<slong>(Nb));
    }
    struct DiagGuard {
        std::vector<fmpq_mat_struct>& v;
        ~DiagGuard() { for (auto& m : v) fmpq_mat_clear(&m); }
    } diag_guard{diag};
    for (std::size_t i = 0; i < Nblk; ++i) {
        RationalMatrix sub = out.B.submatrix(blocks[i], blocks[i]);
        fmpq_mat_from_rational_residue(&diag[i], sub, /*shift=*/1);
    }

    for (std::size_t i = 0; i < Nblk; ++i) {
        for (long jj = static_cast<long>(i) - 1; jj >= 0; --jj) {
            std::size_t j = jj;
            for (;;) {
                const auto& rows = blocks[i];
                const auto& cols = blocks[j];
                RationalMatrix off = out.B.submatrix(rows, cols);
                long p = poincare_rank(off);
                if (p <= 0) break;

                std::size_t Nr = rows.size();
                std::size_t Nc = cols.size();
                fmpq_mat_t b0_q;
                fmpq_mat_init(b0_q, static_cast<slong>(Nr),
                              static_cast<slong>(Nc));
                try {
                    fmpq_mat_from_rational_residue(b0_q, off,
                                                   /*shift=*/p + 1);
                } catch (...) {
                    fmpq_mat_clear(b0_q);
                    throw;
                }

                fmpq_mat_t G_q;
                fmpq_mat_init(G_q, static_cast<slong>(Nr),
                              static_cast<slong>(Nc));
                if (!internal::solve_off_diagonal_fmpq(G_q, &diag[j], b0_q,
                                                       &diag[i], p)) {
                    fmpq_mat_clear(b0_q);
                    fmpq_mat_clear(G_q);
                    throw std::runtime_error(
                        "to_fuchsian_global: exact Sylvester solver failed "
                        "(singular)");
                }

                RationalMatrix G = rational_matrix_from_fmpq_mat(G_q);
                fmpq_mat_clear(b0_q);
                fmpq_mat_clear(G_q);

                RationalFunction eta_neg_p = RationalFunction::monomial(-p);
                RationalFunction eta_neg_p1 = RationalFunction::monomial(-(p + 1));

                auto G_over = [&](std::size_t a, std::size_t b) {
                    return G(a, b) * eta_neg_p;
                };
                auto G_over_p1 = [&](std::size_t a, std::size_t b) {
                    return G(a, b) * eta_neg_p1;
                };

                for (std::size_t l = 0; l <= j; ++l) {
                    const auto& l_cols = blocks[l];
                    for (std::size_t a = 0; a < rows.size(); ++a) {
                        for (std::size_t b = 0; b < l_cols.size(); ++b) {
                            RationalFunction acc;
                            for (std::size_t k = 0; k < cols.size(); ++k) {
                                acc += G_over(a, k) * out.B(blocks[j][k], l_cols[b]);
                            }
                            out.B(rows[a], l_cols[b]) -= acc;
                        }
                    }
                }
                for (std::size_t l = i; l < Nblk; ++l) {
                    const auto& l_rows = blocks[l];
                    for (std::size_t a = 0; a < l_rows.size(); ++a) {
                        for (std::size_t b = 0; b < cols.size(); ++b) {
                            RationalFunction acc;
                            for (std::size_t k = 0; k < rows.size(); ++k) {
                                acc += out.B(l_rows[a], rows[k]) * G_over(k, b);
                            }
                            out.B(l_rows[a], cols[b]) += acc;
                        }
                    }
                }
                for (std::size_t a = 0; a < rows.size(); ++a) {
                    for (std::size_t b = 0; b < cols.size(); ++b) {
                        RationalFunction extra = G_over_p1(a, b);
                        fmpq_t pp;
                        fmpq_init(pp);
                        fmpq_set_si(pp, p, 1);
                        extra.multiply_by_fmpq(pp);
                        fmpq_clear(pp);
                        out.B(rows[a], cols[b]) += extra;
                    }
                }
                for (std::size_t l = i; l < Nblk; ++l) {
                    const auto& l_rows = blocks[l];
                    for (std::size_t a = 0; a < l_rows.size(); ++a) {
                        for (std::size_t b = 0; b < cols.size(); ++b) {
                            RationalFunction acc;
                            for (std::size_t k = 0; k < rows.size(); ++k) {
                                acc += out.T(l_rows[a], rows[k]) * G_over(k, b);
                            }
                            out.T(l_rows[a], cols[b]) += acc;
                        }
                    }
                }
                for (std::size_t l = 0; l <= j; ++l) {
                    const auto& l_cols = blocks[l];
                    for (std::size_t a = 0; a < rows.size(); ++a) {
                        for (std::size_t b = 0; b < l_cols.size(); ++b) {
                            RationalFunction acc;
                            for (std::size_t k = 0; k < cols.size(); ++k) {
                                acc += G_over(a, k) * out.invT(blocks[j][k], l_cols[b]);
                            }
                            out.invT(rows[a], l_cols[b]) -= acc;
                        }
                    }
                }
            }
        }
    }

    return out;
}

}  // namespace

// ===========================================================================
//  NormalizeMat (top level)
// ===========================================================================

NormalizationResult normalize_mat(const RationalMatrix& mat, long prec) {
    if (!silent_mode()) log_line("NormalizeMat: starting");
    auto step1 = normalize_diagonal(mat, prec);
    auto step2 = to_fuchsian_global(step1.B, prec);

    NormalizationResult out;
    out.B = step2.B;

    auto blocks = analyze_block(mat);
    out.T    = step2.T;
    out.invT = step2.invT;
    for (std::size_t i = 0; i < blocks.size(); ++i) {
        for (std::size_t j = 0; j <= i; ++j) {
            const auto& rows = blocks[i];
            const auto& cols = blocks[j];
            RationalMatrix t1 = step1.T.submatrix(rows, rows);
            RationalMatrix t2 = step2.T.submatrix(rows, cols);
            RationalMatrix prod = matmul_rat(t1, t2);
            for (std::size_t a = 0; a < rows.size(); ++a)
                for (std::size_t b = 0; b < cols.size(); ++b)
                    out.T(rows[a], cols[b]) = prod(a, b);

            RationalMatrix inv1 = step1.invT.submatrix(cols, cols);
            RationalMatrix inv2 = step2.invT.submatrix(rows, cols);
            RationalMatrix invprod = matmul_rat(inv2, inv1);
            for (std::size_t a = 0; a < rows.size(); ++a)
                for (std::size_t b = 0; b < cols.size(); ++b)
                    out.invT(rows[a], cols[b]) = invprod(a, b);
        }
    }

    if (!silent_mode()) log_line("NormalizeMat: finished");
    return out;
}

CalcZeroNormalization normalize_mat_for_calc_zero(const RationalMatrix& mat,
                                                  long                  prec) {
    if (!silent_mode()) log_line("NormalizeMatForCalcZero: starting");
    auto step1 = normalize_diagonal_for_calc_zero(mat, prec);
    auto step2 = to_fuchsian_global(step1.base.B, prec);

    CalcZeroNormalization out;
    out.base.B = step2.B;
    out.base.T = step2.T;
    out.base.invT = step2.invT;
    out.rotations = std::move(step1.rotations);

    auto blocks = analyze_block(mat);
    for (std::size_t i = 0; i < blocks.size(); ++i) {
        for (std::size_t j = 0; j <= i; ++j) {
            const auto& rows = blocks[i];
            const auto& cols = blocks[j];
            RationalMatrix t1 = step1.base.T.submatrix(rows, rows);
            RationalMatrix t2 = step2.T.submatrix(rows, cols);
            RationalMatrix prod = matmul_rat(t1, t2);
            for (std::size_t a = 0; a < rows.size(); ++a)
                for (std::size_t b = 0; b < cols.size(); ++b)
                    out.base.T(rows[a], cols[b]) = prod(a, b);

            RationalMatrix inv1 = step1.base.invT.submatrix(cols, cols);
            RationalMatrix inv2 = step2.invT.submatrix(rows, cols);
            RationalMatrix invprod = matmul_rat(inv2, inv1);
            for (std::size_t a = 0; a < rows.size(); ++a)
                for (std::size_t b = 0; b < cols.size(); ++b)
                    out.base.invT(rows[a], cols[b]) = invprod(a, b);
        }
    }

    if (!silent_mode()) log_line("NormalizeMatForCalcZero: finished");
    return out;
}

}  // namespace amflow::ode
