// SPDX-License-Identifier: MIT
// ode::zero — implementation.
//

#include "amflow/ode/zero.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <sstream>
#include <stdexcept>
#include <utility>
#include <vector>

#include <flint/acb.h>
#include <flint/acb_mat.h>
#include <flint/acb_poly.h>
#include <flint/arb.h>
#include <flint/ca.h>
#include <flint/ca_mat.h>
#include <flint/fmpq.h>
#include <flint/fmpq_mat.h>

#include "amflow/numeric/options.hpp"
#include "amflow/ode/acb_rational.hpp"
#include "amflow/ode/asy.hpp"
#include "amflow/ode/blocks.hpp"
#include "amflow/ode/jordan.hpp"
#include "amflow/ode/normalize.hpp"

namespace amflow::ode {

using numeric::AcbValue;
using numeric::FmpqPoly;
using numeric::RationalFunction;
using numeric::RationalMatrix;
using numeric::acb_is_chop_zero;
using numeric::chop_pre;
using numeric::decimal_digits_to_bits;
using numeric::learn_x_order;
using numeric::log_line;
using numeric::silent_mode;
using numeric::test_x_order;
using numeric::working_pre;
using numeric::working_prec_bits;
using numeric::x_order;

// Forward decl of the impl with optional rotations.
std::vector<AsyExpansion>
calcx00_impl(const std::vector<BlockEquation>&            nheq,
             const std::vector<BlockEquationNum>&         nheqn,
             const std::vector<AcbValue>&                 bc,
             acb_srcptr                                   x0,
             const AsymptoticBehaviorList&                behavior,
             const std::vector<BlockJordanRotation>*      rotations,
             long                                         prec);

// ===========================================================================
//  Local helpers
// ===========================================================================

namespace {

inline AcbValue zero_acb() { AcbValue v; return v; }

std::string truncate_for_error(std::string s) {
    constexpr std::size_t kMax = 1600;
    if (s.size() <= kMax) return s;
    s.resize(kMax);
    s += "...<truncated>";
    return s;
}

std::string acb_to_error_string(acb_srcptr x) {
    AcbValue v;
    acb_set(v.raw(), x);
    return v.to_string(40);
}

AcbValue evaluate_with_context(const RationalFunction& rf,
                               acb_srcptr x,
                               long prec,
                               const std::string& where) {
    try {
        return rf.evaluate(x, prec);
    } catch (const std::exception& e) {
        throw std::runtime_error(
            where + ": " + e.what()
            + "; x=" + acb_to_error_string(x)
            + "; rf=" + truncate_for_error(rf.to_string()));
    }
}

std::vector<AcbValue> clone_vec(const std::vector<AcbValue>& v) {
    std::vector<AcbValue> out;
    out.reserve(v.size());
    for (const auto& x : v) out.push_back(x.clone());
    return out;
}

bool acb_close(acb_srcptr a, acb_srcptr b, long prec) {
    // Used to match a stored expansion's exponent mu against the behavior
    // mu -> exponent-equality test -> resonance tolerance, so an
    // eps-deformed sector is not aliased onto a neighbouring one.
    AcbValue d;
    acb_sub(d.raw(), a, b, prec);
    return acb_is_chop_zero(d.raw(), numeric::resonance_chop_pre());
}

BlockBehavior union_behavior_local(const BlockBehavior& in, long prec);
BlockBehavior union_behavior2_local(const BlockBehavior& in, long prec);

double acb_mid_abs(acb_srcptr z) {
    arf_t mr, mi;
    arf_init(mr);
    arf_init(mi);
    arf_abs(mr, arb_midref(acb_realref(z)));
    arf_abs(mi, arb_midref(acb_imagref(z)));
    const double out = arf_get_d(mr, ARF_RND_UP)
                     + arf_get_d(mi, ARF_RND_UP);
    arf_clear(mr);
    arf_clear(mi);
    return out;
}

bool acb_close_to_zero(acb_srcptr z) {
    return acb_is_chop_zero(z, chop_pre());
}

long acb_poly_valuation_chop(const acb_poly_t poly) {
    const slong len = acb_poly_length(poly);
    for (slong i = 0; i < len; ++i) {
        AcbValue c;
        acb_poly_get_coeff_acb(c.raw(), poly, i);
        if (!acb_close_to_zero(c.raw())) return static_cast<long>(i);
    }
    return -1;
}

std::vector<AcbValue> acb_poly_coeffs(const acb_poly_t poly,
                                      long            len,
                                      long            prec) {
    std::vector<AcbValue> out;
    out.reserve(static_cast<std::size_t>(std::max<long>(0, len)));
    for (long i = 0; i < len; ++i) {
        AcbValue c;
        acb_poly_get_coeff_acb(c.raw(), poly, i);
        acb_set_round(c.raw(), c.raw(), prec);
        out.push_back(std::move(c));
    }
    return out;
}

std::vector<AcbValue> acb_rational_taylor_coeffs(
        AcbRationalFunction r,
        long                len,
        long                prec) {
    std::vector<AcbValue> out(static_cast<std::size_t>(std::max<long>(0, len)));
    if (len <= 0 || r.is_zero()) return out;
    r.strip_eta_prefix();

    const long vn = acb_poly_valuation_chop(r.num());
    const long vd = acb_poly_valuation_chop(r.den());
    if (vn < 0) return out;
    if (vd < 0) {
        throw std::runtime_error("acb_rational_taylor_coeffs: zero denominator");
    }
    if (vd != 0) {
        throw std::runtime_error(
            "acb_rational_taylor_coeffs: input has a pole at eta=0");
    }

    auto num = acb_poly_coeffs(r.num(), len, prec);
    auto den = acb_poly_coeffs(r.den(), len, prec);
    if (den.empty() || acb_close_to_zero(den[0].raw())) {
        std::ostringstream oss;
        oss << "acb_rational_taylor_coeffs: denominator still has zero "
            << "constant term after eta-prefix stripping"
            << " (vn=" << vn << ", vd=" << vd
            << ")";
        throw std::runtime_error(oss.str());
    }
    std::vector<AcbValue> unit(static_cast<std::size_t>(len));
    unit[0].set_one();
    return rational_expansion(num, den, unit, prec);
}

PowerSeries acb_rational_to_power_series(AcbRationalFunction r,
                                         long                prec) {
    PowerSeries ps;
    if (r.is_zero()) return ps;
    r.strip_eta_prefix();

    const long vn = acb_poly_valuation_chop(r.num());
    const long vd = acb_poly_valuation_chop(r.den());
    if (vn < 0) return ps;
    if (vd < 0) {
        throw std::runtime_error("acb_rational_to_power_series: zero denominator");
    }
    const long offset = vn - vd;
    r.multiply_by_eta_power(-offset);
    r.strip_eta_prefix();

    const long len = static_cast<long>(x_order()) + 1;
    auto num = acb_poly_coeffs(r.num(), len, prec);
    auto den = acb_poly_coeffs(r.den(), len, prec);
    if (den.empty() || acb_close_to_zero(den[0].raw())) {
        std::ostringstream oss;
        oss << "acb_rational_to_power_series: denominator still has zero "
            << "constant term after eta-prefix stripping"
            << " (vn=" << vn << ", vd=" << vd
            << ", offset=" << offset << ")";
        throw std::runtime_error(oss.str());
    }
    std::vector<AcbValue> unit(static_cast<std::size_t>(len));
    unit[0].set_one();
    ps.leading_offset = offset;
    ps.coeffs = rational_expansion(num, den, unit, prec);
    return ps;
}

std::vector<std::vector<PowerSeries>>
to_power_series_matrix_acb(const AcbRationalMatrix& m, long prec) {
    std::vector<std::vector<PowerSeries>> out(m.rows());
    for (std::size_t i = 0; i < m.rows(); ++i) {
        out[i].reserve(m.cols());
        for (std::size_t j = 0; j < m.cols(); ++j) {
            out[i].push_back(acb_rational_to_power_series(m(i, j), prec));
        }
    }
    return out;
}

AcbValue evaluate_acb_rational(const AcbRationalFunction& f,
                               acb_srcptr                 x,
                               long                       prec) {
    AcbValue num;
    AcbValue den;
    acb_poly_evaluate(num.raw(), f.num(), x, prec);
    acb_poly_evaluate(den.raw(), f.den(), x, prec);
    if (acb_close_to_zero(den.raw())) {
        throw std::runtime_error("evaluate_acb_rational: denominator vanished");
    }
    AcbValue out;
    acb_div(out.raw(), num.raw(), den.raw(), prec);
    return out;
}

std::vector<AcbValue> acb_flat_from_ca_mat_local(ca_mat_t mat,
                                                 ca_ctx_t ctx,
                                                 long     prec) {
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

std::vector<AcbValue> acb_flat_from_fmpq_mat_local(const fmpq_mat_t mat,
                                                   long             prec) {
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

std::vector<AcbValue> acb_flat_matmul_local(const std::vector<AcbValue>& A,
                                            const std::vector<AcbValue>& B,
                                            std::size_t                  n,
                                            long                         prec) {
    std::vector<AcbValue> out(n * n);
    AcbValue tmp;
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            AcbValue acc;
            for (std::size_t k = 0; k < n; ++k) {
                acb_mul(tmp.raw(), A[i * n + k].raw(),
                        B[k * n + j].raw(), prec);
                acb_add(acc.raw(), acc.raw(), tmp.raw(), prec);
            }
            out[i * n + j] = std::move(acc);
        }
    }
    return out;
}

bool acb_flat_close_local(const std::vector<AcbValue>& A,
                          const std::vector<AcbValue>& B) {
    if (A.size() != B.size()) return false;
    for (std::size_t i = 0; i < A.size(); ++i) {
        AcbValue diff;
        acb_sub(diff.raw(), A[i].raw(), B[i].raw(), working_prec_bits());
        if (!acb_close_to_zero(diff.raw())) return false;
    }
    return true;
}

struct AcbJordanData {
    std::vector<AcbValue> u;
    std::vector<AcbValue> invu;
    std::vector<AcbValue> jor;
};

AcbJordanData acb_jordan_data_from_l0(const fmpq_mat_t L0, long prec) {
    const std::size_t n = static_cast<std::size_t>(fmpq_mat_nrows(L0));
    ca_ctx_t ctx;
    ca_ctx_init(ctx);

    ca_mat_t A, J_ca, P_ca, Pinv_ca;
    ca_mat_init(A, static_cast<slong>(n), static_cast<slong>(n), ctx);
    ca_mat_init(J_ca, static_cast<slong>(n), static_cast<slong>(n), ctx);
    ca_mat_init(P_ca, static_cast<slong>(n), static_cast<slong>(n), ctx);
    ca_mat_init(Pinv_ca, static_cast<slong>(n), static_cast<slong>(n), ctx);
    ca_mat_set_fmpq_mat(A, L0, ctx);
    const int jordan_status = ca_mat_jordan_form(J_ca, P_ca, A, ctx);
    const truth_t inv_status = ca_mat_inv(Pinv_ca, P_ca, ctx);
    if (jordan_status == 0 || inv_status != T_TRUE) {
        ca_mat_clear(Pinv_ca, ctx);
        ca_mat_clear(P_ca, ctx);
        ca_mat_clear(J_ca, ctx);
        ca_mat_clear(A, ctx);
        ca_ctx_clear(ctx);
        throw std::runtime_error(
            "acb_jordan_data_from_l0: calcium Jordan decomposition failed");
    }

    AcbJordanData out;
    out.u = acb_flat_from_ca_mat_local(P_ca, ctx, prec);
    out.invu = acb_flat_from_ca_mat_local(Pinv_ca, ctx, prec);
    out.jor = acb_flat_from_ca_mat_local(J_ca, ctx, prec);

    const auto A_flat = acb_flat_from_fmpq_mat_local(L0, prec);
    const auto lhs = acb_flat_matmul_local(
        acb_flat_matmul_local(out.invu, A_flat, n, prec), out.u, n, prec);
    if (!acb_flat_close_local(lhs, out.jor)) {
        const auto alt = acb_flat_matmul_local(
            acb_flat_matmul_local(out.u, A_flat, n, prec), out.invu, n, prec);
        if (acb_flat_close_local(alt, out.jor)) {
            std::swap(out.u, out.invu);
        } else {
            ca_mat_clear(Pinv_ca, ctx);
            ca_mat_clear(P_ca, ctx);
            ca_mat_clear(J_ca, ctx);
            ca_mat_clear(A, ctx);
            ca_ctx_clear(ctx);
            throw std::runtime_error(
                "acb_jordan_data_from_l0: calcium Jordan orientation check failed");
        }
    }

    ca_mat_clear(Pinv_ca, ctx);
    ca_mat_clear(P_ca, ctx);
    ca_mat_clear(J_ca, ctx);
    ca_mat_clear(A, ctx);
    ca_ctx_clear(ctx);
    return out;
}

long acb_floor_real_mid(acb_srcptr z) {
    return static_cast<long>(std::floor(
        arf_get_d(arb_midref(acb_realref(z)), ARF_RND_DOWN)));
}

std::vector<AsyExpansion>
calc_zero_acb_shearing_single_block(const NormalizationResult& norm,
                                    const std::vector<AcbValue>& bc,
                                    acb_srcptr x0,
                                    long prec) {
    const std::size_t n = norm.B.rows();
    if (n == 0 || norm.B.cols() != n || bc.size() != n) {
        throw std::runtime_error(
            "calc_zero_acb_shearing_single_block: inconsistent dimensions");
    }

    BlockPartition partition = build_partition(norm.B);
    AcbRationalMatrix T_alg(n, n);
    AcbRationalMatrix invT_alg(n, n);
    std::vector<long> global_floors(n, 0);
    std::vector<BlockBehavior> essential_behavior(partition.size());
    AsymptoticBehaviorList behavior(partition.size());
    bool any_shift = false;

    for (std::size_t bi = 0; bi < partition.size(); ++bi) {
        const auto& blk = partition.block(bi);
        const std::size_t nb = blk.size();
        RationalMatrix sub = norm.B.submatrix(blk, blk);

        fmpq_mat_t L0;
        fmpq_mat_init(L0, static_cast<slong>(nb), static_cast<slong>(nb));
        fmpq_mat_from_rational_residue(L0, sub, /*shift=*/1);
        AcbJordanData jd = acb_jordan_data_from_l0(L0, prec);
        fmpq_mat_clear(L0);

        std::vector<long> floors(nb, 0);
        for (std::size_t i = 0; i < nb; ++i) {
            floors[i] = acb_floor_real_mid(jd.jor[i * nb + i].raw());
            global_floors[blk[i]] = floors[i];
            if (floors[i] != 0) any_shift = true;
        }

        for (std::size_t a = 0; a < nb; ++a) {
            for (std::size_t b = 0; b < nb; ++b) {
                T_alg(blk[a], blk[b]) = AcbRationalFunction::monomial(
                    floors[b], jd.u[a * nb + b].raw());
                invT_alg(blk[a], blk[b]) = AcbRationalFunction::monomial(
                    -floors[a], jd.invu[a * nb + b].raw());
            }
        }

        BlockBehavior raw_beh;
        std::size_t start = 0;
        while (start < nb) {
            std::size_t block_size = 1;
            while (start + block_size < nb) {
                const auto& super = jd.jor[
                    (start + block_size - 1) * nb + start + block_size];
                if (acb_close_to_zero(super.raw())) break;
                ++block_size;
            }
            BehaviorEntry e;
            e.mu = jd.jor[start * nb + start].clone();
            AcbValue shift;
            shift.set_si(floors[start]);
            acb_sub(e.mu.raw(), e.mu.raw(), shift.raw(), prec);
            e.log_power = static_cast<long>(block_size) - 1;
            raw_beh.push_back(std::move(e));
            start += block_size;
        }
        essential_behavior[bi] = union_behavior_local(raw_beh, prec);
    }

    if (!any_shift) {
        throw std::runtime_error(
            "calc_zero_acb_shearing_single_block: no algebraic shearing needed");
    }

    for (std::size_t bi = 0; bi < partition.size(); ++bi) {
        BlockBehavior subbeh;
        for (auto dep : partition.lower_dependencies(bi)) {
            for (const auto& bh : behavior[dep]) {
                BehaviorEntry e;
                e.mu = bh.mu.clone();
                e.log_power = bh.log_power;
                subbeh.push_back(std::move(e));
            }
        }
        subbeh = union_behavior_local(subbeh, prec);

        BlockBehavior joined;
        for (const auto& e : essential_behavior[bi]) {
            joined.push_back(BehaviorEntry{e.mu.clone(), e.log_power});
        }
        for (const auto& e : subbeh) {
            joined.push_back(BehaviorEntry{e.mu.clone(), e.log_power});
        }
        behavior[bi] = union_behavior2_local(joined, prec);
    }

    AcbRationalMatrix B_alg = invT_alg
        .matmul(AcbRationalMatrix::from_rational(norm.B, prec), prec)
        .matmul(T_alg, prec);
    for (std::size_t i = 0; i < n; ++i) {
        if (global_floors[i] == 0) continue;
        AcbValue c;
        c.set_si(global_floors[i]);
        B_alg(i, i) -= AcbRationalFunction::monomial(-1, c.raw());
    }

    AcbRationalMatrix T_total =
        AcbRationalMatrix::from_rational(norm.T, prec).matmul(T_alg, prec);
    AcbRationalMatrix invT_total =
        invT_alg.matmul(AcbRationalMatrix::from_rational(norm.invT, prec), prec);

    std::vector<AcbValue> bcT(n);
    for (std::size_t i = 0; i < n; ++i) {
        AcbValue acc;
        for (std::size_t j = 0; j < n; ++j) {
            AcbValue tij = evaluate_acb_rational(invT_total(i, j), x0, prec);
            AcbValue tmp;
            acb_mul(tmp.raw(), tij.raw(), bc[j].raw(), prec);
            acb_add(acc.raw(), acc.raw(), tmp.raw(), prec);
        }
        bcT[i] = std::move(acc);
    }

    std::vector<BlockEquationNum> nheqn;
    nheqn.reserve(partition.size());
    const long series_len = static_cast<long>(x_order()) + 1;
    std::vector<AcbValue> den_one(1);
    den_one[0].set_one();
    for (std::size_t bi = 0; bi < partition.size(); ++bi) {
        const auto& blk = partition.block(bi);
        const auto sub = partition.sub_rows(bi);
        const std::size_t nb = blk.size();

        BlockEquationNum eq;
        eq.dxexp.resize(1);
        eq.dxexp[0].set_one();
        eq.axexp.resize(nb);
        eq.bxexpn.resize(nb);
        eq.bxexpd.resize(nb);
        for (std::size_t a = 0; a < nb; ++a) {
            eq.axexp[a].resize(nb);
            for (std::size_t b = 0; b < nb; ++b) {
                AcbRationalFunction r = B_alg(blk[a], blk[b]);
                r *= AcbRationalFunction::monomial(1);
                eq.axexp[a][b] = acb_rational_taylor_coeffs(
                    std::move(r), series_len, prec);
            }

            eq.bxexpn[a].resize(sub.size());
            eq.bxexpd[a].resize(sub.size());
            for (std::size_t s = 0; s < sub.size(); ++s) {
                AcbRationalFunction r = B_alg(blk[a], sub[s]);
                r *= AcbRationalFunction::monomial(1);
                eq.bxexpn[a][s] = acb_rational_taylor_coeffs(
                    std::move(r), series_len, prec);
                eq.bxexpd[a][s] = clone_vec(den_one);
            }
        }
        eq.block = blk;
        eq.sub = sub;
        nheqn.push_back(std::move(eq));
    }

    std::vector<BlockEquation> symbolic_unused;
    auto raw = calcx00_impl(symbolic_unused, nheqn, bcT, x0,
                            behavior, nullptr, prec);
    auto T_ps = to_power_series_matrix_acb(T_total, prec);
    return ps_map_rule_set(T_ps, raw, prec);
}

}  // namespace

// ===========================================================================
//  Calcx00 helpers
// ===========================================================================

namespace {

struct ChainStorage {
    std::vector<std::vector<std::vector<std::vector<AcbValue>>>> bb;
    std::vector<std::vector<std::vector<std::vector<AcbValue>>>> kk;

    long Nblk = 0;

    void resize(long n_behaviors, long max_p_plus_one, long n_eta) {
        bb.clear();
        kk.clear();
        bb.resize(n_behaviors);
        kk.resize(n_behaviors);
        for (long m = 0; m < n_behaviors; ++m) {
            bb[m].resize(max_p_plus_one);
            kk[m].resize(max_p_plus_one);
            for (long p = 0; p < max_p_plus_one; ++p) {
                bb[m][p].resize(n_eta);
                kk[m][p].resize(n_eta);
            }
        }
    }

    std::vector<AcbValue>& bb_at(long m, long p, long n) {
        ensure_size(m, p, n);
        auto& v = bb[m][p][n];
        if (v.empty()) v.resize(Nblk);
        return v;
    }
    std::vector<AcbValue>& kk_at(long m, long p, long n) {
        ensure_size(m, p, n);
        auto& v = kk[m][p][n];
        if (v.empty()) v.resize(Nblk * Nblk);
        return v;
    }

private:
    void ensure_size(long m, long p, long n) {
        if (m >= static_cast<long>(bb.size())) {
            bb.resize(m + 1);
            kk.resize(m + 1);
        }
        if (p >= static_cast<long>(bb[m].size())) {
            bb[m].resize(p + 1);
            kk[m].resize(p + 1);
        }
        if (n >= static_cast<long>(bb[m][p].size())) {
            bb[m][p].resize(n + 1);
            kk[m][p].resize(n + 1);
        }
    }
};

acb_srcptr mat_at(const std::vector<AcbValue>& v, long Nblk, long i, long j) {
    return v[i * Nblk + j].raw();
}

std::vector<AcbValue> matvec(const std::vector<AcbValue>& M,
                             const std::vector<AcbValue>& v,
                             long Nblk, long prec) {
    std::vector<AcbValue> out(Nblk);
    AcbValue tmp;
    for (long i = 0; i < Nblk; ++i) {
        AcbValue acc;
        for (long j = 0; j < Nblk; ++j) {
            acb_mul(tmp.raw(), mat_at(M, Nblk, i, j), v[j].raw(), prec);
            acb_add(acc.raw(), acc.raw(), tmp.raw(), prec);
        }
        out[i] = std::move(acc);
    }
    return out;
}

std::vector<AcbValue> matmul_flat(const std::vector<AcbValue>& M,
                                  const std::vector<AcbValue>& X,
                                  long Nblk, long prec) {
    std::vector<AcbValue> out(Nblk * Nblk);
    AcbValue tmp;
    for (long i = 0; i < Nblk; ++i) {
        for (long j = 0; j < Nblk; ++j) {
            AcbValue acc;
            for (long k = 0; k < Nblk; ++k) {
                acb_mul(tmp.raw(), mat_at(M, Nblk, i, k), mat_at(X, Nblk, k, j), prec);
                acb_add(acc.raw(), acc.raw(), tmp.raw(), prec);
            }
            out[i * Nblk + j] = std::move(acc);
        }
    }
    return out;
}

void vec_add_inplace(std::vector<AcbValue>& v, const std::vector<AcbValue>& w, long prec) {
    for (std::size_t i = 0; i < v.size(); ++i) {
        acb_add(v[i].raw(), v[i].raw(), w[i].raw(), prec);
    }
}

void vec_scale_inplace(std::vector<AcbValue>& v, acb_srcptr s, long prec) {
    AcbValue tmp;
    for (std::size_t i = 0; i < v.size(); ++i) {
        acb_mul(tmp.raw(), v[i].raw(), s, prec);
        acb_set(v[i].raw(), tmp.raw());
    }
}

std::vector<AcbValue> zero_vec(long n) { return std::vector<AcbValue>(n); }
std::vector<AcbValue> identity_mat_flat(long n) {
    std::vector<AcbValue> out(n * n);
    for (long i = 0; i < n; ++i) out[i * n + i].set_one();
    return out;
}

std::vector<AcbValue> compute_a00(
        const std::vector<AcbValue>& dxexp,
        const std::vector<std::vector<std::vector<AcbValue>>>& axexp,
        long Nblk, long prec) {
    if (dxexp.empty() || dxexp[0].is_zero()) {
        throw std::runtime_error(
            "calcx00: dxexp[0] is zero -- the matrix is not in the expected form.");
    }
    AcbValue inv_d0;
    acb_inv(inv_d0.raw(), dxexp[0].raw(), prec);

    std::vector<AcbValue> a00(Nblk * Nblk);
    for (long i = 0; i < Nblk; ++i) {
        for (long j = 0; j < Nblk; ++j) {
            AcbValue v;
            if (!axexp[i][j].empty()) {
                acb_set(v.raw(), axexp[i][j][0].raw());
            }
            acb_mul(a00[i * Nblk + j].raw(), v.raw(), inv_d0.raw(), prec);
        }
    }
    return a00;
}

std::vector<long> find_resonance_positions(acb_srcptr spbeh,
                                           long n,
                                           const std::vector<AcbValue>& a00,
                                           long Nblk,
                                           long prec) {
    std::vector<long> out;
    AcbValue lhs, rhs, diff;
    acb_set(lhs.raw(), spbeh);
    AcbValue n_acb; n_acb.set_si(n);
    acb_add(lhs.raw(), lhs.raw(), n_acb.raw(), prec);

    for (long p = 0; p < Nblk; ++p) {
        acb_set(rhs.raw(), mat_at(a00, Nblk, p, p));
        acb_sub(diff.raw(), lhs.raw(), rhs.raw(), prec);
        // Use the resonance tolerance (≈ working precision), NOT chop_pre:
        // an eps-deformed exponent pair {mu, mu + a*eps} is genuinely
        // distinct for eps != 0 and must NOT be merged just because
        // |a*eps| < 10^-chop_pre -- doing so silently discards the genuine
        // 1/(a*eps) pole and corrupts the eps->0 extrapolation to 0 digits.
        if (acb_is_chop_zero(diff.raw(), numeric::resonance_chop_pre()))
            out.push_back(p);
    }
    return out;
}

std::vector<AcbValue> compute_a00_from_mata(
        const std::vector<AcbValue>& dxexp,
        const std::vector<AcbValue>& mata0,
        long Nblk,
        long prec) {
    if (dxexp.empty() || dxexp[0].is_zero()) {
        throw std::runtime_error(
            "calcx00: dxexp[0] is zero -- the matrix is not in the expected form.");
    }
    AcbValue inv_d0;
    acb_inv(inv_d0.raw(), dxexp[0].raw(), prec);

    std::vector<AcbValue> a00(Nblk * Nblk);
    for (long i = 0; i < Nblk; ++i) {
        for (long j = 0; j < Nblk; ++j) {
            acb_mul(a00[i * Nblk + j].raw(),
                    mata0[i * Nblk + j].raw(), inv_d0.raw(), prec);
        }
    }
    return a00;
}

void rotate_series_rows_left(
        std::vector<std::vector<AcbValue>>& rows_by_n,
        const std::vector<AcbValue>&        left,
        long                                Nblk,
        long                                prec) {
    if (rows_by_n.empty()) return;

    const long n_eta = static_cast<long>(rows_by_n[0].size());
    std::vector<std::vector<AcbValue>> rotated(static_cast<std::size_t>(Nblk));
    for (long i = 0; i < Nblk; ++i) {
        rotated[static_cast<std::size_t>(i)].resize(static_cast<std::size_t>(n_eta));
    }

    for (long n = 0; n < n_eta; ++n) {
        std::vector<AcbValue> column(Nblk);
        for (long i = 0; i < Nblk; ++i) {
            if (n < static_cast<long>(rows_by_n[i].size())) {
                acb_set(column[i].raw(), rows_by_n[i][n].raw());
            }
        }
        std::vector<AcbValue> transformed = matvec(left, column, Nblk, prec);
        for (long i = 0; i < Nblk; ++i) {
            acb_set(rotated[i][n].raw(), transformed[i].raw());
        }
    }

    rows_by_n = std::move(rotated);
}

BlockBehavior union_behavior_local(const BlockBehavior& in, long prec) {
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
            // Exponent-equality test -> resonance tolerance (see
            // resonance_chop_pre): eps-deformed exponents stay distinct.
            if (acb_is_chop_zero(diff.raw(), numeric::resonance_chop_pre())) {
                if (in[j].log_power > combined.log_power) {
                    combined.log_power = in[j].log_power;
                }
                done[j] = true;
            }
        }
        out.push_back(std::move(combined));
    }
    return out;
}

BlockBehavior union_behavior2_local(const BlockBehavior& in, long prec) {
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
            // Exponent-equality test -> resonance tolerance.
            if (acb_is_chop_zero(diff.raw(), numeric::resonance_chop_pre())) {
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

BlockBehavior block_behavior_from_rotation(const BlockJordanRotation& rot,
                                           long                     prec) {
    const long Nblk = static_cast<long>(rot.block.size());
    BlockBehavior raw;
    for (long start = 0; start < Nblk; ) {
        long block_size = 1;
        while (start + block_size < Nblk) {
            AcbValue diff;
            acb_sub(diff.raw(),
                    rot.jor[(start + block_size) * Nblk + (start + block_size)].raw(),
                    rot.jor[start * Nblk + start].raw(), prec);
            // Two Jordan diagonal entries are the SAME exponent only to the
            // resonance tolerance; an eps-deformed pair stays distinct.
            if (!acb_is_chop_zero(diff.raw(), numeric::resonance_chop_pre())) break;

            // Super-diagonal entry: a genuine STRUCTURAL coupling (matrix
            // element), so the coarse coefficient-negligibility chop_pre is
            // the right test here.
            AcbValue super = rot.jor[(start + block_size - 1) * Nblk + start + block_size].clone();
            if (acb_is_chop_zero(super.raw(), chop_pre())) break;
            ++block_size;
        }

        BehaviorEntry e;
        e.mu = rot.jor[start * Nblk + start].clone();
        e.log_power = block_size - 1;
        raw.push_back(std::move(e));
        start += block_size;
    }
    return union_behavior_local(raw, prec);
}

AsymptoticBehaviorList asymptotic_behavior_for_calc_zero(
        const RationalMatrix&                   mat,
        const std::vector<BlockJordanRotation>& rotations,
        long                                    prec) {
    auto partition = build_partition(mat);
    if (partition.size() != rotations.size()) {
        throw std::runtime_error(
            "asymptotic_behavior_for_calc_zero: block/rotation size mismatch");
    }

    AsymptoticBehaviorList behavior(partition.size());
    for (std::size_t i = 0; i < partition.size(); ++i) {
        BlockBehavior essbeh = block_behavior_from_rotation(rotations[i], prec);
        BlockBehavior subbeh;
        for (auto j : partition.lower_dependencies(i)) {
            for (const auto& bh : behavior[j]) {
                BehaviorEntry e;
                e.mu = bh.mu.clone();
                e.log_power = bh.log_power;
                subbeh.push_back(std::move(e));
            }
        }
        subbeh = union_behavior_local(subbeh, prec);

        BlockBehavior joined;
        for (auto& e : essbeh) joined.push_back(BehaviorEntry{e.mu.clone(), e.log_power});
        for (auto& e : subbeh) joined.push_back(BehaviorEntry{e.mu.clone(), e.log_power});
        behavior[i] = union_behavior2_local(joined, prec);
    }
    return behavior;
}

}  // namespace

// ===========================================================================
//  calcx00_impl
// ===========================================================================

std::vector<AsyExpansion>
calcx00_impl(const std::vector<BlockEquation>&            /*nheq*/,
             const std::vector<BlockEquationNum>&         nheqn,
             const std::vector<AcbValue>&                 bc,
             acb_srcptr                                   x0,
             const AsymptoticBehaviorList&                behavior,
             const std::vector<BlockJordanRotation>*      rotations,
             long                                         prec) {
    long xorder = static_cast<long>(x_order());

    std::vector<AsyExpansion> result(bc.size());

    for (std::size_t k = 0; k < nheqn.size(); ++k) {
        const auto& eq        = nheqn[k];
        const auto& dxexp     = eq.dxexp;
        const auto& axexp     = eq.axexp;
        const auto& bxexpn    = eq.bxexpn;
        const auto& bxexpd    = eq.bxexpd;
        const auto& block     = eq.block;
        const auto& sub       = eq.sub;
        const std::size_t Nblk = block.size();
        // Mutable local copy: the Jordan/log auto-rerun may need to bump
        // beh[m].log_power when a resonance forces a log term beyond the
        // supplied power (textbook Frobenius "second solution").
        BlockBehavior beh;
        beh.reserve(behavior[k].size());
        for (const auto& e : behavior[k]) {
            beh.push_back(BehaviorEntry{e.mu.clone(), e.log_power});
        }
        const BlockJordanRotation* rot = nullptr;
        if (rotations != nullptr) {
            if (k >= rotations->size()) {
                throw std::runtime_error(
                    "calcx00: rotation list shorter than block equation list");
            }
            rot = &(*rotations)[k];
            if (rot->block != block) {
                throw std::runtime_error(
                    "calcx00: block/rotation mismatch at block index "
                    + std::to_string(k));
            }
        }

        if (!silent_mode()) {
            std::string block_str = "[";
            for (auto v : block) { block_str += std::to_string(v) + " "; }
            block_str += "]";
            log_line("Calcx00: block " + block_str + " (behaviors=" +
                     std::to_string(beh.size()) + ")");
        }

        long maxd = static_cast<long>(dxexp.size()) - 1;
        long maxa = -1;
        for (const auto& row : axexp)
            for (const auto& vec : row)
                if (static_cast<long>(vec.size()) - 1 > maxa)
                    maxa = static_cast<long>(vec.size()) - 1;

        std::vector<std::vector<AcbValue>> mata(maxa < 0 ? 0 : (maxa + 1));
        for (long m = 0; m <= maxa; ++m) {
            mata[m].clear();
            mata[m].resize(Nblk * Nblk);
            for (long i = 0; i < (long)Nblk; ++i) {
                for (long j = 0; j < (long)Nblk; ++j) {
                    if (m < static_cast<long>(axexp[i][j].size())) {
                        acb_set(mata[m][i * Nblk + j].raw(), axexp[i][j][m].raw());
                    }
                }
            }
            if (rot != nullptr) {
                std::vector<AcbValue> tmp = matmul_flat(rot->invu, mata[m], Nblk, prec);
                mata[m] = matmul_flat(tmp, rot->u, Nblk, prec);
            }
        }

        std::vector<AcbValue> a00 = (rot != nullptr)
            ? clone_vec(rot->jor)
            : (maxa >= 0
                ? compute_a00_from_mata(dxexp, mata[0], Nblk, prec)
                : compute_a00(dxexp, axexp, Nblk, prec));

        auto lookup_f = [&](std::size_t intid, acb_srcptr spbeh, long p) -> std::vector<AcbValue> {
            for (const auto& term : result[intid]) {
                if (acb_close(term.mu.raw(), spbeh, prec)) {
                    if (p < static_cast<long>(term.exp.size())) {
                        return clone_vec(term.exp[p]);
                    } else {
                        return zero_vec(xorder + 1);
                    }
                }
            }
            return zero_vec(xorder + 1);
        };

        // ==================================================================
        //  Jordan/log auto-rerun loop.
        //
        //  The whole Frobenius-chain build + boundary solve + output for
        //  this block is wrapped so that, if the post-solve residue check
        //  finds a nonzero log^{logk+1} term (a resonance has FORCED a log
        //  beyond the supplied beh[m].log_power), we bump that log_power
        //  and redo the block at the higher tower instead of returning a
        //  silently-truncated expansion.  Bounded at +calcx00_max_log_bump
        //  total bumps per block; beyond that we keep the warning so a
        //  pathological resonance ladder cannot loop forever.
        // ==================================================================
        constexpr int calcx00_max_log_bump = 4;
        for (int log_attempt = 0;; ++log_attempt) {

        std::vector<std::vector<std::vector<std::vector<AcbValue>>>> nh_all(beh.size());
        for (std::size_t m = 0; m < beh.size(); ++m) {
            long logk = beh[m].log_power;
            nh_all[m].resize(logk + 1);
            for (long p = 0; p <= logk; ++p) {
                std::vector<std::vector<AcbValue>> sub_expansions(sub.size());
                for (std::size_t s = 0; s < sub.size(); ++s) {
                    sub_expansions[s] = lookup_f(sub[s], beh[m].mu.raw(), p);
                }
                std::vector<std::vector<AcbValue>> nh_per_row =
                    map_rational_expansion(bxexpn, bxexpd, sub_expansions, prec);
                for (auto& v : nh_per_row) {
                    if (v.empty()) v.resize(xorder + 1);
                }
                if (rot != nullptr) {
                    rotate_series_rows_left(nh_per_row, rot->invu, Nblk, prec);
                }
                nh_all[m][p] = std::move(nh_per_row);
            }
        }

        long max_logk = 0;
        for (const auto& b : beh) max_logk = std::max(max_logk, b.log_power);

        ChainStorage chain;
        chain.Nblk = Nblk;
        chain.resize(static_cast<long>(beh.size()), max_logk + 2, xorder + 1);

        std::vector<bool> is_essential(beh.size(), false);
        struct PairData {
            std::vector<AcbValue> kk_logp1;
            std::vector<AcbValue> bb_logp1;
            std::vector<long>     complement_positions;
            // Resonant positions (posi) at this n.  The complement rows give
            // the primary pairing constraints (log^{logk+1} must vanish on
            // the regular subspace).  The resonant-position rows are either
            // identically zero (diagonal A_00) or, for a true Jordan block,
            // genuine extra constraints; they are kept as a SECONDARY pool
            // for the over-determined fallback below so the Jordan/log case
            // contributes the rank it carries.
            std::vector<long>     resonant_positions;
            long                  at_order_n = -1;
            bool                  populated = false;
        };
        // One pairing-constraint snapshot PER RESONANT n, not just the last.
        // A behavior m that resonates at several n (its mu lies an integer
        // below several diagonal entries of A_00) yields one PairData per
        // hit; the boundary system needs ALL of them.  The previous
        // single-slot `pair_data[m]` overwrote earlier hits, which is the
        // root cause of the n_pair undercount that drove
        // compensate > Nblk on the two-mass 5-point top block.
        std::vector<std::vector<PairData>> pair_data(beh.size());
        std::vector<std::vector<std::pair<long, std::vector<long>>>> resonance_log(beh.size());

        AcbValue tmp, tmp2;

        auto rec1 = [&](bool is_bb, long m, long p, long n,
                        acb_srcptr spbeh) -> std::vector<AcbValue> {
            std::vector<AcbValue> acc(Nblk);

            long mm_max1 = std::min<long>(n, maxa);
            for (long mm = 0; mm <= mm_max1; ++mm) {
                std::vector<AcbValue>& src = is_bb
                    ? chain.bb_at(m, p, n - mm)
                    : chain.kk_at(m, p, n - mm);
                if (src.size() != static_cast<std::size_t>(is_bb ? Nblk : Nblk * Nblk)) {
                    src.resize(is_bb ? Nblk : Nblk * Nblk);
                }
                std::vector<AcbValue> contrib;
                if (is_bb) {
                    contrib = matvec(mata[mm], src, Nblk, prec);
                } else {
                    contrib = matmul_flat(mata[mm], src, Nblk, prec);
                }
                if (acc.size() != contrib.size()) { acc.clear(); acc.resize(contrib.size()); }
                vec_add_inplace(acc, contrib, prec);
            }

            if (is_bb) {
                if (p < static_cast<long>(nh_all[m].size())
                    && n < static_cast<long>(nh_all[m][p][0].size())) {
                    if (acc.size() != static_cast<std::size_t>(Nblk)) {
                        acc.clear();
                        acc.resize(Nblk);
                    }
                    for (long i = 0; i < (long)Nblk; ++i) {
                        if (n < static_cast<long>(nh_all[m][p][i].size())) {
                            acb_add(acc[i].raw(), acc[i].raw(),
                                    nh_all[m][p][i][n].raw(), prec);
                        }
                    }
                }
            }

            long mm_max3 = std::min<long>(n, maxd);
            for (long mm = 0; mm <= mm_max3; ++mm) {
                std::vector<AcbValue>& src = is_bb
                    ? chain.bb_at(m, p, n - mm)
                    : chain.kk_at(m, p, n - mm);
                AcbValue factor;
                acb_set(factor.raw(), spbeh);
                AcbValue add; add.set_si(n - mm);
                acb_add(factor.raw(), factor.raw(), add.raw(), prec);
                acb_mul(factor.raw(), factor.raw(), dxexp[mm].raw(), prec);
                acb_neg(factor.raw(), factor.raw());
                std::vector<AcbValue> scaled = clone_vec(src);
                vec_scale_inplace(scaled, factor.raw(), prec);
                if (acc.size() != scaled.size()) { acc.clear(); acc.resize(scaled.size()); }
                vec_add_inplace(acc, scaled, prec);
            }

            long mm_max4 = std::min<long>(n, maxd);
            for (long mm = 1; mm <= mm_max4; ++mm) {
                std::vector<AcbValue>& src = is_bb
                    ? chain.bb_at(m, p + 1, n - mm)
                    : chain.kk_at(m, p + 1, n - mm);
                AcbValue factor;
                AcbValue ppp; ppp.set_si(-(p + 1));
                acb_mul(factor.raw(), dxexp[mm].raw(), ppp.raw(), prec);
                std::vector<AcbValue> scaled = clone_vec(src);
                vec_scale_inplace(scaled, factor.raw(), prec);
                if (acc.size() != scaled.size()) { acc.clear(); acc.resize(scaled.size()); }
                vec_add_inplace(acc, scaled, prec);
            }

            AcbValue denom;
            acb_mul_si(denom.raw(), dxexp[0].raw(), p + 1, prec);
            AcbValue inv_denom;
            acb_inv(inv_denom.raw(), denom.raw(), prec);
            vec_scale_inplace(acc, inv_denom.raw(), prec);
            return acc;
        };

        auto rec2 = [&](bool is_bb, long m, long p, long n,
                        const std::vector<AcbValue>& inv_mat,
                        acb_srcptr spbeh) -> std::vector<AcbValue> {
            std::vector<AcbValue> acc(Nblk);

            long mm_max1 = std::min<long>(n, maxa);
            for (long mm = 1; mm <= mm_max1; ++mm) {
                std::vector<AcbValue>& src = is_bb
                    ? chain.bb_at(m, p, n - mm)
                    : chain.kk_at(m, p, n - mm);
                std::vector<AcbValue> contrib;
                if (is_bb) contrib = matvec(mata[mm], src, Nblk, prec);
                else       contrib = matmul_flat(mata[mm], src, Nblk, prec);
                if (acc.size() != contrib.size()) { acc.clear(); acc.resize(contrib.size()); }
                vec_add_inplace(acc, contrib, prec);
            }

            if (is_bb) {
                if (p < static_cast<long>(nh_all[m].size())
                    && n < static_cast<long>(nh_all[m][p][0].size())) {
                    if (acc.size() != static_cast<std::size_t>(Nblk)) {
                        acc.clear();
                        acc.resize(Nblk);
                    }
                    for (long i = 0; i < (long)Nblk; ++i) {
                        if (n < static_cast<long>(nh_all[m][p][i].size())) {
                            acb_add(acc[i].raw(), acc[i].raw(),
                                    nh_all[m][p][i][n].raw(), prec);
                        }
                    }
                }
            }

            long mm_max3 = std::min<long>(n, maxd);
            for (long mm = 1; mm <= mm_max3; ++mm) {
                std::vector<AcbValue>& src = is_bb
                    ? chain.bb_at(m, p, n - mm)
                    : chain.kk_at(m, p, n - mm);
                AcbValue factor;
                acb_set(factor.raw(), spbeh);
                AcbValue add; add.set_si(n - mm);
                acb_add(factor.raw(), factor.raw(), add.raw(), prec);
                acb_mul(factor.raw(), factor.raw(), dxexp[mm].raw(), prec);
                acb_neg(factor.raw(), factor.raw());
                std::vector<AcbValue> scaled = clone_vec(src);
                vec_scale_inplace(scaled, factor.raw(), prec);
                if (acc.size() != scaled.size()) { acc.clear(); acc.resize(scaled.size()); }
                vec_add_inplace(acc, scaled, prec);
            }

            long mm_max4 = std::min<long>(n, maxd);
            for (long mm = 0; mm <= mm_max4; ++mm) {
                std::vector<AcbValue>& src = is_bb
                    ? chain.bb_at(m, p + 1, n - mm)
                    : chain.kk_at(m, p + 1, n - mm);
                AcbValue factor;
                AcbValue ppp; ppp.set_si(-(p + 1));
                acb_mul(factor.raw(), dxexp[mm].raw(), ppp.raw(), prec);
                std::vector<AcbValue> scaled = clone_vec(src);
                vec_scale_inplace(scaled, factor.raw(), prec);
                if (acc.size() != scaled.size()) { acc.clear(); acc.resize(scaled.size()); }
                vec_add_inplace(acc, scaled, prec);
            }

            std::vector<AcbValue> result_vec;
            if (is_bb) {
                result_vec = matvec(inv_mat, acc, Nblk, prec);
            } else {
                result_vec = matmul_flat(inv_mat, acc, Nblk, prec);
            }
            return result_vec;
        };

        // Residue of the Frobenius recursion at (m, p, n) -- this is the
        // body of rec2 BEFORE the final inv_mat multiply.  Returned as a
        // pair so the rank-deficient (later-resonance) path can both
        // pseudo-solve and extract the solvability constraint from the
        // same evaluation.
        auto rec2_residue = [&](bool is_bb, long m, long p, long n,
                                acb_srcptr spbeh) -> std::vector<AcbValue> {
            std::vector<AcbValue> acc(Nblk);

            long mm_max1 = std::min<long>(n, maxa);
            for (long mm = 1; mm <= mm_max1; ++mm) {
                std::vector<AcbValue>& src = is_bb
                    ? chain.bb_at(m, p, n - mm)
                    : chain.kk_at(m, p, n - mm);
                std::vector<AcbValue> contrib;
                if (is_bb) contrib = matvec(mata[mm], src, Nblk, prec);
                else       contrib = matmul_flat(mata[mm], src, Nblk, prec);
                if (acc.size() != contrib.size()) { acc.clear(); acc.resize(contrib.size()); }
                vec_add_inplace(acc, contrib, prec);
            }

            if (is_bb) {
                if (p < static_cast<long>(nh_all[m].size())
                    && n < static_cast<long>(nh_all[m][p][0].size())) {
                    if (acc.size() != static_cast<std::size_t>(Nblk)) {
                        acc.clear();
                        acc.resize(Nblk);
                    }
                    for (long i = 0; i < (long)Nblk; ++i) {
                        if (n < static_cast<long>(nh_all[m][p][i].size())) {
                            acb_add(acc[i].raw(), acc[i].raw(),
                                    nh_all[m][p][i][n].raw(), prec);
                        }
                    }
                }
            }

            long mm_max3 = std::min<long>(n, maxd);
            for (long mm = 1; mm <= mm_max3; ++mm) {
                std::vector<AcbValue>& src = is_bb
                    ? chain.bb_at(m, p, n - mm)
                    : chain.kk_at(m, p, n - mm);
                AcbValue factor;
                acb_set(factor.raw(), spbeh);
                AcbValue add; add.set_si(n - mm);
                acb_add(factor.raw(), factor.raw(), add.raw(), prec);
                acb_mul(factor.raw(), factor.raw(), dxexp[mm].raw(), prec);
                acb_neg(factor.raw(), factor.raw());
                std::vector<AcbValue> scaled = clone_vec(src);
                vec_scale_inplace(scaled, factor.raw(), prec);
                if (acc.size() != scaled.size()) { acc.clear(); acc.resize(scaled.size()); }
                vec_add_inplace(acc, scaled, prec);
            }

            long mm_max4 = std::min<long>(n, maxd);
            for (long mm = 0; mm <= mm_max4; ++mm) {
                std::vector<AcbValue>& src = is_bb
                    ? chain.bb_at(m, p + 1, n - mm)
                    : chain.kk_at(m, p + 1, n - mm);
                AcbValue factor;
                AcbValue ppp; ppp.set_si(-(p + 1));
                acb_mul(factor.raw(), dxexp[mm].raw(), ppp.raw(), prec);
                std::vector<AcbValue> scaled = clone_vec(src);
                vec_scale_inplace(scaled, factor.raw(), prec);
                if (acc.size() != scaled.size()) { acc.clear(); acc.resize(scaled.size()); }
                vec_add_inplace(acc, scaled, prec);
            }
            return acc;
        };

        auto invert_shifted_a00 = [&](long n, acb_srcptr spbeh,
                                      std::vector<AcbValue>& inv_mat) -> bool {
            std::vector<AcbValue> M(static_cast<std::size_t>(Nblk * Nblk));
            AcbValue n_plus_s;
            acb_set(n_plus_s.raw(), spbeh);
            AcbValue add; add.set_si(n);
            acb_add(n_plus_s.raw(), n_plus_s.raw(), add.raw(), prec);
            for (long i = 0; i < (long)Nblk; ++i) {
                for (long j = 0; j < (long)Nblk; ++j) {
                    AcbValue v;
                    acb_neg(v.raw(), mat_at(a00, Nblk, i, j));
                    if (i == j) acb_add(v.raw(), v.raw(), n_plus_s.raw(), prec);
                    acb_set(M[i * Nblk + j].raw(), v.raw());
                }
            }
            acb_mat_t M_acb, inv_acb;
            acb_mat_init(M_acb, Nblk, Nblk);
            acb_mat_init(inv_acb, Nblk, Nblk);
            for (long i = 0; i < (long)Nblk; ++i)
                for (long j = 0; j < (long)Nblk; ++j)
                    acb_set(acb_mat_entry(M_acb, i, j), M[i * Nblk + j].raw());
            int ok = acb_mat_inv(inv_acb, M_acb, prec);
            if (ok) {
                AcbValue inv_d0;
                acb_inv(inv_d0.raw(), dxexp[0].raw(), prec);
                inv_mat.clear();
                inv_mat.resize(Nblk * Nblk);
                for (long i = 0; i < (long)Nblk; ++i)
                    for (long j = 0; j < (long)Nblk; ++j) {
                        acb_mul(inv_mat[i * Nblk + j].raw(),
                                acb_mat_entry(inv_acb, i, j),
                                inv_d0.raw(), prec);
                    }
            }
            acb_mat_clear(M_acb);
            acb_mat_clear(inv_acb);
            return ok != 0;
        };

        // Deflated (rank-deficient) inverse of M = (mu+n)I - A_00, scaled by
        // 1/dxexp[0] (same scaling as invert_shifted_a00).  Rows AND columns
        // at `posi` are projected to zero before inversion of the complement
        // sub-block, then embedded back.  For diagonal A_00 (the common case
        // after the Jordan rotation) this is the exact Moore-Penrose
        // pseudoinverse.  For Jordan A_00 with off-diagonal coupling between
        // posi and complement the result is approximate; the boolean
        // `coupled` flags that case so the caller can warn / fall back.
        auto deflated_invert_shifted_a00 =
            [&](long n, acb_srcptr spbeh, const std::vector<long>& posi,
                std::vector<AcbValue>& pinv_mat, bool& coupled) -> bool {
            coupled = false;
            const long Nc = (long)Nblk - (long)posi.size();
            pinv_mat.clear();
            pinv_mat.resize(Nblk * Nblk);
            if (Nc <= 0) return true;  // pinv = 0

            std::vector<bool> is_posi(Nblk, false);
            for (long ip : posi) is_posi[ip] = true;
            std::vector<long> comp;
            comp.reserve(Nc);
            for (long i = 0; i < (long)Nblk; ++i)
                if (!is_posi[i]) comp.push_back(i);

            AcbValue n_plus_s;
            acb_set(n_plus_s.raw(), spbeh);
            { AcbValue add; add.set_si(n);
              acb_add(n_plus_s.raw(), n_plus_s.raw(), add.raw(), prec); }

            // Detect off-diagonal coupling between posi and complement.
            for (long ip : posi) {
                for (long jc : comp) {
                    if (!acb_is_chop_zero(mat_at(a00, Nblk, ip, jc),
                                          chop_pre())
                        || !acb_is_chop_zero(mat_at(a00, Nblk, jc, ip),
                                             chop_pre())) {
                        coupled = true;
                    }
                }
            }

            acb_mat_t Mc, Ic;
            acb_mat_init(Mc, Nc, Nc);
            acb_mat_init(Ic, Nc, Nc);
            for (long ri = 0; ri < Nc; ++ri) {
                for (long rj = 0; rj < Nc; ++rj) {
                    AcbValue v;
                    acb_neg(v.raw(),
                            mat_at(a00, Nblk, comp[ri], comp[rj]));
                    if (ri == rj) acb_add(v.raw(), v.raw(),
                                          n_plus_s.raw(), prec);
                    acb_set(acb_mat_entry(Mc, ri, rj), v.raw());
                }
            }
            int ok = acb_mat_inv(Ic, Mc, prec);
            if (ok) {
                AcbValue inv_d0;
                acb_inv(inv_d0.raw(), dxexp[0].raw(), prec);
                for (long ri = 0; ri < Nc; ++ri)
                    for (long rj = 0; rj < Nc; ++rj) {
                        acb_mul(pinv_mat[comp[ri] * Nblk + comp[rj]].raw(),
                                acb_mat_entry(Ic, ri, rj),
                                inv_d0.raw(), prec);
                    }
            }
            acb_mat_clear(Mc);
            acb_mat_clear(Ic);
            return ok != 0;
        };

        bool block_has_later_resonance = false;

        for (long m = 0; m < (long)beh.size(); ++m) {
            const long logk = beh[m].log_power;
            bool seen_first_resonance = false;

            for (long n = 0; n <= xorder; ++n) {
                std::vector<long> posi = find_resonance_positions(
                    beh[m].mu.raw(), n, a00, Nblk, prec);
                std::vector<AcbValue> inv_mat;
                bool inv_ok = false;
                bool full_matrix_resonance = false;
                if (posi.empty()) {
                    inv_ok = invert_shifted_a00(n, beh[m].mu.raw(), inv_mat);
                    if (!inv_ok) {
                        full_matrix_resonance = true;
                        for (long p = 0; p < (long)Nblk; ++p) posi.push_back(p);
                    }
                }

                if (!posi.empty() && !seen_first_resonance) {
                    // -----------------------------------------------------
                    //  FIRST resonance for this behavior: introduce the
                    //  Nblk-vector of free parameters c_m by setting
                    //  kk[m][0][n] = I.  rec1 propagates to log^{logk+1};
                    //  the complement rows of that layer are the primary
                    //  pairing constraints (they pin c_m to the eigenspace
                    //  associated with mu_m at this n).  Unchanged from the
                    //  pre-patch path.
                    // -----------------------------------------------------
                    seen_first_resonance = true;
                    resonance_log[m].push_back({n, posi});
                    is_essential[m] = true;

                    chain.kk_at(m, 0, n) = identity_mat_flat(Nblk);
                    chain.bb_at(m, 0, n) = zero_vec(Nblk);

                    for (long p = 0; p <= logk; ++p) {
                        chain.bb_at(m, p + 1, n) = rec1(true,  m, p, n, beh[m].mu.raw());
                        chain.kk_at(m, p + 1, n) = rec1(false, m, p, n, beh[m].mu.raw());
                    }
                    PairData pd;
                    pd.kk_logp1 = clone_vec(chain.kk_at(m, logk + 1, n));
                    pd.bb_logp1 = clone_vec(chain.bb_at(m, logk + 1, n));
                    std::vector<long> comp;
                    for (long p = 0; p < (long)Nblk; ++p) {
                        if (full_matrix_resonance
                            || std::find(posi.begin(), posi.end(), p) == posi.end()) {
                            comp.push_back(p);
                        }
                    }
                    pd.complement_positions = std::move(comp);
                    pd.resonant_positions   = posi;
                    pd.at_order_n           = n;
                    pd.populated            = true;
                    pair_data[m].push_back(std::move(pd));

                } else if (!posi.empty()) {
                    // -----------------------------------------------------
                    //  LATER resonance for this behavior (mu_m + n hits a
                    //  second/third eigenvalue of A_00 at integer offset).
                    //
                    //  The pre-patch code reset kk[m][0][n] = I here too,
                    //  which (a) discarded the dependence on c_m carried by
                    //  the n'<n recursion, corrupting the chain, and
                    //  (b) generated a pairing row that constrained c_m as
                    //  if it were the n-th coefficient -- generically wrong
                    //  -- and then OVERWROTE the n'=first pairing row.
                    //  That overwrite is the n_pair undercount that drove
                    //  compensate=9 > Nblk=7 on a 5-point production top block.
                    //
                    //  Correct (textbook Frobenius) handling: M = (mu+n)I
                    //  - A_00 is singular with kernel along posi.  Compute
                    //  the recursion residue (rec2 body without inv_mat),
                    //  then
                    //    * complement rows: solve via the deflated inverse
                    //      (exact for diagonal A_00; flagged if Jordan-
                    //      coupled).
                    //    * posi rows: kk[m][p][n] = 0 (the kernel freedom
                    //      at this order is already parameterised by the
                    //      OTHER essential behavior with mu' = mu_m + n at
                    //      its own n=0).
                    //    * SOLVABILITY constraint: residue at posi rows
                    //      must vanish; this is the genuine extra pairing
                    //      row on c_m.  This is also exactly the
                    //      "obstruction" whose nonzero bb-part forces a
                    //      resonance-induced log term -- the Jordan/log
                    //      case.
                    // -----------------------------------------------------
                    resonance_log[m].push_back({n, posi});
                    block_has_later_resonance = true;

                    std::vector<AcbValue> pinv;
                    bool coupled = false;
                    bool pinv_ok = deflated_invert_shifted_a00(
                        n, beh[m].mu.raw(), posi, pinv, coupled);
                    if (!pinv_ok) {
                        throw std::runtime_error(
                            "calcx00: deflated complement sub-block "
                            "is singular at later resonance");
                    }
                    if (coupled) {
                        fprintf(stderr,
                            "[amflow][calcx00] note: behavior m=%ld later "
                            "resonance at n=%ld has Jordan off-diagonal "
                            "posi/complement coupling; deflated inverse is "
                            "approximate (exact only for diagonal A_00).\n",
                            m, n);
                    }

                    // Ensure the (logk+1) layer at this n is zero before
                    // descending p (rec2_residue reads kk[m][p+1][n] at
                    // mm=0).
                    chain.kk_at(m, logk + 1, n).clear();
                    chain.kk_at(m, logk + 1, n).resize(Nblk * Nblk);
                    chain.bb_at(m, logk + 1, n).clear();
                    chain.bb_at(m, logk + 1, n).resize(Nblk);

                    std::vector<AcbValue> solv_kk;  // p=0 residue (kk part)
                    std::vector<AcbValue> solv_bb;  // p=0 residue (bb part)
                    for (long p = logk; p >= 0; --p) {
                        std::vector<AcbValue> res_bb =
                            rec2_residue(true,  m, p, n, beh[m].mu.raw());
                        std::vector<AcbValue> res_kk =
                            rec2_residue(false, m, p, n, beh[m].mu.raw());
                        chain.bb_at(m, p, n) = matvec(pinv, res_bb, Nblk, prec);
                        chain.kk_at(m, p, n) = matmul_flat(pinv, res_kk,
                                                           Nblk, prec);
                        if (p == 0) {
                            solv_bb = std::move(res_bb);
                            solv_kk = std::move(res_kk);
                        }
                    }

                    // Solvability rows on posi: residue_kk[ip,:] . c_m
                    //                          = -residue_bb[ip].
                    // Stored as a PairData with complement_positions = the
                    // NON-TRIVIAL posi rows (so the existing pair_rows
                    // assembly emits them) and kk_logp1/bb_logp1 holding
                    // the RESIDUE, not a log^{logk+1} layer.  Trivially-
                    // satisfied rows (0.c_m = 0) carry no rank and are
                    // dropped here so n_pair_raw stays an honest count.
                    PairData pd;
                    if (solv_kk.size() != (std::size_t)(Nblk * Nblk))
                        solv_kk.resize(Nblk * Nblk);
                    if (solv_bb.size() != (std::size_t)Nblk)
                        solv_bb.resize(Nblk);
                    std::vector<long> solv_emit;
                    for (long ip : posi) {
                        bool nz = !acb_is_chop_zero(solv_bb[ip].raw(),
                                                    chop_pre());
                        for (long k_ = 0; !nz && k_ < (long)Nblk; ++k_) {
                            if (!acb_is_chop_zero(
                                    solv_kk[ip * Nblk + k_].raw(),
                                    chop_pre())) nz = true;
                        }
                        if (nz) solv_emit.push_back(ip);
                    }
                    pd.kk_logp1 = std::move(solv_kk);
                    pd.bb_logp1 = std::move(solv_bb);
                    pd.complement_positions = std::move(solv_emit);
                    pd.resonant_positions   = {};     // none extra
                    pd.at_order_n           = n;
                    pd.populated            = true;
                    pair_data[m].push_back(std::move(pd));

                } else {
                    for (long p = logk; p >= 0; --p) {
                        chain.bb_at(m, p, n) = rec2(true,  m, p, n, inv_mat, beh[m].mu.raw());
                        chain.kk_at(m, p, n) = rec2(false, m, p, n, inv_mat, beh[m].mu.raw());
                    }
                }
            }
        }

        std::vector<long> essentialset;
        for (long m = 0; m < (long)beh.size(); ++m) {
            if (is_essential[m]) essentialset.push_back(m);
        }

        AcbValue x0_log;
        acb_log(x0_log.raw(), x0, prec);

        auto evaluate_chain_to_matrix = [&](long m, bool use_kk) -> std::vector<AcbValue> {
            const long logk = beh[m].log_power;
            std::vector<std::vector<AcbValue>> per_p_value;
            per_p_value.resize(logk + 1);

            for (long p = 0; p <= logk; ++p) {
                long item_size = use_kk ? Nblk * Nblk : Nblk;
                std::vector<AcbValue> val(item_size);
                for (long n = xorder; n >= 0; --n) {
                    for (auto& a : val) {
                        AcbValue tmp_h;
                        acb_mul(tmp_h.raw(), a.raw(), x0, prec);
                        acb_set(a.raw(), tmp_h.raw());
                    }
                    std::vector<AcbValue>& src = use_kk
                        ? chain.kk_at(m, p, n)
                        : chain.bb_at(m, p, n);
                    if (!src.empty()) {
                        if (val.size() != src.size()) { val.clear(); val.resize(src.size()); }
                        vec_add_inplace(val, src, prec);
                    }
                }
                per_p_value[p] = std::move(val);
            }

            std::vector<AcbValue> total(use_kk ? Nblk * Nblk : Nblk);
            AcbValue log_pow;
            log_pow.set_one();
            for (long p = 0; p <= logk; ++p) {
                if (p > 0) acb_mul(log_pow.raw(), log_pow.raw(), x0_log.raw(), prec);
                for (std::size_t i = 0; i < total.size(); ++i) {
                    AcbValue tmp_h;
                    acb_mul(tmp_h.raw(), per_p_value[p][i].raw(), log_pow.raw(), prec);
                    acb_add(total[i].raw(), total[i].raw(), tmp_h.raw(), prec);
                }
            }

            AcbValue x0_spbeh;
            acb_pow(x0_spbeh.raw(), x0, beh[m].mu.raw(), prec);
            for (auto& a : total) {
                AcbValue tmp_h;
                acb_mul(tmp_h.raw(), a.raw(), x0_spbeh.raw(), prec);
                acb_set(a.raw(), tmp_h.raw());
            }
            return total;
        };

        std::vector<std::vector<AcbValue>> finalmat;
        finalmat.reserve(essentialset.size());
        for (long m : essentialset) {
            finalmat.push_back(evaluate_chain_to_matrix(m, /*use_kk=*/true));
        }
        if (std::getenv("AMF_DBG_X00")) {
            // Non-finite scan of the kk-side (bc-INDEPENDENT) matrix.
            for (std::size_t mi = 0; mi < finalmat.size(); ++mi) {
                for (std::size_t i = 0; i < finalmat[mi].size(); ++i) {
                    if (!acb_is_finite(finalmat[mi][i].raw())) {
                        fprintf(stderr,
                                "[X00-kk] block0=%ld: finalmat NON-FINITE at "
                                "ess=%zu entry=%zu (kk chain / x0^mu side)\n",
                                block.empty() ? -1L : (long)block[0], mi, i);
                    }
                }
            }
        }
        std::vector<AcbValue> bc_block(Nblk);
        for (long i = 0; i < (long)Nblk; ++i) {
            acb_set(bc_block[i].raw(), bc[block[i]].raw());
        }
        if (rot != nullptr) {
            bc_block = matvec(rot->invu, bc_block, Nblk, prec);
        }
        std::vector<AcbValue> finalvec = clone_vec(bc_block);
        for (long m = 0; m < (long)beh.size(); ++m) {
            std::vector<AcbValue> bb_eval = evaluate_chain_to_matrix(m, /*use_kk=*/false);
            if (std::getenv("AMF_DBG_X00")) {
                for (long i = 0; i < (long)Nblk; ++i) {
                    if (!acb_is_finite(bb_eval[i].raw())) {
                        char* mu = arb_get_str(acb_realref(beh[m].mu.raw()), 16, 0);
                        char* mui = arb_get_str(acb_imagref(beh[m].mu.raw()), 16, 0);
                        fprintf(stderr,
                                "[X00-nan] block0=%ld m=%ld i=%ld: bb_eval "
                                "NON-FINITE; mu=(%s)+(%s)i log_power=%ld\n",
                                block.empty() ? -1L : (long)block[0], (long)m,
                                (long)i, mu, mui, (long)beh[m].log_power);
                        flint_free(mu);
                        flint_free(mui);
                    }
                }
            }
            for (long i = 0; i < (long)Nblk; ++i) {
                acb_sub(finalvec[i].raw(), finalvec[i].raw(), bb_eval[i].raw(), prec);
            }
        }
        if (std::getenv("AMF_DBG_X00")) {
            for (long i = 0; i < (long)Nblk; ++i) {
                if (!acb_is_finite(finalvec[i].raw())) {
                    fprintf(stderr,
                            "[X00-fv] block0=%ld i=%ld: finalvec NON-FINITE; "
                            "bc_finite=%d beh.size=%zu (rhs side: transported "
                            "bc or bb chain)\n",
                            block.empty() ? -1L : (long)block[0], (long)i,
                            acb_is_finite(bc_block[i].raw()) ? 1 : 0,
                            beh.size());
                }
            }
        }

        long n_ess = static_cast<long>(essentialset.size());
        long n_var = Nblk * n_ess;
        long compensate = 0;

        std::vector<AcbValue> sol_v(n_var);
        if (n_var > 0) {
            acb_mat_t Asys;
            acb_mat_t rhs;
            acb_mat_init(Asys, n_var, n_var);
            acb_mat_init(rhs,  n_var, 1);
            acb_mat_zero(Asys);
            acb_mat_zero(rhs);

            std::vector<std::vector<AcbValue>> match_rows(Nblk);
            std::vector<AcbValue> match_rhs(Nblk);
            for (long row = 0; row < (long)Nblk; ++row) {
                match_rows[row].resize(n_var);
                acb_set(match_rhs[row].raw(), finalvec[row].raw());
                long col_offset = 0;
                for (std::size_t mi = 0; mi < essentialset.size(); ++mi) {
                    for (long c = 0; c < (long)Nblk; ++c) {
                        acb_set(match_rows[row][col_offset + c].raw(),
                                finalmat[mi][row * Nblk + c].raw());
                    }
                    col_offset += Nblk;
                }
            }

            // -----------------------------------------------------------
            //  Pairing rows.
            //
            //  Primary pool: for EACH resonant order n of EACH essential
            //  behavior m, the log^{logk+1} layer must vanish on the
            //  complement (non-resonant) positions -- one row per (m, n,
            //  complement-ip).  Accumulating across ALL resonant n (not
            //  just the last) is what supplies the missing
            //  Nblk*(n_ess-1) constraints in the multi-resonance case.
            //
            //  Secondary (Jordan) pool: the resonant-position rows of the
            //  same log^{logk+1} layer.  For diagonal A_00 these are
            //  identically zero and rank-filtered out, but for a true
            //  Jordan block (off-diagonal coupling) they carry rank and
            //  encode the resonance-induced log term.  They are appended
            //  after the primary pool so that, when the primary pool
            //  alone is rank-deficient, the over-determined fallback
            //  below can draw on them.
            // -----------------------------------------------------------
            std::vector<std::vector<AcbValue>> pair_rows;
            std::vector<AcbValue> pair_rhs;
            std::vector<std::vector<AcbValue>> jordan_rows;
            std::vector<AcbValue> jordan_rhs;
            pair_rows.reserve(n_var);
            pair_rhs.reserve(n_var);
            long col_offset = 0;
            for (std::size_t mi = 0; mi < essentialset.size(); ++mi) {
                long m = essentialset[mi];
                for (const auto& pd : pair_data[m]) {
                    for (long ip : pd.complement_positions) {
                        pair_rows.emplace_back(n_var);
                        pair_rhs.emplace_back();
                        auto& prow = pair_rows.back();
                        for (long k_ = 0; k_ < (long)Nblk; ++k_) {
                            acb_set(prow[col_offset + k_].raw(),
                                    pd.kk_logp1[ip * Nblk + k_].raw());
                        }
                        acb_neg(pair_rhs.back().raw(), pd.bb_logp1[ip].raw());
                    }
                    for (long ip : pd.resonant_positions) {
                        // Row may be all-zero (diagonal A_00) -- rank
                        // selection drops it; nonzero (Jordan) -- it's a
                        // genuine constraint.
                        jordan_rows.emplace_back(n_var);
                        jordan_rhs.emplace_back();
                        auto& jrow = jordan_rows.back();
                        for (long k_ = 0; k_ < (long)Nblk; ++k_) {
                            acb_set(jrow[col_offset + k_].raw(),
                                    pd.kk_logp1[ip * Nblk + k_].raw());
                        }
                        acb_neg(jordan_rhs.back().raw(), pd.bb_logp1[ip].raw());
                    }
                }
                col_offset += Nblk;
            }

            auto rank_of_rows = [&](const std::vector<std::vector<AcbValue>>& rows,
                                    const std::vector<AcbValue>* extra) -> long {
                const long nr = static_cast<long>(rows.size() + (extra ? 1 : 0));
                if (nr == 0 || n_var == 0) return 0;
                acb_mat_t M;
                acb_mat_init(M, nr, n_var);
                for (long r = 0; r < static_cast<long>(rows.size()); ++r) {
                    for (long c = 0; c < n_var; ++c) {
                        acb_set(acb_mat_entry(M, r, c), rows[r][c].raw());
                    }
                }
                if (extra) {
                    const long r = static_cast<long>(rows.size());
                    for (long c = 0; c < n_var; ++c) {
                        acb_set(acb_mat_entry(M, r, c), (*extra)[c].raw());
                    }
                }
                // Rank/pivot tolerance for the boundary-matching solve.  A
                // pivot of order |a*eps| (the genuine coupling produced by an
                // eps-deformed exponent pair) must be KEPT, not declared a
                // zero column: declaring it zero collapses two genuine sectors
                // to one and discards the 1/(a*eps) pole, capping the eps->0
                // extrapolation (calcx00 eps=0 resonance).  The legitimate
                // floor here is the arithmetic noise of acb mid-values at
                // working precision, NOT the coarse coefficient-negligibility
                // chop_pre -- use the resonance tolerance (see
                // resonance_chop_pre).
                const double tol = std::pow(10.0, -numeric::resonance_chop_pre());
                long rank = 0;
                AcbValue tmp_h, inv, factor;
                for (long c = 0; c < n_var && rank < nr; ++c) {
                    long pivot = -1;
                    double best = 0.0;
                    for (long r = rank; r < nr; ++r) {
                        double mag = acb_mid_abs(acb_mat_entry(M, r, c));
                        if (mag > best) {
                            best = mag;
                            pivot = r;
                        }
                    }
                    if (pivot < 0 || !(best > tol)) continue;
                    if (pivot != rank) {
                        for (long j = 0; j < n_var; ++j) {
                            acb_swap(acb_mat_entry(M, rank, j),
                                     acb_mat_entry(M, pivot, j));
                        }
                    }
                    acb_inv(inv.raw(), acb_mat_entry(M, rank, c), prec);
                    for (long j = c; j < n_var; ++j) {
                        acb_mul(acb_mat_entry(M, rank, j),
                                acb_mat_entry(M, rank, j), inv.raw(), prec);
                    }
                    for (long r = 0; r < nr; ++r) {
                        if (r == rank) continue;
                        acb_set(factor.raw(), acb_mat_entry(M, r, c));
                        if (acb_is_chop_zero(factor.raw(), chop_pre())) continue;
                        for (long j = c; j < n_var; ++j) {
                            acb_mul(tmp_h.raw(), factor.raw(),
                                    acb_mat_entry(M, rank, j), prec);
                            acb_sub(acb_mat_entry(M, r, j),
                                    acb_mat_entry(M, r, j), tmp_h.raw(), prec);
                        }
                    }
                    ++rank;
                }
                acb_mat_clear(M);
                return rank;
            };

            // ----------------------------------------------------------------
            //  Resonant-boundary matching count (mirror MMA DESolver.m:1011-1014).
            //
            //  Upstream uses  compensate = Length[allvar] - Length[pairtoequ],
            //  i.e. the RAW COUNT of pairing equations, and then takes the first
            //  `compensate` (highest-LeafCount) match equations so that the total
            //  number of constraints equals exactly Length[allvar] = n_var.
            //
            //  The earlier C++ port instead rank-FILTERED the pairing rows with a
            //  floating `chop_pre` tolerance and set compensate = n_var - RANK.
            //  Because rank <= count, that inflates compensate; near an
            //  integer-spaced (resonant) indicial structure the float rank test
            //  can also drop pairing rows that are exactly independent over the
            //  rationals, which is what pushed a 5-point production top block to
            //  compensate=9 > Nblk=7 and threw.  We restore the MMA count
            //  semantics: keep every raw pairing row and let the exact-rational
            //  Dixon solve below resolve any genuine redundancy.  `n_pair_rank`
            //  is retained purely for diagnostics / the robust fallback.
            const long n_pair_raw = static_cast<long>(pair_rows.size());
            long n_pair_rank = 0;
            {
                std::vector<std::vector<AcbValue>> indep;
                indep.reserve(pair_rows.size());
                long cur = 0;
                for (std::size_t rr = 0; rr < pair_rows.size(); ++rr) {
                    long tr = rank_of_rows(indep, &pair_rows[rr]);
                    if (tr > cur) { indep.push_back(clone_vec(pair_rows[rr])); cur = tr; }
                }
                n_pair_rank = cur;
            }

            auto match_complexity = [&](long row) -> long {
                long count = 0;
                for (long c = 0; c < n_var; ++c) {
                    if (!acb_is_chop_zero(match_rows[row][c].raw(), chop_pre())) {
                        ++count;
                    }
                }
                if (!acb_is_chop_zero(match_rhs[row].raw(), chop_pre())) ++count;
                return count;
            };

            std::vector<long> match_order(Nblk);
            for (long i = 0; i < (long)Nblk; ++i) match_order[i] = i;
            std::sort(match_order.begin(), match_order.end(),
                      [&](long a, long b) {
                          long ca = match_complexity(a);
                          long cb = match_complexity(b);
                          if (ca != cb) return ca > cb;
                          return a < b;
                      });

            auto rank_for_selected = [&](const std::vector<long>& selected) -> long {
                const long rows = static_cast<long>(pair_rows.size() + selected.size());
                if (rows == 0 || n_var == 0) return 0;
                acb_mat_t M;
                acb_mat_init(M, rows, n_var);
                long out_row = 0;
                for (const auto& prow : pair_rows) {
                    for (long c = 0; c < n_var; ++c) {
                        acb_set(acb_mat_entry(M, out_row, c), prow[c].raw());
                    }
                    ++out_row;
                }
                for (long r : selected) {
                    for (long c = 0; c < n_var; ++c) {
                        acb_set(acb_mat_entry(M, out_row, c),
                                match_rows[r][c].raw());
                    }
                    ++out_row;
                }

                // See rank_of_rows above: keep pivots of order |a*eps|, so
                // the rank test sits at the arithmetic noise floor (resonance
                // tolerance), not the coarse chop_pre.
                const double tol = std::pow(10.0, -numeric::resonance_chop_pre());
                long rank = 0;
                AcbValue tmp_h, inv, factor;
                for (long c = 0; c < n_var && rank < rows; ++c) {
                    long pivot = -1;
                    double best = 0.0;
                    for (long r = rank; r < rows; ++r) {
                        double mag = acb_mid_abs(acb_mat_entry(M, r, c));
                        if (mag > best) {
                            best = mag;
                            pivot = r;
                        }
                    }
                    if (pivot < 0 || !(best > tol)) continue;
                    if (pivot != rank) {
                        for (long j = 0; j < n_var; ++j) {
                            acb_swap(acb_mat_entry(M, rank, j),
                                     acb_mat_entry(M, pivot, j));
                        }
                    }
                    acb_inv(inv.raw(), acb_mat_entry(M, rank, c), prec);
                    for (long j = c; j < n_var; ++j) {
                        acb_mul(acb_mat_entry(M, rank, j),
                                acb_mat_entry(M, rank, j), inv.raw(), prec);
                    }
                    for (long r = 0; r < rows; ++r) {
                        if (r == rank) continue;
                        acb_set(factor.raw(), acb_mat_entry(M, r, c));
                        if (acb_is_chop_zero(factor.raw(), chop_pre())) continue;
                        for (long j = c; j < n_var; ++j) {
                            acb_mul(tmp_h.raw(), factor.raw(),
                                    acb_mat_entry(M, rank, j), prec);
                            acb_sub(acb_mat_entry(M, r, j),
                                    acb_mat_entry(M, r, j), tmp_h.raw(), prec);
                        }
                    }
                    ++rank;
                }
                acb_mat_clear(M);
                return rank;
            };

            // MMA count semantics: compensate = n_var - (raw pairing-eqn count).
            // The total constraint set (compensate match rows + all raw pairing
            // rows) then has exactly n_var rows, and the exact-rational Dixon
            // solve below decides solvability (an exactly-redundant raw pairing
            // row makes that square system singular and is reported honestly,
            // rather than being silently dropped by a floating rank test).
            compensate = n_var - n_pair_raw;
            if (compensate < 0) compensate = 0;

            // --------------------------------------------------------------
            //  Row assembly.
            //
            //  Common path (compensate in [0, Nblk] AND n_pair_raw <= n_var):
            //  identical to MMA -- `compensate` highest-LeafCount match rows
            //  + every raw primary pairing row.  Preserves byte-identical
            //  results on every previously-working block.
            //
            //  Resonant fallback (compensate > Nblk, OR n_pair_raw > n_var):
            //  the multi-resonance accumulation may over- or under-shoot
            //  the raw count.  Either way the IVP itself has exactly Nblk
            //  degrees of freedom, so the row-span of {all Nblk match rows,
            //  all primary pairing rows, all Jordan secondary rows} must
            //  have rank n_var for a well-posed block.  Greedily rank-
            //  select n_var independent rows from that pool, preferring
            //  match rows (carry the BC), then primary pairing, then the
            //  Jordan pool (which is nonzero exactly when off-diagonal
            //  A_00 coupling forces a resonance-induced log term).  Only a
            //  genuinely rank-deficient pool now throws.
            // --------------------------------------------------------------
            std::vector<const std::vector<AcbValue>*> sys_rows;
            std::vector<const AcbValue*>              sys_rhs;
            sys_rows.reserve(n_var);
            sys_rhs.reserve(n_var);
            std::vector<long> selected_match;

            const bool use_resonant_fallback =
                (compensate > (long)Nblk) || (n_pair_raw > n_var)
                || block_has_later_resonance;

            if (!use_resonant_fallback) {
                // ---- MMA-count path (unchanged behaviour) ----
                selected_match.reserve(compensate);
                long current_rank = rank_for_selected(selected_match);
                for (long row : match_order) {
                    if ((long)selected_match.size() >= compensate) break;
                    std::vector<long> trial = selected_match;
                    trial.push_back(row);
                    long trial_rank = rank_for_selected(trial);
                    if (trial_rank > current_rank) {
                        selected_match.push_back(row);
                        current_rank = trial_rank;
                    }
                }
                for (long row : match_order) {
                    if ((long)selected_match.size() >= compensate) break;
                    if (std::find(selected_match.begin(), selected_match.end(),
                                  row) == selected_match.end()) {
                        selected_match.push_back(row);
                    }
                }
                for (long src_row : selected_match) {
                    sys_rows.push_back(&match_rows[src_row]);
                    sys_rhs.push_back(&match_rhs[src_row]);
                }
                for (std::size_t r = 0; r < pair_rows.size(); ++r) {
                    sys_rows.push_back(&pair_rows[r]);
                    sys_rhs.push_back(&pair_rhs[r]);
                }
            } else {
                // ---- Resonant over-determined fallback ----
                // Candidate pool, in priority order: primary pairing rows
                // (the structural Frobenius constraints) first, then match
                // rows (the BC), then the Jordan secondary pool.  Match
                // rows fill exactly the rank gap left by the pairing
                // constraints, mirroring the MMA "compensate" semantics
                // but driven by RANK rather than raw count, so redundant
                // / trivially-zero pairing rows cannot poison the square
                // system.
                std::vector<std::pair<const std::vector<AcbValue>*,
                                      const AcbValue*>> pool;
                pool.reserve(Nblk + pair_rows.size() + jordan_rows.size());
                for (std::size_t r = 0; r < pair_rows.size(); ++r) {
                    pool.emplace_back(&pair_rows[r], &pair_rhs[r]);
                }
                for (long row : match_order) {
                    pool.emplace_back(&match_rows[row], &match_rhs[row]);
                }
                for (std::size_t r = 0; r < jordan_rows.size(); ++r) {
                    pool.emplace_back(&jordan_rows[r], &jordan_rhs[r]);
                }

                // Greedy rank selection (incremental Gaussian elimination
                // on a working copy; same tolerance as rank_of_rows).
                std::vector<std::vector<AcbValue>> echelon;
                echelon.reserve(n_var);
                // Same rationale as rank_of_rows: a genuine |a*eps|-scale
                // pivot must survive rank selection, so use the resonance
                // (arithmetic-floor) tolerance, not chop_pre.
                const double tol = std::pow(10.0, -numeric::resonance_chop_pre());
                AcbValue tmp_h, inv, factor;
                for (auto& cand : pool) {
                    if ((long)sys_rows.size() >= n_var) break;
                    // Reduce candidate against current echelon basis.
                    std::vector<AcbValue> red = clone_vec(*cand.first);
                    for (const auto& bas : echelon) {
                        // Find pivot column of bas (first non-chop-zero).
                        long pc = -1;
                        for (long c = 0; c < n_var; ++c) {
                            if (acb_mid_abs(bas[c].raw()) > tol) { pc = c; break; }
                        }
                        if (pc < 0) continue;
                        acb_div(factor.raw(), red[pc].raw(), bas[pc].raw(), prec);
                        for (long c = pc; c < n_var; ++c) {
                            acb_mul(tmp_h.raw(), factor.raw(), bas[c].raw(), prec);
                            acb_sub(red[c].raw(), red[c].raw(), tmp_h.raw(), prec);
                        }
                    }
                    // Is the reduced row non-trivial?
                    bool indep = false;
                    for (long c = 0; c < n_var; ++c) {
                        if (acb_mid_abs(red[c].raw()) > tol) { indep = true; break; }
                    }
                    if (indep) {
                        echelon.push_back(std::move(red));
                        sys_rows.push_back(cand.first);
                        sys_rhs.push_back(cand.second);
                    }
                }
                selected_match.clear();
                for (long row : match_order) {
                    if (std::find(sys_rows.begin(), sys_rows.end(),
                                  &match_rows[row]) != sys_rows.end()) {
                        selected_match.push_back(row);
                    }
                }

                if ((long)sys_rows.size() < n_var) {
                    // Genuinely rank-deficient even with the Jordan pool.
                    std::ostringstream msg;
                    msg << "calcx00: boundary matching rank-deficient after "
                        << "resonant fallback; block=[";
                    for (std::size_t bi = 0; bi < block.size(); ++bi) {
                        if (bi) msg << ",";
                        msg << block[bi];
                    }
                    msg << "]; behaviors=" << beh.size()
                        << "; n_var=" << n_var
                        << "; n_pair_raw=" << n_pair_raw
                        << "; n_pair_rank=" << n_pair_rank
                        << "; n_jordan=" << jordan_rows.size()
                        << "; pool_rank=" << sys_rows.size()
                        << "; Nblk=" << Nblk;
                    acb_mat_clear(Asys);
                    acb_mat_clear(rhs);
                    throw std::runtime_error(msg.str());
                }
            }

            if (std::getenv("AMF_DBG_X00")) {
                std::vector<long> allm; for (long i=0;i<(long)Nblk;++i) allm.push_back(i);
                long rank_pair_allmatch = rank_for_selected(allm);
                fprintf(stderr,
                        "[X00] block0=%ld Nblk=%ld n_ess=%ld n_var=%ld "
                        "n_pair_raw=%ld n_pair_rank=%ld n_jordan=%zu "
                        "rank(pair+allmatch)=%ld compensate=%ld "
                        "selected_match=%zu fallback=%d\n",
                        block.empty()?-1L:(long)block[0], (long)Nblk, n_ess,
                        n_var, n_pair_raw, n_pair_rank, jordan_rows.size(),
                        rank_pair_allmatch, compensate, selected_match.size(),
                        use_resonant_fallback ? 1 : 0);
            }

            long out_row = 0;
            for (std::size_t r = 0; r < sys_rows.size(); ++r) {
                acb_set(acb_mat_entry(rhs, out_row, 0), sys_rhs[r]->raw());
                for (long c = 0; c < n_var; ++c) {
                    acb_set(acb_mat_entry(Asys, out_row, c),
                            (*sys_rows[r])[c].raw());
                }
                ++out_row;
            }

            // Solve via rationalize -> exact rational Dixon -> arb conversion
            // (mirrors MMA DESolver.m:906 Solve, exact over rationals).
            acb_mat_t sol_mat;
            acb_mat_init(sol_mat, n_var, 1);

            // Rationalization tolerance mirrors MMA's RationalizePre
            // (DESolver.m Solve), NOT the full working precision: rounding the
            // boundary-matching matrix at working_pre digits demands that every
            // acb interval be tighter than 10^-working_pre, which a resonant
            // leaf at eps=0 (slowly-converging indicial series) need not reach,
            // producing the spurious "interval too wide" abort.  The rest of the
            // codebase rationalizes at rationalize_pre (see ode/path.cpp,
            // ode/inf.cpp, numeric/matrix.cpp); this block was the lone
            // deviation.  Cap at what the working precision can actually resolve
            // (same guard as acb_real_to_fmpq) so we never fabricate digits.
            int ratdig = numeric::rationalize_pre();
            {
                long resolvable =
                    static_cast<long>(working_prec_bits() * 0.301029995663981195L);
                if (resolvable >= 5) resolvable -= 5;
                if (resolvable < 1) resolvable = 1;
                if (ratdig > resolvable) ratdig = static_cast<int>(resolvable);
                if (ratdig < 1) ratdig = 1;
            }
            long rat_bits = decimal_digits_to_bits(ratdig) + 32;

            fmpq_mat_t Mblock, Bblock, Sblock;
            fmpq_mat_init(Mblock, 2 * n_var, 2 * n_var);
            fmpq_mat_init(Bblock, 2 * n_var, 1);
            fmpq_mat_init(Sblock, 2 * n_var, 1);

            bool rationalize_ok = true;
            {
                arb_t scale;
                arb_init(scale);
                arb_set_si(scale, 10);
                arb_pow_ui(scale, scale,
                           static_cast<unsigned long>(ratdig),
                           rat_bits);

                fmpz_t pow10;
                fmpz_init(pow10);
                fmpz_set_si(pow10, 10);
                fmpz_pow_ui(pow10, pow10,
                            static_cast<unsigned long>(ratdig));

                arb_t scaled, part;
                arb_init(scaled);
                arb_init(part);
                fmpz_t numer;
                fmpz_init(numer);

                auto rationalize_acb = [&](acb_srcptr src,
                                           fmpq_t out_re,
                                           fmpq_t out_im) -> bool {
                    arb_set(part, acb_realref(src));
                    arb_mul(scaled, part, scale, rat_bits);
                    if (!arb_is_finite(scaled)) {
                        if (std::getenv("AMF_DBG_X00")) {
                            char* s = arb_get_str(acb_realref(src), 12,
                                                  ARB_STR_NO_RADIUS);
                            fprintf(stderr,
                                    "[X00-rat] non-finite RE: src.re=%s "
                                    "rat_bits=%ld ratdig=%d\n",
                                    s, rat_bits, ratdig);
                            flint_free(s);
                        }
                        return false;
                    }
                    arf_get_fmpz(numer, arb_midref(scaled), ARF_RND_NEAR);
                    fmpq_set_fmpz_frac(out_re, numer, pow10);

                    arb_set(part, acb_imagref(src));
                    arb_mul(scaled, part, scale, rat_bits);
                    if (!arb_is_finite(scaled)) {
                        if (std::getenv("AMF_DBG_X00")) {
                            char* s = arb_get_str(acb_imagref(src), 12,
                                                  ARB_STR_NO_RADIUS);
                            fprintf(stderr,
                                    "[X00-rat] non-finite IM: src.im=%s "
                                    "rat_bits=%ld ratdig=%d\n",
                                    s, rat_bits, ratdig);
                            flint_free(s);
                        }
                        return false;
                    }
                    arf_get_fmpz(numer, arb_midref(scaled), ARF_RND_NEAR);
                    fmpq_set_fmpz_frac(out_im, numer, pow10);
                    return true;
                };

                fmpq_t a_re, a_im;
                fmpq_init(a_re);
                fmpq_init(a_im);
                for (long r = 0; r < n_var && rationalize_ok; ++r) {
                    for (long c = 0; c < n_var; ++c) {
                        if (!rationalize_acb(acb_mat_entry(Asys, r, c),
                                             a_re, a_im)) {
                            rationalize_ok = false;
                            break;
                        }
                        fmpq_set(fmpq_mat_entry(Mblock, r, c), a_re);
                        fmpq_neg(fmpq_mat_entry(Mblock, r, c + n_var), a_im);
                        fmpq_set(fmpq_mat_entry(Mblock, r + n_var, c), a_im);
                        fmpq_set(fmpq_mat_entry(Mblock, r + n_var, c + n_var),
                                 a_re);
                    }
                    if (!rationalize_ok) break;
                    if (!rationalize_acb(acb_mat_entry(rhs, r, 0),
                                         a_re, a_im)) {
                        rationalize_ok = false;
                        break;
                    }
                    fmpq_set(fmpq_mat_entry(Bblock, r, 0), a_re);
                    fmpq_set(fmpq_mat_entry(Bblock, r + n_var, 0), a_im);
                }
                fmpq_clear(a_re);
                fmpq_clear(a_im);
                fmpz_clear(numer);
                fmpz_clear(pow10);
                arb_clear(scaled);
                arb_clear(part);
                arb_clear(scale);
            }

            int ok = 0;
            if (rationalize_ok) {
                ok = fmpq_mat_solve_dixon(Sblock, Mblock, Bblock);
                if (ok) {
                    acb_t cur;
                    acb_init(cur);
                    for (long i = 0; i < n_var; ++i) {
                        arb_set_fmpq(acb_realref(cur),
                                     fmpq_mat_entry(Sblock, i, 0),
                                     prec);
                        arb_set_fmpq(acb_imagref(cur),
                                     fmpq_mat_entry(Sblock, i + n_var, 0),
                                     prec);
                        acb_set(acb_mat_entry(sol_mat, i, 0), cur);
                    }
                    acb_clear(cur);
                }
            }
            fmpq_mat_clear(Mblock);
            fmpq_mat_clear(Bblock);
            fmpq_mat_clear(Sblock);
            if (!ok) {
                // Distinguish a genuinely non-finite (NaN/Inf) entry -- a
                // DIVERGENT value reached the matching solve; structural, no
                // precision increase can clear it -- from an actual
                // rationalization failure (the calcx00 audit honest-error
                // split).
                bool any_nonfinite = false;
                long nf_r = -1, nf_c = -1;   // -1 col => rhs
                for (long r = 0; r < n_var && !any_nonfinite; ++r) {
                    if (!acb_is_finite(acb_mat_entry(rhs, r, 0))) {
                        any_nonfinite = true; nf_r = r; nf_c = -1;
                    }
                    for (long c = 0; c < n_var && !any_nonfinite; ++c) {
                        if (!acb_is_finite(acb_mat_entry(Asys, r, c))) {
                            any_nonfinite = true; nf_r = r; nf_c = c;
                        }
                    }
                }
                if (std::getenv("AMF_DBG_X00")) {
                    fprintf(stderr,
                            "[X00-fail] block0=%ld rationalize_ok=%d n_var=%ld "
                            "prec=%ld ratdig=%d; entries:\n",
                            block.empty() ? -1L : (long)block[0],
                            rationalize_ok ? 1 : 0, n_var, prec, ratdig);
                    auto dump = [&](acb_srcptr z, const char* tag,
                                    long i, long j) {
                        char* re = arb_get_str(acb_realref(z), 12, 0);
                        char* im = arb_get_str(acb_imagref(z), 12, 0);
                        fprintf(stderr, "  %s[%ld][%ld]= %s + %s*I\n",
                                tag, i, j, re, im);
                        flint_free(re);
                        flint_free(im);
                    };
                    const long rmax = n_var <= 8 ? n_var : 8;
                    for (long r = 0; r < rmax; ++r) {
                        for (long c = 0; c < rmax; ++c) {
                            dump(acb_mat_entry(Asys, r, c), "A", r, c);
                        }
                        dump(acb_mat_entry(rhs, r, 0), "b", r, 0);
                    }
                }
                std::ostringstream msg;
                if (any_nonfinite) {
                    msg << "calcx00: boundary matching system has a non-finite "
                           "(NaN/Inf) entry at "
                        << (nf_c < 0 ? "rhs row " : "matrix row ")
                        << nf_r;
                    if (nf_c >= 0) msg << " col " << nf_c;
                    msg << " (block0="
                        << (block.empty() ? -1L : (long)block[0])
                        << ") -- a divergent value reached the matching "
                           "solve; this is structural and no working-"
                           "precision increase can clear it";
                } else {
                    msg << "calcx00: rationalized boundary matching system "
                        << (rationalize_ok
                            ? "singular under exact rational solve"
                            : "could not be rationalized (interval too wide?)");
                }
                acb_mat_clear(Asys); acb_mat_clear(rhs); acb_mat_clear(sol_mat);
                throw std::runtime_error(msg.str());
            }
            for (long i = 0; i < n_var; ++i) {
                acb_set(sol_v[i].raw(), acb_mat_entry(sol_mat, i, 0));
            }
            acb_mat_clear(Asys); acb_mat_clear(rhs); acb_mat_clear(sol_mat);

            // -----------------------------------------------------------
            //  Jordan/log-power consistency check + auto-rerun trigger.
            //
            //  With sol_v fixed, the log^{logk+1} layer
            //    kk[m][logk+1][n] . sol_v[m] + bb[m][logk+1][n]
            //  must vanish at every (first-)resonant order n -- on the
            //  complement positions by construction (those were the
            //  primary pairing rows), and on the resonant positions iff
            //  the supplied log_power was high enough.  A nonzero residue
            //  at a resonant position means two indicial exponents
            //  differing by an integer have FORCED a log term beyond
            //  beh[m].log_power (the textbook Frobenius "second
            //  solution").  Bump that log_power and re-run the block.
            //
            //  Likewise, a later-resonance solvability row whose kk-part
            //  is zero but whose bb-part is not (so it could not enter
            //  the linear system) is the SAME obstruction at a later n.
            // -----------------------------------------------------------
            {
                std::vector<bool> bump(beh.size(), false);
                AcbValue residue, term;
                for (std::size_t mi = 0; mi < essentialset.size(); ++mi) {
                    const long m = essentialset[mi];
                    for (const auto& pd : pair_data[m]) {
                        // First-resonance log^{logk+1} residue at posi.
                        for (long ip : pd.resonant_positions) {
                            acb_set(residue.raw(), pd.bb_logp1[ip].raw());
                            for (long k_ = 0; k_ < (long)Nblk; ++k_) {
                                acb_mul(term.raw(),
                                        pd.kk_logp1[ip * Nblk + k_].raw(),
                                        sol_v[mi * Nblk + k_].raw(), prec);
                                acb_add(residue.raw(), residue.raw(),
                                        term.raw(), prec);
                            }
                            if (!acb_is_chop_zero(residue.raw(), chop_pre())) {
                                bump[m] = true;
                            }
                        }
                    }
                }
                // The later-resonance "0.c_m = nonzero" obstruction (a
                // solvability row whose kk-part is identically zero but
                // whose bb-part is not) is the same forced-log signal at
                // a higher n.  It is caught on the NEXT attempt by the
                // first-resonance check above once that behavior's
                // log_power has been raised, so no separate scan is run
                // here.

                bool any_bump = false;
                for (std::size_t m = 0; m < beh.size(); ++m) {
                    if (bump[m]) {
                        any_bump = true;
                        ++beh[m].log_power;
                    }
                }
                if (any_bump && log_attempt < calcx00_max_log_bump) {
                    if (!silent_mode() || std::getenv("AMF_DBG_X00")) {
                        std::ostringstream which;
                        for (std::size_t m = 0; m < beh.size(); ++m)
                            if (bump[m]) which << " m=" << m
                                               << "->logk=" << beh[m].log_power;
                        fprintf(stderr,
                            "[amflow][calcx00] block0=%ld attempt %d: "
                            "resonance-forced log term detected; bumping "
                            "log_power and re-running block (%s ).\n",
                            block.empty()?-1L:(long)block[0],
                            log_attempt, which.str().c_str());
                    }
                    continue;  // redo the whole block at higher log_power
                }
                if (any_bump) {
                    fprintf(stderr,
                        "[amflow][calcx00] WARNING: block0=%ld still has "
                        "nonzero log^{logk+1} residue after %d bumps -- "
                        "output may drop a high log power.\n",
                        block.empty()?-1L:(long)block[0],
                        calcx00_max_log_bump);
                }
            }
        }

        // Compute f0[block, beh[m].mu, p, n] = kk[m, p, n] . variables[m] + bb[m, p, n]
        for (long m = 0; m < (long)beh.size(); ++m) {
            std::vector<AcbValue> var_m(Nblk);
            auto it = std::find(essentialset.begin(), essentialset.end(), m);
            if (it != essentialset.end()) {
                long mi = static_cast<long>(it - essentialset.begin());
                for (long k_ = 0; k_ < (long)Nblk; ++k_) {
                    acb_set(var_m[k_].raw(), sol_v[mi * Nblk + k_].raw());
                }
            }

            const long logk = beh[m].log_power;
            std::vector<std::vector<std::vector<AcbValue>>> per_integral(Nblk);
            for (long i = 0; i < (long)Nblk; ++i) {
                per_integral[i].resize(logk + 1);
                for (long p = 0; p <= logk; ++p) {
                    per_integral[i][p].resize(xorder + 1);
                }
            }

            for (long p = 0; p <= logk; ++p) {
                for (long n = 0; n <= xorder; ++n) {
                    auto& kk_pn = chain.kk_at(m, p, n);
                    auto& bb_pn = chain.bb_at(m, p, n);
                    if (kk_pn.empty()) kk_pn.resize(Nblk * Nblk);
                    if (bb_pn.empty()) bb_pn.resize(Nblk);
                    std::vector<AcbValue> kv = matvec(kk_pn, var_m, Nblk, prec);
                    std::vector<AcbValue> coeff_vec(Nblk);
                    for (long i = 0; i < (long)Nblk; ++i) {
                        acb_add(coeff_vec[i].raw(), kv[i].raw(), bb_pn[i].raw(), prec);
                    }
                    if (rot != nullptr) {
                        coeff_vec = matvec(rot->u, coeff_vec, Nblk, prec);
                    }
                    for (long i = 0; i < (long)Nblk; ++i) {
                        acb_set(per_integral[i][p][n].raw(), coeff_vec[i].raw());
                    }
                }
            }

            for (long i = 0; i < (long)Nblk; ++i) {
                AsyTerm term;
                term.mu = beh[m].mu.clone();
                term.exp.resize(logk + 1);
                for (long p = 0; p <= logk; ++p) {
                    term.exp[p] = clone_vec(per_integral[i][p]);
                }
                result[block[i]].push_back(std::move(term));
            }
        }
        break;  // success: exit the Jordan/log auto-rerun loop
        }       // end for(log_attempt)
    }

    return result;
}

std::vector<AsyExpansion>
calcx00(const std::vector<BlockEquation>&     nheq,
        const std::vector<BlockEquationNum>&  nheqn,
        const std::vector<AcbValue>&          bc,
        acb_srcptr                            x0,
        const AsymptoticBehaviorList&         behavior,
        long                                  prec) {
    return calcx00_impl(nheq, nheqn, bc, x0, behavior, nullptr, prec);
}

// ===========================================================================
//  Learn / CalcZero / PickZero
// ===========================================================================

BlockBehavior learn_from_rule_s(const AsyExpansion& asy) {
    int test = test_x_order();
    int learn = learn_x_order();
    long upto = std::min<long>(test, learn) + 1;
    int chop = chop_pre();

    BlockBehavior out;
    for (const auto& term : asy) {
        long k = static_cast<long>(term.exp.size());
        while (k >= 1) {
            bool any = false;
            const auto& row = term.exp[k - 1];
            for (long n = 0; n < (long)row.size() && n < upto; ++n) {
                if (!acb_is_chop_zero(row[n].raw(), chop)) { any = true; break; }
            }
            if (!any) --k;
            else break;
        }
        long log_power = k - 1;
        if (log_power >= 0) {
            BehaviorEntry e;
            e.mu = term.mu.clone();
            e.log_power = log_power;
            out.push_back(std::move(e));
        }
    }
    return out;
}

AsymptoticBehaviorList
learn_from_rule_s_all(const std::vector<AsyExpansion>& asy_list,
                      const std::vector<std::vector<std::size_t>>& blocks) {
    AsymptoticBehaviorList out(blocks.size());
    long prec = working_prec_bits();

    for (std::size_t bi = 0; bi < blocks.size(); ++bi) {
        BlockBehavior joined;
        for (auto idx : blocks[bi]) {
            BlockBehavior contrib = learn_from_rule_s(asy_list[idx]);
            for (auto& e : contrib) joined.push_back(BehaviorEntry{e.mu.clone(), e.log_power});
        }
        BlockBehavior unioned;
        std::vector<bool> done(joined.size(), false);
        for (std::size_t i = 0; i < joined.size(); ++i) {
            if (done[i]) continue;
            BehaviorEntry combined;
            combined.mu = joined[i].mu.clone();
            combined.log_power = joined[i].log_power;
            done[i] = true;
            for (std::size_t j = i + 1; j < joined.size(); ++j) {
                if (done[j]) continue;
                AcbValue diff;
                acb_sub(diff.raw(), joined[i].mu.raw(), joined[j].mu.raw(), prec);
                // Exponent-equality merge -> resonance tolerance.
                if (acb_is_chop_zero(diff.raw(), numeric::resonance_chop_pre())) {
                    if (joined[j].log_power > combined.log_power)
                        combined.log_power = joined[j].log_power;
                    done[j] = true;
                }
            }
            unioned.push_back(std::move(combined));
        }
        out[bi] = std::move(unioned);
    }
    return out;
}

std::vector<AsyExpansion>
calc_zero(const RationalMatrix& de,
          const std::vector<AcbValue>& bc,
          acb_srcptr x0,
          long prec) {
    if (!silent_mode()) log_line("CalcZero: starting");

    NormalizationResult norm;
    std::vector<BlockJordanRotation> rotations;
    bool use_block_rotations = false;
    bool use_acb_shearing = false;
    try {
        norm = normalize_mat(de, prec);
    } catch (const std::exception& e) {
        const std::string msg = e.what();
        if (msg.find("algebraic Jordan rotation") == std::string::npos &&
            msg.find("algebraic leading residue") == std::string::npos &&
            msg.find("algebraic-number shearing") == std::string::npos) {
            throw;
        }
        if (!silent_mode()) {
            log_line("CalcZero: falling back to block-local algebraic Jordan rotations");
        }
        auto fallback = normalize_mat_for_calc_zero(de, prec);
        norm = std::move(fallback.base);
        rotations = std::move(fallback.rotations);
        use_block_rotations = true;
        use_acb_shearing =
            (msg.find("algebraic-number shearing") != std::string::npos);
        if (use_acb_shearing) {
            return calc_zero_acb_shearing_single_block(norm, bc, x0, prec);
        }
    }

    std::vector<AcbValue> bcT(de.rows());
    for (std::size_t i = 0; i < de.rows(); ++i) {
        AcbValue acc;
        for (std::size_t j = 0; j < de.cols(); ++j) {
            std::ostringstream where;
            where << "calc_zero: invT(x0) entry=(" << i << "," << j << ")";
            AcbValue tij = evaluate_with_context(norm.invT(i, j), x0,
                                                  prec, where.str());
            AcbValue tmp;
            acb_mul(tmp.raw(), tij.raw(), bc[j].raw(), prec);
            acb_add(acc.raw(), acc.raw(), tmp.raw(), prec);
        }
        bcT[i] = std::move(acc);
    }

    auto nheq  = nh_equations(norm.B, EquationMode::Singular);
    auto nheqn = nh_equations_num(nheq, prec);

    AsymptoticBehaviorList behavior = use_block_rotations
        ? asymptotic_behavior_for_calc_zero(norm.B, rotations, prec)
        : asymptotic_behavior(norm.B, prec);

    if (learn_x_order() >= 0) {
        numeric::GlobalScope guard;
        guard.expansion.x_order = learn_x_order();
        guard.global.silent_mode = true;
        guard.commit();
        auto trial = use_block_rotations
            ? calcx00_impl(nheq, nheqn, bcT, x0, behavior, &rotations, prec)
            : calcx00_impl(nheq, nheqn, bcT, x0, behavior, nullptr, prec);
        auto blocks = analyze_block(norm.B);
        behavior = learn_from_rule_s_all(trial, blocks);
    }

    auto raw = use_block_rotations
        ? calcx00_impl(nheq, nheqn, bcT, x0, behavior, &rotations, prec)
        : calcx00_impl(nheq, nheqn, bcT, x0, behavior, nullptr, prec);

    auto T_ps = to_power_series_matrix(norm.T, prec);

    auto out = ps_map_rule_set(T_ps, raw, prec);
    if (!silent_mode()) log_line("CalcZero: finished");
    return out;
}

// ===========================================================================
//  calc_zero_multi  (fp-probe column batching)
// ===========================================================================
//
//  Hoists every bc-independent piece of calc_zero out of the per-column
//  loop: normalize_mat (the NormalizeMat / ToFuchsian elephant — one call
//  for ALL columns), the invT(x0) numeric evaluation, nh_equations(+num) of
//  norm.B, asymptotic_behavior, analyze_block(norm.B) and
//  to_power_series_matrix(norm.T).  The per-column body (bcT rotation,
//  learn pass, calcx00_impl, ps_map_rule_set) is the exact operation
//  sequence of calc_zero on the same inputs.

std::vector<std::vector<AsyExpansion>>
calc_zero_multi(const RationalMatrix& de,
                const std::vector<std::vector<AcbValue>>& bc_cols,
                acb_srcptr x0,
                long prec) {
    const std::size_t ncols = bc_cols.size();
    std::vector<std::vector<AsyExpansion>> out_cols(ncols);
    if (ncols == 0) return out_cols;
    for (const auto& b : bc_cols) {
        if (b.size() != de.rows()) {
            throw std::invalid_argument(
                "calc_zero_multi: bc column size does not match de.rows()");
        }
    }

    if (!silent_mode()) {
        log_line("CalcZero: starting (batched "
                 + std::to_string(ncols) + " cols)");
    }

    // ---- calcx00 audit diagnostics, env-gated ----------
    // AMF_DBG_X00: scan the incoming bc columns for non-finite entries --
    // says whether a NaN/Inf ARRIVED with the transported bc (upstream
    // poisoning) or is created below (CalcZero-internal).
    if (std::getenv("AMF_DBG_X00")) {
        long total = 0, shown = 0;
        for (std::size_t c = 0; c < ncols; ++c) {
            for (std::size_t i = 0; i < bc_cols[c].size(); ++i) {
                if (!acb_is_finite(bc_cols[c][i].raw())) {
                    ++total;
                    if (shown < 8) {
                        fprintf(stderr, "[X00-scan] calc_zero_multi entry: "
                                "NON-FINITE bc at col=%zu master=%zu\n", c, i);
                        ++shown;
                    }
                }
            }
        }
        fprintf(stderr, "[X00-scan] calc_zero_multi entry: %ld non-finite "
                "bc entries total (%zu cols x %zu masters)\n",
                total, ncols, bc_cols.empty() ? (std::size_t)0
                                              : bc_cols[0].size());
    }
    // AMFLOW_X00_DUMP_BC=path: exact checkpoint of the transported bc so a
    // later run can re-enter at CalcZero (AMFLOW_X00_BC_FILE) without
    // replaying ~1 h of transport.  arb_dump_str is exact (hex mantissa).
    if (const char* dump_path = std::getenv("AMFLOW_X00_DUMP_BC")) {
        std::FILE* f = std::fopen(dump_path, "w");
        if (f) {
            fprintf(f, "X00BC v1 %zu %zu\n", ncols, bc_cols[0].size());
            for (std::size_t c = 0; c < ncols; ++c) {
                for (std::size_t i = 0; i < bc_cols[c].size(); ++i) {
                    char* re = arb_dump_str(acb_realref(bc_cols[c][i].raw()));
                    char* im = arb_dump_str(acb_imagref(bc_cols[c][i].raw()));
                    fprintf(f, "%s\n%s\n", re, im);
                    flint_free(re);
                    flint_free(im);
                }
            }
            std::fclose(f);
            log_line(std::string("CalcZero: [X00_DUMP_BC] checkpoint written: ")
                     + dump_path);
        } else {
            fprintf(stderr, "[X00-dump] cannot open %s for writing\n",
                    dump_path);
        }
    }
    // ---------------------------------------------------------------------

    NormalizationResult norm;
    std::vector<BlockJordanRotation> rotations;
    bool use_block_rotations = false;
    bool use_acb_shearing = false;
    try {
        norm = normalize_mat(de, prec);
    } catch (const std::exception& e) {
        const std::string msg = e.what();
        if (msg.find("algebraic Jordan rotation") == std::string::npos &&
            msg.find("algebraic leading residue") == std::string::npos &&
            msg.find("algebraic-number shearing") == std::string::npos) {
            throw;
        }
        if (!silent_mode()) {
            log_line("CalcZero: falling back to block-local algebraic Jordan rotations");
        }
        auto fallback = normalize_mat_for_calc_zero(de, prec);
        norm = std::move(fallback.base);
        rotations = std::move(fallback.rotations);
        use_block_rotations = true;
        use_acb_shearing =
            (msg.find("algebraic-number shearing") != std::string::npos);
    }

    if (use_acb_shearing) {
        for (std::size_t c = 0; c < ncols; ++c) {
            out_cols[c] =
                calc_zero_acb_shearing_single_block(norm, bc_cols[c], x0, prec);
        }
        return out_cols;
    }

    // invT(x0): evaluated ONCE (bc-independent).
    const std::size_t N = de.rows();
    std::vector<std::vector<AcbValue>> invT_x0(N);
    for (std::size_t i = 0; i < N; ++i) {
        invT_x0[i].resize(N);
        for (std::size_t j = 0; j < N; ++j) {
            std::ostringstream where;
            where << "calc_zero: invT(x0) entry=(" << i << "," << j << ")";
            invT_x0[i][j] = evaluate_with_context(norm.invT(i, j), x0,
                                                   prec, where.str());
        }
    }

    auto nheq  = nh_equations(norm.B, EquationMode::Singular);
    auto nheqn = nh_equations_num(nheq, prec);

    const AsymptoticBehaviorList behavior0 = use_block_rotations
        ? asymptotic_behavior_for_calc_zero(norm.B, rotations, prec)
        : asymptotic_behavior(norm.B, prec);

    const bool learn = (learn_x_order() >= 0);
    std::vector<std::vector<std::size_t>> learn_blocks;
    if (learn) learn_blocks = analyze_block(norm.B);

    auto T_ps = to_power_series_matrix(norm.T, prec);

    const std::vector<BlockJordanRotation>* rot_ptr =
        use_block_rotations ? &rotations : nullptr;

    for (std::size_t c = 0; c < ncols; ++c) {
        // bcT = invT(x0) . bc  (same product/sum order as calc_zero).
        std::vector<AcbValue> bcT(N);
        for (std::size_t i = 0; i < N; ++i) {
            AcbValue acc;
            for (std::size_t j = 0; j < N; ++j) {
                AcbValue tmp;
                acb_mul(tmp.raw(), invT_x0[i][j].raw(),
                        bc_cols[c][j].raw(), prec);
                acb_add(acc.raw(), acc.raw(), tmp.raw(), prec);
            }
            bcT[i] = std::move(acc);
        }

        AsymptoticBehaviorList behavior_learned;
        const AsymptoticBehaviorList* behavior_use = &behavior0;
        if (learn) {
            numeric::GlobalScope guard;
            guard.expansion.x_order = learn_x_order();
            guard.global.silent_mode = true;
            guard.commit();
            auto trial = calcx00_impl(nheq, nheqn, bcT, x0, behavior0,
                                      rot_ptr, prec);
            behavior_learned = learn_from_rule_s_all(trial, learn_blocks);
            behavior_use = &behavior_learned;
        }

        auto raw = calcx00_impl(nheq, nheqn, bcT, x0, *behavior_use,
                                rot_ptr, prec);

        out_cols[c] = ps_map_rule_set(T_ps, raw, prec);
        if (!silent_mode()) {
            log_line("CalcZero: batched col "
                     + std::to_string(c + 1) + "/" + std::to_string(ncols)
                     + " done");
        }
    }

    if (!silent_mode()) log_line("CalcZero: finished");
    return out_cols;
}

AcbValue pick_zero_rule_s(const AsyExpansion& asy, long prec) {
    std::vector<long> integer_mus;
    std::vector<std::size_t> integer_idx;

    // Integer-exponent detection must use the resonance tolerance, not the
    // coarse chop_pre: an eps-deformed exponent mu = a*eps is NOT an integer
    // for eps != 0, but |a*eps| < 10^-chop_pre would round to the integer 0 and
    // be (wrongly) counted alongside the genuine integer sector, throwing the
    // "multiple integer mu" error and aliasing two distinct sectors.  The
    // resonance tolerance (~0.8*working_pre) sits at the arithmetic noise
    // floor, so only a genuine integer (residual identically 0) is selected.
    int chop = numeric::resonance_chop_pre();
    if (chop <= 0) chop = 20;
    const double eps_d = std::pow(10.0, -chop);

    for (std::size_t i = 0; i < asy.size(); ++i) {
        arf_t mid_im;
        arf_init(mid_im);
        arf_abs(mid_im, arb_midref(acb_imagref(asy[i].mu.raw())));
        bool im_zero = (arf_cmp_d(mid_im, eps_d) < 0);
        arf_clear(mid_im);
        if (!im_zero) continue;

        fmpz_t z;
        fmpz_init(z);
        arf_get_fmpz(z, arb_midref(acb_realref(asy[i].mu.raw())), ARF_RND_NEAR);

        arb_t resid;
        arb_init(resid);
        arb_set_fmpz(resid, z);
        arb_sub(resid, acb_realref(asy[i].mu.raw()), resid, prec);
        arf_t mr;
        arf_init(mr);
        arf_abs(mr, arb_midref(resid));
        bool re_int = (arf_cmp_d(mr, eps_d) < 0);
        arf_clear(mr); arb_clear(resid);

        if (re_int) {
            integer_mus.push_back(fmpz_get_si(z));
            integer_idx.push_back(i);
        }
        fmpz_clear(z);
    }

    if (integer_mus.empty()) return AcbValue();
    if (integer_mus.size() > 1) {
        throw std::runtime_error(
            "pick_zero_rule_s: multiple integer mu values found; "
            "the asymptotic regions have not been fully merged");
    }

    long key = integer_mus[0];
    std::size_t which = integer_idx[0];
    if (key > 0) return AcbValue();

    if (key == 0) {
        const auto& exp = asy[which].exp;
        if (exp.empty() || exp[0].empty()) return AcbValue();
        return exp[0][0].clone();
    }

    const auto& exp = asy[which].exp;
    if (exp.empty()) return AcbValue();
    long want = -key;
    if (want >= (long)exp[0].size()) return AcbValue();
    return exp[0][want].clone();
}

}  // namespace amflow::ode
