// SPDX-License-Identifier: MIT
// qft::complete — ToCompleteExplicit + SPListToDListSymbol +
// SquaredDenominators + ToSquareAll.
//

#include "amflow/qft/region.hpp"

#include <algorithm>
#include <map>
#include <set>
#include <stdexcept>
#include <utility>

#include <flint/flint.h>
#include <flint/fmpq.h>
#include <flint/fmpz.h>
#include <flint/fmpz_mpoly.h>
#include <flint/fmpz_mpoly_q.h>

#include "amflow/qft/dense_q.hpp"
#include "amflow/qft/family_uf.hpp"

namespace amflow::qft {

using algebra::Mfrac;
using algebra::Mpoly;
using algebra::MpolyContext;

DListContext make_dlist_context(const FamilyConfig& fc,
                                  std::size_t n_denominators) {
    std::vector<std::string> names;
    names.reserve((std::size_t)fc.ctx->n_vars() + n_denominators);
    for (long i = 0; i < fc.ctx->n_vars(); ++i) {
        names.push_back(fc.ctx->var_name(i));
    }
    for (std::size_t i = 0; i < n_denominators; ++i) {
        names.push_back("__amf_D" + std::to_string(i + 1));
    }
    auto d_ctx = std::make_shared<MpolyContext>(std::move(names));

    DListContext dctx;
    dctx.ctx          = d_ctx;
    dctx.first_d_var  = fc.ctx->n_vars();
    dctx.n_d          = (long)n_denominators;
    dctx.lift_gens.reserve((std::size_t)fc.ctx->n_vars());
    for (long i = 0; i < fc.ctx->n_vars(); ++i) {
        dctx.lift_gens.push_back(Mpoly::variable(d_ctx, i));
    }
    return dctx;
}

void coeff_over_splist(const Mpoly& den,
                        const std::vector<Mpoly>& sp_list,
                        std::vector<fmpq*>& out) {
    auto ctx = den.ctx();
    long n_vars = ctx->n_vars();
    std::vector<unsigned long> exp((std::size_t)n_vars);
    for (std::size_t k = 0; k < sp_list.size(); ++k) {
        long len_sp = fmpz_mpoly_length(sp_list[k].raw(), ctx->raw());
        if (len_sp != 1) {
            throw std::runtime_error(
                "coeff_over_splist: sp_list entry not a monomial");
        }
        fmpz_mpoly_get_term_exp_ui(exp.data(), sp_list[k].raw(),
                                    0, ctx->raw());
        fmpz_t c;
        fmpz_init(c);
        fmpz_mpoly_get_coeff_fmpz_ui(c, den.raw(), exp.data(),
                                      ctx->raw());
        fmpq_set_fmpz(out[k], c);
        fmpz_clear(c);
    }
}

namespace {

std::map<std::vector<unsigned long>, long>
build_splist_signature_index(const FamilyConfig& fc) {
    std::map<std::vector<unsigned long>, long> out;
    const long n_loop = static_cast<long>(fc.n_loops());
    const long n_red = static_cast<long>(fc.n_red_legs());
    const long n_loop_leg = n_loop + n_red;
    for (long k = 0; k < static_cast<long>(fc.sp_list.size()); ++k) {
        std::vector<unsigned long> sig((std::size_t)n_loop_leg, 0);
        std::vector<unsigned long> exp((std::size_t)fc.ctx->n_vars(), 0);
        fmpz_mpoly_get_term_exp_ui(exp.data(), fc.sp_list[(std::size_t)k].raw(),
                                   0, fc.ctx->raw());
        for (long i = 0; i < n_loop_leg; ++i) {
            sig[(std::size_t)i] = exp[(std::size_t)i];
        }
        out.emplace(std::move(sig), k);
    }
    return out;
}

struct DenominatorSplistDecomposition {
    explicit DenominatorSplistDecomposition(
            const std::shared_ptr<MpolyContext>& ctx)
        : constant(Mfrac::zero(ctx)) {}

    std::vector<Mpoly> sp_coef;
    Mfrac constant;
};

DenominatorSplistDecomposition
decompose_denominator_over_splist(const FamilyConfig& fc, const Mpoly& den) {
    const long n_loop = static_cast<long>(fc.n_loops());
    const long n_red = static_cast<long>(fc.n_red_legs());
    const long n_loop_leg = n_loop + n_red;
    const long n_total = fc.ctx->n_vars();
    auto sp_index = build_splist_signature_index(fc);

    DenominatorSplistDecomposition out(fc.ctx);
    out.sp_coef.reserve(fc.sp_list.size());
    for (std::size_t k = 0; k < fc.sp_list.size(); ++k) {
        out.sp_coef.emplace_back(fc.ctx);
    }

    const long len = fmpz_mpoly_length(den.raw(), fc.ctx->raw());
    std::vector<unsigned long> exp((std::size_t)n_total);
    fmpz_t coeff;
    fmpz_init(coeff);
    for (long t = 0; t < len; ++t) {
        fmpz_mpoly_get_term_coeff_fmpz(coeff, den.raw(), t, fc.ctx->raw());
        fmpz_mpoly_get_term_exp_ui(exp.data(), den.raw(), t, fc.ctx->raw());

        bool has_loop = false;
        std::vector<unsigned long> sig((std::size_t)n_loop_leg, 0);
        for (long i = 0; i < n_loop_leg; ++i) {
            sig[(std::size_t)i] = exp[(std::size_t)i];
            if (i < n_loop && exp[(std::size_t)i] != 0) has_loop = true;
        }

        auto sit = sp_index.find(sig);
        if (sit != sp_index.end()) {
            std::vector<unsigned long> inv_exp = exp;
            for (long i = 0; i < n_loop_leg; ++i) {
                inv_exp[(std::size_t)i] = 0;
            }
            Mpoly inv_mono = Mpoly::zero(fc.ctx);
            fmpz_mpoly_set_coeff_fmpz_ui(inv_mono.raw(), coeff,
                                         inv_exp.data(), fc.ctx->raw());
            out.sp_coef[(std::size_t)sit->second] += inv_mono;
            continue;
        }

        if (has_loop) {
            fmpz_clear(coeff);
            throw std::runtime_error(
                "decompose_denominator_over_splist: loop-dependent monomial "
                "is not in SPList");
        }

        Mpoly mono = Mpoly::zero(fc.ctx);
        fmpz_mpoly_set_coeff_fmpz_ui(mono.raw(), coeff,
                                     exp.data(), fc.ctx->raw());
        out.constant += fc.apply_replacement(mono);
    }
    fmpz_clear(coeff);
    return out;
}

std::vector<std::vector<Mfrac>>
zero_mfrac_matrix(const std::shared_ptr<MpolyContext>& ctx,
                  std::size_t rows,
                  std::size_t cols) {
    std::vector<std::vector<Mfrac>> out;
    out.reserve(rows);
    for (std::size_t r = 0; r < rows; ++r) {
        out.emplace_back();
        out.back().reserve(cols);
        for (std::size_t c = 0; c < cols; ++c) {
            out.back().push_back(Mfrac::zero(ctx));
        }
    }
    return out;
}

void mfrac_rref(std::vector<std::vector<Mfrac>>& mat) {
    if (mat.empty() || mat.front().empty()) return;
    const std::size_t rows = mat.size();
    const std::size_t cols = mat.front().size();
    std::size_t pivot_row = 0;
    for (std::size_t col = 0; col < cols && pivot_row < rows; ++col) {
        std::size_t found = rows;
        for (std::size_t r = pivot_row; r < rows; ++r) {
            if (!mat[r][col].is_zero()) {
                found = r;
                break;
            }
        }
        if (found == rows) continue;
        if (found != pivot_row) std::swap(mat[pivot_row], mat[found]);

        Mfrac pivot = mat[pivot_row][col].clone();
        for (std::size_t c = col; c < cols; ++c) {
            mat[pivot_row][c] /= pivot;
        }

        for (std::size_t r = 0; r < rows; ++r) {
            if (r == pivot_row || mat[r][col].is_zero()) continue;
            Mfrac mult = mat[r][col].clone();
            for (std::size_t c = col; c < cols; ++c) {
                Mfrac prod = mat[pivot_row][c].clone();
                prod *= mult;
                mat[r][c] -= prod;
            }
        }
        ++pivot_row;
    }
}

std::vector<std::vector<Mfrac>>
invert_mfrac_matrix(const std::shared_ptr<MpolyContext>& ctx,
                    const std::vector<std::vector<Mfrac>>& mat) {
    const std::size_t n = mat.size();
    if (n == 0) return {};
    for (const auto& row : mat) {
        if (row.size() != n) {
            throw std::invalid_argument("invert_mfrac_matrix: non-square matrix");
        }
    }

    auto aug = zero_mfrac_matrix(ctx, n, 2 * n);
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            aug[i][j] = mat[i][j].clone();
        }
        aug[i][n + i] = Mfrac::one(ctx);
    }
    mfrac_rref(aug);

    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            if (i == j) {
                if (!aug[i][j].is_one()) {
                    throw std::runtime_error(
                        "invert_mfrac_matrix: matrix is singular");
                }
            } else if (!aug[i][j].is_zero()) {
                throw std::runtime_error(
                    "invert_mfrac_matrix: matrix is singular");
            }
        }
    }

    auto inv = zero_mfrac_matrix(ctx, n, n);
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            inv[i][j] = aug[i][n + j].clone();
        }
    }
    return inv;
}

std::vector<std::size_t>
maximal_group_rows_mfrac(const std::shared_ptr<MpolyContext>& ctx,
                         std::size_t rows,
                         std::size_t cols,
                         const std::vector<std::vector<Mfrac>>& flat_rows) {
    if (flat_rows.size() != rows) {
        throw std::invalid_argument(
            "maximal_group_rows_mfrac: row count mismatch");
    }
    if (rows == 0) return {};

    auto transpose = zero_mfrac_matrix(ctx, cols, rows);
    for (std::size_t r = 0; r < rows; ++r) {
        if (flat_rows[r].size() != cols) {
            throw std::invalid_argument(
                "maximal_group_rows_mfrac: column count mismatch");
        }
        for (std::size_t c = 0; c < cols; ++c) {
            transpose[c][r] = flat_rows[r][c].clone();
        }
    }
    mfrac_rref(transpose);

    std::vector<std::size_t> picks;
    for (std::size_t r = 0; r < cols; ++r) {
        for (std::size_t c = 0; c < rows; ++c) {
            if (!transpose[r][c].is_zero()) {
                if (transpose[r][c].is_one()) picks.push_back(c);
                break;
            }
        }
    }
    std::sort(picks.begin(), picks.end());
    picks.erase(std::unique(picks.begin(), picks.end()), picks.end());
    return picks;
}

Mpoly lift_mpoly_to_dctx(const Mpoly& src, const DListContext& dctx) {
    Mpoly out(dctx.ctx);
    const long src_vars = src.ctx()->n_vars();
    const long len = fmpz_mpoly_length(src.raw(), src.ctx()->raw());
    std::vector<unsigned long> src_exp((std::size_t)src_vars);
    std::vector<unsigned long> dst_exp((std::size_t)dctx.ctx->n_vars(), 0);
    fmpz_t coeff;
    fmpz_init(coeff);
    for (long t = 0; t < len; ++t) {
        fmpz_mpoly_get_term_coeff_fmpz(coeff, src.raw(), t, src.ctx()->raw());
        fmpz_mpoly_get_term_exp_ui(src_exp.data(), src.raw(), t, src.ctx()->raw());
        std::fill(dst_exp.begin(), dst_exp.end(), 0);
        for (long i = 0; i < src_vars; ++i) {
            const long di = dctx.ctx->var_index(src.ctx()->var_name(i));
            if (di < 0) {
                fmpz_clear(coeff);
                throw std::runtime_error(
                    "lift_mpoly_to_dctx: variable missing in destination");
            }
            dst_exp[(std::size_t)di] = src_exp[(std::size_t)i];
        }
        fmpz_mpoly_set_coeff_fmpz_ui(out.raw(), coeff,
                                     dst_exp.data(), dctx.ctx->raw());
    }
    fmpz_clear(coeff);
    return out;
}

Mfrac lift_mfrac_to_dctx(const Mfrac& src, const DListContext& dctx) {
    return Mfrac(lift_mpoly_to_dctx(src.numerator(), dctx),
                 lift_mpoly_to_dctx(src.denominator(), dctx));
}

}  // namespace

std::vector<Mfrac>
sp_list_to_dlist_symbol(const FamilyConfig& fc,
                          const DListContext& dctx,
                          const std::vector<Mpoly>& denominators) {
    long N = (long)denominators.size();
    long M = (long)fc.sp_list.size();

    if (N != M) {
        throw std::invalid_argument(
            "sp_list_to_dlist_symbol: denominator count must equal "
            "#sp_list (use to_complete_explicit first); got "
            + std::to_string(N) + " vs " + std::to_string(M));
    }

    auto C = zero_mfrac_matrix(fc.ctx, (std::size_t)M, (std::size_t)M);
    std::vector<Mfrac> consts;
    consts.reserve((std::size_t)M);
    for (long i = 0; i < N; ++i) {
        auto dec = decompose_denominator_over_splist(
            fc, denominators[(std::size_t)i]);
        for (long k = 0; k < M; ++k) {
            const Mpoly& cf = dec.sp_coef[(std::size_t)k];
            if (!fmpz_mpoly_is_fmpz(cf.raw(), fc.ctx->raw())) {
                throw std::runtime_error(
                    "sp_list_to_dlist_symbol: sp_list bilinear coefficient "
                    "is not a constant integer -- non-physical (e.g. "
                    "linear / Wilson-line) denominator?");
            }
            C[(std::size_t)i][(std::size_t)k] = Mfrac::from_mpoly(cf.clone());
        }
        consts.push_back(std::move(dec.constant));
    }
    auto Cinv = invert_mfrac_matrix(fc.ctx, C);

    std::vector<Mfrac> out;
    out.reserve((std::size_t)M);
    for (long k = 0; k < M; ++k) {
        Mfrac acc = Mfrac::zero(dctx.ctx);
        for (long i = 0; i < M; ++i) {
            const Mfrac& qval = Cinv[(std::size_t)k][(std::size_t)i];
            if (qval.is_zero()) continue;

            Mfrac di_minus_const = Mfrac::from_mpoly(
                Mpoly::variable(dctx.ctx, dctx.first_d_var + i));
            di_minus_const -= lift_mfrac_to_dctx(consts[(std::size_t)i], dctx);
            di_minus_const *= lift_mfrac_to_dctx(qval, dctx);
            acc += di_minus_const;
        }
        out.push_back(std::move(acc));
    }
    return out;
}

std::vector<Mpoly>
to_complete_explicit(const FamilyConfig& fc,
                       const std::vector<Mpoly>& denominators) {
    long n_loop = (long)fc.n_loops();
    long n_red  = (long)fc.n_red_legs();

    std::set<std::string> seen;
    std::vector<Mpoly> pde;

    auto add_unique = [&](Mpoly&& p) {
        std::string k = p.to_string();
        if (seen.insert(k).second) pde.push_back(std::move(p));
    };

    for (long i = 0; i < n_loop; ++i) {
        Mpoly l = Mpoly::variable(fc.ctx, i);
        add_unique(l * l);
    }
    for (long i = 0; i < n_loop; ++i) {
        for (long j = i + 1; j < n_loop; ++j) {
            Mpoly s = Mpoly::variable(fc.ctx, i)
                    + Mpoly::variable(fc.ctx, j);
            add_unique(s * s);
        }
    }
    for (long i = 0; i < n_loop; ++i) {
        for (long k = 0; k < n_red; ++k) {
            Mpoly s = Mpoly::variable(fc.ctx, i)
                    + Mpoly::variable(fc.ctx, n_loop + k);
            add_unique(s * s);
        }
    }

    std::vector<Mpoly> all;
    all.reserve(denominators.size() + pde.size());
    for (const auto& d : denominators) all.push_back(d.clone());
    for (auto& p : pde) all.push_back(std::move(p));

    long M = (long)fc.sp_list.size();
    long rows = (long)all.size();
    std::vector<std::vector<Mfrac>> flat;
    flat.reserve((std::size_t)rows);
    for (long i = 0; i < rows; ++i) {
        auto dec = decompose_denominator_over_splist(fc, all[(std::size_t)i]);
        flat.emplace_back();
        flat.back().reserve((std::size_t)M);
        for (long k = 0; k < M; ++k) {
            flat.back().push_back(
                Mfrac::from_mpoly(dec.sp_coef[(std::size_t)k].clone()));
        }
    }

    auto picks = maximal_group_rows_mfrac(
        fc.ctx, (std::size_t)rows, (std::size_t)M, flat);

    if ((long)picks.size() < M) {
        throw std::runtime_error(
            "to_complete_explicit: cannot complete to "
            + std::to_string(M) + " denominators");
    }

    std::vector<Mpoly> out;
    out.reserve((std::size_t)M);
    for (std::size_t k = 0; k < (std::size_t)M; ++k) {
        out.push_back(all[picks[k]].clone());
    }
    return out;
}

namespace {
Mpoly mfrac_to_mpoly_integer(const Mfrac& f) {
    if (!f.is_polynomial()) {
        throw std::runtime_error(
            "mfrac_to_mpoly_integer: expected a polynomial-valued Mfrac "
            "(mom^2 + mass should clear denominators)");
    }
    auto ctx = f.ctx();
    Mpoly den = f.denominator();
    if (!fmpz_mpoly_is_fmpz(den.raw(), ctx->raw())) {
        throw std::runtime_error(
            "mfrac_to_mpoly_integer: denominator is not an integer");
    }
    fmpz_t c;
    fmpz_init(c);
    fmpz_mpoly_get_fmpz(c, den.raw(), ctx->raw());
    if (fmpz_is_zero(c)) {
        fmpz_clear(c);
        throw std::runtime_error("mfrac_to_mpoly_integer: zero denominator");
    }
    Mpoly num = f.numerator();
    Mpoly out = num.clone();
    fmpz_mpoly_scalar_divexact_fmpz(out.raw(), num.raw(), c, ctx->raw());
    fmpz_clear(c);
    return out;
}
}  // namespace

std::vector<Mpoly>
squared_denominators(const FamilyConfig& fc,
                       const std::vector<Mpoly>& denominators) {
    std::vector<Mpoly> out;
    out.reserve(denominators.size());
    long n_loop = (long)fc.n_loops();

    for (const auto& d : denominators) {
        long chosen = -1;
        for (long j = 0; j < n_loop; ++j) {
            if (!d.coeff_of(j, 2).is_zero()) { chosen = j; break; }
        }
        if (chosen < 0) {
            out.push_back(d.clone());
            continue;
        }

        std::vector<Mpoly> one;
        one.push_back(d.clone());
        auto ts = to_square_all(fc, one);
        Mfrac mom = ts.momenta[0].clone();
        Mfrac comb = mom * mom;
        comb += ts.masses[0];
        out.push_back(mfrac_to_mpoly_integer(comb));
    }
    return out;
}

ToSquareResult
to_square_all(const FamilyConfig& fc,
                const std::vector<Mpoly>& denominators) {
    long n_loop = (long)fc.n_loops();

    ToSquareResult out;
    out.momenta.reserve(denominators.size());
    out.masses.reserve(denominators.size());

    for (const auto& d : denominators) {
        if (d.ctx().get() != fc.ctx.get()) {
            throw std::invalid_argument(
                "to_square_all: denominator ctx mismatch");
        }
        long chosen = -1;
        for (long j = 0; j < n_loop; ++j) {
            if (!d.coeff_of(j, 2).is_zero()) { chosen = j; break; }
        }
        if (chosen < 0) {
            // Bilinear ISP (e.g. `k*a`, `k*p1`): no quadratic loop
            // term to complete-the-square against.  `kira_yaml` only
            // consumes `.masses` for kinematic-invariant detection
            // (loop-free symbols), and renders the propagator from
            // the raw form via the `try_render_closed_form` fallback,
            // so push (0, 0) and continue.  The SingleMass /
            // FactorizeFamily callers (amfsystem.cpp) never reach
            // this path on bilinear ISPs because their per-prop
            // `find_loop_in_prop` guard returns -1 first.
            out.momenta.push_back(Mfrac::zero(fc.ctx));
            out.masses.push_back(Mfrac::zero(fc.ctx));
            continue;
        }

        Mpoly coe1 = d.coeff_of(chosen, 2);
        Mpoly coe2 = d.coeff_of(chosen, 1);
        Mpoly coe3 = d.coeff_of(chosen, 0);

        Mpoly two_coe1(fc.ctx);
        {
            fmpz_t two; fmpz_init_set_ui(two, 2);
            fmpz_mpoly_scalar_mul_fmpz(two_coe1.raw(), coe1.raw(),
                                        two, fc.ctx->raw());
            fmpz_clear(two);
        }
        Mpoly mom_num(fc.ctx);
        {
            Mpoly l = Mpoly::variable(fc.ctx, chosen);
            mom_num = two_coe1 * l + coe2;
        }
        Mfrac momentum(std::move(mom_num), two_coe1.clone());

        Mpoly four_coe1(fc.ctx);
        {
            fmpz_t four; fmpz_init_set_ui(four, 4);
            fmpz_mpoly_scalar_mul_fmpz(four_coe1.raw(), coe1.raw(),
                                        four, fc.ctx->raw());
            fmpz_clear(four);
        }
        Mpoly mass_num = four_coe1 * coe3 - coe2 * coe2;
        Mfrac mass_raw(std::move(mass_num), std::move(four_coe1));

        Mfrac mass = apply_uf_replacement(fc, fc.ctx, mass_raw);

        out.momenta.push_back(std::move(momentum));
        out.masses.push_back(std::move(mass));
    }
    return out;
}

}  // namespace amflow::qft
