// SPDX-License-Identifier: MIT
// ibp::libp_deriv — implementation.
// Provenance (AMFlow Kira interface; LiteIBP by T. Peraro): see libp_deriv.hpp.

#include "amflow/ibp/libp_deriv.hpp"

#include <map>
#include <set>
#include <stdexcept>
#include <utility>

#include <flint/fmpq.h>
#include <flint/fmpz.h>
#include <flint/fmpz_mpoly.h>
#include <flint/fmpz_mpoly_q.h>

#include "amflow/algebra/numeric_subst.hpp"

namespace amflow::ibp {

using algebra::Mfrac;
using algebra::Mpoly;

namespace {

long find_var_index(const qft::FamilyConfig& fc, const std::string& name) {
    long idx = fc.ctx->var_index(name);
    if (idx < 0) {
        throw std::invalid_argument(
            "libp_deriv: invariant '" + name + "' not in fc.ctx");
    }
    return idx;
}

std::map<std::vector<unsigned long>, long>
build_sp_index(const qft::FamilyConfig& fc) {
    std::map<std::vector<unsigned long>, long> out;
    long n_loop = (long)fc.n_loops();
    long n_red  = (long)fc.n_red_legs();
    long n_loop_leg = n_loop + n_red;
    for (long k = 0; k < (long)fc.sp_list.size(); ++k) {
        std::vector<unsigned long> sig((std::size_t)n_loop_leg, 0);
        std::vector<unsigned long> exp((std::size_t)fc.ctx->n_vars(), 0);
        fmpz_mpoly_get_term_exp_ui(exp.data(), fc.sp_list[(std::size_t)k].raw(),
                                    0, fc.ctx->raw());
        for (long i = 0; i < n_loop_leg; ++i) {
            sig[(std::size_t)i] = exp[(std::size_t)i];
        }
        out[sig] = k;
    }
    return out;
}

struct SpDecomposition {
    std::vector<Mpoly> sp_coef;
    Mpoly              constant;
};

SpDecomposition decompose_into_splist(const qft::FamilyConfig& fc,
                                        const Mpoly& P) {
    long n_loop = (long)fc.n_loops();
    long n_red  = (long)fc.n_red_legs();
    long n_loop_leg = n_loop + n_red;
    long n_total = fc.ctx->n_vars();

    SpDecomposition out;
    out.constant = Mpoly::zero(fc.ctx);
    out.sp_coef.reserve(fc.sp_list.size());
    for (std::size_t k = 0; k < fc.sp_list.size(); ++k) {
        out.sp_coef.emplace_back(fc.ctx);
    }

    auto sp_index = build_sp_index(fc);

    long len = fmpz_mpoly_length(P.raw(), fc.ctx->raw());
    std::vector<unsigned long> exp((std::size_t)n_total);
    fmpz_t coeff;
    fmpz_init(coeff);
    for (long t = 0; t < len; ++t) {
        fmpz_mpoly_get_term_coeff_fmpz(coeff, P.raw(), t, fc.ctx->raw());
        fmpz_mpoly_get_term_exp_ui(exp.data(), P.raw(), t, fc.ctx->raw());

        std::vector<unsigned long> sig((std::size_t)n_loop_leg, 0);
        for (long i = 0; i < n_loop_leg; ++i) {
            sig[(std::size_t)i] = exp[(std::size_t)i];
        }
        std::vector<unsigned long> inv_exp = exp;
        for (long i = 0; i < n_loop_leg; ++i) inv_exp[(std::size_t)i] = 0;

        unsigned long deg_loopleg = 0;
        for (long i = 0; i < n_loop_leg; ++i) deg_loopleg += sig[(std::size_t)i];

        if (deg_loopleg == 0) {
            Mpoly mono = Mpoly::zero(fc.ctx);
            fmpz_mpoly_set_coeff_fmpz_ui(mono.raw(), coeff,
                                          inv_exp.data(), fc.ctx->raw());
            out.constant += mono;
        } else if (deg_loopleg == 2) {
            auto it = sp_index.find(sig);
            if (it == sp_index.end()) {
                fmpz_clear(coeff);
                throw std::runtime_error(
                    "decompose_into_splist: monomial not in sp_list");
            }
            Mpoly inv_mono = Mpoly::zero(fc.ctx);
            fmpz_mpoly_set_coeff_fmpz_ui(inv_mono.raw(), coeff,
                                          inv_exp.data(), fc.ctx->raw());
            out.sp_coef[(std::size_t)it->second] += inv_mono;
        } else {
            fmpz_clear(coeff);
            throw std::runtime_error(
                "decompose_into_splist: monomial has loop+leg degree " +
                std::to_string(deg_loopleg) + ", expected 0 or 2");
        }
    }
    fmpz_clear(coeff);
    return out;
}

struct SpToDMap {
    std::vector<Mpoly>                completede;
    std::vector<std::vector<Mfrac>>   alpha;
    std::vector<Mfrac>                beta;
};

// ===========================================================================
//  Kinematic-invariant momentum-derivative chain rule.
//
//  For a Replacement-defined invariant `x` (e.g. `s` defined via
//  `(p1+p2)^2 -> s`) the propagators `D_k` (stored pre-replacement on
//  `fc.ctx`) carry NO explicit `x`-dependence, so the direct
//  `D_k.derivative(x_var)` is identically zero.  The integral's
//  dependence on `x` enters only through the external Gram matrix
//  G_{ij} = p_i·p_j = R_{ij}(x)  (the `fc.reduced_replacement` map).
//  By Lorentz invariance I = I(G), and the standard Gehrmann-Remiddi
//  / LiteRed momentum-derivative chain rule gives
//
//      ∂I/∂x = Σ_{a,i} C_{ai}  (p_a · ∂I/∂p_i)
//
//  with the (n_red × n_red) coefficient matrix
//
//      C  =  ½  G⁻¹ · J ,       J_{ij} = ∂R_{ij}/∂x .
//
//  Under the integral sign  p_a · ∂I/∂p_i  =  Σ_k (-a_k) (p_a · ∂D_k/∂p_i)
//  · J(a + e_k - …),  so the net effect on `libp_denoms_deriv` is that
//  the per-propagator "effective ∂D_k/∂x" becomes
//
//      ∂D_k/∂x|_direct  +  Σ_{a,i} C_{ai} · (p_a · ∂D_k/∂p_i)
//
//  with the second term decomposed into sp_list (loop·leg → D_j via
//  `sp_to_d_map`) plus a leg·leg / constant remainder (→ invariants via
//  `fc.apply_replacement`).  When `x` does not appear in any
//  reduced_replacement value (η, propagator-level masses) `J ≡ 0`,
//  `C ≡ 0`, and the chain term is skipped — byte-identical to the
//  pre-patch path.
// ===========================================================================

std::string legleg_key(const qft::FamilyConfig& fc, long i, long j) {
    const std::string& a = fc.reduced_legs[(std::size_t)i];
    const std::string& b = fc.reduced_legs[(std::size_t)j];
    return (a <= b) ? (a + "*" + b) : (b + "*" + a);
}

// Tiny Gauss-Jordan inverse over Mfrac (n_red ≤ ~5 in practice).  Mirrors
// `qft/complete.cpp::invert_mfrac_matrix` (file-static there).
std::vector<std::vector<Mfrac>>
invert_mfrac_matrix_local(const std::shared_ptr<algebra::MpolyContext>& ctx,
                          const std::vector<std::vector<Mfrac>>& M) {
    const std::size_t n = M.size();
    std::vector<std::vector<Mfrac>> aug;
    aug.reserve(n);
    for (std::size_t i = 0; i < n; ++i) {
        std::vector<Mfrac> row;
        row.reserve(2 * n);
        for (std::size_t j = 0; j < n; ++j) row.push_back(M[i][j].clone());
        for (std::size_t j = 0; j < n; ++j) {
            row.push_back(i == j ? Mfrac::one(ctx) : Mfrac::zero(ctx));
        }
        aug.push_back(std::move(row));
    }
    for (std::size_t col = 0; col < n; ++col) {
        std::size_t piv = col;
        while (piv < n && aug[piv][col].is_zero()) ++piv;
        if (piv == n) {
            throw std::runtime_error(
                "libp_deriv: external Gram matrix is singular; "
                "kinematic-invariant derivative is undefined at this "
                "kinematic configuration (pin a scale or perturb the "
                "numeric point off the degenerate surface).");
        }
        if (piv != col) std::swap(aug[piv], aug[col]);
        Mfrac inv = Mfrac::one(ctx); inv /= aug[col][col];
        for (std::size_t c = col; c < 2 * n; ++c) aug[col][c] *= inv;
        for (std::size_t r = 0; r < n; ++r) {
            if (r == col || aug[r][col].is_zero()) continue;
            Mfrac f = aug[r][col].clone();
            for (std::size_t c = col; c < 2 * n; ++c) {
                Mfrac t = aug[col][c].clone(); t *= f;
                aug[r][c] -= t;
            }
        }
    }
    std::vector<std::vector<Mfrac>> out;
    out.reserve(n);
    for (std::size_t i = 0; i < n; ++i) {
        std::vector<Mfrac> row;
        row.reserve(n);
        for (std::size_t j = 0; j < n; ++j) {
            row.push_back(std::move(aug[i][n + j]));
        }
        out.push_back(std::move(row));
    }
    return out;
}

struct KinematicChain {
    bool active = false;                     // false ⇒ J≡0, skip chain
    std::vector<std::vector<Mfrac>> C;       // [n_red][n_red], on fc.ctx
};

KinematicChain build_kinematic_chain(const qft::FamilyConfig& fc, long s_var) {
    KinematicChain out;
    long n_red = (long)fc.n_red_legs();
    if (n_red == 0) return out;

    // J_{ij} = ∂R_{ij}/∂x.  If all zero (η, mass) → inactive.
    std::vector<std::vector<Mfrac>> J;
    J.reserve((std::size_t)n_red);
    bool any_nonzero = false;
    for (long i = 0; i < n_red; ++i) {
        std::vector<Mfrac> row;
        row.reserve((std::size_t)n_red);
        for (long j = 0; j < n_red; ++j) {
            auto it = fc.reduced_replacement.find(legleg_key(fc, i, j));
            if (it == fc.reduced_replacement.end()) {
                // Missing leg·leg replacement.  Only fatal if the chain is
                // actually needed; defer the throw until we know J≠0.
                row.push_back(Mfrac::zero(fc.ctx));
            } else {
                Mfrac d = it->second.derivative(s_var);
                if (!d.is_zero()) any_nonzero = true;
                row.push_back(std::move(d));
            }
        }
        J.push_back(std::move(row));
    }
    if (!any_nonzero) return out;   // η / mass / non-kinematic → skip chain

    // G_{ij} = R_{ij}(x).  Now every key MUST be present (kinematic case).
    std::vector<std::vector<Mfrac>> G;
    G.reserve((std::size_t)n_red);
    for (long i = 0; i < n_red; ++i) {
        std::vector<Mfrac> row;
        row.reserve((std::size_t)n_red);
        for (long j = 0; j < n_red; ++j) {
            auto it = fc.reduced_replacement.find(legleg_key(fc, i, j));
            if (it == fc.reduced_replacement.end()) {
                throw std::runtime_error(
                    "libp_deriv: kinematic derivative w.r.t. '"
                    + fc.ctx->var_name(s_var)
                    + "' needs a replacement rule for every external "
                      "scalar product, but '" + legleg_key(fc, i, j)
                    + "' is missing.");
            }
            row.push_back(it->second.clone());
        }
        G.push_back(std::move(row));
    }
    auto Ginv = invert_mfrac_matrix_local(fc.ctx, G);

    // C = ½ G⁻¹ J.
    Mfrac half = Mfrac::from_si(fc.ctx, 1);
    half /= Mfrac::from_si(fc.ctx, 2);
    out.C.reserve((std::size_t)n_red);
    for (long a = 0; a < n_red; ++a) {
        std::vector<Mfrac> row;
        row.reserve((std::size_t)n_red);
        for (long i = 0; i < n_red; ++i) {
            Mfrac acc = Mfrac::zero(fc.ctx);
            for (long j = 0; j < n_red; ++j) {
                if (Ginv[(std::size_t)a][(std::size_t)j].is_zero()
                    || J[(std::size_t)j][(std::size_t)i].is_zero()) continue;
                Mfrac t = Ginv[(std::size_t)a][(std::size_t)j].clone();
                t *= J[(std::size_t)j][(std::size_t)i];
                acc += t;
            }
            acc *= half;
            row.push_back(std::move(acc));
        }
        out.C.push_back(std::move(row));
    }
    out.active = true;
    return out;
}

// Decompose `Q = p_a · ∂D_k/∂p_i` (degree-2 in {loops,legs}) into sp_list
// coefficients (loop·loop, loop·leg) plus a constant Mfrac (leg·leg →
// invariants via `fc.apply_replacement`, plus any genuine degree-0 part).
struct ChainDecomposition {
    std::vector<Mpoly> sp_coef;   // length = fc.sp_list.size(), on fc.ctx
    Mfrac              constant;  // on fc.ctx, leg·leg already replaced
};

ChainDecomposition decompose_chain(const qft::FamilyConfig& fc,
                                   const Mpoly& Q,
                                   const std::map<std::vector<unsigned long>,
                                                  long>& sp_index) {
    long n_loop = (long)fc.n_loops();
    long n_red  = (long)fc.n_red_legs();
    long n_loop_leg = n_loop + n_red;
    long n_total = fc.ctx->n_vars();

    ChainDecomposition out;
    out.constant = Mfrac::zero(fc.ctx);
    out.sp_coef.reserve(fc.sp_list.size());
    for (std::size_t k = 0; k < fc.sp_list.size(); ++k) {
        out.sp_coef.emplace_back(fc.ctx);
    }

    long len = fmpz_mpoly_length(Q.raw(), fc.ctx->raw());
    std::vector<unsigned long> exp((std::size_t)n_total);
    fmpz_t coeff;
    fmpz_init(coeff);
    for (long t = 0; t < len; ++t) {
        fmpz_mpoly_get_term_coeff_fmpz(coeff, Q.raw(), t, fc.ctx->raw());
        fmpz_mpoly_get_term_exp_ui(exp.data(), Q.raw(), t, fc.ctx->raw());

        std::vector<unsigned long> sig((std::size_t)n_loop_leg, 0);
        bool has_loop = false;
        for (long i = 0; i < n_loop_leg; ++i) {
            sig[(std::size_t)i] = exp[(std::size_t)i];
            if (i < n_loop && exp[(std::size_t)i] != 0) has_loop = true;
        }
        std::vector<unsigned long> inv_exp = exp;
        for (long i = 0; i < n_loop_leg; ++i) inv_exp[(std::size_t)i] = 0;

        auto it = sp_index.find(sig);
        if (it != sp_index.end()) {
            Mpoly inv_mono = Mpoly::zero(fc.ctx);
            fmpz_mpoly_set_coeff_fmpz_ui(inv_mono.raw(), coeff,
                                          inv_exp.data(), fc.ctx->raw());
            out.sp_coef[(std::size_t)it->second] += inv_mono;
            continue;
        }
        if (has_loop) {
            fmpz_clear(coeff);
            throw std::runtime_error(
                "decompose_chain: loop-dependent monomial not in sp_list");
        }
        // Leg·leg or pure-invariant remainder → apply_replacement.
        Mpoly mono = Mpoly::zero(fc.ctx);
        fmpz_mpoly_set_coeff_fmpz_ui(mono.raw(), coeff,
                                      exp.data(), fc.ctx->raw());
        out.constant += fc.apply_replacement(mono);
    }
    fmpz_clear(coeff);
    return out;
}

// Add Σ_{a,i} C_{ai} · (p_a · ∂D_k/∂p_i), decomposed, into the running
// (sp_coef_acc[l], const_acc).  All on fc.ctx.
void accumulate_chain_for_prop(const qft::FamilyConfig& fc,
                               const Mpoly& D_k,
                               const KinematicChain& chain,
                               const std::map<std::vector<unsigned long>,
                                              long>& sp_index,
                               std::vector<Mfrac>& sp_coef_acc,
                               Mfrac& const_acc) {
    long n_loop = (long)fc.n_loops();
    long n_red  = (long)fc.n_red_legs();
    for (long i = 0; i < n_red; ++i) {
        Mpoly dDk_dpi = D_k.derivative(n_loop + i);
        if (dDk_dpi.is_zero()) continue;
        for (long a = 0; a < n_red; ++a) {
            const Mfrac& Cai = chain.C[(std::size_t)a][(std::size_t)i];
            if (Cai.is_zero()) continue;
            Mpoly Q = Mpoly::variable(fc.ctx, n_loop + a) * dDk_dpi;
            ChainDecomposition cd = decompose_chain(fc, Q, sp_index);
            for (std::size_t l = 0; l < cd.sp_coef.size(); ++l) {
                if (cd.sp_coef[l].is_zero()) continue;
                Mfrac t = Mfrac::from_mpoly(cd.sp_coef[l].clone());
                t *= Cai;
                sp_coef_acc[l] += t;
            }
            if (!cd.constant.is_zero()) {
                Mfrac t = cd.constant.clone();
                t *= Cai;
                const_acc += t;
            }
        }
    }
}

SpToDMap build_sp_to_d_map(const qft::FamilyConfig& fc) {
    auto completede = qft::to_complete_explicit(fc, fc.propagators_after_conservation);
    auto dctx = qft::make_dlist_context(fc, completede.size());
    auto sp_to_d = qft::sp_list_to_dlist_symbol(fc, dctx, completede);

    SpToDMap out;
    out.completede = std::move(completede);
    long n_d = dctx.n_d;
    long first_d = dctx.first_d_var;
    long n_total_d = dctx.ctx->n_vars();

    out.alpha.reserve(sp_to_d.size());
    out.beta.reserve(sp_to_d.size());

    auto extract_in_fc = [&](const Mpoly& p) -> Mpoly {
        Mpoly out_p(fc.ctx);
        long len = fmpz_mpoly_length(p.raw(), dctx.ctx->raw());
        std::vector<unsigned long> dexp((std::size_t)n_total_d);
        std::vector<unsigned long> fexp((std::size_t)fc.ctx->n_vars(), 0);
        fmpz_t cf;
        fmpz_init(cf);
        for (long t = 0; t < len; ++t) {
            fmpz_mpoly_get_term_exp_ui(dexp.data(), p.raw(), t,
                                        dctx.ctx->raw());
            fmpz_mpoly_get_term_coeff_fmpz(cf, p.raw(), t,
                                            dctx.ctx->raw());
            for (long k = 0; k < n_d; ++k) {
                if (dexp[(std::size_t)(first_d + k)] > 0) {
                    fmpz_clear(cf);
                    throw std::runtime_error("extract_in_fc: D dependence found");
                }
            }
            for (long j = 0; j < fc.ctx->n_vars(); ++j) {
                fexp[(std::size_t)j] = dexp[(std::size_t)j];
            }
            fmpz_mpoly_set_coeff_fmpz_ui(out_p.raw(), cf,
                                          fexp.data(), fc.ctx->raw());
        }
        fmpz_clear(cf);
        return out_p;
    };

    for (const auto& m : sp_to_d) {
        Mpoly num = m.numerator();
        Mpoly den = m.denominator();

        std::vector<Mpoly> alpha_num;
        alpha_num.reserve((std::size_t)n_d);
        for (long k_d = 0; k_d < n_d; ++k_d) alpha_num.emplace_back(fc.ctx);
        Mpoly beta_num(fc.ctx);

        long len = fmpz_mpoly_length(num.raw(), dctx.ctx->raw());
        std::vector<unsigned long> dexp((std::size_t)n_total_d);
        std::vector<unsigned long> rest_exp((std::size_t)n_total_d);
        fmpz_t cf;
        fmpz_init(cf);
        for (long t = 0; t < len; ++t) {
            fmpz_mpoly_get_term_exp_ui(dexp.data(), num.raw(), t,
                                        dctx.ctx->raw());
            fmpz_mpoly_get_term_coeff_fmpz(cf, num.raw(), t,
                                            dctx.ctx->raw());

            long which_d = -1;
            unsigned long which_d_exp = 0;
            for (long k = 0; k < n_d; ++k) {
                if (dexp[(std::size_t)(first_d + k)] > 0) {
                    if (which_d != -1) {
                        fmpz_clear(cf);
                        throw std::runtime_error(
                            "build_sp_to_d_map: monomial involves multiple D vars");
                    }
                    which_d = k;
                    which_d_exp = dexp[(std::size_t)(first_d + k)];
                }
            }
            if (which_d != -1 && which_d_exp != 1) {
                fmpz_clear(cf);
                throw std::runtime_error(
                    "build_sp_to_d_map: D variable appears with non-linear exponent");
            }

            rest_exp = dexp;
            for (long k = 0; k < n_d; ++k) {
                rest_exp[(std::size_t)(first_d + k)] = 0;
            }

            std::vector<unsigned long> fexp((std::size_t)fc.ctx->n_vars(), 0);
            for (long j = 0; j < fc.ctx->n_vars(); ++j) {
                fexp[(std::size_t)j] = rest_exp[(std::size_t)j];
            }
            Mpoly rest_mono = Mpoly::zero(fc.ctx);
            fmpz_mpoly_set_coeff_fmpz_ui(rest_mono.raw(), cf,
                                          fexp.data(), fc.ctx->raw());

            if (which_d == -1) {
                beta_num += rest_mono;
            } else {
                alpha_num[(std::size_t)which_d] += rest_mono;
            }
        }
        fmpz_clear(cf);

        Mpoly den_fc = extract_in_fc(den);

        std::vector<Mfrac> alpha_row;
        alpha_row.reserve((std::size_t)n_d);
        for (long k = 0; k < n_d; ++k) {
            alpha_row.push_back(Mfrac(std::move(alpha_num[(std::size_t)k]),
                                        den_fc.clone()));
        }
        out.alpha.push_back(std::move(alpha_row));
        out.beta.push_back(Mfrac(std::move(beta_num), den_fc.clone()));
    }
    return out;
}

}  // namespace

LibpDenomsDerivResult
libp_denoms_deriv(const qft::FamilyConfig& fc, const std::string& s_name) {
    long s_var = find_var_index(fc, s_name);

    SpToDMap sp_map = build_sp_to_d_map(fc);

    LibpDenomsDerivResult out;
    out.completede = std::move(sp_map.completede);
    long n_d = (long)out.completede.size();

    // Kinematic momentum-derivative chain.  Inactive (`!chain.active`) for
    // any variable not appearing in `fc.reduced_replacement` (η, propagator
    // masses) — exactly the pre-patch path.
    KinematicChain chain = build_kinematic_chain(fc, s_var);
    auto sp_index = build_sp_index(fc);

    std::size_t N = fc.propagators_after_conservation.size();
    std::size_t n_sp = fc.sp_list.size();
    out.coef.reserve(N);
    out.constant.reserve(N);
    for (std::size_t k = 0; k < N; ++k) {
        Mpoly d_k = fc.propagators_after_conservation[k]
                        .derivative(s_var);
        SpDecomposition dec = decompose_into_splist(fc, d_k);

        // Per-k sp_list coefficients (Mfrac, fc.ctx) and constant.  Seeded
        // with the direct-derivative decomposition; the chain term (if
        // active) accumulates into the same buffers.
        std::vector<Mfrac> sp_coef;
        sp_coef.reserve(n_sp);
        for (std::size_t l = 0; l < n_sp; ++l) {
            sp_coef.push_back(Mfrac::from_mpoly(dec.sp_coef[l].clone()));
        }
        Mfrac const_term = Mfrac::from_mpoly(dec.constant.clone());

        if (chain.active) {
            accumulate_chain_for_prop(fc,
                fc.propagators_after_conservation[k],
                chain, sp_index, sp_coef, const_term);
        }

        std::vector<Mfrac> coef_row;
        coef_row.reserve((std::size_t)n_d);
        for (long j = 0; j < n_d; ++j) {
            Mfrac acc = Mfrac::zero(fc.ctx);
            for (std::size_t l = 0; l < n_sp; ++l) {
                if (sp_coef[l].is_zero()) continue;
                Mfrac t = sp_coef[l].clone();
                if (t.ctx().get() != sp_map.alpha[l][(std::size_t)j].ctx().get()) {
                    throw std::runtime_error(
                        "libp_denoms_deriv: ctx mismatch t vs alpha at l="
                        + std::to_string(l) + " j=" + std::to_string(j));
                }
                t *= sp_map.alpha[l][(std::size_t)j];
                acc += t;
            }
            coef_row.push_back(std::move(acc));
        }
        for (std::size_t l = 0; l < n_sp; ++l) {
            if (sp_coef[l].is_zero()) continue;
            Mfrac t = sp_coef[l].clone();
            t *= sp_map.beta[l];
            const_term += t;
        }
        // `sp_map.beta[l]` already has leg·leg replaced (see
        // `qft/complete.cpp::decompose_denominator_over_splist`), so
        // `const_term` is invariants-only.
        out.coef.push_back(std::move(coef_row));
        out.constant.push_back(std::move(const_term));
    }
    return out;
}

std::vector<DerivTerm>
libp_deriv(const qft::FamilyConfig& fc,
           const qft::JIntegral& j,
           const std::string& s_name) {
    auto dd = libp_denoms_deriv(fc, s_name);
    long n_d = (long)dd.completede.size();
    long N   = (long)fc.propagators_after_conservation.size();

    if ((long)j.n_indices() != N) {
        throw std::invalid_argument(
            "libp_deriv: J integral indices length != #propagators");
    }

    std::vector<DerivTerm> out;
    for (long k = 0; k < N; ++k) {
        long ak = j.indices()[(std::size_t)k];
        if (ak == 0) continue;

        Mfrac neg_ak = Mfrac::from_si(fc.ctx, -ak);

        for (long jp = 0; jp < n_d; ++jp) {
            const Mfrac& cf = dd.coef[(std::size_t)k][(std::size_t)jp];
            if (cf.is_zero()) continue;
            std::vector<long> new_idx = j.indices();
            new_idx[(std::size_t)k] += 1;
            if (jp < N) {
                new_idx[(std::size_t)jp] -= 1;
            }
            qft::JIntegral shifted(j.family(), new_idx);
            DerivTerm dt;
            dt.integ = shifted;
            dt.coef = neg_ak.clone();
            dt.coef *= cf;
            out.push_back(std::move(dt));
        }

        const Mfrac& ck = dd.constant[(std::size_t)k];
        if (!ck.is_zero()) {
            std::vector<long> new_idx = j.indices();
            new_idx[(std::size_t)k] += 1;
            qft::JIntegral shifted(j.family(), new_idx);
            DerivTerm dt;
            dt.integ = shifted;
            dt.coef = neg_ak.clone();
            dt.coef *= ck;
            out.push_back(std::move(dt));
        }
    }
    return out;
}

std::vector<DerivTerm>
compute_derivative(const qft::FamilyConfig& fc,
                   const std::vector<DerivTerm>& expr,
                   const std::string& x_name) {
    long x_var = find_var_index(fc, x_name);

    std::vector<DerivTerm> out;
    for (const auto& term : expr) {
        Mfrac dc = term.coef.derivative(x_var);
        if (!dc.is_zero()) {
            DerivTerm dt;
            dt.coef = std::move(dc);
            dt.integ = term.integ;
            out.push_back(std::move(dt));
        }
        auto dj = libp_deriv(fc, term.integ, x_name);
        for (auto& d : dj) {
            DerivTerm dt;
            dt.coef = term.coef.clone();
            dt.coef *= d.coef;
            dt.integ = std::move(d.integ);
            out.push_back(std::move(dt));
        }
    }
    return out;
}

std::vector<DerivTerm>
simplify_terms(std::vector<DerivTerm> terms) {
    std::map<std::string, std::vector<std::size_t>> groups;
    for (std::size_t i = 0; i < terms.size(); ++i) {
        std::string key = terms[i].integ.family();
        for (long v : terms[i].integ.indices()) {
            key += "|"; key += std::to_string(v);
        }
        groups[key].push_back(i);
    }
    std::vector<DerivTerm> out;
    for (auto& [_, idxs] : groups) {
        DerivTerm acc;
        acc.integ = terms[idxs.front()].integ;
        acc.coef = Mfrac::zero(terms[idxs.front()].coef.ctx());
        for (std::size_t i : idxs) acc.coef += terms[i].coef;
        if (!acc.coef.is_zero()) out.push_back(std::move(acc));
    }
    return out;
}

// ===========================================================================
// Narrow-context overloads.
//
// Mirror `libp_denoms_deriv` / `libp_deriv` above but substitute the
// user-supplied `numeric_values` into intermediates and reproject
// everything to a caller-provided narrow polynomial-ring context
// (`target_ctx`).  Avoids the FLINT multivariate-GCD overhead that
// dominates the wide-fc.ctx path on multi-distinct-mass 3L+
// topologies (memory footprint can balloon by 10×+ otherwise).
//
// The wide-ctx version is preserved for tests and any caller that
// genuinely needs symbolic dependence on family variables.
// ===========================================================================

LibpDenomsDerivResult
libp_denoms_deriv(const qft::FamilyConfig& fc, const std::string& s_name,
                  const std::map<std::string, std::string>& numeric_values,
                  const std::shared_ptr<algebra::MpolyContext>& target_ctx) {
    (void)find_var_index(fc, s_name);  // sanity check: s_name is in fc.ctx

    // Build sp_map on fc.ctx (this is essentially free — sp_map's
    // Mfracs are small rational coefficients of the inverse of a
    // SP-coefficient matrix).
    SpToDMap sp_map = build_sp_to_d_map(fc);

    LibpDenomsDerivResult out;
    out.completede = std::move(sp_map.completede);
    long n_d = (long)out.completede.size();

    // Pre-substitute + reproject sp_map.alpha[l][j] and sp_map.beta[l]
    // onto `target_ctx`.  Done ONCE up front so the per-k loop below
    // does arithmetic only on the narrow ctx.
    const std::size_t n_sp = sp_map.alpha.size();
    std::vector<std::vector<Mfrac>> alpha_narrow(n_sp);
    std::vector<Mfrac> beta_narrow;
    beta_narrow.reserve(n_sp);
    for (std::size_t l = 0; l < n_sp; ++l) {
        alpha_narrow[l].reserve((std::size_t)n_d);
        for (long j = 0; j < n_d; ++j) {
            alpha_narrow[l].push_back(
                algebra::substitute_and_narrow(
                    sp_map.alpha[l][(std::size_t)j],
                    numeric_values, target_ctx));
        }
        beta_narrow.push_back(
            algebra::substitute_and_narrow(
                sp_map.beta[l], numeric_values, target_ctx));
    }

    // Kinematic momentum-derivative chain (see wide-ctx overload).  Chain
    // coefficients C_{ai} are built on fc.ctx then narrowed once up front
    // so the per-k loop does only narrow-ctx arithmetic.
    long s_var = find_var_index(fc, s_name);
    KinematicChain chain = build_kinematic_chain(fc, s_var);
    auto sp_index = build_sp_index(fc);

    std::size_t N = fc.propagators_after_conservation.size();
    out.coef.reserve(N);
    out.constant.reserve(N);
    for (std::size_t k = 0; k < N; ++k) {
        Mpoly d_k = fc.propagators_after_conservation[k]
                        .derivative(s_var);
        SpDecomposition dec = decompose_into_splist(fc, d_k);

        // Per-k sp_list coefficients and constant on the WIDE fc.ctx
        // (direct + chain), THEN substitute_and_narrow each.  The chain
        // term carries Mfrac coefficients C_{ai}(invariants), so we
        // accumulate on fc.ctx first.
        std::vector<Mfrac> sp_coef_wide;
        sp_coef_wide.reserve(dec.sp_coef.size());
        for (auto& mp : dec.sp_coef) {
            sp_coef_wide.push_back(Mfrac::from_mpoly(std::move(mp)));
        }
        Mfrac const_wide = Mfrac::from_mpoly(std::move(dec.constant));

        if (chain.active) {
            accumulate_chain_for_prop(fc,
                fc.propagators_after_conservation[k],
                chain, sp_index, sp_coef_wide, const_wide);
        }

        // Pre-substitute + reproject sp_coef_wide[l] / const_wide to target_ctx.
        std::vector<Mfrac> sp_coef_narrow;
        sp_coef_narrow.reserve(sp_coef_wide.size());
        for (auto& mf : sp_coef_wide) {
            sp_coef_narrow.push_back(
                algebra::substitute_and_narrow(mf, numeric_values, target_ctx));
        }
        Mfrac const_narrow =
            algebra::substitute_and_narrow(const_wide, numeric_values, target_ctx);

        std::vector<Mfrac> coef_row;
        coef_row.reserve((std::size_t)n_d);
        for (long j = 0; j < n_d; ++j) {
            Mfrac acc = Mfrac::zero(target_ctx);
            for (std::size_t l = 0; l < sp_coef_narrow.size(); ++l) {
                if (sp_coef_narrow[l].is_zero()) continue;
                Mfrac t = sp_coef_narrow[l].clone();
                t *= alpha_narrow[l][(std::size_t)j];
                acc += t;
            }
            coef_row.push_back(std::move(acc));
        }
        Mfrac const_term = const_narrow.clone();
        for (std::size_t l = 0; l < sp_coef_narrow.size(); ++l) {
            if (sp_coef_narrow[l].is_zero()) continue;
            Mfrac t = sp_coef_narrow[l].clone();
            t *= beta_narrow[l];
            const_term += t;
        }
        out.coef.push_back(std::move(coef_row));
        out.constant.push_back(std::move(const_term));
    }
    return out;
}

std::vector<DerivTerm>
libp_deriv(const qft::FamilyConfig& fc,
           const qft::JIntegral& j,
           const std::string& s_name,
           const std::map<std::string, std::string>& numeric_values,
           const std::shared_ptr<algebra::MpolyContext>& target_ctx) {
    auto dd = libp_denoms_deriv(fc, s_name, numeric_values, target_ctx);
    long n_d = (long)dd.completede.size();
    long N   = (long)fc.propagators_after_conservation.size();

    if ((long)j.n_indices() != N) {
        throw std::invalid_argument(
            "libp_deriv: J integral indices length != #propagators");
    }

    std::vector<DerivTerm> out;
    for (long k = 0; k < N; ++k) {
        long ak = j.indices()[(std::size_t)k];
        if (ak == 0) continue;

        // neg_ak is a constant rational on target_ctx (no symbolic
        // dependence — it's just `-ak`).
        Mfrac neg_ak = Mfrac::from_si(target_ctx, -ak);

        for (long jp = 0; jp < n_d; ++jp) {
            const Mfrac& cf = dd.coef[(std::size_t)k][(std::size_t)jp];
            if (cf.is_zero()) continue;
            std::vector<long> new_idx = j.indices();
            new_idx[(std::size_t)k] += 1;
            if (jp < N) {
                new_idx[(std::size_t)jp] -= 1;
            }
            qft::JIntegral shifted(j.family(), new_idx);
            DerivTerm dt;
            dt.integ = shifted;
            dt.coef = neg_ak.clone();
            dt.coef *= cf;
            out.push_back(std::move(dt));
        }

        const Mfrac& ck = dd.constant[(std::size_t)k];
        if (!ck.is_zero()) {
            std::vector<long> new_idx = j.indices();
            new_idx[(std::size_t)k] += 1;
            qft::JIntegral shifted(j.family(), new_idx);
            DerivTerm dt;
            dt.integ = shifted;
            dt.coef = neg_ak.clone();
            dt.coef *= ck;
            out.push_back(std::move(dt));
        }
    }
    return out;
}

}  // namespace amflow::ibp
