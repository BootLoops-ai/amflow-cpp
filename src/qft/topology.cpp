// SPDX-License-Identifier: MIT
// qft::topology — implementation.
//

#include "amflow/qft/topology.hpp"

#include <algorithm>
#include <map>
#include <set>
#include <stdexcept>
#include <utility>

#include <flint/flint.h>
#include <flint/fmpq.h>
#include <flint/fmpz.h>
#include <flint/fmpz_mpoly.h>
#include <flint/fmpz_mpoly_factor.h>
#include <flint/fmpz_mpoly_q.h>

namespace amflow::qft {

using algebra::Mfrac;
using algebra::Mpoly;
using algebra::MpolyContext;

namespace {

const long kInvariantPrimes[] = {3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37,
                                  41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83};

Mfrac numeric_substitute(const Mfrac& f,
                          long first_x_var,
                          long n_x,
                          const std::vector<long>& extra_keep = {}) {
    auto ctx = f.ctx();
    Mfrac out = f.clone();
    long inv_index = 0;
    for (long i = 0; i < ctx->n_vars(); ++i) {
        if (i >= first_x_var && i < first_x_var + n_x) continue;
        if (std::find(extra_keep.begin(), extra_keep.end(), i)
                != extra_keep.end()) continue;
        long val = kInvariantPrimes[inv_index % (long)(sizeof(kInvariantPrimes)
                                                          / sizeof(long))];
        ++inv_index;
        out = out.substitute(i, val);
    }
    return out;
}

}  // namespace

std::vector<long> feynman_vars_in(const Mpoly& p, long first_x_var, long n_x) {
    auto ctx = p.ctx();
    std::set<long> seen;
    long len = fmpz_mpoly_length(p.raw(), ctx->raw());
    std::vector<unsigned long> exp((std::size_t)ctx->n_vars());
    for (long t = 0; t < len; ++t) {
        fmpz_mpoly_get_term_exp_ui(exp.data(), p.raw(), t, ctx->raw());
        for (long k = 0; k < n_x; ++k) {
            if (exp[(std::size_t)(first_x_var + k)] > 0) {
                seen.insert(first_x_var + k);
            }
        }
    }
    return std::vector<long>(seen.begin(), seen.end());
}

std::vector<long> feynman_vars_in(const Mfrac& p, long first_x_var, long n_x) {
    Mpoly num = p.numerator();
    Mpoly den = p.denominator();
    auto a = feynman_vars_in(num, first_x_var, n_x);
    auto b = feynman_vars_in(den, first_x_var, n_x);
    std::set<long> merged(a.begin(), a.end());
    merged.insert(b.begin(), b.end());
    return std::vector<long>(merged.begin(), merged.end());
}

std::vector<Mpoly> factor_mpoly(const Mpoly& p) {
    auto ctx = p.ctx();
    std::vector<Mpoly> out;

    if (p.is_zero()) {
        out.push_back(p.clone());
        return out;
    }

    fmpz_mpoly_factor_t fac;
    fmpz_mpoly_factor_init(fac, ctx->raw());

    int ok = fmpz_mpoly_factor(fac, p.raw(), ctx->raw());
    if (!ok) {
        fmpz_mpoly_factor_clear(fac, ctx->raw());
        out.push_back(p.clone());
        return out;
    }

    fmpz_t constant;
    fmpz_init(constant);
    fmpz_mpoly_factor_get_constant_fmpz(constant, fac, ctx->raw());

    if (!fmpz_is_one(constant)) {
        Mpoly c(ctx);
        fmpz_mpoly_set_fmpz(c.raw(), constant, ctx->raw());
        out.push_back(std::move(c));
    }

    long n_fac = (long)fmpz_mpoly_factor_length(fac, ctx->raw());
    for (long i = 0; i < n_fac; ++i) {
        Mpoly base(ctx);
        fmpz_mpoly_factor_get_base(base.raw(), fac, i, ctx->raw());
        long e = fmpz_mpoly_factor_get_exp_si(fac, i, ctx->raw());
        for (long k = 0; k < e; ++k) {
            out.push_back(base.clone());
        }
    }

    fmpz_clear(constant);
    fmpz_mpoly_factor_clear(fac, ctx->raw());
    return out;
}

bool zero_sector_q(const FamilyConfig& fc,
                   const std::vector<Mpoly>& denominators) {
    UFResult uf = evaluate_uf(fc, denominators);
    auto uf_ctx  = uf.uf_ctx;
    long first_x = uf.first_x_var;
    long n_x     = (long)denominators.size();

    if (uf.degenerate) {
        return true;
    }

    Mfrac u_num   = numeric_substitute(uf.u,  first_x, n_x);
    Mfrac f_num   = numeric_substitute(uf.f,  first_x, n_x);
    Mfrac f0_num  = numeric_substitute(uf.f0, first_x, n_x);

    Mfrac g = u_num.clone();
    g += f_num;
    Mfrac uf0 = u_num.clone();
    uf0 *= f0_num;
    g += uf0;

    Mpoly g_num = g.numerator();

    // Extended ctx with __zsq_var* placeholders.
    std::vector<std::string> ext_names;
    ext_names.reserve((std::size_t)uf_ctx->n_vars() + (std::size_t)n_x);
    for (long i = 0; i < uf_ctx->n_vars(); ++i) {
        ext_names.push_back(uf_ctx->var_name(i));
    }
    for (long i = 0; i < n_x; ++i) {
        ext_names.push_back("__zsq_var" + std::to_string(i + 1));
    }
    auto ext_ctx = std::make_shared<MpolyContext>(std::move(ext_names));
    long var_first = uf_ctx->n_vars();

    std::vector<Mpoly> lift_gens;
    lift_gens.reserve((std::size_t)uf_ctx->n_vars());
    for (long i = 0; i < uf_ctx->n_vars(); ++i) {
        lift_gens.push_back(Mpoly::variable(ext_ctx, i));
    }
    Mpoly g_ext(ext_ctx);
    {
        std::vector<fmpz_mpoly_struct*> ptrs(lift_gens.size());
        for (std::size_t i = 0; i < lift_gens.size(); ++i) {
            ptrs[i] = const_cast<fmpz_mpoly_struct*>(lift_gens[i].raw());
        }
        int rc = fmpz_mpoly_compose_fmpz_mpoly(
            g_ext.raw(), g_num.raw(), ptrs.data(),
            uf_ctx->raw(), ext_ctx->raw());
        if (rc == 0) {
            throw std::runtime_error("zero_sector_q: lift compose failed");
        }
    }

    Mpoly cri = -g_ext;
    for (long i = 0; i < n_x; ++i) {
        Mpoly dgdx = g_ext.derivative(first_x + i);
        Mpoly term = Mpoly::variable(ext_ctx, var_first + i)
                     * Mpoly::variable(ext_ctx, first_x + i)
                     * dgdx;
        cri += term;
    }

    long ext_nvars = ext_ctx->n_vars();
    long len = fmpz_mpoly_length(cri.raw(), ext_ctx->raw());
    std::vector<unsigned long> exp((std::size_t)ext_nvars);
    fmpz_t coeff;
    fmpz_init(coeff);

    std::map<std::vector<unsigned long>, Mpoly> rows;
    for (long t = 0; t < len; ++t) {
        fmpz_mpoly_get_term_coeff_fmpz(coeff, cri.raw(), t, ext_ctx->raw());
        fmpz_mpoly_get_term_exp_ui(exp.data(), cri.raw(), t, ext_ctx->raw());

        std::vector<unsigned long> key(exp);
        for (long i = 0; i < n_x; ++i) {
            key[(std::size_t)(var_first + i)] = 0;
        }
        std::vector<unsigned long> rest_exp((std::size_t)ext_nvars, 0);
        for (long i = 0; i < n_x; ++i) {
            rest_exp[(std::size_t)(var_first + i)] = exp[(std::size_t)(var_first + i)];
        }

        Mpoly mono = Mpoly::zero(ext_ctx);
        fmpz_mpoly_set_coeff_fmpz_ui(mono.raw(), coeff,
                                      rest_exp.data(), ext_ctx->raw());
        auto it = rows.find(key);
        if (it == rows.end()) {
            rows.emplace(key, std::move(mono));
        } else {
            it->second += mono;
        }
    }
    fmpz_clear(coeff);

    long n_rows = (long)rows.size();
    long n_cols = n_x + 1;
    if (n_rows == 0) {
        return true;
    }

    std::vector<fmpq*> Mptr((std::size_t)(n_rows * n_cols), nullptr);
    auto Midx = [&](long r, long c) -> fmpq* {
        return Mptr[(std::size_t)(r * n_cols + c)];
    };
    auto Mset = [&](long r, long c, fmpq* p) {
        Mptr[(std::size_t)(r * n_cols + c)] = p;
    };
    for (long r = 0; r < n_rows; ++r) {
        for (long c = 0; c < n_cols; ++c) {
            fmpq* p = (fmpq*)flint_malloc(sizeof(fmpq));
            fmpq_init(p);
            Mset(r, c, p);
        }
    }
    long row_idx = 0;
    for (auto& [key, row] : rows) {
        for (long i = 0; i < n_x; ++i) {
            Mpoly co = row.coeff_of(var_first + i, 1);
            if (!co.is_zero()) {
                fmpz_t c;
                fmpz_init(c);
                std::vector<unsigned long> z((std::size_t)ext_nvars, 0);
                fmpz_mpoly_get_coeff_fmpz_ui(c, co.raw(),
                                              z.data(), ext_ctx->raw());
                fmpq_set_fmpz(Midx(row_idx, i), c);
                fmpz_clear(c);
            }
        }
        Mpoly co0 = row.clone();
        fmpz_t zero;
        fmpz_init(zero);
        for (long i = 0; i < n_x; ++i) {
            Mpoly tmp(ext_ctx);
            fmpz_mpoly_evaluate_one_fmpz(tmp.raw(), co0.raw(),
                                          var_first + i, zero,
                                          ext_ctx->raw());
            co0 = std::move(tmp);
        }
        fmpz_clear(zero);
        if (!co0.is_zero()) {
            fmpz_t c;
            fmpz_init(c);
            std::vector<unsigned long> z((std::size_t)ext_nvars, 0);
            fmpz_mpoly_get_coeff_fmpz_ui(c, co0.raw(),
                                          z.data(), ext_ctx->raw());
            fmpq_set_fmpz(Midx(row_idx, n_cols - 1), c);
            fmpz_clear(c);
        }
        ++row_idx;
    }

    auto swap_rows = [&](long a, long b) {
        for (long c = 0; c < n_cols; ++c) {
            std::swap(Mptr[(std::size_t)(a * n_cols + c)],
                      Mptr[(std::size_t)(b * n_cols + c)]);
        }
    };
    for (long col = 0; col < n_cols; ++col) {
        long pivot_row = -1;
        for (long r = col; r < n_rows; ++r) {
            if (!fmpq_is_zero(Midx(r, col))) {
                pivot_row = r;
                break;
            }
        }
        if (pivot_row < 0) continue;
        if (pivot_row != col) swap_rows(col, pivot_row);

        fmpq_t pivot;
        fmpq_init(pivot);
        fmpq_set(pivot, Midx(col, col));
        for (long c2 = col; c2 < n_cols; ++c2) {
            fmpq_div(Midx(col, c2), Midx(col, c2), pivot);
        }
        for (long r = 0; r < n_rows; ++r) {
            if (r == col) continue;
            if (fmpq_is_zero(Midx(r, col))) continue;
            fmpq_t mult;
            fmpq_init(mult);
            fmpq_set(mult, Midx(r, col));
            for (long c2 = col; c2 < n_cols; ++c2) {
                fmpq_t prod;
                fmpq_init(prod);
                fmpq_mul(prod, Midx(col, c2), mult);
                fmpq_sub(Midx(r, c2), Midx(r, c2), prod);
                fmpq_clear(prod);
            }
            fmpq_clear(mult);
        }
        fmpq_clear(pivot);
    }

    bool found_unit_vec = false;
    for (long r = 0; r < n_rows; ++r) {
        bool all_zero_lhs = true;
        for (long c = 0; c < n_x; ++c) {
            if (!fmpq_is_zero(Midx(r, c))) {
                all_zero_lhs = false;
                break;
            }
        }
        if (all_zero_lhs && !fmpq_is_zero(Midx(r, n_cols - 1))) {
            found_unit_vec = true;
            break;
        }
    }

    for (auto* p : Mptr) {
        if (p) {
            fmpq_clear(p);
            flint_free(p);
        }
    }

    return !found_unit_vec;
}

namespace {

ComponentInfo make_component(const Mpoly& u0_factor,
                              const Mfrac& f,
                              const Mfrac& f0,
                              long first_x_var,
                              long n_x) {
    ComponentInfo ci;
    ci.u0  = Mfrac::from_mpoly(u0_factor.clone());
    ci.var = feynman_vars_in(u0_factor, first_x_var, n_x);

    auto ctx = u0_factor.ctx();
    long len = fmpz_mpoly_length(u0_factor.raw(), ctx->raw());
    if (len == 0) {
        ci.loopnum = 0;
    } else {
        std::vector<unsigned long> exp((std::size_t)ctx->n_vars());
        fmpz_mpoly_get_term_exp_ui(exp.data(), u0_factor.raw(), 0,
                                    ctx->raw());
        long count = 0;
        for (long k = 0; k < n_x; ++k) {
            if (exp[(std::size_t)(first_x_var + k)] > 0) {
                ++count;
            }
        }
        ci.loopnum = count;
    }

    Mpoly f0_num = f0.numerator();
    Mpoly f0_den = f0.denominator();
    auto uf_ctx = f0.ctx();
    for (long v : ci.var) {
        Mpoly num_co = f0_num.coeff_of(v, 1);
        Mpoly den_co = f0_den.clone();
        ci.mass.push_back(Mfrac(std::move(num_co), std::move(den_co)));
    }

    Mpoly f_num = f.numerator();
    long fn_len = fmpz_mpoly_length(f_num.raw(), uf_ctx->raw());
    bool vacQ = true;
    {
        std::vector<unsigned long> exp((std::size_t)uf_ctx->n_vars());
        for (long t = 0; t < fn_len; ++t) {
            fmpz_mpoly_get_term_exp_ui(exp.data(), f_num.raw(), t,
                                        uf_ctx->raw());
            long inter = 0;
            for (long v : ci.var) {
                if (exp[(std::size_t)v] > 0) ++inter;
            }
            if (inter > ci.loopnum) {
                vacQ = false;
                break;
            }
        }
    }
    ci.vacQ = vacQ;
    return ci;
}

}  // namespace

std::vector<ComponentInfo>
analyze_topology(const FamilyConfig& fc,
                 const std::vector<Mpoly>& denominators) {
    UFResult uf = evaluate_uf(fc, denominators);
    long first_x = uf.first_x_var;
    long n_x     = (long)denominators.size();

    if (uf.degenerate) {
        return {};
    }

    Mpoly u_num = uf.u.numerator();
    auto factors = factor_mpoly(u_num);

    std::vector<ComponentInfo> result;
    for (auto& fac : factors) {
        auto vars = feynman_vars_in(fac, first_x, n_x);
        if (vars.empty()) continue;

        ComponentInfo ci = make_component(fac, uf.f, uf.f0, first_x, n_x);
        if (ci.loopnum > 0) {
            result.push_back(std::move(ci));
        }
    }

    // Linear (worldline) propagator fix: a propagator linear in the loop
    // momenta (worldline 1/(2u·k), HQET, ...) contributes to F but NOT to U,
    // so its Feynman parameter x_k never appears in any U-factor and the
    // prop is silently dropped from every component.  Downstream, the
    // eta-candidate selection (amf_position) and the boundary recursion
    // never see it, so eta is never placed on a linear prop and the
    // boundary leaf integrals retain unresolved linear denominators.
    //
    // Attach each orphan x_k to the component that owns the loop momenta
    // appearing in prop_k.  Its mass entry is the f0 coefficient (zero for
    // a pure 2u·k term).  This is a strict no-op when every propagator is
    // quadratic in some loop (then no x_k is orphaned).
    {
        std::set<long> assigned;
        for (const auto& ci : result)
            for (long v : ci.var) assigned.insert(v);

        Mpoly f0_num = uf.f0.numerator();
        Mpoly f0_den = uf.f0.denominator();

        long n_loop = (long)fc.n_loops();
        for (long k = 0; k < n_x; ++k) {
            long xk = first_x + k;
            if (assigned.count(xk)) continue;
            // Which loops appear (at any degree) in this propagator?
            const Mpoly& d = denominators[(std::size_t)k];
            std::set<long> loops_in;
            {
                long len = fmpz_mpoly_length(d.raw(), fc.ctx->raw());
                std::vector<unsigned long> exp((std::size_t)fc.ctx->n_vars());
                for (long t = 0; t < len; ++t) {
                    fmpz_mpoly_get_term_exp_ui(exp.data(), d.raw(), t,
                                                fc.ctx->raw());
                    for (long li = 0; li < n_loop; ++li) {
                        if (exp[(std::size_t)li] > 0) loops_in.insert(li);
                    }
                }
            }
            if (loops_in.empty()) continue;   // pure constant prop (ignored)
            // Find the component whose var-set contains a Feynman param of
            // a quadratic prop sharing one of these loops.
            long target = -1;
            for (std::size_t ci_idx = 0; ci_idx < result.size() && target < 0;
                 ++ci_idx) {
                for (long v : result[ci_idx].var) {
                    long pk = v - first_x;
                    const Mpoly& dp = denominators[(std::size_t)pk];
                    for (long li : loops_in) {
                        if (!dp.coeff_of(li, 2).is_zero()
                                || !dp.coeff_of(li, 1).is_zero()) {
                            target = (long)ci_idx;
                            break;
                        }
                    }
                    if (target >= 0) break;
                }
            }
            if (target < 0 && !result.empty()) target = 0;
            if (target < 0) continue;
            result[(std::size_t)target].var.push_back(xk);
            Mpoly num_co = f0_num.coeff_of(xk, 1);
            result[(std::size_t)target].mass.push_back(
                Mfrac(std::move(num_co), f0_den.clone()));
            // Re-evaluate vacQ: the linear prop's x_k brings in F-terms
            // x_k * x_j (cross with quadratic params), which exceeds
            // loopnum and makes the component non-vacuum (correctly --
            // it now carries an external direction u).
            {
                Mpoly f_num = uf.f.numerator();
                auto uf_ctx = uf.f.ctx();
                long fn_len = fmpz_mpoly_length(f_num.raw(), uf_ctx->raw());
                std::vector<unsigned long> exp((std::size_t)uf_ctx->n_vars());
                bool vacQ = true;
                for (long t = 0; t < fn_len && vacQ; ++t) {
                    fmpz_mpoly_get_term_exp_ui(exp.data(), f_num.raw(), t,
                                                uf_ctx->raw());
                    long inter = 0;
                    for (long v : result[(std::size_t)target].var) {
                        if (exp[(std::size_t)v] > 0) ++inter;
                    }
                    if (inter > result[(std::size_t)target].loopnum) vacQ = false;
                }
                result[(std::size_t)target].vacQ = vacQ;
            }
        }
    }

    std::stable_sort(result.begin(), result.end(),
                     [](const ComponentInfo& a, const ComponentInfo& b) {
                         return a.loopnum > b.loopnum;
                     });

    return result;
}

}  // namespace amflow::qft
