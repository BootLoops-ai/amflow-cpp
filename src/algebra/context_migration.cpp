// SPDX-License-Identifier: MIT
// algebra::context_migration — implementation.
//
// New primitives.  Per
// notes/refactor_design/algebra.md §6, these consolidate the inline
// monomial-walking code that earlier ports duplicated across amflow / blackbox /
// complete / family_uf source files.
//
// Implementation strategy: walk monomials, build src_to_dst index map by
// name lookup, re-emit on dst.  Fast-path when src and dst share identity.

#include "amflow/algebra/context_migration.hpp"

#include <algorithm>
#include <map>
#include <stdexcept>
#include <utility>
#include <vector>

#include <flint/fmpz.h>
#include <flint/fmpz_mpoly.h>
#include <flint/fmpz_mpoly_q.h>

namespace amflow::algebra {

namespace {

// Build src_to_dst[i] : dst-ctx index for src var i, or -2 if missing in
// dst (caller decides: throw vs drop), -1 if dropped (sentinel-prefix
// match).  drop_prefixes is a list of prefixes that mark a missing src
// var as "drop instead of throw".
//
// Returns nullopt-equivalent:  if any var is -2 (missing AND not in
// drop_prefixes), throws std::invalid_argument unless allow_missing is
// true.
std::vector<long>
build_var_map(const MpolyContext& src,
              const MpolyContext& dst,
              const std::vector<std::string>& drop_prefixes,
              bool drop_unknown_with_prefixes) {
    std::map<std::string, long> name_to_dst;
    for (long i = 0; i < dst.n_vars(); ++i) {
        name_to_dst[dst.var_name(i)] = i;
    }

    std::vector<long> src_to_dst((std::size_t)src.n_vars(), -2);
    for (long i = 0; i < src.n_vars(); ++i) {
        const std::string& name = src.var_name(i);
        auto it = name_to_dst.find(name);
        if (it != name_to_dst.end()) {
            src_to_dst[(std::size_t)i] = it->second;
            continue;
        }
        if (drop_unknown_with_prefixes) {
            bool prefix_match = false;
            for (const auto& prefix : drop_prefixes) {
                if (name.compare(0, prefix.size(), prefix) == 0) {
                    prefix_match = true;
                    break;
                }
            }
            if (prefix_match) {
                src_to_dst[(std::size_t)i] = -1;
                continue;
            }
        }
        throw std::invalid_argument(
            std::string("algebra::context_migration: variable '") + name
            + "' missing in destination context");
    }
    return src_to_dst;
}

// Walk monomials of `p` (on src.ctx) and re-emit on dst, applying the
// supplied src_to_dst map.  Monomials that touch any var with map[v] == -1
// (drop) are skipped entirely.
Mpoly walk_mpoly(const Mpoly& p,
                 const std::shared_ptr<MpolyContext>& dst_ctx,
                 const std::vector<long>& src_to_dst) {
    auto src_ctx = p.ctx();
    Mpoly out(dst_ctx);

    long len = fmpz_mpoly_length(p.raw(), src_ctx->raw());
    std::vector<unsigned long> sexp((std::size_t)src_ctx->n_vars());
    std::vector<unsigned long> dexp((std::size_t)dst_ctx->n_vars(), 0);
    fmpz_t coeff;
    fmpz_init(coeff);

    for (long t = 0; t < len; ++t) {
        fmpz_mpoly_get_term_exp_ui(sexp.data(), p.raw(), t, src_ctx->raw());
        fmpz_mpoly_get_term_coeff_fmpz(coeff, p.raw(), t, src_ctx->raw());

        std::fill(dexp.begin(), dexp.end(), 0);
        bool drop = false;
        for (long v = 0; v < src_ctx->n_vars(); ++v) {
            if (sexp[(std::size_t)v] == 0) continue;
            long pi = src_to_dst[(std::size_t)v];
            if (pi == -1) { drop = true; break; }
            // pi == -2 already triggered a throw in build_var_map.
            dexp[(std::size_t)pi] = sexp[(std::size_t)v];
        }
        if (drop) continue;
        fmpz_mpoly_set_coeff_fmpz_ui(out.raw(), coeff,
                                      dexp.data(), dst_ctx->raw());
    }
    fmpz_clear(coeff);
    return out;
}

}  // namespace

// ---------------------------------------------------------------------------
//  Strict lift
// ---------------------------------------------------------------------------

Mpoly mpoly_to_ctx(const Mpoly& src,
                    const std::shared_ptr<MpolyContext>& dst) {
    if (src.ctx().get() == dst.get()) return src.clone();
    auto map = build_var_map(*src.ctx(), *dst, {}, false);
    return walk_mpoly(src, dst, map);
}

Mfrac mfrac_to_ctx(const Mfrac& src,
                    const std::shared_ptr<MpolyContext>& dst) {
    if (src.ctx().get() == dst.get()) return src.clone();
    auto map = build_var_map(*src.ctx(), *dst, {}, false);
    Mpoly num = walk_mpoly(src.numerator(), dst, map);
    Mpoly den = walk_mpoly(src.denominator(), dst, map);
    if (den.is_zero()) {
        throw std::runtime_error(
            "mfrac_to_ctx: denominator vanished after migration");
    }
    return Mfrac(std::move(num), std::move(den));
}

// ---------------------------------------------------------------------------
//  Lossy projection
// ---------------------------------------------------------------------------

Mpoly project_mpoly_dropping(const Mpoly& src,
                              const std::shared_ptr<MpolyContext>& dst,
                              std::vector<std::string> drop_prefixes) {
    if (src.ctx().get() == dst.get()) return src.clone();
    auto map = build_var_map(*src.ctx(), *dst, drop_prefixes, true);
    return walk_mpoly(src, dst, map);
}

Mfrac project_mfrac_dropping(const Mfrac& src,
                              const std::shared_ptr<MpolyContext>& dst,
                              std::vector<std::string> drop_prefixes) {
    if (src.ctx().get() == dst.get()) return src.clone();
    auto map = build_var_map(*src.ctx(), *dst, drop_prefixes, true);
    Mpoly num = walk_mpoly(src.numerator(), dst, map);
    Mpoly den = walk_mpoly(src.denominator(), dst, map);
    if (den.is_zero()) {
        throw std::runtime_error(
            "project_mfrac_dropping: denominator vanished after projection");
    }
    return Mfrac(std::move(num), std::move(den));
}

}  // namespace amflow::algebra
