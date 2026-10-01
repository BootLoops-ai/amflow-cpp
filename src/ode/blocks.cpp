// SPDX-License-Identifier: MIT
// ode::blocks — implementation.
//

#include "amflow/ode/blocks.hpp"

#include <algorithm>
#include <iterator>
#include <set>
#include <stdexcept>
#include <utility>

#include <flint/fmpq.h>
#include <flint/fmpq_poly.h>

#include "amflow/numeric/log.hpp"

namespace amflow::ode {

using numeric::AcbValue;
using numeric::FmpqPoly;
using numeric::RationalFunction;
using numeric::RationalMatrix;

// ===========================================================================
//  BlockPartition
// ===========================================================================

std::vector<std::size_t> BlockPartition::sub_rows(std::size_t i) const {
    std::vector<std::size_t> rows;
    for (auto j : sub_block_ids_.at(i)) {
        const auto& bj = blocks_.at(j);
        rows.insert(rows.end(), bj.begin(), bj.end());
    }
    return rows;
}

// ===========================================================================
//  AnalyzeBlockOld  (.m line 193-201)
// ===========================================================================

namespace {

template <class IsNonZero>
std::vector<std::vector<std::size_t>>
analyze_block_old_impl(std::size_t n, IsNonZero is_nonzero) {
    if (n == 0) return {};

    // Right-most nonzero column index in each row (0 = none, otherwise 1-based).
    std::vector<std::size_t> right(n);
    for (std::size_t i = 0; i < n; ++i) {
        std::size_t r = 0;
        for (std::size_t j = n; j-- > 0; ) {
            if (is_nonzero(i, j)) { r = j + 1; break; }
        }
        right[i] = r;
    }

    // Running max.
    std::size_t running = 0;
    for (std::size_t i = 0; i < n; ++i) {
        if (right[i] > running) running = right[i];
        right[i] = running;
    }

    // Markers: i+1 >= right[i]  =>  block ends at i (inclusive).
    std::vector<std::size_t> markers;
    for (std::size_t i = 0; i < n; ++i) {
        if (i + 1 >= right[i]) markers.push_back(i);
    }

    // Build inclusive ranges between consecutive markers.
    std::vector<std::vector<std::size_t>> blocks;
    std::size_t prev_end = 0;
    for (std::size_t k = 0; k < markers.size(); ++k) {
        std::size_t start = (k == 0) ? 0 : prev_end + 1;
        std::size_t end   = markers[k];
        std::vector<std::size_t> rng;
        for (std::size_t v = start; v <= end; ++v) rng.push_back(v);
        blocks.push_back(std::move(rng));
        prev_end = end;
    }
    return blocks;
}

}  // namespace

std::vector<std::vector<std::size_t>>
analyze_block_old(const RationalMatrix& mat) {
    auto pred = [&mat](std::size_t i, std::size_t j) { return !mat(i, j).is_zero(); };
    return analyze_block_old_impl(mat.rows(), pred);
}

std::vector<std::vector<std::size_t>>
analyze_block_old_pattern(std::size_t n,
                          std::function<bool(std::size_t, std::size_t)> is_nonzero) {
    return analyze_block_old_impl(n, is_nonzero);
}

// ===========================================================================
//  AnalyzeBlock  (.m line 204-237) — closure-based partition
// ===========================================================================

namespace {

bool same_set(const std::vector<std::size_t>& a, const std::vector<std::size_t>& b) {
    if (a.size() != b.size()) return false;
    for (auto v : a) if (std::find(b.begin(), b.end(), v) == b.end()) return false;
    return true;
}

bool subset_of(const std::vector<std::size_t>& s,
               const std::vector<std::size_t>& t) {
    if (s.size() > t.size()) return false;
    for (auto v : s) {
        if (!std::binary_search(t.begin(), t.end(), v)) return false;
    }
    return true;
}

std::vector<std::size_t> sorted_unique(std::vector<std::size_t> v) {
    std::sort(v.begin(), v.end());
    v.erase(std::unique(v.begin(), v.end()), v.end());
    return v;
}

}  // namespace

std::vector<std::vector<std::size_t>>
analyze_block(const RationalMatrix& mat) {
    const std::size_t n = mat.rows();
    if (n == 0) return {};
    if (mat.cols() != n) {
        throw std::invalid_argument("analyze_block: matrix must be square");
    }

    // 1. subint[i] = { j : mat(i,j) != 0 } ∪ {i}.
    std::vector<std::vector<std::size_t>> subint(n);
    for (std::size_t i = 0; i < n; ++i) {
        std::vector<std::size_t> tmp;
        for (std::size_t j = 0; j < n; ++j) {
            if (!mat(i, j).is_zero()) tmp.push_back(j);
        }
        tmp.push_back(i);
        subint[i] = sorted_unique(std::move(tmp));
    }

    // Forward closure: iterate `bl ← ∪ subint[i] for i in bl` until fixpoint.
    // Matches MMA `extend` (DESolver.m:240) verbatim.  Do NOT additionally
    // filter `dep` by "subint[j] ∩ bl ≠ ∅" — that discards forward-only
    // edges (j has no back-edge into bl) and shrinks closures to
    // fragments of an SCC, which lets oversized blocks reach
    // `determine_block_boundary_order` and over-estimate boundary orders.
    auto extend = [&](std::vector<std::size_t> bl) -> std::vector<std::size_t> {
        bl = sorted_unique(std::move(bl));
        for (;;) {
            std::vector<std::size_t> next_;
            for (auto i : bl) {
                next_.insert(next_.end(), subint[i].begin(), subint[i].end());
            }
            next_ = sorted_unique(std::move(next_));
            if (next_ == bl) return bl;
            bl = std::move(next_);
        }
    };

    // 2. Seed each row, take closure.
    std::vector<std::vector<std::size_t>> blocks(n);
    for (std::size_t i = 0; i < n; ++i) blocks[i] = extend({i});

    // 3. Dedupe closures by SAME-SET (.m line 245:
    //    `First/@Gather[blocks, samesetQ]`).  Do NOT merge by
    //    subset/superset — that collapses nested-closure chains
    //    ({4} ⊂ {4,5} ⊂ {4,5,6} ⊂ …) into a single big block instead
    //    of n singletons, which over-couples masters in
    //    `determine_block_boundary_order` and produces incorrect orders.
    std::vector<std::vector<std::size_t>> unique_closures;
    for (std::size_t i = 0; i < blocks.size(); ++i) {
        bool dup = false;
        for (const auto& uc : unique_closures) {
            if (same_set(blocks[i], uc)) { dup = true; break; }
        }
        if (!dup) unique_closures.push_back(blocks[i]);
    }

    // 4. MMA `sortblocks` (.m line 242):  iteratively peel off "top blocks",
    //    i.e. closures that have no strict superset among the remaining set.
    //    Output order: largest closures first, smaller ones last.
    auto is_topblock = [&](const std::vector<std::vector<std::size_t>>& bls,
                            const std::vector<std::size_t>& block) -> bool {
        for (const auto& b : bls) {
            if (b.size() == block.size() && same_set(b, block)) continue;
            if (subset_of(block, b)) return false;  // strict superset exists
        }
        return true;
    };

    std::vector<std::vector<std::size_t>> sorted;
    std::vector<std::vector<std::size_t>> remaining_bls = std::move(unique_closures);
    while (!remaining_bls.empty()) {
        std::vector<std::vector<std::size_t>> tops;
        std::vector<std::vector<std::size_t>> rest;
        for (auto& b : remaining_bls) {
            if (is_topblock(remaining_bls, b)) tops.push_back(std::move(b));
            else                                rest.push_back(std::move(b));
        }
        if (tops.empty()) {
            numeric::log_line(
                "analyze_block: sortblocks stalled; emitting remaining as-is");
            for (auto& r : rest) sorted.push_back(std::move(r));
            break;
        }
        for (auto& t : tops) sorted.push_back(std::move(t));
        remaining_bls = std::move(rest);
    }

    // 5. .m line 247: `Reverse[Table[Complement[Sequence@@blocks[[i;;]]], ...]]`.
    //    Each block becomes its members MINUS all subsequent (smaller) blocks'
    //    members, then the whole list is reversed (smallest first).
    //    This is the key step that turns nested closures into disjoint SCCs.
    std::vector<std::vector<std::size_t>> trimmed(sorted.size());
    for (std::size_t i = 0; i < sorted.size(); ++i) {
        std::vector<std::size_t> downstream;
        for (std::size_t j = i + 1; j < sorted.size(); ++j) {
            downstream.insert(downstream.end(),
                              sorted[j].begin(), sorted[j].end());
        }
        downstream = sorted_unique(std::move(downstream));
        for (auto v : sorted[i]) {
            if (!std::binary_search(downstream.begin(), downstream.end(), v)) {
                trimmed[i].push_back(v);
            }
        }
    }
    std::reverse(trimmed.begin(), trimmed.end());

    // 6. Cover/duplication sanity check (matches .m guard at line 248).
    {
        std::vector<std::size_t> all;
        for (const auto& r : trimmed) all.insert(all.end(), r.begin(), r.end());
        all = sorted_unique(std::move(all));
        std::vector<std::size_t> expected(n);
        for (std::size_t i = 0; i < n; ++i) expected[i] = i;
        if (all != expected) {
            throw std::runtime_error(
                "analyze_block: bad blocks (cover or duplication mismatch)");
        }
    }

    return trimmed;
}

// ===========================================================================
//  SubBlockID
// ===========================================================================

std::vector<std::vector<std::size_t>>
sub_block_ids(const RationalMatrix& mat,
              const std::vector<std::vector<std::size_t>>& blocks) {
    std::vector<std::vector<std::size_t>> result(blocks.size());
    for (std::size_t i = 0; i < blocks.size(); ++i) {
        for (std::size_t j = 0; j < i; ++j) {
            bool any = false;
            for (auto r : blocks[i]) {
                for (auto c : blocks[j]) {
                    if (!mat(r, c).is_zero()) { any = true; break; }
                }
                if (any) break;
            }
            if (any) result[i].push_back(j);
        }
    }
    return result;
}

BlockPartition build_partition(const RationalMatrix& mat) {
    BlockPartition p;
    p.set_blocks(analyze_block(mat));
    p.set_sub_block_ids(sub_block_ids(mat, p.blocks()));
    return p;
}

// ===========================================================================
//  PoincareRank
// ===========================================================================

long poincare_rank(const RationalMatrix& mat) {
    long max_v = 0;
    bool any = false;
    for (std::size_t i = 0; i < mat.rows(); ++i) {
        for (std::size_t j = 0; j < mat.cols(); ++j) {
            if (mat(i, j).is_zero()) continue;
            long v = mat(i, j).denominator().valuation();
            if (v < 0) v = 0;
            if (!any || v > max_v) { max_v = v; any = true; }
        }
    }
    return any ? (max_v - 1) : -1;
}

// ===========================================================================
//  NHEquations  (.m line 271-288)
// ===========================================================================

namespace {

FmpqPoly multiply_poly(const FmpqPoly& a, const FmpqPoly& b) {
    FmpqPoly out;
    fmpq_poly_mul(out.raw(), a.raw(), b.raw());
    return out;
}

// PolynomialLCM over a list, ignoring zero entries.  Result is monic
// (FLINT's fmpq_poly_lcm normalises).
FmpqPoly poly_lcm(const std::vector<const FmpqPoly*>& polys) {
    FmpqPoly out;
    bool init = false;
    for (auto* p : polys) {
        if (p == nullptr || p->is_zero()) continue;
        if (!init) {
            fmpq_poly_set(out.raw(), p->raw());
            init = true;
        } else {
            fmpq_poly_lcm(out.raw(), out.raw(), p->raw());
        }
    }
    if (!init) out.set_one();
    return out;
}

FmpqPoly factor_poly(EquationMode mode, long rank) {
    FmpqPoly f;
    switch (mode) {
        case EquationMode::Singular:
            f.set_coeff_si(1, 1);          // eta
            break;
        case EquationMode::Regular:
            f.set_one();
            break;
        case EquationMode::Taylor: {
            // factor = eta^(PoincareRank + 1) — mirrors upstream
            // DESolver.m line 279.  For rank = -1 (regular matrix) this
            // is eta^0 = 1, NOT eta^1.
            long power = rank + 1;
            if (power <= 0) {
                f.set_one();
            } else {
                f.set_coeff_si(power, 1);
            }
            break;
        }
    }
    return f;
}

}  // namespace

std::vector<BlockEquation>
nh_equations(const RationalMatrix& mat,
             const BlockPartition& partition,
             EquationMode mode) {
    if (mat.rows() != mat.cols()) {
        throw std::invalid_argument("nh_equations: matrix must be square");
    }
    long rank = (mode == EquationMode::Taylor) ? poincare_rank(mat) : 0;

    FmpqPoly factor = factor_poly(mode, rank);
    RationalFunction factor_rf = RationalFunction::from_polynomial(factor.clone());

    std::vector<BlockEquation> eqs;
    eqs.reserve(partition.size());

    for (std::size_t i = 0; i < partition.size(); ++i) {
        const auto& blk = partition.block(i);
        auto sub_rows = partition.sub_rows(i);

        // Compute factor * mat[blk, blk] symbolically.
        std::vector<RationalFunction> on_diag(blk.size() * blk.size());
        for (std::size_t a = 0; a < blk.size(); ++a) {
            for (std::size_t b = 0; b < blk.size(); ++b) {
                on_diag[a * blk.size() + b] = mat(blk[a], blk[b]) * factor_rf;
            }
        }

        // dx = LCM of denominators of on_diag entries.
        std::vector<const FmpqPoly*> denoms;
        denoms.reserve(on_diag.size());
        for (const auto& e : on_diag) {
            if (!e.is_zero()) denoms.push_back(&e.denominator());
        }
        FmpqPoly dx = poly_lcm(denoms);

        // ax = factor * dx * mat[blk, blk].  Result is polynomial because
        // dx contains every denominator.
        std::vector<std::vector<FmpqPoly>> ax(blk.size());
        for (auto& row : ax) row.resize(blk.size());

        RationalFunction dx_rf = RationalFunction::from_polynomial(dx.clone());
        for (std::size_t a = 0; a < blk.size(); ++a) {
            for (std::size_t b = 0; b < blk.size(); ++b) {
                RationalFunction val = on_diag[a * blk.size() + b] * dx_rf;
                if (!val.is_polynomial()) {
                    throw std::runtime_error(
                        "nh_equations: ax entry is not polynomial — "
                        "block analysis or LCM step is inconsistent");
                }
                // Numerator may have a constant denominator scalar baked
                // in, but is_polynomial() means denominator is monic
                // constant 1, so the numerator IS the value.
                ax[a][b] = val.numerator().clone();
            }
        }

        // bx = factor * dx * mat[blk, sub_rows].  Still rational (sub_rows
        // entries are off-block).
        std::vector<std::vector<RationalFunction>> bx(blk.size());
        for (auto& row : bx) row.resize(sub_rows.size());
        for (std::size_t a = 0; a < blk.size(); ++a) {
            for (std::size_t k = 0; k < sub_rows.size(); ++k) {
                bx[a][k] = mat(blk[a], sub_rows[k]) * factor_rf * dx_rf;
            }
        }

        // For non-Singular modes, multiply dx by factor at the end (after
        // ax/bx are computed using the original dx).  Mirrors .m line 286.
        if (mode != EquationMode::Singular) {
            FmpqPoly dx_new = multiply_poly(dx, factor);
            dx = std::move(dx_new);
        }

        BlockEquation eq;
        eq.dx    = std::move(dx);
        eq.ax    = std::move(ax);
        eq.bx    = std::move(bx);
        eq.block = blk;
        eq.sub   = std::move(sub_rows);
        eqs.push_back(std::move(eq));
    }
    return eqs;
}

std::vector<BlockEquation>
nh_equations(const RationalMatrix& mat, EquationMode mode) {
    BlockPartition p = build_partition(mat);
    return nh_equations(mat, p, mode);
}

// ===========================================================================
//  ToNum / NHEquationsNum
// ===========================================================================

std::vector<AcbValue> to_num(const FmpqPoly& poly, long prec) {
    if (poly.is_zero()) return {};
    long len = poly.length();
    std::vector<AcbValue> out;
    out.reserve(len);
    fmpq_t c;
    fmpq_init(c);
    for (long k = 0; k < len; ++k) {
        AcbValue v;
        poly.coeff(k, c);
        v.set_fmpq(c, prec);
        out.push_back(std::move(v));
    }
    fmpq_clear(c);
    return out;
}

std::vector<BlockEquationNum>
nh_equations_num(const std::vector<BlockEquation>& eqs, long prec) {
    std::vector<BlockEquationNum> out;
    out.reserve(eqs.size());

    auto poly_mat_to_num = [&](const std::vector<std::vector<FmpqPoly>>& m) {
        std::vector<std::vector<std::vector<AcbValue>>> r(m.size());
        for (std::size_t i = 0; i < m.size(); ++i) {
            r[i].resize(m[i].size());
            for (std::size_t j = 0; j < m[i].size(); ++j) {
                r[i][j] = to_num(m[i][j], prec);
            }
        }
        return r;
    };

    auto rat_mat_num_part = [&](const std::vector<std::vector<RationalFunction>>& m,
                                bool numerator_part) {
        std::vector<std::vector<std::vector<AcbValue>>> r(m.size());
        for (std::size_t i = 0; i < m.size(); ++i) {
            r[i].resize(m[i].size());
            for (std::size_t j = 0; j < m[i].size(); ++j) {
                const FmpqPoly& p = numerator_part ? m[i][j].numerator()
                                                   : m[i][j].denominator();
                r[i][j] = to_num(p, prec);
            }
        }
        return r;
    };

    for (const auto& eq : eqs) {
        BlockEquationNum n;
        n.dxexp  = to_num(eq.dx, prec);
        n.axexp  = poly_mat_to_num(eq.ax);
        n.bxexpn = rat_mat_num_part(eq.bx, true);
        n.bxexpd = rat_mat_num_part(eq.bx, false);
        n.block  = eq.block;
        n.sub    = eq.sub;
        out.push_back(std::move(n));
    }
    return out;
}

// ===========================================================================
//  PickElement / PickList / PickMat
// ===========================================================================

AcbValue pick_element(const std::vector<AcbValue>& list, std::size_t order) {
    if (order >= list.size()) return AcbValue();
    return list[order].clone();
}

std::vector<AcbValue>
pick_list(const std::vector<std::vector<AcbValue>>& lists, std::size_t order) {
    std::vector<AcbValue> out;
    out.reserve(lists.size());
    for (const auto& l : lists) out.push_back(pick_element(l, order));
    return out;
}

std::vector<std::vector<AcbValue>>
pick_mat(const std::vector<std::vector<std::vector<AcbValue>>>& mat,
         std::size_t order) {
    std::vector<std::vector<AcbValue>> out;
    out.reserve(mat.size());
    for (const auto& row : mat) out.push_back(pick_list(row, order));
    return out;
}

}  // namespace amflow::ode
