// SPDX-License-Identifier: MIT
// Tests for amflow::pipeline::factorize.

#include <gtest/gtest.h>

#include <flint/fmpz.h>
#include <flint/fmpz_mpoly.h>

#include "amflow/algebra/mpoly.hpp"
#include "amflow/pipeline/factorize.hpp"
#include "amflow/qft/family_config.hpp"

namespace alg     = amflow::algebra;
namespace pipeline = amflow::pipeline;
namespace qft     = amflow::qft;

namespace {

std::vector<alg::Mpoly>
zero_legs(const qft::FamilyConfig& fc) {
    std::vector<alg::Mpoly> props;
    long n_loop = (long)fc.n_loops();
    long n_red  = (long)fc.n_red_legs();
    fmpz_t zero;
    fmpz_init(zero);
    for (const auto& p : fc.propagators_after_conservation) {
        alg::Mpoly tmp = p.clone();
        for (long k = 0; k < n_red; ++k) {
            alg::Mpoly out(fc.ctx);
            fmpz_mpoly_evaluate_one_fmpz(out.raw(), tmp.raw(),
                                          n_loop + k, zero, fc.ctx->raw());
            tmp = std::move(out);
        }
        props.push_back(std::move(tmp));
    }
    fmpz_clear(zero);
    return props;
}

}  // namespace

TEST(FactorizeTest, OneLoopBubble_Identity) {
    auto fc = qft::FamilyConfig::build(
        "bubble", {"l"}, {"p"}, {}, {{"p^2", "s"}},
        {"l^2 - msq", "(l - p)^2 - msq"});

    auto props = zero_legs(fc);
    std::vector<std::vector<long>> patts = {{1, 1}};
    auto fams = pipeline::factorize_family(fc, props, patts);

    ASSERT_EQ(fams.size(), 1u);
    EXPECT_EQ(fams[0].loops.size(), 1u);
    EXPECT_EQ(fams[0].loops[0], "l");
    EXPECT_EQ(fams[0].propagators.size(), 2u);
    EXPECT_EQ(fams[0].patterns.size(), 1u);
    EXPECT_EQ(fams[0].patterns[0].size(), 2u);
    EXPECT_EQ(fams[0].patterns[0][0], 1);
    EXPECT_EQ(fams[0].patterns[0][1], 1);
}

TEST(FactorizeTest, SortBy_PicksMonomialContainingMassOne) {
    auto fc = qft::FamilyConfig::build(
        "two_indep_loops", {"l1", "l2"}, {}, {}, {},
        {"l1^2 - 1", "l2^2"});
    std::vector<alg::Mpoly> props;
    for (const auto& p : fc.propagators_after_conservation) {
        props.push_back(p.clone());
    }
    std::vector<std::vector<long>> patts = {{1, 1}};
    auto fams = pipeline::factorize_family(fc, props, patts);
    ASSERT_EQ(fams.size(), 2u);
    for (auto& f : fams) {
        EXPECT_EQ(f.loops.size(), 1u);
        EXPECT_EQ(f.propagators.size(), 1u);
    }
}

TEST(FactorizeTest, TwoLoopSunrise_OneComponent) {
    auto fc = qft::FamilyConfig::build(
        "sunrise", {"l1", "l2"}, {"p"}, {}, {{"p^2", "s"}},
        {"l1^2 - m1sq", "l2^2 - m2sq", "(l1 + l2 - p)^2 - m3sq"});

    auto props = zero_legs(fc);
    std::vector<std::vector<long>> patts = {{1, 1, 1}};
    auto fams = pipeline::factorize_family(fc, props, patts);

    ASSERT_EQ(fams.size(), 1u);
    EXPECT_EQ(fams[0].loops.size(), 2u);
    EXPECT_EQ(fams[0].propagators.size(), 3u);
    EXPECT_EQ(fams[0].patterns.size(), 1u);
    EXPECT_EQ(fams[0].patterns[0].size(), 3u);
}

// cross-loop bilinear mass=-1 propagator
// (`(l1-l2)^2 - 1`) inside a 4-loop vacuum top sector.  Mirrors the
// SingleMass input that arises on the region-1 boundary of
// `sunset_bubble_4L_1mass_l4`.  After factorize_family applies its
// loop redefinition, the mass=-1 propagator MUST be a single-loop
// square — otherwise `find_loop_in_prop` later picks the wrong loop
// to promote to leg in SingleMass, and the corner integral collapses
// to a sub-sector value.  Pre-fix this assertion failed because the
// redef matrix was built from `coeff_of(j,1).constant_term()` (which
// silently drops cross-loop bilinear contributions), producing a
// GL(4,Z) automorphism that left the mass propagator cross-loop.
TEST(FactorizeTest, CrossLoopMassProp_RedefIsSingleLoopSquared) {
    auto fc = qft::FamilyConfig::build(
        "cross_loop_mass",
        {"l1", "l2", "l3", "l4"}, {}, {}, {},
        {"l2^2", "l3^2", "l4^2", "(l1 + l3 + l4)^2", "(l1 - l2)^2 - 1"});

    std::vector<alg::Mpoly> props;
    for (const auto& p : fc.propagators_after_conservation) {
        props.push_back(p.clone());
    }
    std::vector<std::vector<long>> patts = {{1, 1, 1, 1, 1}};
    auto fams = pipeline::factorize_family(fc, props, patts);

    ASSERT_EQ(fams.size(), 1u);
    auto& f = fams[0];
    ASSERT_EQ(f.propagators.size(), 5u);

    auto count_quadratic_loops = [&](const alg::Mpoly& p) {
        long n = 0;
        for (long j = 0; j < (long)f.loops.size(); ++j) {
            if (!p.coeff_of(j, 2).is_zero()) ++n;
        }
        return n;
    };

    long mass_minus_one_idx = -1;
    for (long i = 0; i < (long)f.propagators.size(); ++i) {
        const alg::Mpoly& p = f.propagators[i];
        std::vector<unsigned long> zero_exp(
            (std::size_t)p.ctx()->n_vars(), 0);
        fmpz_t c;
        fmpz_init(c);
        fmpz_mpoly_get_coeff_fmpz_ui(c, p.raw(),
                                      zero_exp.data(),
                                      p.ctx()->raw());
        bool const_is_minus_one = (fmpz_cmp_si(c, -1) == 0);
        fmpz_clear(c);
        if (const_is_minus_one && count_quadratic_loops(p) >= 1) {
            mass_minus_one_idx = i;
            break;
        }
    }
    ASSERT_GE(mass_minus_one_idx, 0)
        << "Did not find a mass=-1 propagator after factorize_family.";

    EXPECT_EQ(count_quadratic_loops(f.propagators[(std::size_t)mass_minus_one_idx]), 1)
        << "After factorize_family, the mass=-1 propagator must be "
           "single-loop² (one loop appearing quadratically).  Cross-loop "
           "form means SingleMass will pick the wrong droploop.";
}
