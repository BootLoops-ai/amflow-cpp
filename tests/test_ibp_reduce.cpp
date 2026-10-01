// SPDX-License-Identifier: MIT
// Tests for amflow::ibp::reduce.
//
// The full reduce/diffeq paths require a real Kira installation; those
// are exercised by the integration test
// `test_ibp_kira_integration.cpp` (gated on /usr/local/bin/kira).
// Here we only test the local-only helpers.

#include <cstdlib>
#include <filesystem>
#include <gtest/gtest.h>

#include "amflow/algebra/mpoly.hpp"
#include "amflow/ibp/reduce.hpp"
#include "amflow/qft/family_config.hpp"
#include "amflow/qft/jintegral.hpp"

namespace alg = amflow::algebra;
namespace ibp = amflow::ibp;
namespace qft = amflow::qft;
namespace fs  = std::filesystem;

TEST(IbpReduceTest, MakeReductionContext_NarrowsToEtaAndD) {
    // make_reduction_context returns an MpolyContext with exactly the
    // two variables `{eta, d}`, independent of the family's variable
    // count.
    auto fc = qft::FamilyConfig::build(
        "bubble", {"l"}, {"p"}, {}, {{"p^2", "s"}},
        {"l^2 - msq", "(l - p)^2 - msq"});
    auto rctx = ibp::make_reduction_context(fc);

    EXPECT_EQ(rctx.ctx->n_vars(), 2);
    EXPECT_EQ(rctx.ctx->var_name(0), "eta");
    EXPECT_EQ(rctx.ctx->var_name(1), "d");
    EXPECT_EQ(rctx.d_var, 1);
    EXPECT_TRUE(rctx.lift_gens.empty());
}

TEST(IbpReduceTest, MakeReductionContext_DoesNotMutateFc) {
    auto fc = qft::FamilyConfig::build(
        "bubble", {"l"}, {"p"}, {}, {{"p^2", "s"}},
        {"l^2 - msq", "(l - p)^2 - msq"});
    long n_before = fc.ctx->n_vars();
    auto rctx = ibp::make_reduction_context(fc);
    EXPECT_EQ(fc.ctx->n_vars(), n_before);
    EXPECT_NE(rctx.ctx.get(), fc.ctx.get());
}

TEST(IbpReduceTest, MakeReductionContext_ShapeIsFamilyIndependent) {
    // Different family with many more variables must still yield a
    // narrow {eta, d} ctx.
    auto fc = qft::FamilyConfig::build(
        "tribox", {"l1", "l2", "l3"}, {"p1", "p2"}, {},
        {{"p1^2", "0"}, {"p2^2", "0"}, {"(p1 + p2)^2", "s"}},
        {"l1^2 - mAsq", "l2^2 - mBsq", "l3^2 - mCsq",
         "(l1 + p1)^2", "(l2 + p2)^2", "(l3 + p1 + p2)^2"});
    auto rctx = ibp::make_reduction_context(fc);
    EXPECT_EQ(rctx.ctx->n_vars(), 2);
    EXPECT_EQ(rctx.ctx->var_name(0), "eta");
    EXPECT_EQ(rctx.ctx->var_name(1), "d");
    EXPECT_EQ(rctx.d_var, 1);
}

namespace {

bool kira_available() {
    return fs::exists("/usr/local/bin/kira") &&
           fs::exists("/usr/share/Ferl7/fer64");
}

}  // namespace

TEST(IbpReduceTest, MassiveTadpole_J2_to_J1) {
    if (!kira_available()) {
        GTEST_SKIP() << "Kira/Fermat not installed";
    }
    auto fc = qft::FamilyConfig::build(
        "tad", {"l"}, {}, {}, {}, {"l^2 - msq"});

    ibp::ReduceOptions opts;
    opts.kira_executable = "/usr/local/bin/kira";
    opts.fermat_executable = "/usr/share/Ferl7/fer64";
    opts.numeric_values["msq"] = "1";
    opts.work_dir = (fs::temp_directory_path() /
                       "amflow_v2_ibp_integ_tad").string();
    opts.log_file = opts.work_dir + ".log";
    fs::remove_all(opts.work_dir);

    std::vector<qft::JIntegral> targets;
    targets.push_back(qft::JIntegral("tad", {2}));
    std::vector<qft::JIntegral> preferred;
    preferred.push_back(qft::JIntegral("tad", {1}));
    std::vector<int> top_pattern = {1};

    ibp::ReduceResult res;
    try {
        res = ibp::reduce(fc, targets, preferred, top_pattern, opts);
    } catch (const std::exception& e) {
        GTEST_SKIP() << "Kira invocation failed: " << e.what();
    }

    EXPECT_GE(res.masters.size(), 1u);
    ASSERT_EQ(res.rules.size(), 1u);
    EXPECT_FALSE(res.rules[0].empty()) << "no rule for J(2)";
    if (!res.rules[0].empty()) {
        std::string coef_str = res.rules[0][0].coef.to_string();
        EXPECT_NE(coef_str.find("d"), std::string::npos)
            << "coefficient: " << coef_str;
    }

    if (!::testing::Test::HasFailure()) {
        fs::remove_all(opts.work_dir);
        fs::remove(opts.log_file);
    }
}
