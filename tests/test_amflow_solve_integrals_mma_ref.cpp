// SPDX-License-Identifier: MIT
// Layer 17 MMA-reference differential tests, ported from
// tests/test_layer17_solve_integrals.cpp.

#include <algorithm>
#include <cmath>
#include <complex>
#include <filesystem>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include <flint/acb.h>
#include <flint/arb.h>
#include <flint/fmpq.h>

#include "amflow/pipeline/amfsystem.hpp"
#include "amflow/pipeline/solve_integrals.hpp"
#include "amflow/numeric/acb_wrap.hpp"
#include "amflow/numeric/options.hpp"
#include "amflow/qft/family_config.hpp"
#include "amflow/qft/jintegral.hpp"

#include "util/math_ref_util.hpp"

namespace pipeline = amflow::pipeline;
namespace numeric = amflow::numeric;
namespace qft     = amflow::qft;
namespace fs      = std::filesystem;
using amflow::refactor_test::acb_imag_mid;
using amflow::refactor_test::acb_real_mid;
using amflow::refactor_test::j_to_math_ref_key;
using amflow::refactor_test::kira_available;
using amflow::refactor_test::load_math_ref_laurent;
using amflow::refactor_test::load_math_ref_values;
using amflow::refactor_test::source_root;

namespace {

numeric::AcbValue make_rational(long num, long den) {
    fmpq_t q;
    fmpq_init(q);
    fmpq_set_si(q, num, den);
    numeric::AcbValue out;
    out.set_fmpq(q);
    fmpq_clear(q);
    return out;
}

::testing::AssertionResult acb_close_relative(const numeric::AcbValue& value,
                                                double re, double im,
                                                double rel_tol = 1e-12,
                                                double abs_floor = 1e-20) {
    const double got_re = acb_real_mid(value);
    const double got_im = acb_imag_mid(value);
    const double err_re = std::abs(got_re - re);
    const double err_im = std::abs(got_im - im);
    const double scale_re = std::max(std::abs(re), abs_floor);
    const double scale_im = std::max(std::abs(im), abs_floor);
    if (err_re <= rel_tol * scale_re && err_im <= rel_tol * scale_im) {
        return ::testing::AssertionSuccess();
    }
    return ::testing::AssertionFailure()
        << "got (" << got_re << ", " << got_im << ")"
        << " vs expected (" << re << ", " << im << ")"
        << " with relative tol " << rel_tol;
}

qft::FamilyConfig make_bubble_family() {
    return qft::FamilyConfig::build(
        "bubblefam", {"l"}, {"p"}, {}, {{"p^2", "s"}},
        {"l^2 - msq", "(l - p)^2 - msq"});
}

pipeline::AMFSystemOptions make_bubble_opts(const std::string& cache_root) {
    pipeline::AMFSystemOptions opts;
    opts.bb.kira_executable = "/usr/local/bin/kira";
    opts.bb.fermat_executable = "/usr/share/Ferl7/fer64";
    opts.bb.numeric_values["s"] = "5";
    opts.bb.numeric_values["msq"] = "1";
    opts.cache_root = cache_root;
    return opts;
}

qft::FamilyConfig make_box1_family() {
    return qft::FamilyConfig::build(
        "box1",
        {"l"}, {"p1", "p2", "p3", "p4"},
        {{"p4", "-p1 - p2 - p3"}},
        {{"p1^2", "0"}, {"p2^2", "0"}, {"p3^2", "0"}, {"p4^2", "0"},
         {"(p1 + p2)^2", "s"}, {"(p1 + p3)^2", "t"}},
        {"l^2", "(l + p1)^2", "(l + p1 + p2)^2", "(l + p1 + p2 + p4)^2"});
}

pipeline::AMFSystemOptions make_box1_opts(const std::string& cache_root) {
    pipeline::AMFSystemOptions opts;
    opts.bb.kira_executable = "/usr/local/bin/kira";
    opts.bb.fermat_executable = "/usr/share/Ferl7/fer64";
    opts.bb.numeric_values["s"] = "100";
    opts.bb.numeric_values["t"] = "-1";
    opts.bb.log_file = cache_root + "_kira.log";
    opts.cache_root = cache_root;
    return opts;
}

const std::vector<numeric::AcbValue>&
find_sampled_values(const std::vector<pipeline::SampledIntegralSolution>& sols,
                    const qft::JIntegral& target) {
    for (const auto& row : sols) {
        if (row.integral == target) return row.values;
    }
    throw std::runtime_error("sampled target not found");
}

const pipeline::LaurentIntegralSolution&
find_laurent_solution(const std::vector<pipeline::LaurentIntegralSolution>& sols,
                       const qft::JIntegral& target) {
    for (const auto& row : sols) {
        if (row.integral == target) return row;
    }
    throw std::runtime_error("laurent target not found");
}

}  // namespace

TEST(SolveIntegralsMmaRefTest, BlackBoxAMFlow_SplitTargetBubbleTadpolesMatchMathReferenceJson) {
    if (!kira_available()) {
        GTEST_SKIP() << "Kira/Fermat not installed";
    }

    numeric::GlobalScope s;
    s.global.working_pre = 120;
    s.global.chop_pre = 20;
    s.global.rationalize_pre = 20;
    s.global.silent_mode = true;
    s.expansion.x_order = 120;
    s.expansion.extra_x_order = 160;
    s.expansion.learn_x_order = -1;
    s.expansion.test_x_order = 5;
    s.commit();

    auto fc = make_bubble_family();
    auto opts = make_bubble_opts(
        (fs::temp_directory_path() / "amflow_v2_l17_bb_split").string());
    fs::remove_all(opts.cache_root);

    const std::vector<qft::JIntegral> targets = {
        qft::JIntegral("bubblefam", {1, 0}),
        qft::JIntegral("bubblefam", {0, 1}),
    };
    std::vector<numeric::AcbValue> eps_samples;
    eps_samples.push_back(make_rational(1, 100));

    const auto got = pipeline::black_box_amflow(
        fc, targets, eps_samples, opts, opts.cache_root + "_driver");
    ASSERT_EQ(got.size(), targets.size());

    const auto ref = load_math_ref_values(
        source_root() / "tests/data/math_ref/bubble_1L.json", 1.0 / 100.0);
    for (const auto& target : targets) {
        auto it = ref.find(j_to_math_ref_key(target));
        ASSERT_NE(it, ref.end());
        const auto& row = find_sampled_values(got, target);
        ASSERT_EQ(row.size(), 1u);
        EXPECT_NEAR(acb_real_mid(row[0]), it->second.real(), 1e-8)
            << target.to_string();
        EXPECT_NEAR(acb_imag_mid(row[0]), it->second.imag(), 1e-10)
            << target.to_string();
    }
}

TEST(SolveIntegralsMmaRefTest, SolveIntegrals_BubbleLaurentCoefficientsMatchMathReferenceJson) {
    if (!kira_available()) {
        GTEST_SKIP() << "Kira/Fermat not installed";
    }

    numeric::GlobalScope s;
    s.global.working_pre = 40;
    s.global.chop_pre = 20;
    s.global.rationalize_pre = 20;
    s.global.silent_mode = true;
    s.expansion.x_order = 20;
    s.expansion.extra_x_order = 20;
    s.expansion.learn_x_order = -1;
    s.expansion.test_x_order = 5;
    s.commit();

    auto fc = make_bubble_family();
    auto opts = make_bubble_opts(
        (fs::temp_directory_path() / "amflow_v2_l17_solve").string());
    fs::remove_all(opts.cache_root);

    const qft::JIntegral target("bubblefam", {1, 1});
    const auto got = pipeline::solve_integrals(
        fc, {target}, /*goal_digits=*/25, /*eps_order=*/2,
        opts, opts.cache_root + "_driver");

    ASSERT_EQ(got.size(), 1u);
    const auto& row = find_laurent_solution(got, target);
    EXPECT_EQ(row.leading_order, -1);
    ASSERT_EQ(row.coefficients.size(), 2u);

    const auto ref = load_math_ref_laurent(
        source_root() / "tests/data/math_ref/bubble_1L.json",
        j_to_math_ref_key(target));
    ASSERT_NE(ref.find(-1), ref.end());
    ASSERT_NE(ref.find(0), ref.end());

    EXPECT_NEAR(acb_real_mid(row.coefficients[0]), ref.at(-1).real(), 1e-7);
    EXPECT_NEAR(acb_imag_mid(row.coefficients[0]), ref.at(-1).imag(), 1e-10);
    EXPECT_NEAR(acb_real_mid(row.coefficients[1]), ref.at(0).real(), 1e-6);
    EXPECT_NEAR(acb_imag_mid(row.coefficients[1]), ref.at(0).imag(), 1e-6);
}

TEST(SolveIntegralsMmaRefTest, BlackBoxAMFlowSingle_Box1AutomaticLoopTargetsMatchMathReferenceJson) {
    if (!kira_available()) {
        GTEST_SKIP() << "Kira/Fermat not installed";
    }

    numeric::GlobalScope s;
    s.global.working_pre = 120;
    s.global.chop_pre = 20;
    s.global.rationalize_pre = 20;
    s.global.silent_mode = true;
    s.expansion.x_order = 120;
    s.expansion.extra_x_order = 160;
    s.expansion.learn_x_order = -1;
    s.expansion.test_x_order = 5;
    s.commit();

    auto fc = make_box1_family();
    auto opts = make_box1_opts(
        (fs::temp_directory_path() / "amflow_v2_l17_box1_auto_loop").string());
    fs::remove_all(opts.cache_root);

    const std::vector<qft::JIntegral> targets = {
        qft::JIntegral("box1", {1, 0, 1, 0}),
        qft::JIntegral("box1", {-2, 1, 1, 1}),
        qft::JIntegral("box1", {1, 1, 1, 1}),
    };
    std::vector<numeric::AcbValue> eps_samples;
    eps_samples.push_back(make_rational(1, 100));

    const auto got = pipeline::black_box_amflow_single(
        fc, targets, eps_samples, opts, opts.cache_root + "_driver");
    ASSERT_EQ(got.size(), targets.size());

    const auto ref = load_math_ref_values(
        source_root() / "tests/data/math_ref/box1_automatic_loop.json",
        1.0 / 100.0);
    for (const auto& target : targets) {
        auto it = ref.find(j_to_math_ref_key(target));
        ASSERT_NE(it, ref.end());
        const auto& row = find_sampled_values(got, target);
        ASSERT_EQ(row.size(), 1u);
        EXPECT_TRUE(acb_close_relative(
            row[0], it->second.real(), it->second.imag(), 5e-4, 1e-8))
            << target.to_string();
    }
}

TEST(SolveIntegralsMmaRefTest, SolveIntegrals_Box1SingleTargetLaurentCoefficientsMatchMathReferenceJson) {
    if (!kira_available()) {
        GTEST_SKIP() << "Kira/Fermat not installed";
    }

    numeric::GlobalScope s;
    s.global.chop_pre = 20;
    s.global.rationalize_pre = 20;
    s.global.silent_mode = true;
    s.expansion.learn_x_order = -1;
    s.expansion.test_x_order = 5;
    s.commit();

    auto fc = make_box1_family();
    auto opts = make_box1_opts(
        (fs::temp_directory_path() / "amflow_v2_l17_solve_box1_single").string());
    fs::remove_all(opts.cache_root);

    const qft::JIntegral target("box1", {1, 0, 1, 0});
    const auto got = pipeline::solve_integrals(
        fc, {target}, /*goal_digits=*/10, /*eps_order=*/3,
        opts, opts.cache_root + "_driver");

    ASSERT_EQ(got.size(), 1u);
    const auto& row = find_laurent_solution(got, target);
    const auto ref = load_math_ref_laurent(
        source_root() / "tests/data/math_ref/box1_automatic_loop.json",
        j_to_math_ref_key(target));
    ASSERT_FALSE(ref.empty()) << target.to_string();
    EXPECT_EQ(row.leading_order, ref.begin()->first) << target.to_string();
    ASSERT_EQ(row.coefficients.size(), ref.size()) << target.to_string();

    for (const auto& [order, want] : ref) {
        const auto idx = static_cast<std::size_t>(order - row.leading_order);
        EXPECT_TRUE(acb_close_relative(
            row.coefficients[idx], want.real(), want.imag(), 5e-3, 1e-8))
            << "order " << order << ", target " << target.to_string();
    }
}
