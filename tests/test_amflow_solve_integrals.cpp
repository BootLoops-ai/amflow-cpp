// SPDX-License-Identifier: MIT
// Tests for amflow::pipeline::solve_integrals.
//
// BlackBoxAMFlow / SolveIntegrals end-to-end paths require ode::amflow
// (Layer 9) which is currently a forward-declaration stub.  Tests below
// only exercise the local-only helpers:
//   - generate_numerical_config
//   - fit_eps
// Full integration coverage will land once Layer 9 is ported in Phase 6.

#include <complex>

#include <gtest/gtest.h>

#include <flint/acb.h>
#include <flint/arb.h>
#include <flint/arf.h>
#include <flint/fmpq.h>

#include "amflow/pipeline/solve_integrals.hpp"

namespace pipeline = amflow::pipeline;
namespace numeric = amflow::numeric;

namespace {

double acb_real_mid(const numeric::AcbValue& value) {
    return arf_get_d(arb_midref(acb_realref(value.raw())), ARF_RND_NEAR);
}

double acb_imag_mid(const numeric::AcbValue& value) {
    return arf_get_d(arb_midref(acb_imagref(value.raw())), ARF_RND_NEAR);
}

numeric::AcbValue make_real(double value) {
    numeric::AcbValue out;
    out.set_d_d(value, 0.0);
    return out;
}

numeric::AcbValue make_complex(double re, double im) {
    numeric::AcbValue out;
    out.set_d_d(re, im);
    return out;
}

}  // namespace

TEST(SolveIntegralsTest, GenerateNumericalConfig_SunriseMatchesMathSampling) {
    // AMFLOW_TUNE: this test pins the *legacy* AMFlow.m sampling
    // (n_eps = ceil(2.5*eps_order+2L), equispaced (100+i)/100 grid,
    // x_order = 4*single_pre).  The tuned default changes all three; gate on
    // AMFLOW_NEPS_LEGACY which restores bit-exact legacy behaviour.
    if (std::getenv("AMFLOW_NEPS_LEGACY") == nullptr) {
        GTEST_SKIP() << "legacy-sampling test; set AMFLOW_NEPS_LEGACY=1";
    }
    const auto cfg = pipeline::generate_numerical_config(/*loop_count=*/2,
                                                         /*goal_digits=*/25,
                                                         /*eps_order=*/3);

    EXPECT_EQ(cfg.working_pre, 232);
    EXPECT_EQ(cfg.x_order, 464);
    ASSERT_EQ(cfg.eps_samples.size(), 12u);

    const std::vector<std::pair<long, long>> expected = {
        {101, 1778279400},
        {17, 296379900},
        {103, 1778279400},
        {13, 222284925},
        {7, 118551960},
        {53, 889139700},
        {107, 1778279400},
        {1, 16465550},
        {109, 1778279400},
        {11, 177827940},
        {37, 592759800},
        {14, 222284925},
    };

    for (std::size_t i = 0; i < expected.size(); ++i) {
        const double want =
            static_cast<double>(expected[i].first) /
            static_cast<double>(expected[i].second);
        EXPECT_NEAR(acb_real_mid(cfg.eps_samples[i]), want, 1e-20)
            << "eps sample " << i;
        EXPECT_NEAR(acb_imag_mid(cfg.eps_samples[i]), 0.0, 1e-30)
            << "eps sample " << i << " imag part";
    }
}

TEST(SolveIntegralsTest, GenerateNumericalConfig_RejectsTooLargeOrder) {
    // AMFLOW_TUNE: tuned n_eps caps at the legacy value, so
    // eps_order=40 no longer overflows (n_eps=69<100).  The throw is a legacy
    // guard; verify it under AMFLOW_NEPS_LEGACY only.
    if (std::getenv("AMFLOW_NEPS_LEGACY") == nullptr) {
        EXPECT_NO_THROW(pipeline::generate_numerical_config(1, 25, 40));
        return;
    }
    EXPECT_THROW(pipeline::generate_numerical_config(/*loop_count=*/1,
                                                    /*goal_digits=*/25,
                                                    /*eps_order=*/40),
                 std::runtime_error);
}

TEST(SolveIntegralsTest, FitEps_RecoversThreeRealLaurentCoefficients) {
    std::vector<numeric::AcbValue> eps_samples;
    eps_samples.push_back(make_real(0.1));
    eps_samples.push_back(make_real(0.2));
    eps_samples.push_back(make_real(0.4));

    std::vector<numeric::AcbValue> values;
    for (double eps : {0.1, 0.2, 0.4}) {
        values.push_back(make_real(2.0 / (eps * eps) - 3.0 / eps + 5.0));
    }

    const auto coeffs =
        pipeline::fit_eps(eps_samples, values, /*leading_order=*/-2, 256);
    ASSERT_EQ(coeffs.size(), 3u);
    EXPECT_NEAR(acb_real_mid(coeffs[0]), 2.0, 1e-20);
    EXPECT_NEAR(acb_real_mid(coeffs[1]), -3.0, 1e-20);
    EXPECT_NEAR(acb_real_mid(coeffs[2]), 5.0, 1e-20);
    EXPECT_NEAR(acb_imag_mid(coeffs[0]), 0.0, 1e-20);
    EXPECT_NEAR(acb_imag_mid(coeffs[1]), 0.0, 1e-20);
    EXPECT_NEAR(acb_imag_mid(coeffs[2]), 0.0, 1e-20);
}

TEST(SolveIntegralsTest, FitEps_RecoversComplexLaurentCoefficients) {
    std::vector<numeric::AcbValue> eps_samples;
    eps_samples.push_back(make_real(0.25));
    eps_samples.push_back(make_real(0.5));

    const std::complex<double> c0(1.5, -0.5);
    const std::complex<double> c1(-2.0, 3.0);

    std::vector<numeric::AcbValue> values;
    for (double eps : {0.25, 0.5}) {
        const auto v = c0 / eps + c1;
        values.push_back(make_complex(v.real(), v.imag()));
    }

    const auto coeffs =
        pipeline::fit_eps(eps_samples, values, /*leading_order=*/-1, 256);
    ASSERT_EQ(coeffs.size(), 2u);
    EXPECT_NEAR(acb_real_mid(coeffs[0]), c0.real(), 1e-18);
    EXPECT_NEAR(acb_imag_mid(coeffs[0]), c0.imag(), 1e-18);
    EXPECT_NEAR(acb_real_mid(coeffs[1]), c1.real(), 1e-18);
    EXPECT_NEAR(acb_imag_mid(coeffs[1]), c1.imag(), 1e-18);
}
