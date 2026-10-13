// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Anthropic, PBC. Added by the BootLoops fork of AMFlow.cpp; see NOTICE and PATCHES.md.
// Tests for the per-eps rational frame of amflow::pipeline (the "eps"
// entry of the numeric substitution map built for each eps sample).
//
// With AMFLOW_DIFFEQ_NUMD=P/Q the differential-equation tables carry
// d = P/Q exactly, so the region exponents and boundary coefficients
// must be evaluated at eps = (4 - P/Q)/2 exactly as well; a sample
// rationalized to `rationalize_pre` decimal digits (a truncation) would
// sit in a different frame and the mismatch is pivoted on by the sparse
// solver once it exceeds the chop threshold.  Without the env var the
// entry is the sample's midpoint recovered as the simplest rational within
// the rationalization tolerance, which for a sample of 1/113 is 1/113
// itself.

#include <gtest/gtest.h>

#include <flint/acb.h>
#include <flint/arb.h>
#include <flint/fmpq.h>
#include <flint/fmpz.h>

#include <cstdlib>
#include <string>

#include "amflow/numeric/acb_wrap.hpp"
#include "amflow/numeric/options.hpp"
#include "amflow/pipeline/amfsystem.hpp"

namespace pipeline = amflow::pipeline;
namespace dsn      = amflow::numeric;

namespace {

// Save / restore one environment variable around a test body.
class EnvGuard {
public:
    explicit EnvGuard(const char* name) : name_(name) {
        const char* v = std::getenv(name);
        had_ = (v != nullptr);
        if (had_) old_ = v;
    }
    ~EnvGuard() {
        if (had_) ::setenv(name_.c_str(), old_.c_str(), 1);
        else      ::unsetenv(name_.c_str());
    }
    void set(const char* v) { ::setenv(name_.c_str(), v, 1); }
    void unset()            { ::unsetenv(name_.c_str()); }
private:
    std::string name_;
    bool        had_ = false;
    std::string old_;
};

class NumdFrameFixture : public ::testing::Test {
protected:
    void SetUp() override {
        dsn::set_default_options();
        auto g = dsn::global_options();
        g.working_pre     = 100;
        g.rationalize_pre = 60;
        dsn::set_global_options(g);
    }
    void TearDown() override { dsn::set_default_options(); }
};

struct Fmpq {
    fmpq_t q;
    Fmpq()  { fmpq_init(q); }
    ~Fmpq() { fmpq_clear(q); }
};

// The sample 1/113 exactly as the JSON front end builds it: a ball at
// working precision around the exact rational.
dsn::AcbValue sample_1_over_113() {
    Fmpq e;
    fmpq_set_si(e.q, 1, 113);
    dsn::AcbValue v;
    v.set_fmpq(e.q, dsn::working_prec_bits());
    return v;
}

}  // namespace

TEST_F(NumdFrameFixture, NumdSet_EpsEntryIsExactRational) {
    EnvGuard env("AMFLOW_DIFFEQ_NUMD");
    env.set("450/113");   // d = 450/113  <=>  eps = (4 - 450/113)/2 = 1/113

    const auto eps = sample_1_over_113();
    Fmpq got;
    pipeline::detail::eps_sample_to_fmpq(got.q, eps, dsn::working_prec_bits());

    Fmpq want;
    fmpq_set_si(want.q, 1, 113);
    EXPECT_TRUE(fmpq_equal(got.q, want.q))
        << "expected exactly 1/113, got numerator/denominator "
        << fmpz_get_str(nullptr, 10, fmpq_numref(got.q)) << "/"
        << fmpz_get_str(nullptr, 10, fmpq_denref(got.q));
}

TEST_F(NumdFrameFixture, NumdSet_HelperParsesEnv) {
    EnvGuard env("AMFLOW_DIFFEQ_NUMD");
    env.set("450/113");
    Fmpq got;
    ASSERT_TRUE(pipeline::detail::numd_exact_eps(got.q));
    Fmpq want;
    fmpq_set_si(want.q, 1, 113);
    EXPECT_TRUE(fmpq_equal(got.q, want.q));

    env.unset();
    Fmpq untouched;
    fmpq_set_si(untouched.q, 7, 3);
    EXPECT_FALSE(pipeline::detail::numd_exact_eps(untouched.q));
    Fmpq seven_thirds;
    fmpq_set_si(seven_thirds.q, 7, 3);
    EXPECT_TRUE(fmpq_equal(untouched.q, seven_thirds.q));
}

TEST_F(NumdFrameFixture, NumdUnset_EpsEntryIsSimplestRational) {
    EnvGuard env("AMFLOW_DIFFEQ_NUMD");
    env.unset();

    const auto eps = sample_1_over_113();
    Fmpq got;
    pipeline::detail::eps_sample_to_fmpq(got.q, eps, dsn::working_prec_bits());

    // The sample's midpoint is 1/113 to working precision, so the simplest
    // rational within 10^-60 of it is 1/113 itself: the same frame as the
    // exact numeric-d value.
    Fmpq exact;
    fmpq_set_si(exact.q, 1, 113);
    EXPECT_TRUE(fmpq_equal(got.q, exact.q))
        << "expected exactly 1/113, got numerator/denominator "
        << fmpz_get_str(nullptr, 10, fmpq_numref(got.q)) << "/"
        << fmpz_get_str(nullptr, 10, fmpq_denref(got.q));
}

TEST_F(NumdFrameFixture, NumdUnset_GenericSampleStaysWithinTolerance) {
    // A sample that is not a short rational (a 200-digit-irrational-looking
    // value: the ball around sqrt(2)/100) rationalizes to a value within
    // 10^-60 whose denominator is below 10^60.
    EnvGuard env("AMFLOW_DIFFEQ_NUMD");
    env.unset();

    dsn::AcbValue eps;
    arb_sqrt_ui(acb_realref(eps.raw()), 2, dsn::working_prec_bits());
    arb_div_ui(acb_realref(eps.raw()), acb_realref(eps.raw()), 100,
               dsn::working_prec_bits());
    Fmpq got;
    pipeline::detail::eps_sample_to_fmpq(got.q, eps, dsn::working_prec_bits());

    arb_t diff, bound;
    arb_init(diff);
    arb_init(bound);
    arb_set_fmpq(diff, got.q, 400);
    arb_sub(diff, diff, acb_realref(eps.raw()), 400);
    arb_abs(diff, diff);
    arb_set_si(bound, 10);
    arb_pow_ui(bound, bound, 60, 400);
    arb_inv(bound, bound, 400);
    EXPECT_TRUE(arb_lt(diff, bound));
    arb_clear(diff);
    arb_clear(bound);

    fmpz_t pow10;
    fmpz_init(pow10);
    fmpz_set_ui(pow10, 10);
    fmpz_pow_ui(pow10, pow10, 60);
    EXPECT_LT(fmpz_cmp(fmpq_denref(got.q), pow10), 0);
    fmpz_clear(pow10);
}

TEST_F(NumdFrameFixture, NumdSet_AgreesWithRecoveredSample) {
    // For a sample whose midpoint is the exact value the two paths agree;
    // the env path remains authoritative (it does not depend on the sample
    // being exact to working precision).
    EnvGuard env("AMFLOW_DIFFEQ_NUMD");
    const auto eps = sample_1_over_113();
    Fmpq exact_path, sample_path;

    env.set("450/113");
    pipeline::detail::eps_sample_to_fmpq(exact_path.q, eps, dsn::working_prec_bits());
    env.unset();
    pipeline::detail::eps_sample_to_fmpq(sample_path.q, eps, dsn::working_prec_bits());

    EXPECT_TRUE(fmpq_equal(exact_path.q, sample_path.q));
}
