// SPDX-License-Identifier: MIT
// Tests for amflow::numeric::options + log.
//
// Domain contract: docs/REFACTOR_DESIGN.md §3.1, §5.

#include <gtest/gtest.h>

#include <flint/fmpq.h>

#include "amflow/numeric/log.hpp"
#include "amflow/numeric/options.hpp"

#include <cstdlib>
#include <sstream>
#include <string>

namespace dsn = amflow::numeric;

class NumericOptionsResetFixture : public ::testing::Test {
protected:
    void SetUp()    override { dsn::set_default_options(); }
    void TearDown() override { dsn::set_default_options(); }
};

TEST_F(NumericOptionsResetFixture, DefaultsMatchSpec) {
    EXPECT_EQ(dsn::working_pre(), 100);
    EXPECT_EQ(dsn::chop_pre(), 20);
    EXPECT_EQ(dsn::rationalize_pre(), 100);
    EXPECT_FALSE(dsn::silent_mode());
    EXPECT_EQ(dsn::d0(), "4");

    EXPECT_EQ(dsn::x_order(), 100);
    EXPECT_EQ(dsn::extra_x_order(), 20);
    EXPECT_EQ(dsn::learn_x_order(), -1);
    EXPECT_EQ(dsn::test_x_order(), 5);

    EXPECT_EQ(dsn::run_radius(), 2);
    EXPECT_EQ(dsn::run_length(), 1000);
    EXPECT_EQ(dsn::run_candidate(), 10);
    EXPECT_EQ(dsn::run_direction(), dsn::RunningOptions::Direction::NegIm);
}

TEST_F(NumericOptionsResetFixture, SetGlobalOptionsRoundTrips) {
    auto g = dsn::global_options();
    g.working_pre  = 50;
    g.silent_mode  = true;
    g.d0           = "7/3";
    dsn::set_global_options(g);

    EXPECT_EQ(dsn::working_pre(), 50);
    EXPECT_TRUE(dsn::silent_mode());
    EXPECT_EQ(dsn::d0(), "7/3");
}

TEST_F(NumericOptionsResetFixture, DecimalDigitsToBitsClampsAtFiftyThree) {
    EXPECT_EQ(dsn::decimal_digits_to_bits(0), 53);
    EXPECT_EQ(dsn::decimal_digits_to_bits(-5), 53);
    EXPECT_EQ(dsn::decimal_digits_to_bits(15), 53);
}

TEST_F(NumericOptionsResetFixture, DecimalDigitsToBitsCeils) {
    // ceil(100 * log2(10)) = ceil(332.19...) = 333
    EXPECT_EQ(dsn::decimal_digits_to_bits(100), 333);
    EXPECT_EQ(dsn::decimal_digits_to_bits(30), 100);   // ceil(99.66) = 100
}

TEST_F(NumericOptionsResetFixture, WorkingPrecBitsTracksGlobal) {
    auto g = dsn::global_options();
    g.working_pre = 50;
    dsn::set_global_options(g);
    EXPECT_EQ(dsn::working_prec_bits(), dsn::decimal_digits_to_bits(50));
}

TEST_F(NumericOptionsResetFixture, D0ShiftDefaultIsZero) {
    fmpq_t shift;
    fmpq_init(shift);
    dsn::d0_eps_shift_fmpq(shift);
    EXPECT_NE(fmpq_is_zero(shift), 0);
    fmpq_clear(shift);
}

TEST_F(NumericOptionsResetFixture, D0ShiftRationalIsExact) {
    auto g = dsn::global_options();
    g.d0 = "7/3";
    dsn::set_global_options(g);

    // (4 - 7/3) / 2 = (12/3 - 7/3) / 2 = (5/3) / 2 = 5/6
    fmpq_t shift, expected;
    fmpq_init(shift);
    fmpq_init(expected);
    dsn::d0_eps_shift_fmpq(shift);
    fmpq_set_si(expected, 5, 6);
    EXPECT_NE(fmpq_equal(shift, expected), 0);
    fmpq_clear(shift);
    fmpq_clear(expected);
}

TEST_F(NumericOptionsResetFixture, D0ShiftThrowsOnGarbage) {
    auto g = dsn::global_options();
    g.d0 = "not-a-rational";
    dsn::set_global_options(g);

    fmpq_t shift;
    fmpq_init(shift);
    EXPECT_THROW(dsn::d0_eps_shift_fmpq(shift), std::runtime_error);
    fmpq_clear(shift);
}

TEST_F(NumericOptionsResetFixture, GlobalScopeRestoresOnDestruction) {
    {
        dsn::GlobalScope guard;
        guard.global.working_pre  = 50;
        guard.expansion.x_order   = 30;
        guard.commit();
        EXPECT_EQ(dsn::working_pre(), 50);
        EXPECT_EQ(dsn::x_order(), 30);
    }
    // After dtor: defaults restored.
    EXPECT_EQ(dsn::working_pre(), 100);
    EXPECT_EQ(dsn::x_order(), 100);
}

TEST_F(NumericOptionsResetFixture, GlobalScopeWithoutCommitDoesNothing) {
    {
        dsn::GlobalScope guard;
        guard.global.working_pre = 50;
        // No commit() — value should not propagate.
    }
    EXPECT_EQ(dsn::working_pre(), 100);
}

TEST_F(NumericOptionsResetFixture, GlobalScopeCommitIsIdempotent) {
    {
        dsn::GlobalScope guard;
        guard.global.working_pre = 50;
        guard.commit();
        // Mutate after commit — should not propagate (commit is one-shot).
        guard.global.working_pre = 75;
        guard.commit();
        EXPECT_EQ(dsn::working_pre(), 50);
    }
    EXPECT_EQ(dsn::working_pre(), 100);
}

TEST_F(NumericOptionsResetFixture, LogLineRespectsSilentMode) {
    std::stringstream capture;
    auto* old_buf = dsn::log_stream().rdbuf(capture.rdbuf());

    dsn::log_line("noisy");
    EXPECT_NE(capture.str().find("noisy"), std::string::npos);

    capture.str("");
    auto g = dsn::global_options();
    g.silent_mode = true;
    dsn::set_global_options(g);
    dsn::log_line("silenced");
    EXPECT_TRUE(capture.str().empty());

    dsn::log_stream().rdbuf(old_buf);
}

// ---------------------------------------------------------------------------
//  Trace gate
// ---------------------------------------------------------------------------

TEST(NumericLog, TraceDisabledByDefault) {
    ::unsetenv("AMFLOW_TEST_TRACE_X");
    EXPECT_FALSE(::amflow::numeric::log::trace_enabled(
        "AMFLOW_TEST_TRACE_X"));
}

TEST(NumericLog, TraceEnabledWhenEnvNonEmpty) {
    ::setenv("AMFLOW_TEST_TRACE_X", "1", /*overwrite=*/1);
    EXPECT_TRUE(::amflow::numeric::log::trace_enabled(
        "AMFLOW_TEST_TRACE_X"));
    ::unsetenv("AMFLOW_TEST_TRACE_X");
}

TEST(NumericLog, TraceDisabledWhenEnvEmpty) {
    ::setenv("AMFLOW_TEST_TRACE_X", "", /*overwrite=*/1);
    EXPECT_FALSE(::amflow::numeric::log::trace_enabled(
        "AMFLOW_TEST_TRACE_X"));
    ::unsetenv("AMFLOW_TEST_TRACE_X");
}

TEST(NumericLog, TraceMacroFiresGuardedBlock) {
    ::setenv("AMFLOW_TEST_TRACE_Y", "1", /*overwrite=*/1);
    int hits = 0;
    AMFLOW_TRACE("AMFLOW_TEST_TRACE_Y") {
        ++hits;
    }
    EXPECT_EQ(hits, 1);
    ::unsetenv("AMFLOW_TEST_TRACE_Y");

    AMFLOW_TRACE("AMFLOW_TEST_TRACE_Y") {
        ++hits;
    }
    EXPECT_EQ(hits, 1);   // still 1
}
