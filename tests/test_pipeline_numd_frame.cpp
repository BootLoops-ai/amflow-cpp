// SPDX-License-Identifier: MIT
// Tests for the per-eps rational frame of amflow::pipeline (the "eps"
// entry of the numeric substitution map built for each eps sample).
//
// With AMFLOW_DIFFEQ_NUMD=P/Q the differential-equation tables carry
// d = P/Q exactly, so the region exponents and boundary coefficients
// must be evaluated at eps = (4 - P/Q)/2 exactly as well; a sample
// rationalized to `rationalize_pre` decimal digits (a truncation) would
// sit in a different frame and the mismatch is pivoted on by the sparse
// solver once it exceeds the chop threshold.  Without the env var the
// entry is the truncation, unchanged.

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

TEST_F(NumdFrameFixture, NumdUnset_EpsEntryIsDecimalTruncation) {
    EnvGuard env("AMFLOW_DIFFEQ_NUMD");
    env.unset();

    const auto eps = sample_1_over_113();
    Fmpq got;
    pipeline::detail::eps_sample_to_fmpq(got.q, eps, dsn::working_prec_bits());

    // Existing semantics: numer / 10^rationalize_pre, i.e. NOT 1/113 ...
    Fmpq exact;
    fmpq_set_si(exact.q, 1, 113);
    EXPECT_FALSE(fmpq_equal(got.q, exact.q));

    // ... the denominator divides 10^60 ...
    Fmpq pow10;
    {
        fmpz_t p;
        fmpz_init(p);
        fmpz_set_ui(p, 10);
        fmpz_pow_ui(p, p, 60);
        fmpq_set_fmpz(pow10.q, p);
        fmpz_clear(p);
    }
    Fmpq scaled;
    fmpq_mul(scaled.q, got.q, pow10.q);
    EXPECT_TRUE(fmpz_is_one(fmpq_denref(scaled.q)))
        << "10^60 * eps_q should be an integer";

    // ... and |eps_q - 1/113| < 10^-59 (in fact <= 0.5 * 10^-60).
    Fmpq diff;
    fmpq_sub(diff.q, got.q, exact.q);
    fmpq_abs(diff.q, diff.q);
    Fmpq bound;
    {
        fmpz_t p;
        fmpz_init(p);
        fmpz_set_ui(p, 10);
        fmpz_pow_ui(p, p, 59);
        fmpq_set_si(bound.q, 1, 1);
        fmpq_div_fmpz(bound.q, bound.q, p);
        fmpz_clear(p);
    }
    EXPECT_LT(fmpq_cmp(diff.q, bound.q), 0);
    EXPECT_GT(fmpq_sgn(diff.q), 0);
}

TEST_F(NumdFrameFixture, NumdSet_WrongFrameWouldDifferFromTruncation) {
    // The two paths disagree for a denominator that does not divide a
    // power of ten — this is the frame difference the fix removes.
    EnvGuard env("AMFLOW_DIFFEQ_NUMD");
    const auto eps = sample_1_over_113();
    Fmpq exact_path, trunc_path;

    env.set("450/113");
    pipeline::detail::eps_sample_to_fmpq(exact_path.q, eps, dsn::working_prec_bits());
    env.unset();
    pipeline::detail::eps_sample_to_fmpq(trunc_path.q, eps, dsn::working_prec_bits());

    EXPECT_FALSE(fmpq_equal(exact_path.q, trunc_path.q));
}
