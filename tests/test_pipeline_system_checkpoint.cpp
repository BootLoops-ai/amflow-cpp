// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Anthropic, PBC. Added by the BootLoops fork of AMFlow.cpp; see NOTICE and PATCHES.md.
// Tests for amflow::pipeline::system_checkpoint — key coverage.
//
// A per-system solve checkpoint is replayed whenever its key matches.
// The key must therefore change with every numeric option and every
// environment override that changes the solve; otherwise a re-run with
// (say) a different rationalize_pre, or with AMFLOW_DIFFEQ_NUMD cleared,
// would silently return the solution of the previous configuration.

#include <gtest/gtest.h>

#include <flint/fmpq.h>

#include <cstdlib>
#include <string>
#include <vector>

#include "amflow/numeric/acb_wrap.hpp"
#include "amflow/numeric/options.hpp"
#include "amflow/pipeline/system_checkpoint.hpp"

namespace pipeline = amflow::pipeline;
namespace dsn      = amflow::numeric;

namespace {

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

class SystemCheckpointKeyFixture : public ::testing::Test {
protected:
    void SetUp() override {
        dsn::set_default_options();
        // Pin the env knobs the key folds in so the baseline is fixed
        // regardless of the ambient environment.
        for (auto* g : guards_) g->unset();
    }
    void TearDown() override { dsn::set_default_options(); }

    EnvGuard numd_{"AMFLOW_DIFFEQ_NUMD"};
    EnvGuard chop_{"AMFLOW_SPARSE_CHOP_DIGITS"};
    EnvGuard nochop_{"AMFLOW_NO_SPARSE_CHOP"};
    EnvGuard reson_{"AMFLOW_RESONANCE_CHOP"};
    EnvGuard guard_{"AMFLOW_X00_PIVOT_GUARD"};
    EnvGuard ending_{"AMFLOW_FORCE_ENDING_DEPTH"};
    std::vector<EnvGuard*> guards_{&numd_, &chop_, &nochop_, &reson_, &guard_,
                                   &ending_};

    std::string key() const {
        std::vector<dsn::AcbValue> eps;
        {
            fmpq_t q;
            fmpq_init(q);
            fmpq_set_si(q, 1, 113);
            dsn::AcbValue v;
            v.set_fmpq(q, dsn::working_prec_bits());
            eps.push_back(std::move(v));
            fmpq_clear(q);
        }
        return pipeline::system_checkpoint_key(
            "fam", {"fam|1|1|0", "fam|1|1|1"}, {1, 0, 0}, eps,
            dsn::working_prec_bits(), /*direction=*/3, /*is_ending=*/false);
    }
};

}  // namespace

TEST_F(SystemCheckpointKeyFixture, StableWhenNothingChanges) {
    const std::string a = key();
    const std::string b = key();
    EXPECT_EQ(a, b);
    EXPECT_EQ(a.size(), 16u);
}

TEST_F(SystemCheckpointKeyFixture, ChangesWithRationalizePre) {
    const std::string base = key();
    auto g = dsn::global_options();
    g.rationalize_pre = 60;
    dsn::set_global_options(g);
    EXPECT_NE(base, key());
}

TEST_F(SystemCheckpointKeyFixture, ChangesWithChopPre) {
    const std::string base = key();
    auto g = dsn::global_options();
    g.chop_pre = 40;
    dsn::set_global_options(g);
    EXPECT_NE(base, key());
}

TEST_F(SystemCheckpointKeyFixture, ChangesWithXOrder) {
    const std::string base = key();
    auto x = dsn::expansion_options();
    x.x_order = 164;
    dsn::set_expansion_options(x);
    const std::string xo = key();
    EXPECT_NE(base, xo);

    x.extra_x_order = 164;
    dsn::set_expansion_options(x);
    EXPECT_NE(xo, key());
}

TEST_F(SystemCheckpointKeyFixture, ChangesWithDiffeqNumdEnv) {
    const std::string base = key();
    numd_.set("450/113");
    const std::string on = key();
    EXPECT_NE(base, on);
    numd_.set("502/127");
    EXPECT_NE(on, key());
    numd_.unset();
    EXPECT_EQ(base, key());
}

TEST_F(SystemCheckpointKeyFixture, ChangesWithSolveEnvOverrides) {
    const std::string base = key();
    chop_.set("55");
    const std::string k1 = key();
    EXPECT_NE(base, k1);
    chop_.unset();

    nochop_.set("1");
    const std::string k2 = key();
    EXPECT_NE(base, k2);
    EXPECT_NE(k1, k2);
    nochop_.unset();

    reson_.set("70");
    const std::string k3 = key();
    EXPECT_NE(base, k3);
    reson_.unset();

    guard_.set("1");
    const std::string k4 = key();
    EXPECT_NE(base, k4);
    guard_.unset();

    EXPECT_EQ(base, key());
}

TEST_F(SystemCheckpointKeyFixture, ChangesWithForceEndingDepthEnv) {
    // The ending scheme of every node below the depth changes, and an
    // ancestor's stored solution subsumes its descendants, so the value
    // must re-key.
    const std::string base = key();
    ending_.set("2");
    const std::string k1 = key();
    EXPECT_NE(base, k1);
    ending_.set("3");
    const std::string k2 = key();
    EXPECT_NE(base, k2);
    EXPECT_NE(k1, k2);
    ending_.unset();
    EXPECT_EQ(base, key());
}
