// SPDX-License-Identifier: MIT
// Tests for amflow::ibp::kira (yaml writers + parser + readers).
//
// kira_run is exercised separately in a (currently absent) integration test.

#include <cstdio>
#include <filesystem>
#include <fstream>
#include <gtest/gtest.h>
#include <sstream>

#include "amflow/algebra/mpoly.hpp"
#include "amflow/ibp/kira.hpp"
#include "amflow/qft/family_config.hpp"
#include "amflow/qft/jintegral.hpp"

namespace alg = amflow::algebra;
namespace ibp = amflow::ibp;
namespace qft = amflow::qft;
namespace fs  = std::filesystem;

namespace {

std::string read_file(const std::string& path) {
    std::ifstream f(path);
    if (!f) return "";
    std::ostringstream oss;
    oss << f.rdbuf();
    return oss.str();
}

}  // namespace

TEST(KiraTest, ParseExpression_IntegerLiteral) {
    auto ctx = std::make_shared<alg::MpolyContext>(
        std::vector<std::string>{"d", "s"});
    EXPECT_EQ(ibp::kira_parse_expression(ctx, "5").to_string(), "5");
    EXPECT_EQ(ibp::kira_parse_expression(ctx, "-7").to_string(), "-7");
    EXPECT_EQ(ibp::kira_parse_expression(ctx, "(3)").to_string(), "3");
}

TEST(KiraTest, ParseExpression_SingleVariable) {
    auto ctx = std::make_shared<alg::MpolyContext>(
        std::vector<std::string>{"d", "s"});
    EXPECT_EQ(ibp::kira_parse_expression(ctx, "d").to_string(), "d");
    EXPECT_EQ(ibp::kira_parse_expression(ctx, "-d").to_string(), "-d");
}

TEST(KiraTest, ParseExpression_Arithmetic) {
    auto ctx = std::make_shared<alg::MpolyContext>(
        std::vector<std::string>{"d", "s"});
    EXPECT_EQ(ibp::kira_parse_expression(ctx, "d + 2").to_string(), "d+2");
    EXPECT_EQ(ibp::kira_parse_expression(ctx, "d - 4").to_string(), "d-4");
    EXPECT_EQ(ibp::kira_parse_expression(ctx, "2 * d").to_string(), "2*d");
}

TEST(KiraTest, ParseExpression_Power) {
    auto ctx = std::make_shared<alg::MpolyContext>(
        std::vector<std::string>{"d", "s"});
    EXPECT_EQ(ibp::kira_parse_expression(ctx, "d^2").to_string(), "d^2");
    EXPECT_EQ(ibp::kira_parse_expression(ctx, "(d+1)^3").to_string(),
              "d^3+3*d^2+3*d+1");
}

TEST(KiraTest, ParseExpression_NestedParens) {
    auto ctx = std::make_shared<alg::MpolyContext>(
        std::vector<std::string>{"d", "s"});
    auto r = ibp::kira_parse_expression(ctx, "(d-4)*(d-2)/(s+1)");
    alg::Mpoly num = r.numerator();
    alg::Mpoly den = r.denominator();
    EXPECT_FALSE(num.is_zero());
    EXPECT_FALSE(den.is_one());
}

TEST(KiraTest, ParseExpression_DivisionAndUnaryMinus) {
    auto ctx = std::make_shared<alg::MpolyContext>(
        std::vector<std::string>{"d"});
    auto r = ibp::kira_parse_expression(ctx, "-3/(d-4)");
    alg::Mfrac expected = alg::Mfrac::from_si(ctx, -3) /
                          (alg::Mfrac::from_mpoly(alg::Mpoly::variable(ctx, 0))
                             - alg::Mfrac::from_si(ctx, 4));
    EXPECT_EQ(r.to_string(), expected.to_string());
}

TEST(KiraTest, ParseExpression_UnknownIdentifierThrows) {
    auto ctx = std::make_shared<alg::MpolyContext>(
        std::vector<std::string>{"d"});
    EXPECT_THROW(ibp::kira_parse_expression(ctx, "x"), std::runtime_error);
}

TEST(KiraTest, ParseExpression_TrailingInputThrows) {
    auto ctx = std::make_shared<alg::MpolyContext>(
        std::vector<std::string>{"d"});
    EXPECT_THROW(ibp::kira_parse_expression(ctx, "d ) extra"),
                 std::runtime_error);
}

// --- Coefficient parser grammar lock ---------
//
// `kira_parse_expression` (`src/ibp/kira_parse.cpp`) supports a
// strict subset of Mathematica algebraic syntax — Kira's coefficient
// output uses only integer literals, identifiers (registered in
// the MpolyContext), unary `-`, binary `+ - * /`, and `^` followed
// by an integer.  No decimals, no fractional/symbolic exponents.
// This parser is not
// fuzzed; production correctness comes from the fact that Kira's
// normal output never trips the unsupported forms.
//
// These negative-acceptance tests lock the contract: any input
// that uses an unsupported lexical form throws (no silent partial
// parse, no swallow of the remainder).  Positive-acceptance is
// covered by the seven `ParseExpression_*` tests above.

TEST(KiraTest, ParseExpression_DecimalLiteralThrows) {
    auto ctx = std::make_shared<alg::MpolyContext>(
        std::vector<std::string>{"d"});
    // "1.5" — decimal point in a numeric literal.
    EXPECT_THROW(ibp::kira_parse_expression(ctx, "1.5"), std::runtime_error);
}

TEST(KiraTest, ParseExpression_NonIntegerExponentThrows) {
    auto ctx = std::make_shared<alg::MpolyContext>(
        std::vector<std::string>{"d"});
    // "d^d" — RHS of `^` must be an integer literal, not an identifier.
    EXPECT_THROW(ibp::kira_parse_expression(ctx, "d^d"), std::runtime_error);
    // "d^(2+1)" — RHS of `^` must be a *literal* integer, not even
    // a parenthesised expression that evaluates to one.
    EXPECT_THROW(ibp::kira_parse_expression(ctx, "d^(2+1)"),
                 std::runtime_error);
}

TEST(KiraTest, ParseExpression_GarbageCharacterThrows) {
    auto ctx = std::make_shared<alg::MpolyContext>(
        std::vector<std::string>{"d"});
    // "#" — outside the supported lexical alphabet.
    EXPECT_THROW(ibp::kira_parse_expression(ctx, "d # 2"),
                 std::runtime_error);
    // "@" — same.
    EXPECT_THROW(ibp::kira_parse_expression(ctx, "@d"),
                 std::runtime_error);
}

TEST(KiraTest, ParseExpression_UnclosedParenThrows) {
    auto ctx = std::make_shared<alg::MpolyContext>(
        std::vector<std::string>{"d"});
    EXPECT_THROW(ibp::kira_parse_expression(ctx, "(d + 1"),
                 std::runtime_error);
}

TEST(KiraTest, WriteConfig_OneLoopBubble_HasExpectedYaml) {
    auto fc = qft::FamilyConfig::build(
        "bubble", {"l"}, {"p"}, {}, {{"p^2", "s"}},
        {"l^2 - msq", "(l - p)^2 - msq"});

    ibp::KiraConfig cfg;
    cfg.fc = &fc;
    cfg.top_pattern = {1, 1};
    cfg.numeric_values["s"] = "100";
    cfg.numeric_values["msq"] = "1";
    cfg.kira_executable = "/usr/local/bin/kira";
    cfg.fermat_executable = "/usr/share/Ferl7/fer64";

    std::string dir =
        (fs::temp_directory_path() / "amflow_v2_kira_t1").string();
    fs::remove_all(dir);
    fs::create_directories(dir);

    ibp::kira_write_config(cfg, dir);

    auto fam_yaml = read_file(dir + "/config/integralfamilies.yaml");
    auto kin_yaml = read_file(dir + "/config/kinematics.yaml");

    EXPECT_NE(fam_yaml.find("name: \"bubble\""), std::string::npos)
        << fam_yaml;
    EXPECT_NE(fam_yaml.find("loop_momenta: [l]"), std::string::npos);
    EXPECT_NE(fam_yaml.find("top_level_sectors: [3]"), std::string::npos)
        << fam_yaml;
    EXPECT_NE(kin_yaml.find("incoming_momenta: [p, bubbleAuxLeg]"),
              std::string::npos) << kin_yaml;
    EXPECT_NE(kin_yaml.find("scalarproduct_rules:"), std::string::npos);
    EXPECT_NE(kin_yaml.find("[p,p]"), std::string::npos) << kin_yaml;

    // Pin propagator yaml format to MMA-aligned `[Propagator, 0]` form
    // (matches MMA Kira/interface.m:126).  Both propagators must appear
    // as a full denominator expression with mass-field=0, NOT the older
    // `[Momentum, -Mass]` form (commented out at MMA :125).  The
    // closed-form `(l - p)^2` (not the expanded `l^2-2*l*p+p^2`) is
    // produced by try_render_closed_form's rank-1 factoring.
    EXPECT_NE(fam_yaml.find("\"l^2 - 1\", 0"), std::string::npos)
        << "expected propagator l^2 - msq (numeric-subbed msq=1) emitted "
           "as `[\"l^2 - 1\", 0]`; got:\n" << fam_yaml;
    EXPECT_NE(fam_yaml.find("\"(l - p)^2 - 1\", 0"), std::string::npos)
        << "expected propagator (l - p)^2 - msq emitted as closed form "
           "`[\"(l - p)^2 - 1\", 0]` (not expanded `l^2-2*l*p+p^2-1`); "
           "got:\n" << fam_yaml;
    EXPECT_EQ(fam_yaml.find("\", l^2 - 1"), std::string::npos)
        << "saw signs of the old `[Momentum, -Mass]` writer form; the "
           "writer should emit `[FullDenominator, 0]` per MMA "
           "interface.m:126.  yaml:\n" << fam_yaml;

    fs::remove_all(dir);
}

TEST(KiraTest, WriteConfig_NoLegs_NoConservation) {
    auto fc = qft::FamilyConfig::build(
        "tad", {"l"}, {}, {}, {}, {"l^2 - msq"});

    ibp::KiraConfig cfg;
    cfg.fc = &fc;
    cfg.top_pattern = {1};
    cfg.kira_executable = "/usr/local/bin/kira";

    std::string dir =
        (fs::temp_directory_path() / "amflow_v2_kira_t2").string();
    fs::remove_all(dir);
    fs::create_directories(dir);

    ibp::kira_write_config(cfg, dir);

    auto kin_yaml = read_file(dir + "/config/kinematics.yaml");
    EXPECT_NE(kin_yaml.find("incoming_momenta: []"), std::string::npos);
    EXPECT_NE(kin_yaml.find("momentum_conservation: []"), std::string::npos);

    fs::remove_all(dir);
}

TEST(KiraTest, WriteJobs_MastersMode) {
    auto fc = qft::FamilyConfig::build(
        "bubble", {"l"}, {"p"}, {}, {{"p^2", "s"}},
        {"l^2 - msq", "(l - p)^2 - msq"});

    ibp::KiraConfig cfg;
    cfg.fc = &fc;
    cfg.top_pattern = {1, 1};
    cfg.ibp_rank = 5;
    cfg.ibp_dot = 0;

    std::string dir =
        (fs::temp_directory_path() / "amflow_v2_kira_t3").string();
    fs::remove_all(dir);
    fs::create_directories(dir);

    ibp::kira_write_jobs(cfg, dir, ibp::KiraReductionMode::Masters);

    auto jobs = read_file(dir + "/jobs.yaml");
    EXPECT_NE(jobs.find("run_initiate: masters"), std::string::npos) << jobs;
    EXPECT_EQ(jobs.find("kira2math"), std::string::npos)
        << "Masters mode should not produce kira2math";

    fs::remove_all(dir);
}

// --- r = nonzero(top_pattern) + IBPDot arithmetic -
//
// Lock the formula `r = Length[TopSector] - Count[TopSector, 0] +
// IBPDot` (Kira/interface.m yaml emit) against the C++ implementation
// at `src/ibp/kira_yaml.cpp` (post-fix: `count_if(... != 0)` so the
// formula matches upstream literally).  These tests cover:
//
//   1. zero IBPDot (production default) — most common case;
//   2. non-trivial IBPDot — exercises the `+ ibp_dot` term;
//   3. mixed-presence top_pattern with several zeros and ones —
//      exercises the count-of-nonzero semantic;
//   4. defensive: a non-binary top_pattern entry (e.g. `2`) is
//      counted as "present" matching upstream's nonzero semantics
//      (production `get_top_sector` only emits 0/1 but the writer
//      should remain parity-correct under any future relaxation).
//
// Each test reads the emitted jobs.yaml and asserts the `r:` field
// has the expected integer value.

namespace {

// Extract the first `r: <int>` value from a jobs.yaml string.
long extract_r_from_jobs(const std::string& jobs) {
    auto pos = jobs.find("r: ");
    if (pos == std::string::npos) return -1;
    auto start = pos + 3;
    auto end = jobs.find_first_of(",}\n", start);
    if (end == std::string::npos) return -1;
    return std::stol(jobs.substr(start, end - start));
}

}  // namespace

TEST(KiraTest, WriteJobs_R_ZeroIbpDot_AllOnes) {
    auto fc = qft::FamilyConfig::build(
        "bubble", {"l"}, {"p"}, {}, {{"p^2", "s"}},
        {"l^2 - msq", "(l - p)^2 - msq"});

    ibp::KiraConfig cfg;
    cfg.fc = &fc;
    cfg.top_pattern = {1, 1};
    cfg.ibp_dot = 0;

    std::string dir =
        (fs::temp_directory_path() / "amflow_kira_r_t1").string();
    fs::remove_all(dir);
    fs::create_directories(dir);
    ibp::kira_write_jobs(cfg, dir, ibp::KiraReductionMode::Masters);

    auto jobs = read_file(dir + "/jobs.yaml");
    EXPECT_EQ(extract_r_from_jobs(jobs), 2)
        << "r should equal nonzero(top_pattern) + ibp_dot = 2 + 0";

    fs::remove_all(dir);
}

TEST(KiraTest, WriteJobs_R_NontrivialIbpDot) {
    // Larger family so we can vary the dot meaningfully.
    auto fc = qft::FamilyConfig::build(
        "box", {"l"}, {"p1", "p2", "p3", "p4"},
        {{"p4", "-p1 - p2 - p3"}},
        {{"p1^2", "0"}, {"p2^2", "0"}, {"p3^2", "0"},
         {"(p1 + p2)^2", "s"}, {"(p2 + p3)^2", "t"}},
        {"l^2", "(l + p1)^2", "(l + p1 + p2)^2",
         "(l + p1 + p2 + p3)^2"});

    ibp::KiraConfig cfg;
    cfg.fc = &fc;
    cfg.top_pattern = {1, 1, 1, 1};
    cfg.ibp_dot = 3;

    std::string dir =
        (fs::temp_directory_path() / "amflow_kira_r_t2").string();
    fs::remove_all(dir);
    fs::create_directories(dir);
    ibp::kira_write_jobs(cfg, dir, ibp::KiraReductionMode::Masters);

    auto jobs = read_file(dir + "/jobs.yaml");
    EXPECT_EQ(extract_r_from_jobs(jobs), 7)
        << "r should equal nonzero(top_pattern) + ibp_dot = 4 + 3";

    fs::remove_all(dir);
}

TEST(KiraTest, WriteJobs_R_MixedZerosAndOnes) {
    // Sub-sector top_pattern: only a strict subset of propagators
    // active.  Confirms the "count nonzero" semantics matches upstream
    // for a non-trivial sparse pattern.
    auto fc = qft::FamilyConfig::build(
        "box", {"l"}, {"p1", "p2", "p3", "p4"},
        {{"p4", "-p1 - p2 - p3"}},
        {{"p1^2", "0"}, {"p2^2", "0"}, {"p3^2", "0"},
         {"(p1 + p2)^2", "s"}, {"(p2 + p3)^2", "t"}},
        {"l^2", "(l + p1)^2", "(l + p1 + p2)^2",
         "(l + p1 + p2 + p3)^2"});

    ibp::KiraConfig cfg;
    cfg.fc = &fc;
    cfg.top_pattern = {1, 0, 1, 0};   // 2 nonzeros
    cfg.ibp_dot = 2;

    std::string dir =
        (fs::temp_directory_path() / "amflow_kira_r_t3").string();
    fs::remove_all(dir);
    fs::create_directories(dir);
    ibp::kira_write_jobs(cfg, dir, ibp::KiraReductionMode::Masters);

    auto jobs = read_file(dir + "/jobs.yaml");
    EXPECT_EQ(extract_r_from_jobs(jobs), 4)
        << "r should equal nonzero(top_pattern) + ibp_dot = 2 + 2";

    fs::remove_all(dir);
}

TEST(KiraTest, WriteJobs_R_NonBinaryEntryCountedAsNonzero) {
    // Defensive: production `get_top_sector` always emits 0/1, but
    // the yaml writer should remain parity-correct (matching
    // upstream's `Length - Count(0)` formula) under any relaxation
    // of that guarantee.  Lock the count-nonzero semantic by
    // injecting a `2` into top_pattern and verifying it is treated
    // the same as a `1`.
    auto fc = qft::FamilyConfig::build(
        "bubble", {"l"}, {"p"}, {}, {{"p^2", "s"}},
        {"l^2 - msq", "(l - p)^2 - msq"});

    ibp::KiraConfig cfg;
    cfg.fc = &fc;
    cfg.top_pattern = {2, 1};   // 2 nonzeros (one of which is non-1)
    cfg.ibp_dot = 1;

    std::string dir =
        (fs::temp_directory_path() / "amflow_kira_r_t4").string();
    fs::remove_all(dir);
    fs::create_directories(dir);
    ibp::kira_write_jobs(cfg, dir, ibp::KiraReductionMode::Masters);

    auto jobs = read_file(dir + "/jobs.yaml");
    EXPECT_EQ(extract_r_from_jobs(jobs), 3)
        << "r should equal nonzero(top_pattern) + ibp_dot = 2 + 1";

    fs::remove_all(dir);
}

TEST(KiraTest, WriteJobs_ReduceMode) {
    auto fc = qft::FamilyConfig::build(
        "bubble", {"l"}, {"p"}, {}, {{"p^2", "s"}},
        {"l^2 - msq", "(l - p)^2 - msq"});

    ibp::KiraConfig cfg;
    cfg.fc = &fc;
    cfg.top_pattern = {1, 1};

    std::string dir =
        (fs::temp_directory_path() / "amflow_v2_kira_t4").string();
    fs::remove_all(dir);
    fs::create_directories(dir);

    ibp::kira_write_jobs(cfg, dir, ibp::KiraReductionMode::Reduce);

    auto jobs = read_file(dir + "/jobs.yaml");
    EXPECT_NE(jobs.find("run_initiate: true"), std::string::npos);
    EXPECT_NE(jobs.find("run_triangular: true"), std::string::npos);
    EXPECT_NE(jobs.find("run_back_substitution: true"), std::string::npos);
    EXPECT_NE(jobs.find("kira2math"), std::string::npos);

    // No split requested: only jobs.yaml exists (byte-identical to the
    // pre-split combined reduce+kira2math job; no jobs_kira2math.yaml).
    EXPECT_FALSE(fs::exists(dir + "/jobs_kira2math.yaml"))
        << "Default (split off) must not emit jobs_kira2math.yaml";

    fs::remove_all(dir);
}

// Kira 3.1 FireFly export split (the export fault fix): with
// `split_kira2math` set, the Reduce-mode jobs.yaml must carry the reduce
// job WITHOUT the in-process kira2math job, and a separate
// jobs_kira2math.yaml must carry ONLY the kira2math export.  This is what
// lets the caller run the export in a fresh kira process (sharing kira.db)
// and dodge the in-process kira2math-after-firefly SIG6/bad_alloc.
TEST(KiraTest, WriteJobs_ReduceMode_FireflySplit) {
    auto fc = qft::FamilyConfig::build(
        "bubble", {"l"}, {"p"}, {}, {{"p^2", "s"}},
        {"l^2 - msq", "(l - p)^2 - msq"});

    ibp::KiraConfig cfg;
    cfg.fc = &fc;
    cfg.top_pattern = {1, 1};
    cfg.run_firefly = true;
    cfg.split_kira2math = true;

    std::string dir =
        (fs::temp_directory_path() / "amflow_kira_ff_split").string();
    fs::remove_all(dir);
    fs::create_directories(dir);

    ibp::kira_write_jobs(cfg, dir, ibp::KiraReductionMode::Reduce);

    auto jobs = read_file(dir + "/jobs.yaml");
    // Reduce job present (FireFly path), but kira2math is NOT in jobs.yaml.
    EXPECT_NE(jobs.find("run_initiate: true"), std::string::npos) << jobs;
    EXPECT_NE(jobs.find("run_firefly: true"), std::string::npos) << jobs;
    EXPECT_EQ(jobs.find("kira2math"), std::string::npos)
        << "Split: jobs.yaml must NOT contain the in-process kira2math job";

    // The standalone export file exists and contains exactly the kira2math
    // job (and no reduce job).
    ASSERT_TRUE(fs::exists(dir + "/jobs_kira2math.yaml"))
        << "Split: jobs_kira2math.yaml must be emitted";
    auto k2m = read_file(dir + "/jobs_kira2math.yaml");
    EXPECT_NE(k2m.find("kira2math"), std::string::npos) << k2m;
    EXPECT_EQ(k2m.find("run_firefly"), std::string::npos)
        << "Export file must not re-run the reduce";
    EXPECT_EQ(k2m.find("reduce_sectors"), std::string::npos)
        << "Export file must not re-run the reduce";

    fs::remove_all(dir);
}

// Regression: kira_run must SURFACE the real child failure text, not just a
// bare signal number.  A self-capped std::bad_alloc (OOM under RLIMIT_AS /
// AMFLOW_KIRA_MEM_CAP_GB) reaches the parent as SIG6 (std::terminate->abort),
// identical at the wait4 level to a Kira internal assert.  We install a fake
// "kira" that prints the libstdc++ terminate-after-bad_alloc banner to stderr
// and then abort()s, and require the thrown error to (a) quote that banner
// verbatim and (b) correctly flag the SIG6 as a possible self-capped OOM
// rather than only a Kira assert.  Before the fix, the std::bad_alloc text was
// dropped and SIG6 was mislabeled "Kira abort/assert", sending diagnosis the
// wrong way and costing a standalone re-run to recover.
TEST(KiraTest, RunSurfacesChildBadAllocOnSig6) {
    std::string dir =
        (fs::temp_directory_path() / "amflow_kira_badalloc_surface").string();
    fs::remove_all(dir);
    fs::create_directories(dir);

    // Fake kira: mimic the libstdc++ uncaught-exception path that a
    // std::bad_alloc takes under a process-self-capped allocator, then abort
    // (SIGABRT == SIG6), exactly as the real terminate handler would.
    std::string fake_kira = dir + "/fake_kira.sh";
    {
        std::ofstream f(fake_kira);
        f << "#!/bin/sh\n"
          << "echo 'reading jobs.yaml' >&2\n"
          << "echo 'terminate called after throwing an instance of "
             "'\\''std::bad_alloc'\\''' >&2\n"
          << "echo '  what():  std::bad_alloc' >&2\n"
          << "kill -ABRT $$\n";
    }
    fs::permissions(fake_kira,
                    fs::perms::owner_all | fs::perms::group_read |
                        fs::perms::group_exec,
                    fs::perm_options::add);

    ibp::KiraConfig cfg;
    cfg.kira_executable = fake_kira;
    cfg.n_thread = 1;

    std::string log = dir + "/kira.log";
    std::string thrown;
    try {
        ibp::kira_run(cfg, dir, log, "jobs.yaml");
        FAIL() << "kira_run should have thrown on the aborting fake kira";
    } catch (const std::exception& e) {
        thrown = e.what();
    }

    // (1) The real C++ exception text is surfaced verbatim.
    EXPECT_NE(thrown.find("std::bad_alloc"), std::string::npos)
        << "thrown error must quote the child's std::bad_alloc text. Got:\n"
        << thrown;
    EXPECT_NE(thrown.find("terminate called"), std::string::npos)
        << "thrown error must include the child terminate banner. Got:\n"
        << thrown;
    // (2) SIG6 hint is corrected: it must mention the self-capped OOM /
    // RLIMIT_AS possibility, not assert-only.
    EXPECT_NE(thrown.find("SIGNAL 6"), std::string::npos) << thrown;
    EXPECT_NE(thrown.find("RLIMIT_AS"), std::string::npos)
        << "SIG6 hint must mention the self-capped OOM (RLIMIT_AS) path. Got:\n"
        << thrown;

    fs::remove_all(dir);
}

TEST(KiraTest, WritePreferredAndTargets) {
    auto fc = qft::FamilyConfig::build(
        "bubble", {"l"}, {"p"}, {}, {{"p^2", "s"}},
        {"l^2 - msq", "(l - p)^2 - msq"});

    std::string dir =
        (fs::temp_directory_path() / "amflow_v2_kira_t5").string();
    fs::remove_all(dir);
    fs::create_directories(dir);

    std::vector<qft::JIntegral> prefs;
    prefs.push_back(qft::JIntegral("bubble", {1, 1}));
    prefs.push_back(qft::JIntegral("bubble", {2, 1}));
    ibp::kira_write_preferred(prefs, fc, dir);

    auto pref_str = read_file(dir + "/preferred");
    EXPECT_NE(pref_str.find("bubble[1, 1]"), std::string::npos);
    EXPECT_NE(pref_str.find("bubble[2, 1]"), std::string::npos);

    std::vector<qft::JIntegral> targets;
    targets.push_back(qft::JIntegral("bubble", {1, 1}));
    ibp::kira_write_targets(targets, fc, dir);

    auto tgt_str = read_file(dir + "/target");
    EXPECT_NE(tgt_str.find("bubble[1, 1]"), std::string::npos);

    fs::remove_all(dir);
}

TEST(KiraTest, ReadMasters_Synthetic) {
    auto fc = qft::FamilyConfig::build(
        "bubble", {"l"}, {"p"}, {}, {{"p^2", "s"}},
        {"l^2 - msq", "(l - p)^2 - msq"});

    std::string dir =
        (fs::temp_directory_path() / "amflow_v2_kira_t6").string();
    fs::remove_all(dir);
    fs::create_directories(dir + "/results/bubble");

    {
        std::ofstream f(dir + "/results/bubble/masters");
        f << "bubble(1,1)  some_extra_tag\n"
          << "\n"
          << "bubble(0,1)  another_tag\n";
    }

    auto masters = ibp::kira_read_masters(fc, dir);
    ASSERT_EQ(masters.size(), 2u);
    EXPECT_EQ(masters[0].family(), "bubble");
    EXPECT_EQ(masters[0].indices(), (std::vector<long>{1, 1}));
    EXPECT_EQ(masters[1].indices(), (std::vector<long>{0, 1}));

    fs::remove_all(dir);
}

TEST(KiraTest, ReadMasters_MissingFileReturnsEmpty) {
    auto fc = qft::FamilyConfig::build(
        "bubble", {"l"}, {"p"}, {}, {{"p^2", "s"}},
        {"l^2 - msq", "(l - p)^2 - msq"});
    std::string dir =
        (fs::temp_directory_path() / "amflow_v2_kira_t7").string();
    fs::remove_all(dir);
    auto masters = ibp::kira_read_masters(fc, dir);
    EXPECT_TRUE(masters.empty());
}

TEST(KiraTest, ReadTargetTable_Synthetic) {
    auto fc = qft::FamilyConfig::build(
        "bubble", {"l"}, {"p"}, {}, {{"p^2", "s"}},
        {"l^2 - msq", "(l - p)^2 - msq"});

    std::string dir =
        (fs::temp_directory_path() / "amflow_v2_kira_t8").string();
    fs::remove_all(dir);
    fs::create_directories(dir + "/results/bubble");

    {
        std::ofstream f(dir + "/results/bubble/kira_target.m");
        f << "{\n"
          << "  bubble(2,1) -> (d-3)/(2*s) * bubble(1,1) "
          << "- 1/(2*s) * bubble(0,1)\n"
          << "}\n";
    }

    auto rules = ibp::kira_read_target_table(fc, dir);
    ASSERT_EQ(rules.size(), 1u);
    EXPECT_EQ(rules[0].lhs.indices(), (std::vector<long>{2, 1}));
    ASSERT_EQ(rules[0].rhs.size(), 2u);
    EXPECT_EQ(rules[0].rhs[0].second.indices(),
              (std::vector<long>{1, 1}));
    EXPECT_EQ(rules[0].rhs[1].second.indices(),
              (std::vector<long>{0, 1}));

    auto cctx = std::make_shared<alg::MpolyContext>(
        std::vector<std::string>{"d", "s"});
    auto c0 = ibp::kira_parse_expression(cctx, rules[0].rhs[0].first);
    auto c1 = ibp::kira_parse_expression(cctx, rules[0].rhs[1].first);
    EXPECT_FALSE(c0.is_zero());
    EXPECT_FALSE(c1.is_zero());

    fs::remove_all(dir);
}

TEST(KiraTest, ReadTargetTable_HandlesMultipleRules) {
    auto fc = qft::FamilyConfig::build(
        "bubble", {"l"}, {"p"}, {}, {{"p^2", "s"}},
        {"l^2 - msq", "(l - p)^2 - msq"});

    std::string dir =
        (fs::temp_directory_path() / "amflow_v2_kira_t9").string();
    fs::remove_all(dir);
    fs::create_directories(dir + "/results/bubble");

    {
        std::ofstream f(dir + "/results/bubble/kira_target.m");
        f << "{\n"
          << "  bubble(2,1) -> 3 * bubble(1,1),\n"
          << "  bubble(1,2) -> -2 * bubble(0,1) + d * bubble(1,1)\n"
          << "}\n";
    }

    auto rules = ibp::kira_read_target_table(fc, dir);
    ASSERT_EQ(rules.size(), 2u);
    EXPECT_EQ(rules[1].rhs.size(), 2u);
    EXPECT_EQ(rules[1].rhs[1].first, "d");
}

// --- kira_target.m parser strictness ----------
//
// `kira_read_target_table` parses Kira's kira2math-emitted reduction
// table.  The tokenizer is custom and was flagged by the audit as
// "not fuzzed".  Production Kira output never produces malformed
// rules, but the parser should still reject malformed input cleanly
// rather than silently producing a partial or wrong rule set.
//
// These tests cover the canonical malformed shapes:
//   1. LHS not a valid JIntegral (missing parens, garbage chars)
//   2. RHS term without a J integral (bare coefficient)
//   3. Missing file is *accepted* — returns empty (mirrors the
//      "no targets requested" Kira mode); existing
//      ReadMasters_MissingFileReturnsEmpty tests the analogous
//      behavior for masters.  This is the only "missing input
//      becomes empty output, not an error" exception in the
//      parser; lock it explicitly so future refactors don't
//      change the contract.

namespace {

// Write a synthetic kira_target.m to a fresh temp dir.  Returns
// the dir path so the test can call kira_read_target_table on it.
std::string write_synthetic_target_table(const std::string& subdir,
                                          const std::string& fam,
                                          const std::string& content) {
    std::string dir = (fs::temp_directory_path() / subdir).string();
    fs::remove_all(dir);
    fs::create_directories(dir + "/results/" + fam);
    {
        std::ofstream f(dir + "/results/" + fam + "/kira_target.m");
        f << content;
    }
    return dir;
}

}  // namespace

TEST(KiraTest, ReadTargetTable_MissingFileReturnsEmpty) {
    auto fc = qft::FamilyConfig::build(
        "bubble", {"l"}, {"p"}, {}, {{"p^2", "s"}},
        {"l^2 - msq", "(l - p)^2 - msq"});
    std::string dir = (fs::temp_directory_path()
                       / "amflow_kira_target_t_missing").string();
    fs::remove_all(dir);
    fs::create_directories(dir);
    // No results/bubble/kira_target.m written.
    auto rules = ibp::kira_read_target_table(fc, dir);
    EXPECT_TRUE(rules.empty())
        << "missing kira_target.m should produce empty rules, not throw";
    fs::remove_all(dir);
}

TEST(KiraTest, ReadTargetTable_MalformedLHSThrows) {
    auto fc = qft::FamilyConfig::build(
        "bubble", {"l"}, {"p"}, {}, {{"p^2", "s"}},
        {"l^2 - msq", "(l - p)^2 - msq"});
    // LHS is garbage that doesn't parse as a JIntegral.
    std::string dir = write_synthetic_target_table(
        "amflow_kira_target_t_lhs", "bubble",
        "{\n"
        "  garbage_no_parens -> 3 * bubble(1,1)\n"
        "}\n");
    EXPECT_THROW(ibp::kira_read_target_table(fc, dir),
                 std::runtime_error)
        << "malformed LHS should raise";
    fs::remove_all(dir);
}

TEST(KiraTest, ReadTargetTable_RHSTermWithoutJIntegralThrows) {
    auto fc = qft::FamilyConfig::build(
        "bubble", {"l"}, {"p"}, {}, {{"p^2", "s"}},
        {"l^2 - msq", "(l - p)^2 - msq"});
    // RHS contains a coefficient with no `* J(...)` part.  Note `0`
    // is allowed (filtered out as a zero rule), but a non-zero bare
    // coefficient must raise.
    std::string dir = write_synthetic_target_table(
        "amflow_kira_target_t_rhs", "bubble",
        "{\n"
        "  bubble(2,1) -> 3\n"
        "}\n");
    EXPECT_THROW(ibp::kira_read_target_table(fc, dir),
                 std::runtime_error)
        << "RHS term without J integral should raise";
    fs::remove_all(dir);
}

TEST(KiraTest, ReadTargetTable_NoOuterBraceReturnsEmpty) {
    auto fc = qft::FamilyConfig::build(
        "bubble", {"l"}, {"p"}, {}, {{"p^2", "s"}},
        {"l^2 - msq", "(l - p)^2 - msq"});
    // File present but has no `{...}` block.  This shape would
    // result from an upstream Kira misconfiguration.  Locked
    // behavior: parse silently produces empty (callers detect
    // empty masters and abort with a clear message).
    std::string dir = write_synthetic_target_table(
        "amflow_kira_target_t_nobrace", "bubble",
        "totally garbage content with no braces\n");
    auto rules = ibp::kira_read_target_table(fc, dir);
    EXPECT_TRUE(rules.empty())
        << "no-outer-brace input should yield empty rules";
    fs::remove_all(dir);
}

TEST(KiraTest, ReadTargetTable_RHSZeroFiltered) {
    // Positive sanity: a rule with `-> 0` produces a valid rule
    // with empty rhs.  This is the "zero reduction" pattern (the
    // target reduces to nothing), and the parser must filter the
    // literal `0` term without throwing.
    auto fc = qft::FamilyConfig::build(
        "bubble", {"l"}, {"p"}, {}, {{"p^2", "s"}},
        {"l^2 - msq", "(l - p)^2 - msq"});
    std::string dir = write_synthetic_target_table(
        "amflow_kira_target_t_zero", "bubble",
        "{\n"
        "  bubble(2,1) -> 0\n"
        "}\n");
    auto rules = ibp::kira_read_target_table(fc, dir);
    ASSERT_EQ(rules.size(), 1u);
    EXPECT_TRUE(rules[0].rhs.empty())
        << "literal `0` RHS should yield a rule with empty rhs vector";
    fs::remove_all(dir);
}
