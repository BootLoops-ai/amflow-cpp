// SPDX-License-Identifier: MIT
// Tests for amflow::api::run_json.

#include <cstdlib>

#include <gtest/gtest.h>
#include <nlohmann/json.hpp>

#include "amflow/api/run_json.hpp"

namespace api = amflow::api;
using nlohmann::json;

TEST(ApiRunJsonTest, UnknownModeThrows) {
    json input = {{"mode", "no_such_mode"}};
    EXPECT_THROW(api::run_json(input), std::runtime_error);
}

TEST(ApiRunJsonTest, OptionsAreEchoedBack) {
    // Setting a recognised option in input should be visible in output.options.
    json input = {
        {"mode", "amflow"},
        {"options", {{"silent_mode", true}, {"working_pre", 60}}},
        // Provide a tiny 1x1 DE that solves the boundary mu=0,value=1 trivially.
        {"matrix", json::array({json::array({json{
            {"num", json::array({"0"})},
            {"den", json::array({"1"})},
        }})})},
        {"boundaries", json::array({json::array({json{
            {"mu",    json{{"re", "0"}, {"im", "0"}}},
            {"value", json{{"re", "1"}, {"im", "0"}}},
        }})})},
    };

    json out = api::run_json(input);
    EXPECT_EQ(out.at("mode"), "amflow");
    ASSERT_TRUE(out.contains("options"));
    EXPECT_EQ(out.at("options").at("silent_mode"), true);
    EXPECT_EQ(out.at("options").at("working_pre"), 60);
    ASSERT_TRUE(out.contains("result"));
    EXPECT_TRUE(out.at("result").is_array());
    EXPECT_EQ(out.at("result").size(), 1u);
}

TEST(ApiRunJsonTest, AmflowMissingMatrixThrows) {
    json input = {
        {"mode", "amflow"},
        {"boundaries", json::array()},
    };
    EXPECT_THROW(api::run_json(input), std::exception);
}

TEST(ApiRunJsonTest, AmflowBoundariesSizeMismatchThrows) {
    json input = {
        {"mode", "amflow"},
        {"matrix", json::array({json::array({json{
            {"num", json::array({"0"})},
            {"den", json::array({"1"})},
        }})})},
        // 0 boundaries for a 1x1 system.
        {"boundaries", json::array()},
    };
    EXPECT_THROW(api::run_json(input), std::runtime_error);
}

// Complex-valued numeric kinematics is intentionally not ported
// (upstream `IBPRule` / `CompensateRule` at `Kira/interface.m:50-57`).
// The dispatcher must reject the object form `{"re":..,"im":..}` in
// `amf_options.blackbox.numeric_values` with a clear error rather
// than silently truncating to the real part.
TEST(ApiRunJsonTest, BlackboxComplexNumericValueIsRejected) {
    json input = {
        {"mode", "black_box_amflow"},
        {"family", {
            {"family", "fam"},
            {"loops", json::array({"l1"})},
            {"legs", json::array()},
            {"propagators", json::array({"l1*l1"})},
        }},
        {"targets", json::array({json::array({1})})},
        {"eps_samples", json::array({json{{"re", "1/1000"}, {"im", "0"}}})},
        {"amf_options", {{"blackbox", {{"numeric_values", {
            {"s", {{"re", "100"}, {"im", "1"}}},
        }}}}}},
    };
    try {
        api::run_json(input);
        FAIL() << "expected complex-form numeric_values to be rejected";
    } catch (const std::runtime_error& e) {
        const std::string what(e.what());
        EXPECT_NE(what.find("complex-numeric"), std::string::npos)
            << "error should name the unsupported feature, got: " << what;
        EXPECT_NE(what.find("AUDIT_MMA_PARITY.md"), std::string::npos)
            << "error should point at the parity audit doc, got: " << what;
    }
}

// -------------------------------------------------------------------------
// Config hygiene: a knob in the wrong JSON object must warn
// loudly on stderr (once per key), with a hint naming where the key WOULD
// be read.  The canonical disease: `eps_parallel` under `options` instead
// of `amf_options` parsed cleanly, did nothing, and cost a serial night.
// AMFLOW_STRICT_CONFIG=1 upgrades the warning to an abort.

namespace {

json trivial_amflow_input(json options) {
    return json{
        {"mode", "amflow"},
        {"options", std::move(options)},
        {"matrix", json::array({json::array({json{
            {"num", json::array({"0"})},
            {"den", json::array({"1"})},
        }})})},
        {"boundaries", json::array({json::array({json{
            {"mu",    json{{"re", "0"}, {"im", "0"}}},
            {"value", json{{"re", "1"}, {"im", "0"}}},
        }})})},
    };
}

}  // namespace

TEST(ApiRunJsonTest, MisplacedEpsParallelWarnsOnStderrOnce) {
    json input = trivial_amflow_input(
        {{"silent_mode", true}, {"eps_parallel", 4}});

    testing::internal::CaptureStderr();
    json out = api::run_json(input);
    const std::string err = testing::internal::GetCapturedStderr();
    EXPECT_NE(err.find("'eps_parallel' found in options — not read; "
                       "use amf_options.eps_parallel or AMFLOW_EPS_PARALLEL"),
              std::string::npos)
        << "stderr was: " << err;
    ASSERT_TRUE(out.contains("result"));

    // Once per (object, key) per process: a second identical run is silent.
    testing::internal::CaptureStderr();
    api::run_json(input);
    const std::string err2 = testing::internal::GetCapturedStderr();
    EXPECT_EQ(err2.find("eps_parallel"), std::string::npos)
        << "second stderr was: " << err2;
}

TEST(ApiRunJsonTest, MisplacedNodeParallelWarnsOnStderrOnce) {
    json input = trivial_amflow_input(
        {{"silent_mode", true}, {"node_parallel", 4}});

    testing::internal::CaptureStderr();
    json out = api::run_json(input);
    const std::string err = testing::internal::GetCapturedStderr();
    EXPECT_NE(err.find("'node_parallel' found in options — not read; "
                       "use amf_options.node_parallel or AMFLOW_NODE_PARALLEL"),
              std::string::npos)
        << "stderr was: " << err;
    ASSERT_TRUE(out.contains("result"));

    // Once per (object, key) per process: a second identical run is silent.
    testing::internal::CaptureStderr();
    api::run_json(input);
    const std::string err2 = testing::internal::GetCapturedStderr();
    EXPECT_EQ(err2.find("node_parallel"), std::string::npos)
        << "second stderr was: " << err2;
}

TEST(ApiRunJsonTest, UnknownEverywhereKeyListsAcceptedSet) {
    json input = trivial_amflow_input(
        {{"silent_mode", true}, {"working_prec_typo", 60}});

    testing::internal::CaptureStderr();
    api::run_json(input);
    const std::string err = testing::internal::GetCapturedStderr();
    EXPECT_NE(err.find("'working_prec_typo' found in options"),
              std::string::npos) << "stderr was: " << err;
    EXPECT_NE(err.find("accepted keys in options"), std::string::npos)
        << "stderr was: " << err;
    EXPECT_NE(err.find("working_pre"), std::string::npos)
        << "stderr was: " << err;
}

TEST(ApiRunJsonTest, UnderscoreCommentKeysNeverWarn) {
    json input = trivial_amflow_input(
        {{"silent_mode", true}, {"_comment", "comment-key convention"}});

    testing::internal::CaptureStderr();
    api::run_json(input);
    const std::string err = testing::internal::GetCapturedStderr();
    EXPECT_EQ(err.find("_comment"), std::string::npos)
        << "stderr was: " << err;
}

TEST(ApiRunJsonTest, StrictConfigUpgradesWarningToAbort) {
    ::setenv("AMFLOW_STRICT_CONFIG", "1", 1);
    json input = trivial_amflow_input(
        {{"silent_mode", true}, {"eps_paralel_strict_typo", 4}});
    try {
        api::run_json(input);
        ::unsetenv("AMFLOW_STRICT_CONFIG");
        FAIL() << "expected AMFLOW_STRICT_CONFIG=1 to abort on unknown key";
    } catch (const std::runtime_error& e) {
        ::unsetenv("AMFLOW_STRICT_CONFIG");
        const std::string what(e.what());
        EXPECT_NE(what.find("AMFLOW_STRICT_CONFIG"), std::string::npos)
            << "error was: " << what;
        EXPECT_NE(what.find("eps_paralel_strict_typo"), std::string::npos)
            << "error was: " << what;
    }
}
