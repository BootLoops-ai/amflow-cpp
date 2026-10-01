// SPDX-License-Identifier: MIT
// Integration tests for amflow_cli.
//
// We invoke the binary as a subprocess with a trivial 1x1 DE and
// inspect the JSON output.

#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <gtest/gtest.h>
#include <nlohmann/json.hpp>
#include <sstream>
#include <string>

namespace fs = std::filesystem;
using nlohmann::json;

namespace {

#ifndef AMFLOW_CLI_PATH
#define AMFLOW_CLI_PATH ""
#endif

constexpr const char* kCliPath = AMFLOW_CLI_PATH;

bool cli_available() {
    return *kCliPath && fs::exists(kCliPath);
}

}  // namespace

TEST(CliIntegration, TrivialAmflowEchoesBoundary) {
    if (!cli_available()) {
        GTEST_SKIP() << "amflow_cli not built";
    }

    fs::path input_path =
        fs::temp_directory_path() / "amflow_cli_int_in.json";
    fs::path output_path =
        fs::temp_directory_path() / "amflow_cli_int_out.json";

    {
        std::ofstream of(input_path);
        of << R"({
            "mode": "amflow",
            "options": { "silent_mode": true },
            "matrix": [[ {"num": ["0"], "den": ["1"]} ]],
            "boundaries": [[
                {"mu": {"re": "0", "im": "0"},
                 "value": {"re": "7", "im": "0"}}
            ]]
        })";
    }

    std::string cmd = std::string(kCliPath) + " " + input_path.string()
                       + " " + output_path.string();
    int rc = std::system(cmd.c_str());
    ASSERT_EQ(rc, 0) << "cli exited with " << rc;

    std::ifstream in(output_path);
    json out; in >> out;
    ASSERT_EQ(out.at("mode"), "amflow");
    ASSERT_TRUE(out.contains("result"));
    ASSERT_EQ(out.at("result").size(), 1u);
    EXPECT_EQ(out.at("result")[0].at("re").get<std::string>().substr(0, 1), "7");

    fs::remove(input_path);
    fs::remove(output_path);
}

TEST(CliIntegration, UnknownModeExitsNonZero) {
    if (!cli_available()) {
        GTEST_SKIP() << "amflow_cli not built";
    }

    fs::path input_path =
        fs::temp_directory_path() / "amflow_cli_bad_in.json";
    {
        std::ofstream of(input_path);
        of << R"({"mode": "no_such_mode"})";
    }

    std::string cmd = std::string(kCliPath) + " " + input_path.string()
                       + " 2>/dev/null > /dev/null";
    int rc = std::system(cmd.c_str());
    EXPECT_NE(WEXITSTATUS(rc), 0);

    fs::remove(input_path);
}
