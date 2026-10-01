// SPDX-License-Identifier: MIT
// Regression for the monolithic-eta FireFly "Regenerate" 0-CPU deadlock
// (the aux-mass-DE stall) and the kira_run stall watchdog that prevents it.
//
// The deadlock signature is process-level: the kira tree burns ZERO CPU
// and blocks on a pipe read forever (main kira on anon_pipe_read, FireFly
// workers idle).  Nothing dies, so the pre-fix blocking wait4 in kira_run
// hangs INDEFINITELY and the SIG6->FireFly auto-staging never triggers.
//
// We do not need a multi-GB Kira reduce to exercise this: we point
// KiraConfig::kira_executable at a tiny fake-kira that reproduces exactly
// the signature (idle worker children + a blocking FIFO read at 0 CPU),
// and assert the watchdog detects the stall, kills the tree, and throws a
// recognizable "STALLED" error PROMPTLY -- and that it does NOT fire on a
// legitimately busy (CPU-burning) child.

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <string>

#include <gtest/gtest.h>

#include "amflow/ibp/kira.hpp"

namespace ibp = amflow::ibp;

namespace {

std::string make_dir() {
    char tmpl[] = "/tmp/eta_stall_XXXXXX";
    char* d = mkdtemp(tmpl);
    return d ? std::string(d) : std::string();
}

std::string write_script(const std::string& dir, const std::string& name,
                         const std::string& body) {
    std::string path = dir + "/" + name;
    std::ofstream f(path);
    f << body;
    f.close();
    [[maybe_unused]] int rc = std::system(("chmod +x " + path).c_str());
    return path;
}

}  // namespace

// A 0-CPU stall must be detected, the tree killed, and a STALLED error
// thrown promptly -- NOT hang forever.
TEST(KiraStallWatchdog, ZeroCpuStallDetectedAndKilled) {
    std::string dir = make_dir();
    ASSERT_FALSE(dir.empty());
    // Idle 'FireFly worker' children + main 'kira' blocking on a FIFO read
    // with no writer: the whole tree sits at 0 CPU indefinitely.
    std::string exe = write_script(dir, "fake_kira_stall.sh",
        "#!/bin/bash\n"
        "sleep 100000 &\n"
        "sleep 100000 &\n"
        "fifo=$(mktemp -u)\n"
        "mkfifo \"$fifo\"\n"
        "exec 3<\"$fifo\"\n"
        "read -u 3 _ignored\n");

    ibp::KiraConfig cfg;
    cfg.kira_executable = exe;   // kira_run never dereferences cfg.fc here
    cfg.n_thread = 1;
    setenv("AMFLOW_KIRA_STALL_SECS", "4", 1);  // short window for the test

    auto t0 = std::chrono::steady_clock::now();
    bool threw_stall = false;
    try {
        ibp::kira_run(cfg, dir, dir + "/stall.log");
    } catch (const std::exception& e) {
        threw_stall =
            std::string(e.what()).find("STALLED") != std::string::npos;
    }
    double secs = std::chrono::duration<double>(
                      std::chrono::steady_clock::now() - t0).count();

    EXPECT_TRUE(threw_stall) << "watchdog did not throw a STALLED error";
    EXPECT_LT(secs, 30.0) << "watchdog too slow / nearly hung (" << secs
                          << " s)";
}

// A legitimately busy child (CPU-burning) longer than the stall window
// must NOT be falsely killed: runtime is never the trigger, only 0-CPU is.
TEST(KiraStallWatchdog, BusyChildNotFalselyStalled) {
    std::string dir = make_dir();
    ASSERT_FALSE(dir.empty());
    // Burn CPU for ~6 s with a 4 s stall window: only continuous CPU
    // progress keeps it alive past the window.
    std::string exe = write_script(dir, "fake_kira_busy.sh",
        "#!/bin/bash\n"
        "end=$((SECONDS+6))\n"
        "x=0\n"
        "while [ $SECONDS -lt $end ]; do x=$((x+1)); done\n"
        "exit 0\n");

    ibp::KiraConfig cfg;
    cfg.kira_executable = exe;
    cfg.n_thread = 1;
    setenv("AMFLOW_KIRA_STALL_SECS", "4", 1);

    EXPECT_NO_THROW(ibp::kira_run(cfg, dir, dir + "/busy.log"));
}

// Watchdog disabled (AMFLOW_KIRA_STALL_SECS=0) restores the original
// blocking wait4: a quick busy child still completes normally.
TEST(KiraStallWatchdog, DisabledRestoresBlockingWait) {
    std::string dir = make_dir();
    ASSERT_FALSE(dir.empty());
    std::string exe = write_script(dir, "fake_kira_quick.sh",
        "#!/bin/bash\n"
        "end=$((SECONDS+2))\n"
        "x=0\n"
        "while [ $SECONDS -lt $end ]; do x=$((x+1)); done\n"
        "exit 0\n");

    ibp::KiraConfig cfg;
    cfg.kira_executable = exe;
    cfg.n_thread = 1;
    setenv("AMFLOW_KIRA_STALL_SECS", "0", 1);

    EXPECT_NO_THROW(ibp::kira_run(cfg, dir, dir + "/quick.log"));
    unsetenv("AMFLOW_KIRA_STALL_SECS");
}
