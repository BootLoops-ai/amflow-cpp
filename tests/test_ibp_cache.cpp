// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Anthropic, PBC. Added by the BootLoops fork of AMFlow.cpp; see NOTICE and PATCHES.md.
// Tests for amflow::ibp::ibp_cache — key canonicalization.
//
// The cache key is FNV-1a over the Kira input files.  Env knobs that
// change jobs.yaml content but NOT the Kira result (DCAP search-bound
// cap, n_thread parallelism) used to perturb the key → ~2h cold-cache
// kira on relaunch (measured on a cold-start farm).  These tests pin the
// canonicalization that fixes that pitfall.

#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <string>
#include <gtest/gtest.h>
#include <sys/wait.h>
#include <unistd.h>

#include "amflow/ibp/ibp_cache.hpp"
#include "amflow/ibp/load_adaptive.hpp"

namespace ibp = amflow::ibp;
namespace fs  = std::filesystem;

namespace {

// Reduce-mode jobs.yaml WITHOUT the DCAP `, d: 3` (AMFLOW_REDUCE_DCAP unset).
const char* kReduceNoDcap =
    "jobs:\n"
    " - reduce_sectors:\n"
    "    reduce:\n"
    "     - {topologies: [cbox3l], sectors: [63], r: 9, s: 1}\n"
    "    truncate_sp:\n"
    "     - {topologies: [cbox3l], l: 3}\n"
    "    select_integrals:\n"
    "      select_mandatory_list:\n"
    "        - [cbox3l, target]\n"
    "    preferred_masters: preferred\n"
    "    integral_ordering: 5\n"
    "    run_initiate: true\n"
    "    run_triangular: true\n"
    "    run_back_substitution: true\n"
    " - kira2math:\n"
    "    target:\n"
    "     - [cbox3l, target]\n";

// Same Reduce-mode jobs.yaml WITH `, d: 3` (AMFLOW_REDUCE_DCAP=1).
const char* kReduceDcap =
    "jobs:\n"
    " - reduce_sectors:\n"
    "    reduce:\n"
    "     - {topologies: [cbox3l], sectors: [63], r: 9, s: 1, d: 3}\n"
    "    truncate_sp:\n"
    "     - {topologies: [cbox3l], l: 3}\n"
    "    select_integrals:\n"
    "      select_mandatory_list:\n"
    "        - [cbox3l, target]\n"
    "    preferred_masters: preferred\n"
    "    integral_ordering: 5\n"
    "    run_initiate: true\n"
    "    run_triangular: true\n"
    "    run_back_substitution: true\n"
    " - kira2math:\n"
    "    target:\n"
    "     - [cbox3l, target]\n";

// Masters-mode jobs.yaml: `d:` lives in select_mandatory_recursively
// (RESULT-determining) and must be PRESERVED.
const char* kMastersD4 =
    "jobs:\n"
    " - reduce_sectors:\n"
    "    reduce:\n"
    "     - {topologies: [cbox3l], sectors: [63], r: 10, s: 1}\n"
    "    select_integrals:\n"
    "      select_mandatory_recursively:\n"
    "        - {topologies: [cbox3l], sectors: [63], r: 10, s: 1, d: 4}\n"
    "    preferred_masters: preferred\n"
    "    integral_ordering: 5\n"
    "    run_initiate: masters\n";

const char* kMastersD5 =
    "jobs:\n"
    " - reduce_sectors:\n"
    "    reduce:\n"
    "     - {topologies: [cbox3l], sectors: [63], r: 10, s: 1}\n"
    "    select_integrals:\n"
    "      select_mandatory_recursively:\n"
    "        - {topologies: [cbox3l], sectors: [63], r: 10, s: 1, d: 5}\n"
    "    preferred_masters: preferred\n"
    "    integral_ordering: 5\n"
    "    run_initiate: masters\n";

void write(const fs::path& p, const std::string& s) {
    fs::create_directories(p.parent_path());
    std::ofstream(p) << s;
}

struct TmpDir {
    fs::path p;
    TmpDir() {
        p = fs::temp_directory_path() /
            ("amflow_ibp_cache_test_" + std::to_string(::getpid()) + "_" +
             std::to_string(reinterpret_cast<std::uintptr_t>(this)));
        fs::create_directories(p / "config");
    }
    ~TmpDir() { std::error_code ec; fs::remove_all(p, ec); }
};

}  // namespace

TEST(IbpCacheCanon, ReduceDcapStripped) {
    // The DCAP `, d: <int>` in the reduce: seed line must canonicalize away.
    EXPECT_EQ(ibp::canonicalize_jobs_yaml(kReduceDcap),
              ibp::canonicalize_jobs_yaml(kReduceNoDcap));
    // And the canonical form is exactly the no-DCAP yaml (idempotent).
    EXPECT_EQ(ibp::canonicalize_jobs_yaml(kReduceNoDcap),
              std::string(kReduceNoDcap));
}

TEST(IbpCacheCanon, MastersDPreserved) {
    // The `d:` inside select_mandatory_recursively is RESULT-determining
    // (bounds which masters are enumerated) and must NOT be stripped.
    EXPECT_NE(ibp::canonicalize_jobs_yaml(kMastersD4),
              ibp::canonicalize_jobs_yaml(kMastersD5));
    // No-op on Masters yaml (no d: in the reduce: seed line).
    EXPECT_EQ(ibp::canonicalize_jobs_yaml(kMastersD4), std::string(kMastersD4));
}

TEST(IbpCacheCanon, ThreadsLineDropped) {
    std::string with =
        "jobs:\n"
        " - reduce_sectors:\n"
        "    n_thread: 6\n"
        "    reduce:\n"
        "     - {topologies: [f], sectors: [1], r: 2, s: 0}\n";
    std::string without =
        "jobs:\n"
        " - reduce_sectors:\n"
        "    reduce:\n"
        "     - {topologies: [f], sectors: [1], r: 2, s: 0}\n";
    EXPECT_EQ(ibp::canonicalize_jobs_yaml(with),
              ibp::canonicalize_jobs_yaml(without));
    // Variant spellings.
    EXPECT_EQ(ibp::canonicalize_jobs_yaml("  threads: 12\nfoo\n"),  "foo\n");
    EXPECT_EQ(ibp::canonicalize_jobs_yaml("n_threads: 1\nfoo\n"),   "foo\n");
}

TEST(IbpCacheCanon, KeyInvariantUnderDcap) {
    // End-to-end: two work dirs differing ONLY in the DCAP `, d: 3`
    // must produce the SAME ibp_cache_key().
    TmpDir a, b;
    const std::string kin =
        "kinematics:\n  incoming_momenta: [p1]\n  kinematic_invariants:\n";
    const std::string fam =
        "integralfamilies:\n  - name: \"cbox3l\"\n";
    write(a.p / "config/kinematics.yaml",       kin);
    write(b.p / "config/kinematics.yaml",       kin);
    write(a.p / "config/integralfamilies.yaml", fam);
    write(b.p / "config/integralfamilies.yaml", fam);
    write(a.p / "preferred", "cbox3l[1,1,1,1,1,1,0,0,0,0,0,0,0,0,0]\n\n");
    write(b.p / "preferred", "cbox3l[1,1,1,1,1,1,0,0,0,0,0,0,0,0,0]\n\n");
    write(a.p / "target",    "cbox3l[1,1,2,1,1,1,0,0,0,0,0,0,0,0,0]\n\n");
    write(b.p / "target",    "cbox3l[1,1,2,1,1,1,0,0,0,0,0,0,0,0,0]\n\n");
    write(a.p / "jobs.yaml", kReduceNoDcap);
    write(b.p / "jobs.yaml", kReduceDcap);

    EXPECT_EQ(ibp::ibp_cache_key(a.p.string()), ibp::ibp_cache_key(b.p.string()));

    // Negative control: a genuine difference (different target file) must
    // still key differently.
    write(b.p / "target", "cbox3l[1,1,3,1,1,1,0,0,0,0,0,0,0,0,0]\n\n");
    EXPECT_NE(ibp::ibp_cache_key(a.p.string()), ibp::ibp_cache_key(b.p.string()));
}

// ==========================================================================
//  IbpCacheLock — stampede lock
// ==========================================================================

namespace {

// One child's worth of the stampede race.  Simulates the
// MISS→lock→re-check→(compute|HIT) path of `kira_run_cached`, with a
// fake 1 s "kira".  Touches `<sentinel_dir>/ran.<pid>` iff it actually
// ran the fake kira.  Writes "WAIT" / "HIT" / "RAN" to
// `<sentinel_dir>/log.<pid>` so the parent can assert the protocol.
// Returns an exit code (0 ok, 1 error) — never returns to gtest.
int stampede_child(const fs::path& cache_root,
                   const std::string& key,
                   const std::string& fam,
                   const fs::path& sentinel_dir) {
    std::ofstream log(sentinel_dir / ("log." + std::to_string(::getpid())));

    // Per-child work_dir so concurrent ibp_cache_lookup() copies don't
    // collide on the same destination path.
    fs::path work = sentinel_dir / ("work." + std::to_string(::getpid()));
    std::error_code ec;
    fs::create_directories(work, ec);

    // Initial lookup (the one before the lock in kira_run_cached).
    if (ibp::ibp_cache_lookup(cache_root.string(), key, fam, work.string())) {
        log << "HIT\n";
        return 0;
    }

    ibp::IbpCacheLock lock(cache_root.string(), key);
    if (lock.waited()) {
        log << "WAIT\n";
        if (ibp::ibp_cache_lookup(cache_root.string(), key, fam,
                                  work.string())) {
            log << "HIT\n";
            return 0;
        }
    }

    // We are the (only) computer.  Fake "kira": sleep 1 s, write a
    // results/<fam>/masters file, store, touch sentinel.
    log << "RAN\n";
    ::sleep(1);
    fs::create_directories(work / "results" / fam, ec);
    std::ofstream(work / "results" / fam / "masters")
        << "cbox3l[1,1,1,1,1,1,0,0,0,0,0,0,0,0,0]\n";
    ibp::ibp_cache_store(cache_root.string(), key, fam, work.string());
    std::ofstream(sentinel_dir / ("ran." + std::to_string(::getpid()))) << "1\n";
    return 0;
}

}  // namespace

TEST(IbpCacheLock, ThreeWayStampedeOneRunsKira) {
    // Three children race on the same cache key.  Exactly ONE should run
    // the fake kira (touch a `ran.*` sentinel); the other two should log
    // WAIT then HIT.
    TmpDir root;       // cache root
    TmpDir sentinels;  // ran.* / log.* / work.* live here
    const std::string key = "deadbeefcafebabe";
    const std::string fam = "cbox3l";

    ::setenv("AMFLOW_IBP_CACHE_LOCK", "1", 1);
    ::setenv("AMFLOW_IBP_LOCK_TIMEOUT_S", "30", 1);

    constexpr int N = 3;
    pid_t pids[N];
    for (int i = 0; i < N; ++i) {
        pid_t p = ::fork();
        ASSERT_GE(p, 0) << "fork failed";
        if (p == 0) {
            // Child.  _exit so gtest's atexit hooks don't run here.
            ::_exit(stampede_child(root.p, key, fam, sentinels.p));
        }
        pids[i] = p;
    }
    for (int i = 0; i < N; ++i) {
        int status = 0;
        ASSERT_EQ(::waitpid(pids[i], &status, 0), pids[i]);
        EXPECT_TRUE(WIFEXITED(status));
        EXPECT_EQ(WEXITSTATUS(status), 0);
    }

    // Exactly one `ran.*` sentinel.
    int n_ran = 0;
    for (const auto& e : fs::directory_iterator(sentinels.p)) {
        if (e.path().filename().string().rfind("ran.", 0) == 0) ++n_ran;
    }
    EXPECT_EQ(n_ran, 1) << "expected exactly one child to run kira";

    // The two non-runners must have logged WAIT then HIT.
    int n_wait_hit = 0;
    for (const auto& e : fs::directory_iterator(sentinels.p)) {
        auto name = e.path().filename().string();
        if (name.rfind("log.", 0) != 0) continue;
        std::ifstream f(e.path());
        std::string a, b;
        std::getline(f, a);
        std::getline(f, b);
        if (a == "WAIT" && b == "HIT") ++n_wait_hit;
    }
    EXPECT_EQ(n_wait_hit, N - 1)
        << "expected " << (N - 1) << " children to log WAIT then HIT";

    // The cache entry exists and is non-empty.
    fs::path stored = root.p / "ibp" / key / "results" / fam / "masters";
    EXPECT_TRUE(fs::exists(stored));
    EXPECT_GT(fs::file_size(stored), 0u);

    ::unsetenv("AMFLOW_IBP_LOCK_TIMEOUT_S");
    ::unsetenv("AMFLOW_IBP_CACHE_LOCK");
}

TEST(IbpCacheLock, DisabledIsNoop) {
    TmpDir root;
    ::setenv("AMFLOW_IBP_CACHE_LOCK", "0", 1);
    {
        ibp::IbpCacheLock lock(root.p.string(), "0123456789abcdef");
        EXPECT_FALSE(lock.held());
        EXPECT_FALSE(lock.waited());
    }
    ::unsetenv("AMFLOW_IBP_CACHE_LOCK");
    // Default ON.
    EXPECT_TRUE(ibp::IbpCacheLock::enabled());
}

TEST(IbpCacheLock, TimeoutFallthrough) {
    // Holder grabs the lock and sleeps past the (short) timeout; the
    // waiter must log WAIT, time out, and proceed (waited()==true,
    // held()==false).
    TmpDir root;
    const std::string key = "feedfacefeedface";
    ::setenv("AMFLOW_IBP_CACHE_LOCK", "1", 1);
    ::setenv("AMFLOW_IBP_LOCK_TIMEOUT_S", "1", 1);

    pid_t holder = ::fork();
    ASSERT_GE(holder, 0);
    if (holder == 0) {
        ibp::IbpCacheLock lk(root.p.string(), key);
        ::sleep(3);  // hold past the 1 s timeout
        ::_exit(0);
    }
    // Give the holder a head start.
    ::usleep(200 * 1000);
    {
        ibp::IbpCacheLock lk(root.p.string(), key);
        EXPECT_TRUE(lk.waited());
        EXPECT_FALSE(lk.held());  // timed out, proceeding anyway
    }
    int st = 0;
    ::waitpid(holder, &st, 0);

    ::unsetenv("AMFLOW_IBP_LOCK_TIMEOUT_S");
    ::unsetenv("AMFLOW_IBP_CACHE_LOCK");
}

// ==========================================================================
//  LoadAdaptive — load-adaptive parallelism formula
// ==========================================================================

TEST(LoadAdaptive, OversubscribedClampsToOne) {
    // 96 cores, load 250, 14 farm jobs → (96-250)/14 < 0 → clamp to 1.
    ::setenv("AMFLOW_NPROC_OVERRIDE",   "96",    1);
    ::setenv("AMFLOW_LOADAVG_OVERRIDE", "250.0", 1);
    ::setenv("AMFLOW_FARM_N",           "14",    1);
    EXPECT_EQ(ibp::load_adaptive_threads(/*cap=*/6), 1);
    EXPECT_EQ(ibp::load_adaptive_threads(/*cap=*/0), 1);  // eps_parallel: no cap
    ::unsetenv("AMFLOW_NPROC_OVERRIDE");
    ::unsetenv("AMFLOW_LOADAVG_OVERRIDE");
    ::unsetenv("AMFLOW_FARM_N");
}

TEST(LoadAdaptive, SpareCapacityCappedByOptsNThread) {
    // 96 cores, load 10, 4 farm jobs → floor((96-10)/4)=21 → min(6,21)=6.
    ::setenv("AMFLOW_NPROC_OVERRIDE",   "96",   1);
    ::setenv("AMFLOW_LOADAVG_OVERRIDE", "10.0", 1);
    ::setenv("AMFLOW_FARM_N",           "4",    1);
    EXPECT_EQ(ibp::load_adaptive_threads(/*cap=*/6),  6);
    EXPECT_EQ(ibp::load_adaptive_threads(/*cap=*/40), 21);
    EXPECT_EQ(ibp::load_adaptive_threads(/*cap=*/0),  21);  // uncapped
    ::unsetenv("AMFLOW_NPROC_OVERRIDE");
    ::unsetenv("AMFLOW_LOADAVG_OVERRIDE");
    ::unsetenv("AMFLOW_FARM_N");
}

TEST(LoadAdaptive, FarmNDefaultsToOne) {
    ::unsetenv("AMFLOW_FARM_N");
    ::setenv("AMFLOW_NPROC_OVERRIDE",   "8",   1);
    ::setenv("AMFLOW_LOADAVG_OVERRIDE", "3.0", 1);
    // (8-3)/1 = 5, cap 12 → 5.
    EXPECT_EQ(ibp::load_adaptive_threads(/*cap=*/12), 5);
    ::unsetenv("AMFLOW_NPROC_OVERRIDE");
    ::unsetenv("AMFLOW_LOADAVG_OVERRIDE");
}
