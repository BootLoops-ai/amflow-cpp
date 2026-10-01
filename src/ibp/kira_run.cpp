// SPDX-License-Identifier: MIT
// ibp::kira_run + kira_read_masters + kira_read_target_table.
//

#include "amflow/ibp/kira.hpp"

#include "amflow/ibp/reduce.hpp"  // kira_ffsave_enabled

#include <algorithm>
#include <cctype>
#include <chrono>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <errno.h>
#include <fcntl.h>
#include <signal.h>
#include <sys/resource.h>
#include <sys/wait.h>
#include <unistd.h>

extern "C" char** environ;

namespace amflow::ibp {

namespace fs = std::filesystem;

namespace {

// === Aux-mass-DE monolithic-eta FireFly STALL watchdog ===========
//
// The monolithic eta-deformed top-sector reduce (dot-cap OFF / large
// mandatory list, e.g. a 15-prop 3-loop family at r:17,s:3, mandatory ~499) can DEADLOCK
// rather than abort: the main kira process blocks on anon_pipe_read in
// the FireFly "Regenerate:" scan phase while its FireFly black-box
// workers sit idle -- the whole tree burns 0 CPU and consumes ~32 MB
// RSS indefinitely.  Because nothing dies, a plain blocking wait4 hangs
// FOREVER (no SIG6/SIG11/SIG9 ever arrives), so the existing
// SIG6->FireFly auto-staging in kira_run_cached never triggers.
//
// The design here makes the hang IMPOSSIBLE: we wait NON-blocking and
// watch the aggregate CPU time of the child process TREE.  As long as
// the tree keeps burning CPU (a legitimately long reduce) we wait
// indefinitely -- runtime is never the trigger.  Only when the tree
// makes ZERO measurable CPU progress for a continuous window
// (AMFLOW_KIRA_STALL_SECS, default 180 s) do we declare a deadlock,
// kill the tree, and throw a recognizable "STALLED" error that
// kira_run_cached catches and re-stages via the dot-capped / FireFly
// path (so it auto-recovers instead of hanging).  Disable with
// AMFLOW_KIRA_STALL_SECS=0 (restores the old blocking wait4).

double stall_window_secs() {
    if (const char* v = std::getenv("AMFLOW_KIRA_STALL_SECS"); v && *v) {
        char* end = nullptr;
        double s = std::strtod(v, &end);
        if (end != v && s >= 0.0) return s;
    }
    return 180.0;  // default: 3 min of zero CPU progress => deadlock
}

// Parse field 14 (utime) + field 15 (stime), in clock ticks, from
// /proc/<pid>/stat.  Returns false if the process is gone/unreadable.
// The comm field (2) is parenthesized and may contain spaces, so we
// scan past the LAST ')' before tokenizing.
bool read_proc_cpu_ticks(pid_t pid, unsigned long long& out_ticks) {
    std::string path = "/proc/" + std::to_string(pid) + "/stat";
    std::ifstream f(path);
    if (!f) return false;
    std::string content;
    std::getline(f, content);
    std::size_t rp = content.rfind(')');
    if (rp == std::string::npos) return false;
    std::istringstream is(content.substr(rp + 1));
    // After ')': field 3 (state) onward.  utime=14, stime=15 overall,
    // i.e. tokens 12 and 13 of this remainder (1-based: state is 1).
    std::string tok;
    unsigned long long utime = 0, stime = 0;
    for (int idx = 3; is >> tok; ++idx) {
        if (idx == 14) { try { utime = std::stoull(tok); } catch (...) {} }
        else if (idx == 15) {
            try { stime = std::stoull(tok); } catch (...) {}
            break;
        }
    }
    out_ticks = utime + stime;
    return true;
}

// Sum utime+stime over the whole tree rooted at `root` by walking
// /proc and following ppid (field 4).  One /proc sweep per poll tick
// (~every few seconds) is cheap relative to a Kira reduce.  If /proc
// is unreadable we return the root's own ticks so the watchdog degrades
// to single-process monitoring rather than spuriously firing.
unsigned long long tree_cpu_ticks(pid_t root) {
    // Build pid -> ppid and pid -> ticks from /proc in one pass.
    std::vector<std::pair<pid_t, pid_t>> ppid_of;     // (pid, ppid)
    std::map<pid_t, unsigned long long> ticks_of;
    std::error_code ec;
    bool any = false;
    for (const auto& de : fs::directory_iterator("/proc", ec)) {
        if (ec) break;
        const std::string name = de.path().filename().string();
        if (name.empty() || !std::isdigit((unsigned char)name[0])) continue;
        pid_t pid = 0;
        try { pid = (pid_t)std::stol(name); } catch (...) { continue; }
        std::ifstream f(de.path() / "stat");
        if (!f) continue;
        std::string content;
        std::getline(f, content);
        std::size_t rp = content.rfind(')');
        if (rp == std::string::npos) continue;
        std::istringstream is(content.substr(rp + 1));
        std::string tok;
        pid_t ppid = 0;
        unsigned long long utime = 0, stime = 0;
        for (int idx = 3; is >> tok; ++idx) {
            if (idx == 4) { try { ppid = (pid_t)std::stol(tok); } catch (...) {} }
            else if (idx == 14) { try { utime = std::stoull(tok); } catch (...) {} }
            else if (idx == 15) {
                try { stime = std::stoull(tok); } catch (...) {}
                break;
            }
        }
        ppid_of.emplace_back(pid, ppid);
        ticks_of[pid] = utime + stime;
        any = true;
    }
    if (!any) {
        unsigned long long t = 0;
        return read_proc_cpu_ticks(root, t) ? t : 0;
    }
    // BFS from root over the ppid graph.
    std::vector<pid_t> stack{root};
    std::map<pid_t, bool> seen;
    unsigned long long total = 0;
    while (!stack.empty()) {
        pid_t cur = stack.back();
        stack.pop_back();
        if (seen[cur]) continue;
        seen[cur] = true;
        auto it = ticks_of.find(cur);
        if (it != ticks_of.end()) total += it->second;
        for (const auto& [pid, ppid] : ppid_of) {
            if (ppid == cur && !seen[pid]) stack.push_back(pid);
        }
    }
    return total;
}

// SIGTERM then (after a grace period) SIGKILL the whole process group.
// kira_run launches the child in its own group (setpgid) so one signal
// to -pgid reaches the FireFly/Fermat workers too.
void kill_tree(pid_t pgid) {
    ::kill(-pgid, SIGTERM);
    for (int i = 0; i < 20; ++i) {       // up to ~2 s grace
        if (::kill(-pgid, 0) != 0 && errno == ESRCH) return;
        struct timespec ts{0, 100 * 1000 * 1000};  // 100 ms
        nanosleep(&ts, nullptr);
    }
    ::kill(-pgid, SIGKILL);
}

// Read the last ~`max_bytes` of a (log) file as a string.  The Kira/Fermat/
// FireFly child's stdout+stderr were dup2'd onto this file, so its tail holds
// the *real* failure text: an uncaught C++ exception ("terminate called after
// throwing an instance of 'std::bad_alloc'"), a Fermat/FLINT error, a Kira
// assert message, etc.  Without surfacing this verbatim the bare
// "KILLED BY SIGNAL n" signal number is ambiguous (a self-capped OOM under
// RLIMIT_AS aborts with SIG6 via std::terminate, exactly like a Kira assert).
std::string read_log_tail(const std::string& log_path,
                          std::size_t max_bytes = 4096) {
    if (log_path.empty()) return {};
    std::error_code ec;
    if (!fs::exists(log_path, ec)) return {};
    std::ifstream f(log_path, std::ios::binary | std::ios::ate);
    if (!f) return {};
    std::streamoff size = f.tellg();
    if (size <= 0) return {};
    std::streamoff start =
        size > (std::streamoff)max_bytes ? size - (std::streamoff)max_bytes : 0;
    f.seekg(start, std::ios::beg);
    std::ostringstream oss;
    oss << f.rdbuf();
    std::string tail = oss.str();
    // Trim trailing whitespace/newlines for a tidier one-block report.
    while (!tail.empty() &&
           (tail.back() == '\n' || tail.back() == '\r' || tail.back() == ' ' ||
            tail.back() == '\t')) {
        tail.pop_back();
    }
    return tail;
}

// True if the child log tail shows the std::terminate -> abort(SIG6) path that
// a self-capped std::bad_alloc (RLIMIT_AS / AMFLOW_KIRA_MEM_CAP_GB) produces,
// as opposed to an internal Kira assert.  Used only to sharpen the diagnostic.
bool log_shows_bad_alloc(const std::string& tail) {
    return tail.find("std::bad_alloc") != std::string::npos ||
           tail.find("bad_alloc") != std::string::npos ||
           tail.find("std::bad_array_new_length") != std::string::npos ||
           tail.find("Cannot allocate memory") != std::string::npos ||
           tail.find("out of memory") != std::string::npos ||
           tail.find("OutOfMemory") != std::string::npos;
}

}  // namespace

double kira_run(const KiraConfig& cfg,
                const std::string& dir,
                const std::string& log_path,
                const std::string& jobs_file) {
    if (cfg.kira_executable.empty()) {
        throw std::invalid_argument("kira_run: kira_executable not set");
    }
    if (!fs::exists(cfg.kira_executable)) {
        throw std::runtime_error(
            "kira_run: kira executable not found at " + cfg.kira_executable);
    }
    if (!cfg.fermat_executable.empty() && !fs::exists(cfg.fermat_executable)) {
        throw std::runtime_error(
            "kira_run: fermat executable not found at " + cfg.fermat_executable);
    }
    if (!fs::is_directory(dir)) {
        throw std::runtime_error("kira_run: dir not found " + dir);
    }

    std::vector<std::string> args;
    args.push_back(cfg.kira_executable);
    args.push_back("-p" + std::to_string(cfg.n_thread));
    if (cfg.run_firefly && cfg.firefly_bunch > 0) {
        args.push_back("--bunch_size=" + std::to_string(cfg.firefly_bunch));
    }
    args.push_back(jobs_file.empty() ? std::string("jobs.yaml") : jobs_file);
    // Mirrors upstream `FilterRules[IBPRule, Prepend[MassScale, ep]]`
    // (Kira/interface.m:274).  After `ibp::kira_write_config` baked the
    // numeric_values into kinematics.yaml directly, the SP variables no
    // longer appear in `kinematic_invariants` -- passing `-s<var>=<val>`
    // for them confuses Kira (it locked up `set: s12 = -2` waiting on a
    // symbol that no longer exists in the YAML).  We only forward `-sd`
    // for the dimensional regulator `eps`, which is never substituted
    // into the YAML.
    for (const auto& [name, value] : cfg.numeric_values) {
        if (name == "eps") {
            args.push_back("-sd=4-2*" + value);
        }
    }

    std::vector<char*> argv_c;
    argv_c.reserve(args.size() + 1);
    for (auto& s : args) argv_c.push_back(const_cast<char*>(s.c_str()));
    argv_c.push_back(nullptr);

    std::vector<std::string> env_strings;
    bool fermat_set = false;
    bool ldpreload_set = false;
    bool mallocconf_set = false;
    for (char** e = environ; *e; ++e) {
        std::string entry(*e);
        if (entry.rfind("FERMATPATH=", 0) == 0) {
            if (!cfg.fermat_executable.empty()) {
                entry = "FERMATPATH=" + cfg.fermat_executable;
                fermat_set = true;
            }
        }
        if (entry.rfind("LD_PRELOAD=", 0) == 0)  ldpreload_set  = true;
        if (entry.rfind("MALLOC_CONF=", 0) == 0) mallocconf_set = true;
        env_strings.push_back(entry);
    }
    if (!fermat_set && !cfg.fermat_executable.empty()) {
        env_strings.push_back("FERMATPATH=" + cfg.fermat_executable);
    }
    // Kira 3.1 built with the glibc allocator fragments badly on the huge
    // "Select equations recursively" phase: a single contiguous allocation
    // can throw std::bad_alloc with hundreds of GB still aggregate-free.
    // jemalloc both fragments far less and (with the decay knobs below)
    // returns dirty pages to the OS, lowering peak RSS.  We now build Kira
    // with jemalloc *linked in*, but as a belt-and-suspenders the launcher
    // also LD_PRELOADs the system libjemalloc into the child if present so
    // older binaries / Fermat / FireFly children inherit it too.  Set
    // AMFLOW_KIRA_NO_JEMALLOC=1 to disable, or set LD_PRELOAD/MALLOC_CONF
    // yourself to override.
    if (!std::getenv("AMFLOW_KIRA_NO_JEMALLOC")) {
        const char* jepath = std::getenv("AMFLOW_JEMALLOC_PATH");
        std::string je = jepath ? jepath
                                : "/lib/x86_64-linux-gnu/libjemalloc.so.2";
        if (::access(je.c_str(), R_OK) == 0) {
            if (!ldpreload_set) {
                env_strings.push_back("LD_PRELOAD=" + je);
            }
            if (!mallocconf_set) {
                env_strings.push_back(
                    "MALLOC_CONF=background_thread:true,"
                    "dirty_decay_ms:5000,muzzy_decay_ms:5000,narenas:32");
            }
        }
    }
    std::vector<char*> envp_c;
    envp_c.reserve(env_strings.size() + 1);
    for (auto& s : env_strings) envp_c.push_back(const_cast<char*>(s.c_str()));
    envp_c.push_back(nullptr);

    // Kira stdout/stderr capture.  When the caller leaves `log_path` empty
    // (the CLI default), default to a *per-working-dir*
    // `kira_stdout.log` next to Kira's own fd-3 `kira.log`.  Previously the
    // CLI hardcoded a SHARED `/tmp/amflow_cli_kira.log`, which (a) was
    // O_TRUNC'ed by every concurrent amflow run and (b) interleaved events
    // from unrelated jobs, so monitors read the wrong job's progress.  The
    // working dir is already job-unique, so this is collision-free.
    const std::string effective_log_path =
        log_path.empty() ? (fs::path(dir) / "kira_stdout.log").string()
                         : log_path;
    // AMFLOW_KIRA_FFSAVE=1: APPEND instead of truncating.  Under the
    // resumable-dir regime a relaunch re-runs Kira in the SAME preserved
    // working dir; O_TRUNC would destroy the previous attempt's streamed
    // FireFly `Done: N / M` progress lines — exactly the salvage/rate
    // datum an interrupted run leaves behind (see the resume notes in BUILD.md).
    // Knob OFF (default): stock O_TRUNC, byte-identical.
    const int log_trunc_mode =
        kira_ffsave_enabled() ? O_APPEND : O_TRUNC;
    int log_fd = ::open(effective_log_path.c_str(),
                        O_WRONLY | O_CREAT | log_trunc_mode, 0644);
    if (log_fd < 0) {
        throw std::runtime_error(
            "kira_run: cannot open log file " + effective_log_path);
    }

    auto t0 = std::chrono::steady_clock::now();

    pid_t pid = fork();
    if (pid < 0) {
        if (log_fd >= 0) close(log_fd);
        throw std::runtime_error(
            "kira_run: fork failed (" + std::string(strerror(errno)) + ")");
    }
    if (pid == 0) {
        // Put the child (and the FireFly/Fermat workers it forks) into
        // its own process group so the stall watchdog can signal the
        // whole tree with kill(-pgid, ...) without touching unrelated
        // amflow threads/processes.  Best-effort: harmless if it fails.
        setpgid(0, 0);
        if (chdir(dir.c_str()) != 0) {
            _exit(127);
        }
        // Optional per-child address-space cap (Route C: memory-capped
        // reduce).  AMFLOW_KIRA_MEM_CAP_GB=<N> applies RLIMIT_AS=N GiB to
        // the Kira child (and, by inheritance, the Fermat/FireFly
        // black-box workers it forks).  This bounds a runaway *symbolic*
        // (Fermat `run_back_substitution`) coefficient swell — e.g. the
        // 3-loop crossed-box aux-mass eta-power tower — so it self-aborts
        // at the cap with a clean per-process allocation failure instead
        // of OOM-killing the whole host (or, worse, driving it into swap
        // death).  The FireFly finite-field path
        // (`AMFLOW_RUN_FIREFLY=1`) is the actual memory-light *fix*;
        // this cap is a fail-safe guard so the Fermat fallback can never
        // take down the box.  Unset (default) => no cap, byte-identical
        // to the prior behaviour.  RLIMIT_AS (not RLIMIT_DATA) is used so
        // mmap-backed FLINT/GMP allocations are counted.
        if (const char* cap = std::getenv("AMFLOW_KIRA_MEM_CAP_GB");
            cap && *cap) {
            char* end = nullptr;
            double gib = std::strtod(cap, &end);
            if (end != cap && gib > 0.0) {
                rlim_t bytes =
                    (rlim_t)(gib * 1024.0 * 1024.0 * 1024.0);
                struct rlimit rl;
                rl.rlim_cur = bytes;
                rl.rlim_max = bytes;
                // Best-effort: if the cap can't be set, proceed uncapped
                // rather than aborting the run.
                setrlimit(RLIMIT_AS, &rl);
            }
        }
        if (log_fd >= 0) {
            dup2(log_fd, STDOUT_FILENO);
            dup2(log_fd, STDERR_FILENO);
            close(log_fd);
        }
        execve(cfg.kira_executable.c_str(), argv_c.data(), envp_c.data());
        _exit(127);
    }
    if (log_fd >= 0) close(log_fd);

    // Parent-side setpgid too: whichever of parent/child wins the race,
    // the child lands in a group whose id == child pid (so pgid == pid).
    setpgid(pid, pid);

    int status = 0;
    struct rusage usage{};
    const double stall_secs = stall_window_secs();

    if (stall_secs <= 0.0) {
        // Watchdog disabled -> original blocking wait4 (byte-identical).
        while (wait4(pid, &status, 0, &usage) < 0) {
            if (errno == EINTR) continue;
            throw std::runtime_error(
                "kira_run: wait4 failed (" + std::string(strerror(errno)) +
                ")");
        }
    } else {
        // Non-blocking poll + CPU-progress watchdog.  We only fire when
        // the child TREE makes zero CPU progress for `stall_secs`; a
        // busy long reduce never trips it (runtime is not the trigger).
        unsigned long long last_ticks = tree_cpu_ticks(pid);
        auto last_progress = std::chrono::steady_clock::now();
        const auto poll = std::chrono::milliseconds(2000);
        bool stalled = false;
        while (true) {
            pid_t r = wait4(pid, &status, WNOHANG, &usage);
            if (r == pid) break;            // child reaped
            if (r < 0) {
                if (errno == EINTR) continue;
                throw std::runtime_error(
                    "kira_run: wait4 failed (" +
                    std::string(strerror(errno)) + ")");
            }
            // r == 0: still running.  Sleep, then check CPU progress.
            struct timespec ts{(time_t)(poll.count() / 1000),
                               (long)((poll.count() % 1000) * 1000000L)};
            nanosleep(&ts, nullptr);
            unsigned long long now_ticks = tree_cpu_ticks(pid);
            auto now = std::chrono::steady_clock::now();
            if (now_ticks != last_ticks) {
                last_ticks = now_ticks;
                last_progress = now;
            } else {
                double idle =
                    std::chrono::duration<double>(now - last_progress).count();
                if (idle >= stall_secs) {
                    stalled = true;
                    break;
                }
            }
        }
        if (stalled) {
            std::cerr << "[amflow][stall-watchdog] kira tree pid=" << pid
                      << " made 0 CPU progress for >= " << stall_secs
                      << " s (stall; observed with monolithic-eta FireFly "
                         "reduces); killing tree and re-staging\n";
            kill_tree(pid);
            // Reap the now-dead child so we don't leak a zombie.
            while (wait4(pid, &status, 0, &usage) < 0) {
                if (errno == EINTR) continue;
                break;
            }
            throw std::runtime_error(
                "kira_run: kira STALLED (no CPU progress for " +
                std::to_string((long)stall_secs) +
                " s; stall, observed with monolithic-eta FireFly reduces) in " + dir +
                " -- killed tree (pgid=" + std::to_string(pid) + ")");
        }
    }

    auto t1 = std::chrono::steady_clock::now();
    double elapsed =
        std::chrono::duration<double>(t1 - t0).count();

    // ru_maxrss is in KiB on Linux but BYTES on macOS; report peak child RSS
    // so a memory-driven kill (an external memory cap or the kernel OOM
    // killer, SIG9) is distinguishable
    // from an internal Kira/Fermat fault (SIG11) without guessing.
#ifdef __APPLE__
    long peak_rss_mb = usage.ru_maxrss / (1024 * 1024);
#else
    long peak_rss_mb = usage.ru_maxrss / 1024;
#endif

    // Pull the real child failure text out of the captured log so the error
    // report names the actual C++ exception / Kira message instead of leaving
    // the caller to guess from a bare signal number.
    std::string log_tail = read_log_tail(log_path);
    bool bad_alloc = log_shows_bad_alloc(log_tail);
    auto append_tail = [&](std::string msg) -> std::string {
        if (!log_tail.empty()) {
            msg += "\n  --- child stderr/stdout tail (" + log_path + ") ---\n";
            msg += log_tail;
            msg += "\n  --- end child tail ---";
        } else if (!log_path.empty()) {
            msg += " (child log " + log_path + " was empty/unreadable)";
        }
        return msg;
    };

    if (WIFSIGNALED(status)) {
        int sig = WTERMSIG(status);
        const char* signame = strsignal(sig);
        std::string msg =
            "kira_run: kira KILLED BY SIGNAL " + std::to_string(sig) + " (" +
            std::string(signame ? signame : "?") + "); peak child RSS = " +
            std::to_string(peak_rss_mb) + " MB. [SIG9=external SIGKILL "
            "(an external memory cap, e.g. cgroup/ulimit) or kernel OOM; SIG11=internal Kira/Fermat "
            "segfault; SIG6=abort -- EITHER a self-capped std::bad_alloc "
            "(OOM under RLIMIT_AS / AMFLOW_KIRA_MEM_CAP_GB, raised as "
            "std::terminate->abort) OR a Kira/Fermat internal assert -- "
            "DISAMBIGUATE via the child tail below]";
        if (sig == SIGABRT && bad_alloc) {
            msg +=
                " >>> child tail shows an out-of-memory exception: this SIG6 "
                "is a SELF-CAPPED std::bad_alloc (RLIMIT_AS / "
                "AMFLOW_KIRA_MEM_CAP_GB), NOT a Kira assert -- raise the cap "
                "or use the FireFly finite-field path (AMFLOW_RUN_FIREFLY=1).";
        }
        throw std::runtime_error(append_tail(std::move(msg)));
    }
    if (!WIFEXITED(status) || WEXITSTATUS(status) != 0) {
        std::string msg =
            "kira_run: kira exited with non-zero status (rc=" +
            std::to_string(WEXITSTATUS(status)) + "); peak child RSS = " +
            std::to_string(peak_rss_mb) + " MB";
        if (bad_alloc) {
            msg +=
                " [child tail shows an out-of-memory exception -- likely "
                "OOM under RLIMIT_AS / AMFLOW_KIRA_MEM_CAP_GB]";
        }
        throw std::runtime_error(append_tail(std::move(msg)));
    }
    // Stage-time instrumentation: emit a marker line for external timing
    // analysis (sum these to get total Kira wallclock; subtract from
    // amflow_cli's total to get non-Kira amflow time).  Enabled when
    // AMFLOW_STAGE_TIMING env var is set to any non-empty value.
    if (const char* env = std::getenv("AMFLOW_STAGE_TIMING");
        env && *env) {
        std::cerr << "[kira_time] " << elapsed << " s (dir=" << dir << ")\n";
    }
    return elapsed;
}

namespace {

void skip_ws(const std::string& s, std::size_t& i) {
    while (i < s.size() && std::isspace((unsigned char)s[i])) ++i;
}

std::string parse_ident(const std::string& s, std::size_t& i) {
    std::string out;
    while (i < s.size() &&
            (std::isalnum((unsigned char)s[i]) || s[i] == '_')) {
        out += s[i];
        ++i;
    }
    return out;
}

std::optional<qft::JIntegral>
parse_jintegral_inplace(const std::string& s, std::size_t& i) {
    skip_ws(s, i);
    std::string fam = parse_ident(s, i);
    if (fam.empty()) return std::nullopt;
    skip_ws(s, i);
    if (i >= s.size() || (s[i] != '(' && s[i] != '[')) return std::nullopt;
    char closer = (s[i] == '(') ? ')' : ']';
    ++i;
    std::vector<long> indices;
    while (true) {
        skip_ws(s, i);
        long sign = 1;
        if (i < s.size() && s[i] == '-') { sign = -1; ++i; }
        else if (i < s.size() && s[i] == '+') { ++i; }
        std::string num;
        while (i < s.size() && std::isdigit((unsigned char)s[i])) {
            num += s[i]; ++i;
        }
        if (num.empty()) return std::nullopt;
        indices.push_back(sign * std::stol(num));
        skip_ws(s, i);
        if (i < s.size() && s[i] == ',') {
            ++i; continue;
        }
        if (i < s.size() && s[i] == closer) {
            ++i; break;
        }
        return std::nullopt;
    }
    return qft::JIntegral(fam, indices);
}

}  // namespace

std::vector<qft::JIntegral>
kira_read_masters(const qft::FamilyConfig& fc, const std::string& dir) {
    std::string path = dir + "/results/" + fc.family + "/masters";
    if (!fs::exists(path)) {
        return {};
    }
    std::ifstream f(path);
    if (!f) throw std::runtime_error("kira_read_masters: cannot open " + path);

    std::vector<qft::JIntegral> out;
    std::string line;
    while (std::getline(f, line)) {
        std::size_t i = 0;
        skip_ws(line, i);
        if (i >= line.size()) continue;
        if (line[i] == '#') continue;
        auto j = parse_jintegral_inplace(line, i);
        if (!j.has_value()) {
            throw std::runtime_error(
                "kira_read_masters: cannot parse line: " + line);
        }
        out.push_back(std::move(*j));
    }
    return out;
}

namespace {

std::string slurp(const std::string& path) {
    std::ifstream f(path, std::ios::binary);
    if (!f) throw std::runtime_error("slurp: cannot open " + path);
    std::ostringstream ss;
    ss << f.rdbuf();
    return ss.str();
}

std::size_t find_matching_close(const std::string& s, std::size_t i,
                                  char open, char close) {
    long depth = 0;
    for (; i < s.size(); ++i) {
        if (s[i] == open) ++depth;
        else if (s[i] == close) {
            if (--depth == 0) return i;
        }
    }
    throw std::runtime_error("find_matching_close: unmatched");
}

std::vector<std::string>
split_top_level(const std::string& s, char sep) {
    std::vector<std::string> out;
    long depth = 0;
    std::size_t start = 0;
    for (std::size_t i = 0; i < s.size(); ++i) {
        char c = s[i];
        if (c == '(' || c == '{' || c == '[') ++depth;
        else if (c == ')' || c == '}' || c == ']') --depth;
        else if (depth == 0 && c == sep) {
            out.push_back(s.substr(start, i - start));
            start = i + 1;
        }
    }
    out.push_back(s.substr(start));
    return out;
}

std::string trim(const std::string& s) {
    std::size_t a = 0;
    while (a < s.size() && std::isspace((unsigned char)s[a])) ++a;
    std::size_t b = s.size();
    while (b > a && std::isspace((unsigned char)s[b - 1])) --b;
    return s.substr(a, b - a);
}

std::size_t find_top_arrow(const std::string& s) {
    long depth = 0;
    for (std::size_t i = 0; i + 1 < s.size(); ++i) {
        char c = s[i];
        if (c == '(' || c == '{' || c == '[') ++depth;
        else if (c == ')' || c == '}' || c == ']') --depth;
        else if (depth == 0 && c == '-' && s[i + 1] == '>') {
            return i;
        }
    }
    return std::string::npos;
}

std::vector<std::string> split_terms(const std::string& s) {
    std::vector<std::string> out;
    long depth = 0;
    std::size_t start = 0;
    bool first = true;
    for (std::size_t i = 0; i < s.size(); ++i) {
        char c = s[i];
        if (c == '(' || c == '{' || c == '[') ++depth;
        else if (c == ')' || c == '}' || c == ']') --depth;
        else if (depth == 0 && (c == '+' || c == '-')) {
            if (i == 0 || first) {
                first = false;
                continue;
            }
            std::size_t k = i;
            while (k > 0 && std::isspace((unsigned char)s[k - 1])) --k;
            if (k == 0) continue;
            char prev = s[k - 1];
            if (prev == '*' || prev == '/' || prev == '^' || prev == '(' ||
                prev == ',' || prev == '=' || prev == '+' || prev == '-') {
                continue;
            }
            out.push_back(s.substr(start, i - start));
            start = i;
            first = false;
        } else if (!std::isspace((unsigned char)c)) {
            first = false;
        }
    }
    out.push_back(s.substr(start));
    return out;
}

struct TermParts {
    std::string coef;
    std::optional<qft::JIntegral> j;
};

TermParts split_term_into_coef_and_jintegral(
    const std::string& term, const std::string& family) {
    TermParts r;
    std::size_t pos = std::string::npos;
    char open_char = '(';
    char close_char = ')';
    {
        std::size_t p1 = term.find(family + "(");
        std::size_t p2 = term.find(family + "[");
        std::size_t cand;
        if (p1 == std::string::npos && p2 == std::string::npos) {
            r.coef = trim(term);
            return r;
        }
        if (p1 == std::string::npos) cand = p2;
        else if (p2 == std::string::npos) cand = p1;
        else cand = std::min(p1, p2);
        pos = cand;
        if (pos < term.size() && term[pos + family.size()] == '[') {
            open_char = '['; close_char = ']';
        }
        if (pos > 0) {
            char p = term[pos - 1];
            if (std::isalnum((unsigned char)p) || p == '_') {
                r.coef = trim(term);
                return r;
            }
        }
    }
    std::size_t open = pos + family.size();
    std::size_t close = find_matching_close(term, open, open_char, close_char);
    std::size_t i = pos;
    auto j = parse_jintegral_inplace(term, i);
    if (!j.has_value() || i != close + 1) {
        r.coef = trim(term);
        return r;
    }
    r.j = std::move(*j);

    auto strip_mul_space = [](std::string& s, bool from_back) {
        if (from_back) {
            while (!s.empty() && (s.back() == '*' ||
                    std::isspace((unsigned char)s.back()))) s.pop_back();
        } else {
            while (!s.empty() && (s.front() == '*' ||
                    std::isspace((unsigned char)s.front())))
                s.erase(s.begin());
        }
    };
    std::string prefix = trim(term.substr(0, pos));
    std::string suffix = trim(term.substr(close + 1));
    strip_mul_space(prefix, true);
    strip_mul_space(suffix, false);
    if (!prefix.empty() && prefix.front() == '+') {
        prefix = trim(prefix.substr(1));
    }

    std::string coef;
    if (prefix.empty() && suffix.empty()) {
        coef = "1";
    } else if (prefix.empty()) {
        coef = suffix;
    } else if (suffix.empty()) {
        coef = prefix;
    } else {
        coef = "(" + prefix + ")*(" + suffix + ")";
    }
    if (coef == "-") coef = "-1";
    r.coef = std::move(coef);
    return r;
}

}  // namespace

std::vector<KiraReductionRule>
kira_read_target_table(const qft::FamilyConfig& fc, const std::string& dir) {
    std::string path = dir + "/results/" + fc.family + "/kira_target.m";
    if (!fs::exists(path)) return {};
    std::string content = slurp(path);

    std::size_t open = content.find('{');
    if (open == std::string::npos) return {};
    std::size_t close = find_matching_close(content, open, '{', '}');
    std::string body = content.substr(open + 1, close - open - 1);

    auto rules_strs = split_top_level(body, ',');

    std::vector<KiraReductionRule> out;
    for (auto& rs : rules_strs) {
        std::string r = trim(rs);
        if (r.empty()) continue;
        std::size_t arrow = find_top_arrow(r);
        if (arrow == std::string::npos) continue;
        std::string lhs_str = trim(r.substr(0, arrow));
        std::string rhs_str = trim(r.substr(arrow + 2));

        std::size_t k = 0;
        auto lhs = parse_jintegral_inplace(lhs_str, k);
        if (!lhs.has_value()) {
            throw std::runtime_error(
                "kira_read_target_table: cannot parse LHS '" + lhs_str + "'");
        }

        KiraReductionRule kr;
        kr.lhs = std::move(*lhs);

        auto terms = split_terms(rhs_str);
        for (auto& t : terms) {
            std::string ts = trim(t);
            if (ts.empty()) continue;
            std::string ts_no_sign = ts;
            if (!ts_no_sign.empty() && (ts_no_sign[0] == '+' || ts_no_sign[0] == '-')) {
                ts_no_sign = trim(ts_no_sign.substr(1));
            }
            if (ts_no_sign == "0") continue;
            TermParts tp = split_term_into_coef_and_jintegral(ts, fc.family);
            if (!tp.j.has_value()) {
                throw std::runtime_error(
                    "kira_read_target_table: term has no J integral: '"
                    + ts + "'");
            }
            kr.rhs.emplace_back(std::move(tp.coef), std::move(*tp.j));
        }
        out.push_back(std::move(kr));
    }
    return out;
}

}  // namespace amflow::ibp
