// SPDX-License-Identifier: MIT
// ibp::ibp_cache — implementation.  See header for the design.

#include "amflow/ibp/ibp_cache.hpp"

#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <system_error>
#include <vector>

#include <cerrno>
#include <chrono>
#include <ctime>

#include <fcntl.h>
#include <sys/file.h>
#include <unistd.h>

namespace amflow::ibp {

namespace fs = std::filesystem;

namespace {

// Canonical list of Kira INPUT files inside a working directory.  These
// are everything `kira_write_*` emits before `kira_run`; the cache key
// is the FNV-1a hash over their concatenated (name, contents) pairs.
// Any file not present is skipped (its name is still NOT folded in, so
// e.g. a Masters-mode dir with no `target` keys differently from a
// Reduce-mode dir only via the divergent `jobs.yaml`, which is
// sufficient — Masters and Reduce always emit different jobs.yaml).
constexpr const char* kInputFiles[] = {
    "config/integralfamilies.yaml",
    "config/kinematics.yaml",
    "jobs.yaml",
    "jobs_kira2math.yaml",
    "preferred",
    "target",
};

// Kira OUTPUT files we cache.  `masters` is the master-integral list
// (read by `kira_read_masters`); `kira_target.m` is the Mathematica-
// format reduction table (read by `kira_read_target_table`).  Both
// live under `results/<family>/`.  We do NOT cache `kira.db` — it is
// large, family-temperature-dependent, and only ever read by the
// in-dir `jobs_kira2math.yaml` export, which the caching wrapper
// arranges to run BEFORE storing so the .m is always present.
constexpr const char* kOutputFiles[] = {
    "masters",
    "kira_target.m",
};

inline void fnv1a(std::uint64_t& h, const void* data, std::size_t n) {
    const auto* p = static_cast<const unsigned char*>(data);
    for (std::size_t i = 0; i < n; ++i) {
        h ^= static_cast<std::uint64_t>(p[i]);
        h *= 0x100000001b3ULL;
    }
}

inline void fnv1a(std::uint64_t& h, const std::string& s) {
    fnv1a(h, s.data(), s.size());
}

std::string hex16(std::uint64_t v) {
    static const char* d = "0123456789abcdef";
    std::string out(16, '0');
    for (int i = 15; i >= 0; --i) {
        out[(std::size_t)i] = d[v & 0xfULL];
        v >>= 4;
    }
    return out;
}

std::string dir_tag(const std::string& dir) {
    // Short human-readable tag for log lines: last path component.
    auto p = fs::path(dir);
    return p.filename().empty() ? p.parent_path().filename().string()
                                : p.filename().string();
}

}  // namespace

// Canonicalize jobs.yaml content for hashing.  The cache key must be
// invariant under env knobs that change the SEARCH BOUND or PARALLELISM
// but not the RESULT (the `results/<fam>/{masters,kira_target.m}` files
// the cache stores).  The known pitfalls:
//
//   * `AMFLOW_REDUCE_DCAP=1` appends `, d: <ibp_dot>` to the `reduce:`
//     seed-spec flow-mapping line (kira_yaml.cpp).  The d-cap shrinks the
//     IBP seed-gen box but `select_mandatory_list: [fam, target]` still
//     determines the actual reduced set, so the kira_target.m output is
//     identical.  A relaunch that forgets DCAP=1 hashes to a different
//     key → ~2h cold kira instead of seconds (measured on a
//     14-job cold-start farm).
//
//   * `n_thread` / `threads:` — Kira parallelism.  Currently passed via
//     `-p<n>` on the command line (kira_run.cpp), NOT written into
//     jobs.yaml, so this is a no-op today; kept as a defensive strip in
//     case a future yaml writer change emits it.
//
// We strip `, d: <int>` ONLY from seed-spec lines (the `{topologies: …,
// sectors: …, r: …, s: …}` flow mappings) that appear BEFORE the
// `select_integrals:` block — i.e. the `reduce:` seeding line, where the
// d-cap is a search-bound knob.  The `d:` inside
// `select_mandatory_recursively:` (Masters mode, AFTER `select_integrals:`)
// is RESULT-determining (it bounds which masters are enumerated) and is
// left untouched.  Any whole line of the form `<ws>n_thread[s]: <int>` or
// `<ws>threads: <int>` is dropped.
std::string canonicalize_jobs_yaml(const std::string& raw) {
    std::string out;
    out.reserve(raw.size());
    bool past_select_integrals = false;
    std::size_t i = 0;
    while (i < raw.size()) {
        std::size_t eol = raw.find('\n', i);
        std::size_t end = (eol == std::string::npos) ? raw.size() : eol;
        std::string line = raw.substr(i, end - i);

        // Track whether we've entered the select_integrals block.
        if (line.find("select_integrals:") != std::string::npos) {
            past_select_integrals = true;
        }

        // Drop n_thread/threads lines wholesale (parallelism knob).
        {
            std::size_t p = line.find_first_not_of(" \t");
            if (p != std::string::npos) {
                auto starts = [&](const char* k) {
                    return line.compare(p, std::char_traits<char>::length(k), k) == 0;
                };
                if (starts("n_thread:") || starts("n_threads:") ||
                    starts("threads:")) {
                    i = (eol == std::string::npos) ? raw.size() : eol + 1;
                    continue;  // skip line entirely (incl. newline)
                }
            }
        }

        // Strip `, d: <int>` from the reduce: seed-spec flow mapping
        // (before select_integrals only — see header comment).
        if (!past_select_integrals &&
            line.find("topologies:") != std::string::npos &&
            line.find("sectors:") != std::string::npos) {
            std::size_t dp = line.find(", d: ");
            if (dp != std::string::npos) {
                std::size_t q = dp + 5;
                while (q < line.size() &&
                       (line[q] == '-' || (line[q] >= '0' && line[q] <= '9')))
                    ++q;
                line.erase(dp, q - dp);
            }
        }

        out.append(line);
        if (eol != std::string::npos) out.push_back('\n');
        i = (eol == std::string::npos) ? raw.size() : eol + 1;
    }
    return out;
}

std::string ibp_cache_resolve_root(const std::string& opt_root) {
    if (!opt_root.empty()) return opt_root;
    if (const char* e = std::getenv("AMFLOW_IBP_CACHE"); e && *e) return e;
    return {};
}

std::string ibp_cache_key(const std::string& dir,
                          const std::string& extra_salt) {
    std::uint64_t h = 0xcbf29ce484222325ULL;  // FNV offset basis
    for (const char* rel : kInputFiles) {
        fs::path p = fs::path(dir) / rel;
        std::error_code ec;
        if (!fs::exists(p, ec)) continue;
        // Fold the relative name in so two dirs that happen to share
        // identical bytes across DIFFERENTLY-named inputs still key
        // distinctly (defensive; in practice jobs.yaml always differs).
        fnv1a(h, std::string("\0@", 2));
        fnv1a(h, std::string(rel));
        fnv1a(h, std::string("\0=", 2));
        std::ifstream f(p, std::ios::binary);
        if (!f) continue;
        // jobs.yaml gets canonicalized (DCAP / n_thread stripped) so env
        // knobs that don't change Kira's RESULT don't perturb the key.
        // All other inputs are hashed byte-exact.  See
        // `canonicalize_jobs_yaml` for the rationale and the precise
        // strip rules.  jobs_kira2math.yaml is NOT canonicalized — it
        // never carries d:/threads.
        if (std::string(rel) == "jobs.yaml") {
            std::ostringstream ss;
            ss << f.rdbuf();
            fnv1a(h, canonicalize_jobs_yaml(ss.str()));
            continue;
        }
        char buf[8192];
        while (f) {
            f.read(buf, sizeof buf);
            fnv1a(h, buf, static_cast<std::size_t>(f.gcount()));
        }
    }
    // Result-determining parameters that are NOT captured by the input
    // files (currently only the AMFLOW_DIFFEQ_NUMD numeric-d
    // substitution — see the header comment).  Empty salt (the default
    // and every pre-existing call site) folds NOTHING in, so legacy
    // keys are byte-identical.
    if (!extra_salt.empty()) {
        fnv1a(h, std::string("\0@salt\0=", 8));
        fnv1a(h, extra_salt);
    }
    return hex16(h);
}

bool ibp_cache_lookup(const std::string& cache_root,
                      const std::string& key,
                      const std::string& family,
                      const std::string& dir) {
    if (cache_root.empty()) return false;
    fs::path src = fs::path(cache_root) / "ibp" / key / "results" / family;
    std::error_code ec;
    if (!fs::is_directory(src, ec)) return false;

    // A HIT requires at least one non-empty cached output; we restore
    // every one that exists.  An entry written by a Masters-mode store
    // has only `masters`; a Reduce-mode store has both.
    bool any = false;
    fs::path dst = fs::path(dir) / "results" / family;
    for (const char* rel : kOutputFiles) {
        fs::path sp = src / rel;
        if (!fs::exists(sp, ec) || fs::file_size(sp, ec) == 0) continue;
        if (!any) fs::create_directories(dst, ec);
        fs::copy_file(sp, dst / rel,
                      fs::copy_options::overwrite_existing, ec);
        if (ec) {
            std::cerr << "[amflow][ibp-cache] WARN copy failed (" << rel
                      << "): " << ec.message() << "\n";
            return false;
        }
        any = true;
    }
    if (any) {
        std::cerr << "[amflow][ibp-cache] HIT  key=" << key
                  << " fam=" << family
                  << " tag=" << dir_tag(dir)
                  << " (skipping Kira)\n";
    }
    return any;
}

void ibp_cache_store(const std::string& cache_root,
                     const std::string& key,
                     const std::string& family,
                     const std::string& dir) {
    if (cache_root.empty()) return;
    fs::path src = fs::path(dir) / "results" / family;
    std::error_code ec;
    if (!fs::is_directory(src, ec)) {
        std::cerr << "[amflow][ibp-cache] WARN nothing to store for key="
                  << key << " (no results/" << family << ")\n";
        return;
    }

    // Write-then-rename: build under a pid-unique temp sibling so a
    // concurrent lookup on the same key never sees a half-written
    // directory.  If the final rename fails because another process
    // beat us to it, that's fine — the entry is content-addressed.
    fs::path final_dir = fs::path(cache_root) / "ibp" / key;
    fs::path tmp_dir   = fs::path(cache_root) / "ibp"
                       / (key + ".tmp." + std::to_string(::getpid()));
    fs::create_directories(tmp_dir / "results" / family, ec);
    if (ec) {
        std::cerr << "[amflow][ibp-cache] WARN mkdir failed: "
                  << ec.message() << "\n";
        return;
    }

    bool any = false;
    for (const char* rel : kOutputFiles) {
        fs::path sp = src / rel;
        if (!fs::exists(sp, ec) || fs::file_size(sp, ec) == 0) continue;
        fs::copy_file(sp, tmp_dir / "results" / family / rel,
                      fs::copy_options::overwrite_existing, ec);
        if (ec) {
            std::cerr << "[amflow][ibp-cache] WARN copy-to-cache failed ("
                      << rel << "): " << ec.message() << "\n";
            fs::remove_all(tmp_dir, ec);
            return;
        }
        any = true;
    }
    if (!any) {
        fs::remove_all(tmp_dir, ec);
        return;
    }

    // Atomic publish.  If `final_dir` already exists (raced), drop our
    // temp — the existing entry is by construction identical.
    fs::rename(tmp_dir, final_dir, ec);
    if (ec) {
        if (fs::exists(final_dir)) {
            fs::remove_all(tmp_dir, ec);
        } else {
            std::cerr << "[amflow][ibp-cache] WARN rename failed: "
                      << ec.message() << "\n";
            fs::remove_all(tmp_dir, ec);
            return;
        }
    }
    std::cerr << "[amflow][ibp-cache] STORED key=" << key
              << " fam=" << family
              << " tag=" << dir_tag(dir) << "\n";
}

// ===========================================================================
//  IbpCacheLock — stampede lock
// ===========================================================================

bool IbpCacheLock::enabled() {
    // Default ON.  AMFLOW_IBP_CACHE_LOCK=0 disables.
    if (const char* v = std::getenv("AMFLOW_IBP_CACHE_LOCK"); v && *v) {
        return !(v[0] == '0' && v[1] == '\0');
    }
    return true;
}

namespace {
long ibp_lock_timeout_s() {
    if (const char* v = std::getenv("AMFLOW_IBP_LOCK_TIMEOUT_S"); v && *v) {
        char* end = nullptr;
        long t = std::strtol(v, &end, 10);
        if (end != v && t > 0) return t;
    }
    return 7200L;  // 2 h default — bounds the worst Kira reduce we cache.
}
}  // namespace

IbpCacheLock::IbpCacheLock(const std::string& cache_root,
                           const std::string& key) {
    if (cache_root.empty() || !enabled()) return;

    // Ensure the lock directory exists (the cache entry dir may not yet).
    std::error_code ec;
    fs::path lock_dir = fs::path(cache_root) / "ibp";
    fs::create_directories(lock_dir, ec);
    fs::path lock_path = lock_dir / (key + ".lock");

    fd_ = ::open(lock_path.c_str(), O_RDWR | O_CREAT | O_CLOEXEC, 0644);
    if (fd_ < 0) {
        std::cerr << "[amflow][ibp-cache] WARN lock open failed for key="
                  << key << " (errno=" << errno << "); proceeding unlocked\n";
        return;
    }

    // Non-blocking probe first so the log distinguishes "first job" from
    // "waited on another job".
    if (::flock(fd_, LOCK_EX | LOCK_NB) == 0) {
        held_ = true;
        return;
    }
    if (errno != EWOULDBLOCK) {
        std::cerr << "[amflow][ibp-cache] WARN flock(NB) failed for key="
                  << key << " (errno=" << errno << "); proceeding unlocked\n";
        ::close(fd_);
        fd_ = -1;
        return;
    }

    // Another job holds it.  Log WAIT and poll LOCK_NB until acquired or
    // timeout.  Polling (200 ms) rather than a blocking flock() so a
    // crashed holder (kernel releases the lock on process death, but a
    // SIGSTOP'd or wedged holder does not) cannot deadlock the farm
    // beyond the timeout.
    waited_ = true;
    std::cerr << "[amflow][ibp-cache] WAIT key=" << key
              << " (another job computing)\n";

    const long timeout_s = ibp_lock_timeout_s();
    const auto t0 = std::chrono::steady_clock::now();
    struct timespec slp { 0, 200L * 1000 * 1000 };  // 200 ms
    for (;;) {
        ::nanosleep(&slp, nullptr);
        if (::flock(fd_, LOCK_EX | LOCK_NB) == 0) {
            held_ = true;
            return;
        }
        if (errno != EWOULDBLOCK) {
            std::cerr << "[amflow][ibp-cache] WARN flock(NB) failed for key="
                      << key << " (errno=" << errno
                      << "); proceeding unlocked\n";
            ::close(fd_);
            fd_ = -1;
            return;
        }
        auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(
                           std::chrono::steady_clock::now() - t0).count();
        if (elapsed >= timeout_s) {
            std::cerr << "[amflow][ibp-cache] WARN lock wait timed out after "
                      << elapsed << "s for key=" << key
                      << "; proceeding anyway (crashed holder?)\n";
            // Leave fd_ open (so a later store still has a handle to
            // release-noop on), but held_ stays false.
            return;
        }
    }
}

IbpCacheLock::~IbpCacheLock() {
    if (fd_ >= 0) {
        if (held_) ::flock(fd_, LOCK_UN);
        ::close(fd_);
    }
}

}  // namespace amflow::ibp
