// SPDX-License-Identifier: MIT
// ibp::diffeq_numd_* — AMFLOW_DIFFEQ_NUMD numeric-d substitution inside
// the diffeq target_reduce leg (the etanumd_convert.sh sed INTERNALIZED).
//
// See KiraConfig::diffeq_numd (kira.hpp) for the full rationale.  This
// file implements the two-stage pre-FireFly hook:
//
//   1. `jobs_initiate.yaml` — a copy of the already-written jobs.yaml
//      with every solve step (run_firefly / run_triangular /
//      run_back_substitution / iterative_reduction) and any kira2math
//      block stripped, so the kira invocation ONLY builds the IBP
//      system (tmp/<fam>/SYSTEM_*.gz + sectormappings + masters).
//   2. token-safe in-place substitution d -> (P/Q) in every
//      tmp/*/SYSTEM_*.gz, porting the proven etanumd_convert.sh checks:
//        * identifier alphabet must be ⊆ {Eq, d, eta} (any other
//          symbol in a SYSTEM file means this is NOT the aux-mass DE
//          system shape the substitution was proven on — abort);
//        * word-boundary replacement with a parenthesized value;
//        * zero residual standalone `d` tokens verified by re-reading
//          every rewritten file;
//        * std::stol range guard on P and Q (Kira parses coefficient
//          tokens with std::stol);
//      Idempotent: files with zero `d` tokens (ffsave resume of an
//      already-substituted dir) are left byte-untouched.
//
// The caller (`kira_run_cached`, reduce.cpp) then runs the normal
// jobs.yaml in the SAME dir; kira resumes past run_initiate (the
// documented same-dir resume behavior — the exact mechanism the external
// etanumd rerun used) and FireFly reconstructs 1-var (eta).

#include "amflow/ibp/kira.hpp"

#include <cctype>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <zlib.h>

namespace amflow::ibp {

namespace fs = std::filesystem;

namespace {

bool is_ident_start(char c) {
    return std::isalpha(static_cast<unsigned char>(c)) || c == '_';
}
bool is_ident_char(char c) {
    return std::isalnum(static_cast<unsigned char>(c)) || c == '_';
}

// Validate "P" or "P/Q" with the etanumd std::stol guard.  Throws on
// anything else.  Returns the canonical value string (as given).
std::string validate_numd(const std::string& v) {
    const std::string::size_type slash = v.find('/');
    const std::string p_str = (slash == std::string::npos) ? v : v.substr(0, slash);
    const std::string q_str = (slash == std::string::npos) ? std::string("1")
                                                           : v.substr(slash + 1);
    auto is_int = [](const std::string& s, bool allow_sign) {
        if (s.empty()) return false;
        std::size_t i = 0;
        if (allow_sign && (s[0] == '-' || s[0] == '+')) i = 1;
        if (i >= s.size()) return false;
        for (; i < s.size(); ++i) {
            if (!std::isdigit(static_cast<unsigned char>(s[i]))) return false;
        }
        return true;
    };
    if (!is_int(p_str, /*allow_sign=*/true) || !is_int(q_str, /*allow_sign=*/false)) {
        throw std::invalid_argument(
            "AMFLOW_DIFFEQ_NUMD: value must be an integer or rational "
            "\"P/Q\" (got '" + v + "')");
    }
    // std::stol guard: Kira parses coefficient tokens with std::stol, so
    // each of |P|, Q must stay below ~9.2e18 (the etanumd guard).
    try {
        const long long p = std::stoll(p_str);
        const long long q = std::stoll(q_str);
        (void)p;
        if (q <= 0) {
            throw std::invalid_argument("denominator must be positive");
        }
    } catch (const std::exception&) {
        throw std::invalid_argument(
            "AMFLOW_DIFFEQ_NUMD: P and Q must each fit in a signed 64-bit "
            "integer (std::stol guard) with Q > 0 (got '" + v + "')");
    }
    return v;
}

// gz streaming helpers -------------------------------------------------------

struct GzIn {
    gzFile f = nullptr;
    explicit GzIn(const std::string& path) {
        f = gzopen(path.c_str(), "rb");
        if (!f) throw std::runtime_error("diffeq_numd: cannot open " + path);
    }
    ~GzIn() { if (f) gzclose(f); }
    // Returns bytes read; 0 at EOF; throws on error.
    int read(char* buf, unsigned n) {
        int r = gzread(f, buf, n);
        if (r < 0) throw std::runtime_error("diffeq_numd: gzread failed");
        return r;
    }
};

struct GzOut {
    gzFile f = nullptr;
    explicit GzOut(const std::string& path) {
        // Level 1 mirrors etanumd_convert.sh's `gzip -1` (Kira re-reads
        // these once; compression speed dominates).
        f = gzopen(path.c_str(), "wb1");
        if (!f) throw std::runtime_error("diffeq_numd: cannot open " + path);
    }
    ~GzOut() { if (f) gzclose(f); }
    void write(const std::string& s) {
        if (s.empty()) return;
        if (gzwrite(f, s.data(), static_cast<unsigned>(s.size())) !=
                static_cast<int>(s.size())) {
            throw std::runtime_error("diffeq_numd: gzwrite failed");
        }
    }
    void close() {
        if (f) { gzclose(f); f = nullptr; }
    }
};

// One streaming pass over a SYSTEM_*.gz: verify the identifier alphabet
// is ⊆ {Eq, d, eta}, count standalone `d` tokens, and (when `out` is
// non-null) emit the content with every standalone `d` replaced by
// `repl`.  Returns the number of `d` tokens seen.
std::size_t scan_system_gz(const std::string& path,
                           const std::string& repl,
                           GzOut* out) {
    GzIn in(path);
    constexpr unsigned kChunk = 1u << 22;  // 4 MiB
    std::vector<char> buf(kChunk);
    std::string carry;        // trailing partial identifier from last chunk
    std::string outbuf;
    outbuf.reserve(kChunk + (kChunk >> 2));
    std::size_t n_d = 0;

    auto emit_ident = [&](const std::string& ident) {
        if (ident == "d") {
            ++n_d;
            if (out) outbuf.append(repl);
            return;
        }
        if (ident != "Eq" && ident != "eta") {
            throw std::runtime_error(
                "diffeq_numd: SYSTEM token alphabet violation in " + path +
                ": identifier '" + ident + "' not in {Eq, d, eta} — "
                "refusing to substitute (etanumd alphabet guard)");
        }
        if (out) outbuf.append(ident);
    };

    bool eof = false;
    while (!eof) {
        int r = in.read(buf.data(), kChunk);
        if (r == 0) eof = true;
        std::string data;
        data.reserve(carry.size() + static_cast<std::size_t>(r));
        data.append(carry);
        carry.clear();
        data.append(buf.data(), static_cast<std::size_t>(r));

        std::size_t i = 0;
        const std::size_t n = data.size();
        while (i < n) {
            char c = data[i];
            if (is_ident_start(c)) {
                std::size_t j = i + 1;
                while (j < n && is_ident_char(data[j])) ++j;
                if (j == n && !eof) {
                    // Identifier may continue in the next chunk — carry it.
                    carry.assign(data, i, n - i);
                    break;
                }
                emit_ident(data.substr(i, j - i));
                i = j;
            } else {
                if (out) outbuf.push_back(c);
                ++i;
            }
        }
        if (out && !outbuf.empty()) {
            out->write(outbuf);
            outbuf.clear();
        }
    }
    if (!carry.empty()) {
        // File ended mid-identifier (eof handling above flushes it as a
        // complete token on the final loop pass; reaching here means the
        // final pass still carried — process it now).
        emit_ident(carry);
        if (out && !outbuf.empty()) out->write(outbuf);
    }
    return n_d;
}

}  // namespace

std::string diffeq_numd_env() {
    const char* v = std::getenv("AMFLOW_DIFFEQ_NUMD");
    if (!v || !*v) return {};
    return validate_numd(v);
}

void diffeq_numd_pre_stage(const KiraConfig& cfg,
                           const std::string& dir,
                           const std::string& log_path) {
    const std::string& numd = cfg.diffeq_numd;
    if (numd.empty()) return;  // defensive; caller gates on this too
    const std::string repl = "(" + numd + ")";

    // --- 1. jobs_initiate.yaml: strip solve steps from jobs.yaml -----------
    const fs::path jobs = fs::path(dir) / "jobs.yaml";
    std::ifstream jf(jobs);
    if (!jf) {
        throw std::runtime_error(
            "diffeq_numd: cannot open " + jobs.string());
    }
    std::ostringstream stripped;
    std::string line;
    bool in_kira2math = false;
    while (std::getline(jf, line)) {
        const std::size_t p = line.find_first_not_of(" \t");
        const std::string t =
            (p == std::string::npos) ? std::string() : line.substr(p);
        // A kira2math job block ends at the next top-level job (` - `)
        // or EOF.
        if (in_kira2math) {
            if (t.rfind("- ", 0) == 0 && p <= 1 &&
                t.find("kira2math") == std::string::npos) {
                in_kira2math = false;  // fall through: keep this line
            } else {
                continue;
            }
        }
        if (t.rfind("- kira2math:", 0) == 0) { in_kira2math = true; continue; }
        if (t.rfind("run_firefly:", 0) == 0)            continue;
        if (t.rfind("run_triangular:", 0) == 0)         continue;
        if (t.rfind("run_back_substitution:", 0) == 0)  continue;
        if (t.rfind("iterative_reduction:", 0) == 0)    continue;
        stripped << line << "\n";
    }
    jf.close();
    {
        std::ofstream out(fs::path(dir) / "jobs_initiate.yaml");
        if (!out) {
            throw std::runtime_error(
                "diffeq_numd: cannot write jobs_initiate.yaml in " + dir);
        }
        out << stripped.str();
    }

    std::cerr << "[amflow][diffeq-numd] ACTIVE d -> " << repl
              << " in " << dir << " (stage 1: run_initiate only)\n";

    // --- 2. kira run_initiate only -----------------------------------------
    kira_run(cfg, dir, log_path, "jobs_initiate.yaml");

    // --- 3. substitute d -> (P/Q) in tmp/*/SYSTEM_*.gz ----------------------
    const fs::path tmp = fs::path(dir) / "tmp";
    std::error_code ec;
    if (!fs::is_directory(tmp, ec)) {
        throw std::runtime_error(
            "diffeq_numd: kira run_initiate produced no tmp/ in " + dir);
    }
    std::size_t n_files_total = 0;
    std::size_t n_subst_total = 0;
    for (const auto& fam_entry : fs::directory_iterator(tmp)) {
        if (!fam_entry.is_directory()) continue;
        const std::string fam = fam_entry.path().filename().string();
        std::vector<fs::path> sys_files;
        for (const auto& f : fs::directory_iterator(fam_entry.path())) {
            const std::string name = f.path().filename().string();
            if (name.rfind("SYSTEM_", 0) == 0 &&
                name.size() > 3 &&
                name.compare(name.size() - 3, 3, ".gz") == 0) {
                sys_files.push_back(f.path());
            }
        }
        if (sys_files.empty()) {
            throw std::runtime_error(
                "diffeq_numd: no SYSTEM_*.gz under " +
                fam_entry.path().string() +
                " after run_initiate (unexpected dir shape)");
        }
        // Pass A over ALL files first: alphabet check + d count, before
        // ANY rewrite — an alphabet violation in file k must not leave
        // files 1..k-1 substituted (partial-state pitfall).
        std::vector<std::size_t> d_counts;
        d_counts.reserve(sys_files.size());
        for (const auto& sp : sys_files) {
            d_counts.push_back(scan_system_gz(sp.string(), repl, nullptr));
        }
        std::size_t n_subst = 0;
        for (std::size_t fi = 0; fi < sys_files.size(); ++fi) {
            const auto& sp = sys_files[fi];
            if (d_counts[fi] == 0) continue;  // already numeric-d (resume)
            // Pass B: rewrite via temp-then-rename so a crash mid-write
            // never leaves a truncated SYSTEM file under the real name.
            const fs::path tmp_out = sp.string() + ".numd_tmp";
            {
                GzOut out(tmp_out.string());
                scan_system_gz(sp.string(), repl, &out);
                out.close();
            }
            // Pass C: verify zero residual standalone `d` tokens.
            const std::size_t residual =
                scan_system_gz(tmp_out.string(), repl, nullptr);
            if (residual != 0) {
                fs::remove(tmp_out, ec);
                throw std::runtime_error(
                    "diffeq_numd: " + std::to_string(residual) +
                    " residual 'd' tokens after substitution in " +
                    sp.string());
            }
            fs::rename(tmp_out, sp);
            ++n_subst;
        }
        // File-count sanity (mirrors etanumd's count check): every file
        // listed above still exists.
        for (const auto& sp : sys_files) {
            if (!fs::exists(sp, ec)) {
                throw std::runtime_error(
                    "diffeq_numd: SYSTEM file vanished during "
                    "substitution: " + sp.string());
            }
        }
        std::cerr << "[amflow][diffeq-numd] fam=" << fam << ": "
                  << n_subst << "/" << sys_files.size()
                  << " SYSTEM files substituted d -> " << repl
                  << (n_subst == 0 ? " (already numeric-d; resume)" : "")
                  << "\n";
        n_files_total += sys_files.size();
        n_subst_total += n_subst;
    }
    std::cerr << "[amflow][diffeq-numd] stage 2 done: " << n_subst_total
              << "/" << n_files_total << " SYSTEM files rewritten in "
              << dir << "; handing back to kira jobs.yaml (FireFly 1-var)\n";
}

}  // namespace amflow::ibp
