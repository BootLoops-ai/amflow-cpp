// SPDX-License-Identifier: MIT
// pipeline::system_checkpoint — implementation.  See header for the design.

#include "amflow/pipeline/system_checkpoint.hpp"

#include "amflow/numeric/options.hpp"
#include "amflow/pipeline/amfsystem.hpp"
#include "amflow/qft/jintegral.hpp"

#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <system_error>
#include <vector>

#include <unistd.h>

#include <flint/acb.h>
#include <flint/arb.h>

namespace amflow::pipeline {

namespace fs = std::filesystem;

namespace {

// Serialization format version.  Bumped if the on-disk layout changes, or
// if the set of inputs folded into the key grows, so a stale file from an
// older binary is rejected as a MISS rather than mis-parsed or replayed.
//   v1: family / preferred / etac / eps grid / prec / direction / ending.
//   v2: + every numeric option and environment override that changes the
//       solve (see numeric_option_fingerprint below).
constexpr int kCkptVersion = 2;

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

// Environment overrides that change the numeric solve.  Their VALUES are
// folded into the key ("" when unset), so a re-run with any of them
// changed, set or cleared re-keys and re-solves instead of replaying a
// checkpoint produced under different numerics.
const char* const kSolveEnvKnobs[] = {
    "AMFLOW_DIFFEQ_NUMD",          // numeric-d substitution in the DE tables
    "AMFLOW_SPARSE_CHOP_DIGITS",   // sparse-layer chop threshold
    "AMFLOW_NO_SPARSE_CHOP",       // sparse-layer chop disabled
    "AMFLOW_RESONANCE_CHOP",       // resonance-detection tolerance
    "AMFLOW_X00_PIVOT_GUARD",      // uncertifiable-pivot guard
    "AMFLOW_FORCE_ENDING_DEPTH",   // trivial ending scheme below a depth
    "AMFLOW_FIXEDPOINT_PROBE",     // fixed-point probe route at the root
    "AMFLOW_ALLOW_FIXEDPOINT_RECURSION",  // recursion guard off
    "AMFLOW_ODE_HALF_COSETS",
    "AMFLOW_BC_PARTIAL_ZEROFILL",
    "AMFLOW_TOFUCHSIAN_KERNEL_RESCUE",
    "AMFLOW_XORDER_TAIL_BOUND",
    "AMFLOW_X00_ONECOL",
    "AMFLOW_X00_SYNTH_BC",
    "AMFLOW_X00_BC_FILE",
    "AMFLOW_ALLOW_CUT_FALSE_ZERO",
    "AMFLOW_EIKONAL_REGIONS",
    "AMFLOW_MERCEDES_CUT",
    "AMFLOW_MERCEDES_CONV_RE",
    "AMFLOW_MERCEDES_CONV_IM",
};

// Fingerprint of everything that selects the numerics of a solve but is
// not part of the system's identity: the global / expansion / running
// option groups (working_pre is folded separately as bits by the caller;
// silent_mode does not change the result) and the env knobs above.
// Order and spelling are fixed — this string is hashed, never parsed.
std::string numeric_option_fingerprint() {
    const auto& g = numeric::global_options();
    const auto& x = numeric::expansion_options();
    const auto& r = numeric::running_options();
    std::ostringstream oss;
    oss << "wp=" << g.working_pre
        << ";cp=" << g.chop_pre
        << ";rp=" << g.rationalize_pre
        << ";d0=" << g.d0
        << ";xo=" << x.x_order
        << ";exo=" << x.extra_x_order
        << ";lxo=" << x.learn_x_order
        << ";txo=" << x.test_x_order
        << ";rr=" << r.run_radius
        << ";rl=" << r.run_length
        << ";rc=" << r.run_candidate
        << ";rd=" << static_cast<int>(r.run_direction)
        << ";cre=" << r.custom_direction_re
        << ";cim=" << r.custom_direction_im;
    for (const char* name : kSolveEnvKnobs) {
        const char* v = std::getenv(name);
        oss << ";" << name << "=" << (v ? v : "");
    }
    return oss.str();
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

// Lossless acb -> "<real_dump>\t<imag_dump>".  arb_dump_str is base-16,
// exact (midpoint + radius); free the C strings it returns.
std::string acb_dump(const numeric::AcbValue& v) {
    char* re = arb_dump_str(acb_realref(v.raw()));
    char* im = arb_dump_str(acb_imagref(v.raw()));
    std::string out = std::string(re ? re : "") + "\t"
                    + std::string(im ? im : "");
    flint_free(re);
    flint_free(im);
    return out;
}

// Inverse of acb_dump.  Returns false on a malformed token.
bool acb_undump(const std::string& line, numeric::AcbValue& out) {
    const auto tab = line.find('\t');
    if (tab == std::string::npos) return false;
    const std::string re = line.substr(0, tab);
    const std::string im = line.substr(tab + 1);
    if (arb_load_str(acb_realref(out.raw()), re.c_str()) != 0) return false;
    if (arb_load_str(acb_imagref(out.raw()), im.c_str()) != 0) return false;
    return true;
}

std::string jint_serialize(const qft::JIntegral& j) {
    std::ostringstream oss;
    oss << j.family();
    for (long idx : j.indices()) oss << "|" << idx;
    return oss.str();
}

qft::JIntegral jint_deserialize(const std::string& s) {
    std::vector<std::string> parts;
    std::string cur;
    std::istringstream iss(s);
    while (std::getline(iss, cur, '|')) parts.push_back(cur);
    if (parts.empty()) return {};
    std::vector<long> idx;
    idx.reserve(parts.size() - 1);
    for (std::size_t i = 1; i < parts.size(); ++i) {
        idx.push_back(std::stol(parts[i]));
    }
    return qft::JIntegral(parts[0], std::move(idx));
}

}  // namespace

std::string system_checkpoint_resolve_root(const std::string& ibp_cache_dir) {
    if (const char* e = std::getenv("AMFLOW_SYSTEM_CHECKPOINT"); e && *e) {
        return e;
    }
    if (!ibp_cache_dir.empty()) {
        return (fs::path(ibp_cache_dir) / "syssolve").string();
    }
    return {};
}

std::string system_checkpoint_key(const std::string& family,
                                  const std::vector<std::string>& preferred_keys,
                                  const std::vector<int>& etac,
                                  const std::vector<numeric::AcbValue>& epslist,
                                  long working_prec_bits,
                                  int direction,
                                  bool is_ending) {
    std::uint64_t h = 0xcbf29ce484222325ULL;  // FNV offset basis
    fnv1a(h, std::string("ckpt-v"));
    fnv1a(h, std::to_string(kCkptVersion));
    fnv1a(h, std::string("\0fam\0", 5));
    fnv1a(h, family);
    fnv1a(h, std::string("\0pref\0", 6));
    for (const auto& k : preferred_keys) {
        fnv1a(h, k);
        fnv1a(h, std::string("\0", 1));
    }
    fnv1a(h, std::string("\0eta\0", 5));
    for (int e : etac) {
        fnv1a(h, std::to_string(e));
        fnv1a(h, std::string(",", 1));
    }
    fnv1a(h, std::string("\0eps\0", 5));
    // The ε-grid (lossless dump) — a different grid is a different solve.
    for (const auto& e : epslist) {
        const std::string d = acb_dump(e);
        fnv1a(h, d);
        fnv1a(h, std::string("\0", 1));
    }
    fnv1a(h, std::string("\0prec\0", 6));
    fnv1a(h, std::to_string(working_prec_bits));
    fnv1a(h, std::string("\0dir\0", 5));
    fnv1a(h, std::to_string(direction));
    fnv1a(h, std::string("\0end\0", 5));
    fnv1a(h, std::to_string(is_ending ? 1 : 0));
    // Numeric options + env overrides: a different rationalize_pre,
    // chop_pre, x_order, numeric-d value, chop override, ... is a
    // different solve and must not replay this key's solution.
    fnv1a(h, std::string("\0opt\0", 5));
    fnv1a(h, numeric_option_fingerprint());
    return hex16(h);
}

// On-disk layout (text, one record):
//   AMFLOW_SYSCKPT v<version>
//   nsol <N>
//   for each of N AMFSystemSolution:
//     sol
//     mv <count>
//       <acb_dump>            (count lines)
//     gp <count>
//       <jint_serialize>      (count lines)
//     gv <count>
//       <acb_dump>            (count lines)
//   end
bool system_checkpoint_load(const std::string& root,
                            const std::string& key,
                            std::vector<AMFSystemSolution>& out) {
    if (root.empty()) return false;
    fs::path path = fs::path(root) / ("system_" + key + ".ckpt");
    std::error_code ec;
    if (!fs::exists(path, ec)) return false;
    std::ifstream f(path);
    if (!f) return false;

    auto fail = [&](const char* /*why*/) -> bool {
        // Treat any malformed content as a clean MISS.
        out.clear();
        return false;
    };

    std::string tag, ver;
    long nsol = 0;
    if (!(f >> tag >> ver)) return fail("hdr");
    if (tag != "AMFLOW_SYSCKPT") return fail("magic");
    if (ver != "v" + std::to_string(kCkptVersion)) return fail("version");
    {
        std::string kw;
        if (!(f >> kw >> nsol) || kw != "nsol" || nsol < 0) return fail("nsol");
    }

    std::vector<AMFSystemSolution> parsed;
    parsed.reserve((std::size_t)nsol);
    // Consume to end of the "nsol N" line before line-oriented reads.
    std::string discard;
    std::getline(f, discard);

    auto read_count = [&](const char* want, long& n) -> bool {
        std::string line;
        if (!std::getline(f, line)) return false;
        std::istringstream iss(line);
        std::string kw;
        if (!(iss >> kw >> n) || kw != want || n < 0) return false;
        return true;
    };

    for (long s = 0; s < nsol; ++s) {
        std::string line;
        if (!std::getline(f, line) || line != "sol") return fail("sol");
        AMFSystemSolution cur;

        long n_mv = 0;
        if (!read_count("mv", n_mv)) return fail("mv");
        for (long i = 0; i < n_mv; ++i) {
            if (!std::getline(f, line)) return fail("mv-line");
            numeric::AcbValue v;
            if (!acb_undump(line, v)) return fail("mv-parse");
            cur.master_values.push_back(std::move(v));
        }

        long n_gp = 0;
        if (!read_count("gp", n_gp)) return fail("gp");
        for (long i = 0; i < n_gp; ++i) {
            if (!std::getline(f, line)) return fail("gp-line");
            cur.global_preferred.push_back(jint_deserialize(line));
        }

        long n_gv = 0;
        if (!read_count("gv", n_gv)) return fail("gv");
        for (long i = 0; i < n_gv; ++i) {
            if (!std::getline(f, line)) return fail("gv-line");
            numeric::AcbValue v;
            if (!acb_undump(line, v)) return fail("gv-parse");
            cur.global_values.push_back(std::move(v));
        }
        parsed.push_back(std::move(cur));
    }

    out = std::move(parsed);
    std::cerr << "[amflow][syssolve-ckpt] HIT  key=" << key
              << " (" << out.size() << " eps-samples restored, skipping solve)\n";
    return true;
}

void system_checkpoint_store(const std::string& root,
                             const std::string& key,
                             const std::vector<AMFSystemSolution>& sol) {
    if (root.empty()) return;
    std::error_code ec;
    fs::create_directories(root, ec);
    if (ec) {
        std::cerr << "[amflow][syssolve-ckpt] WARN mkdir failed: "
                  << ec.message() << "\n";
        return;
    }

    fs::path final_path = fs::path(root) / ("system_" + key + ".ckpt");
    fs::path tmp_path   = fs::path(root)
                        / ("system_" + key + ".tmp." + std::to_string(::getpid()));

    {
        std::ofstream f(tmp_path, std::ios::trunc);
        if (!f) {
            std::cerr << "[amflow][syssolve-ckpt] WARN open temp failed\n";
            return;
        }
        f << "AMFLOW_SYSCKPT v" << kCkptVersion << "\n";
        f << "nsol " << sol.size() << "\n";
        for (const auto& s : sol) {
            f << "sol\n";
            f << "mv " << s.master_values.size() << "\n";
            for (const auto& v : s.master_values) f << acb_dump(v) << "\n";
            f << "gp " << s.global_preferred.size() << "\n";
            for (const auto& j : s.global_preferred) f << jint_serialize(j) << "\n";
            f << "gv " << s.global_values.size() << "\n";
            for (const auto& v : s.global_values) f << acb_dump(v) << "\n";
        }
        f.flush();
        if (!f) {
            std::cerr << "[amflow][syssolve-ckpt] WARN write failed\n";
            f.close();
            fs::remove(tmp_path, ec);
            return;
        }
    }  // close before rename

    // Atomic publish.  If another process beat us to the same key, the
    // content is by construction identical — drop our temp.
    fs::rename(tmp_path, final_path, ec);
    if (ec) {
        if (fs::exists(final_path, ec)) {
            fs::remove(tmp_path, ec);
        } else {
            std::cerr << "[amflow][syssolve-ckpt] WARN rename failed: "
                      << ec.message() << "\n";
            fs::remove(tmp_path, ec);
            return;
        }
    }
    std::cerr << "[amflow][syssolve-ckpt] STORED key=" << key
              << " (" << sol.size() << " eps-samples)\n";
}

}  // namespace amflow::pipeline
