// SPDX-License-Identifier: MIT
// ibp::load_adaptive — farm-mode load-adaptive parallelism.
//
// PROBLEM: N farm jobs launched together each spawn kira `-p<n_thread>`
// (and optionally an ε-parallel dsolve pool) with a *fixed* per-job
// thread count.  On a 96-core box with 14 jobs × kira -p6 the 1-min
// loadavg climbs to ~250 (2.6× oversubscribed), and contention costs
// roughly 2× wall on every job (measured on a 14-job reduction farm).
//
// FIX: when `AMFLOW_KIRA_P_AUTO=1` (kira -p) or
// `AMFLOW_EPS_PARALLEL=auto` (dsolve pool), pick the per-job thread
// count from the *current* spare capacity:
//
//     n = clamp(floor((nproc − loadavg1) / FARM_N), 1, cap)
//
// where `nproc` = std::thread::hardware_concurrency(), `loadavg1` is
// the 1-minute load average from /proc/loadavg, `FARM_N` is the number
// of concurrent farm jobs (set by the launcher via env `AMFLOW_FARM_N`,
// default 1), and `cap` is the configured upper bound (opts.n_thread
// for kira, n_eps for the ε-pool — applied at the call site).
//
// `nproc` and `loadavg1` are overridable via `AMFLOW_NPROC_OVERRIDE`
// and `AMFLOW_LOADAVG_OVERRIDE` for unit testing (no /proc on the test
// harness's mock).

#ifndef AMFLOW_IBP_LOAD_ADAPTIVE_HPP
#define AMFLOW_IBP_LOAD_ADAPTIVE_HPP

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <thread>

namespace amflow::ibp {

// Hardware thread count, overridable via AMFLOW_NPROC_OVERRIDE.
inline long load_adaptive_nproc() {
    if (const char* v = std::getenv("AMFLOW_NPROC_OVERRIDE"); v && *v) {
        char* end = nullptr;
        long n = std::strtol(v, &end, 10);
        if (end != v && n > 0) return n;
    }
    unsigned hc = std::thread::hardware_concurrency();
    return hc > 0 ? static_cast<long>(hc) : 1L;
}

// 1-minute load average, overridable via AMFLOW_LOADAVG_OVERRIDE.
// Returns 0.0 if /proc/loadavg cannot be read (so the formula
// degrades to nproc/FARM_N — safe upper bound).
inline double load_adaptive_loadavg1() {
    if (const char* v = std::getenv("AMFLOW_LOADAVG_OVERRIDE"); v && *v) {
        char* end = nullptr;
        double d = std::strtod(v, &end);
        if (end != v) return d;
    }
    std::ifstream f("/proc/loadavg");
    double l1 = 0.0;
    if (f) f >> l1;
    return l1;
}

// AMFLOW_FARM_N — number of concurrent farm jobs (caller-set; default 1).
inline long load_adaptive_farm_n() {
    if (const char* v = std::getenv("AMFLOW_FARM_N"); v && *v) {
        char* end = nullptr;
        long n = std::strtol(v, &end, 10);
        if (end != v && n > 0) return n;
    }
    return 1L;
}

// clamp(floor((nproc − loadavg1)/FARM_N), 1, cap).  `cap <= 0` is
// treated as "no upper bound" (returns the raw clamped-from-below value).
inline long load_adaptive_threads(long cap) {
    const long   nproc  = load_adaptive_nproc();
    const double load   = load_adaptive_loadavg1();
    const long   farm_n = load_adaptive_farm_n();
    double spare = (static_cast<double>(nproc) - load)
                 / static_cast<double>(std::max(1L, farm_n));
    long n = static_cast<long>(std::floor(spare));
    if (n < 1) n = 1;
    if (cap > 0 && n > cap) n = cap;
    return n;
}

}  // namespace amflow::ibp

#endif  // AMFLOW_IBP_LOAD_ADAPTIVE_HPP
