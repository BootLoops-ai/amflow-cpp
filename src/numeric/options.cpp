// SPDX-License-Identifier: MIT
// numeric::Options + GlobalScope + logging — implementation.
//
// amflow::numeric.

#include "amflow/numeric/options.hpp"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <ostream>
#include <stdexcept>

#include <flint/fmpq.h>

namespace amflow::numeric {

namespace {

// Thread-local singletons via function-local statics.
//
// THREADING NOTE (ε-grid parallelism, Lever 1): the AMF solve loop
// dispatches independent ε-samples to worker threads.  Both
// `AMFSystem::solve_one_eps` (per-node `run_direction`) and
// `ode::calc_zero` (LearnXOrder) wrap their work in `GlobalScope`
// guards.  With process-global singletons concurrent guards would
// race on commit/restore.  Making the singletons thread_local turns
// `GlobalScope` into a per-thread Block[].  Worker threads start
// with package defaults; the spawner must hand them an
// `OptionsSnapshot` (see `capture_options()` / `install_options()`).
// Single-threaded callers see no behaviour change.
GlobalOptions& mutable_global() {
    static thread_local GlobalOptions g{};
    return g;
}
ExpansionOptions& mutable_expansion() {
    static thread_local ExpansionOptions e{};
    return e;
}
RunningOptions& mutable_running() {
    static thread_local RunningOptions r{};
    return r;
}

constexpr double kLog2_10 = 3.3219280948873623478703194294894;

}  // namespace

// ---------------------------------------------------------------------------
GlobalOptions&    global_options()    { return mutable_global();    }
ExpansionOptions& expansion_options() { return mutable_expansion(); }
RunningOptions&   running_options()   { return mutable_running();   }

void set_global_options(const GlobalOptions& opt)       { mutable_global()    = opt; }
void set_expansion_options(const ExpansionOptions& opt) { mutable_expansion() = opt; }
void set_running_options(const RunningOptions& opt)     { mutable_running()   = opt; }

void set_default_options() {
    mutable_global()    = GlobalOptions{};
    mutable_expansion() = ExpansionOptions{};
    mutable_running()   = RunningOptions{};
}

OptionsSnapshot capture_options() {
    return OptionsSnapshot{mutable_global(), mutable_expansion(),
                           mutable_running()};
}

void install_options(const OptionsSnapshot& s) {
    mutable_global()    = s.global;
    mutable_expansion() = s.expansion;
    mutable_running()   = s.running;
}

// ---------------------------------------------------------------------------
long decimal_digits_to_bits(int digits) {
    if (digits <= 0) return 53;
    long bits = static_cast<long>(std::ceil(static_cast<double>(digits) * kLog2_10));
    if (bits < 53) bits = 53;
    return bits;
}

long working_prec_bits() {
    return decimal_digits_to_bits(working_pre());
}

// ---------------------------------------------------------------------------
int resonance_chop_pre() {
    // Audit override.
    if (const char* env = std::getenv("AMFLOW_RESONANCE_CHOP")) {
        char* end = nullptr;
        long v = std::strtol(env, &end, 10);
        if (end != env && v > 0) return static_cast<int>(v);
    }
    // Default: a guard band below the genuine arithmetic floor.  We use
    // ~80% of the working precision, leaving ~20% of working_pre digits of
    // slack so that round-off in the residue matrix cannot fake a NON-zero
    // separation (which would MISS a true resonance), while any separation
    // of order |a*eps| with |a*eps| > 10^-(0.8*working_pre) is correctly
    // kept distinct.  Floored at chop_pre so the test is never COARSER than
    // the legacy behaviour (a true integer/eps-independent resonance has
    // separation identically 0 and is caught at any tolerance).
    const int wp = working_pre();
    int tol = static_cast<int>(std::floor(0.8 * static_cast<double>(wp)));
    const int cp = chop_pre();
    if (tol < cp) tol = cp;
    // CONTRACT: 10^-resonance_chop_pre() must be representable as a POSITIVE
    // double, because several consumers (ode/asy.cpp mu_differ_by_integer,
    // the calcx00 rank/pivot tolerances in ode/zero.cpp) form the threshold
    // via std::pow(10.0, -tol).  For tol >~ 323 that pow underflows to 0.0
    // and every "|x| < threshold" test goes FALSE -- even for x == 0 exactly.
    // Concretely: a legs=[] vacuum solve_integrals job escalates the working
    // precision past 400 digits, tol became 0.8*wp = 342, eps_double
    // underflowed to 0, mu_differ_by_integer stopped merging ANY pair of
    // indicial exponents (including exactly-equal ones), the eta->infinity
    // Frobenius regions never recombined, and every vacuum master silently
    // transported to a hard 0.  Cap at 300 (the same threshold below which
    // numeric::acb_is_chop_zero trusts the double path); a separation below
    // 10^-300 is indistinguishable in any double-mediated comparison anyway,
    // and the calcx00 eps=0-resonance purpose of this tolerance (keeping
    // |a*eps| ~ 10^-6..10^-18 separations distinct) is unaffected.
    if (tol > 300) tol = 300;
    return tol;
}

// ---------------------------------------------------------------------------
void d0_eps_shift_fmpq(fmpq_t out) {
    const std::string& s = global_options().d0;

    fmpq_t d;
    fmpq_init(d);
    if (fmpq_set_str(d, s.c_str(), 10) != 0) {
        fmpq_clear(d);
        throw std::runtime_error(
            "d0_eps_shift_fmpq: failed to parse D0 rational from '" + s + "'");
    }

    // out = (4 - D0) / 2
    fmpq_t four;
    fmpq_init(four);
    fmpq_set_si(four, 4, 1);
    fmpq_sub(out, four, d);

    fmpq_t half;
    fmpq_init(half);
    fmpq_set_si(half, 1, 2);
    fmpq_mul(out, out, half);

    fmpq_clear(half);
    fmpq_clear(four);
    fmpq_clear(d);
}

// ---------------------------------------------------------------------------
GlobalScope::GlobalScope()
    : global(global_options()),
      expansion(expansion_options()),
      running(running_options()),
      saved_global_(global_options()),
      saved_expansion_(expansion_options()),
      saved_running_(running_options()) {}

GlobalScope::~GlobalScope() {
    if (committed_) {
        mutable_global()    = saved_global_;
        mutable_expansion() = saved_expansion_;
        mutable_running()   = saved_running_;
    }
}

void GlobalScope::commit() {
    if (committed_) return;
    committed_         = true;
    mutable_global()    = global;
    mutable_expansion() = expansion;
    mutable_running()   = running;
}

// ---------------------------------------------------------------------------
std::ostream& log_stream() {
    return std::clog;
}

void log_line(const std::string& msg) {
    if (silent_mode()) return;
    log_stream() << msg << '\n';
}

}  // namespace amflow::numeric
