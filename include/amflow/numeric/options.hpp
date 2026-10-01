// SPDX-License-Identifier: MIT
// numeric::Options — three option groups + GlobalScope RAII override.
//
//
// WorkingPre, ChopPre, RationalizePre are *decimal* digit counts (matching
// MMA semantics).  Internally Arb / Acb work in *binary* precision; conversion
// is provided by working_prec_bits().
//
// Options are global per-process state.  This is intentional: AMFlow's
// algorithms read these values dozens of times in inner loops; threading
// them through every signature would clutter without value (the algorithm
// is single-runtime).  Use GlobalScope (RAII) for temporary overrides — the
// analogue of MMA's `Block[{XOrder = ...}, ...]` idiom.

#ifndef AMFLOW_NUMERIC_OPTIONS_HPP
#define AMFLOW_NUMERIC_OPTIONS_HPP

#include <iosfwd>
#include <string>

#include <flint/fmpq.h>

namespace amflow::numeric {

// ---------------------------------------------------------------------------
//  Option groups (POD)
// ---------------------------------------------------------------------------

struct GlobalOptions {
    int  working_pre     = 100;     // decimal digits, global working precision
    int  chop_pre        = 20;      // |x| < 10^-chop_pre is treated as zero
    bool silent_mode     = false;   // suppress informational logging
    int  rationalize_pre = 100;     // tolerance when rationalising floats (mirrors upstream `RationalizePre = 100`)

    // Spacetime dimension D0 in d = D0 - 2*eps; eps -> 0 corresponds to d = D0.
    // Mirrors AMFlow.m's "D0" option (default 4).  Encoded as a rational
    // string ("4", "7/3", "1/3"); parsed by d0_eps_shift_fmpq().
    //
    // Effect: Laurent-mode solve_integrals samples the family at the internal
    // eps grid `user_eps + (4 - D0)/2` (mirrors AMFlow.m:1342/1351) but fits
    // the Laurent expansion against the *original* user-facing grid (mirrors
    // AMFlow.m:1356).  Sampled black-box mode does NOT apply the shift.
    std::string d0 = "4";
};

struct ExpansionOptions {
    int x_order       = 100;        // main truncation order of power series
    int extra_x_order = 20;         // extra orders for boundary determination / CalcTaylor
    int learn_x_order = -1;         // negative => skip learning phase
    int test_x_order  = 5;          // truncation when probing log powers after learning
};

struct RunningOptions {
    enum class Direction {
        Re,         // along +Re axis with small Im perturbation
        Im,         // along +Im axis with small Re perturbation
        NegRe,      // along -Re axis
        NegIm,      // along -Im axis  (the default)
        Custom      // user-supplied direction (custom_direction_re/im below)
    };

    int       run_radius     = 2;
    int       run_length     = 1000;     // mirrors upstream `RunLength = 1000`
    int       run_candidate  = 10;
    Direction run_direction  = Direction::NegIm;

    std::string custom_direction_re = "0";
    std::string custom_direction_im = "0";
};

// ---------------------------------------------------------------------------
//  Accessors
// ---------------------------------------------------------------------------

GlobalOptions&    global_options();
ExpansionOptions& expansion_options();
RunningOptions&   running_options();

// Inline getters used in inner loops (avoid struct copy).
inline int  working_pre()     { return global_options().working_pre; }
inline int  chop_pre()        { return global_options().chop_pre; }
inline bool silent_mode()     { return global_options().silent_mode; }
inline int  rationalize_pre() { return global_options().rationalize_pre; }

// Tolerance (decimal digits) for deciding whether two indicial exponents are
// EQUAL / integer-separated, i.e. for resonance detection in the eta->0
// Frobenius solve.  This is a fundamentally different question from "is this
// coefficient numerically negligible" (chop_pre): an eps-deformed exponent
// pair {mu, mu + a*eps} is GENUINELY distinct for any eps != 0, and merging
// it as a resonance silently destroys the genuine 1/(a*eps) pole structure
// (the master then loses its eps-Laurent pole and the eps->0 extrapolation is
// corrupted to 0 digits).  At a single numeric eps we cannot see the
// eps-dependence, so we must NOT merge unless the separation is below the
// numerical noise floor -- ~ the working precision, not the coarse chop_pre.
//
// Concretely:  chop_pre is typically O(20-40) digits while the eps grid sits
// at |eps| ~ 10^-O(goal/eps_order); for high-precision 3-loop targets
// |a*eps| can be far below 10^-chop_pre yet far above 10^-working_pre, which
// is exactly the regime where the old chop_pre test spuriously resonates.
// Using ~80% of working_pre keeps a safe guard band below the genuine
// arithmetic floor while still catching true (eps-independent) resonances,
// whose separation is identically 0.  Floored at chop_pre so behaviour is
// never coarser than before; env AMFLOW_RESONANCE_CHOP overrides for audits.
int resonance_chop_pre();
inline const std::string& d0() { return global_options().d0; }

inline int x_order()       { return expansion_options().x_order; }
inline int extra_x_order() { return expansion_options().extra_x_order; }
inline int learn_x_order() { return expansion_options().learn_x_order; }
inline int test_x_order()  { return expansion_options().test_x_order; }

inline int                       run_radius()    { return running_options().run_radius; }
inline int                       run_length()    { return running_options().run_length; }
inline int                       run_candidate() { return running_options().run_candidate; }
inline RunningOptions::Direction run_direction() { return running_options().run_direction; }

// ---------------------------------------------------------------------------
//  Mutators
// ---------------------------------------------------------------------------

void set_global_options(const GlobalOptions& opt);
void set_expansion_options(const ExpansionOptions& opt);
void set_running_options(const RunningOptions& opt);

// Re-initialise all three groups to package defaults (= MMA defaults).
void set_default_options();

// ---------------------------------------------------------------------------
//  Precision helpers
// ---------------------------------------------------------------------------

// Convert the current working_pre (decimal digits) to bit precision for
// Arb/Acb.  Always >= 53 (so Arb does not get a non-positive prec).
//
// Formula: ceil(working_pre * log2(10)).
long working_prec_bits();

// Same conversion for an arbitrary digit count.
long decimal_digits_to_bits(int digits);

// Compute the eps shift `(4 - D0) / 2` from GlobalOptions::d0, writing it
// into `out` (an initialised fmpq_t).  Throws std::runtime_error if d0 cannot
// be parsed as a rational.
//
// Mirrors AMFlow.m:1342 / 1351:  epslist = epslist0 + (4-$D0)*1/2.
void d0_eps_shift_fmpq(fmpq_t out);

// ---------------------------------------------------------------------------
//  Thread-local option snapshot (for ε-parallel workers)
// ---------------------------------------------------------------------------
//
// The three option singletons are *thread-local* (see options.cpp).  A
// worker thread therefore starts with package defaults, NOT the parent
// thread's options.  Capture in the parent → install in the worker.
struct OptionsSnapshot {
    GlobalOptions    global;
    ExpansionOptions expansion;
    RunningOptions   running;
};
OptionsSnapshot capture_options();
void            install_options(const OptionsSnapshot& s);

// ---------------------------------------------------------------------------
//  GlobalScope  (RAII override)
// ---------------------------------------------------------------------------
//
//   {
//     GlobalScope guard;
//     guard.expansion.x_order = 30;
//     guard.global.silent_mode = true;
//     guard.commit();          // installs the override
//     ... do work ...
//   }   // dtor restores the previous values
//
// Mirrors MMA's `Block[{XOrder = 30, SilentMode = True}, ...]`.

class GlobalScope {
public:
    GlobalScope();
    ~GlobalScope();

    GlobalScope(const GlobalScope&)            = delete;
    GlobalScope& operator=(const GlobalScope&) = delete;

    GlobalOptions    global;
    ExpansionOptions expansion;
    RunningOptions   running;

    // Push the current values of `global`, `expansion`, `running` into the
    // global state.  Idempotent (subsequent calls are no-ops).
    void commit();

private:
    GlobalOptions    saved_global_;
    ExpansionOptions saved_expansion_;
    RunningOptions   saved_running_;
    bool             committed_ = false;
};

// ---------------------------------------------------------------------------
//  Logging  (mirrors AMFPrint in DESolver.m)
// ---------------------------------------------------------------------------
//
// Routes to std::clog so callers can redirect.  log_line() respects
// silent_mode (suppresses output when true).

std::ostream& log_stream();
void          log_line(const std::string& msg);

}  // namespace amflow::numeric

#endif  // AMFLOW_NUMERIC_OPTIONS_HPP
