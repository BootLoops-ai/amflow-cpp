// SPDX-License-Identifier: MIT
// ode::amflow — implementation.
//

#include "amflow/ode/amflow.hpp"

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <flint/acb.h>
#include <flint/arb.h>
#include <flint/arf.h>
#include <flint/fmpq.h>
#include <flint/mag.h>

#include "amflow/numeric/log.hpp"
#include "amflow/numeric/options.hpp"
#include "amflow/numeric/rational.hpp"
#include "amflow/ode/asy.hpp"
#include "amflow/ode/blocks.hpp"
#include "amflow/ode/inf.hpp"
#include "amflow/ode/path.hpp"
#include "amflow/ode/regular.hpp"
#include "amflow/ode/zero.hpp"

namespace amflow::ode {

using numeric::AcbValue;
using numeric::RationalComplex;
using numeric::RationalMatrix;

namespace {

bool debug_stage_enabled() {
    return numeric::log::trace_enabled("AMFLOW_DEBUG_STAGES");
}

void dump_acb_vector(const char* label,
                     const std::vector<AcbValue>& values) {
    if (!debug_stage_enabled()) return;
    std::cerr << "[system_stage] " << label
              << " size=" << values.size() << std::endl;
    for (std::size_t i = 0; i < values.size(); ++i) {
        std::cerr << "  [" << i << "] = "
                  << values[i].to_string(20) << std::endl;
    }
}

void dump_run(const std::vector<RationalComplex>& run) {
    if (!debug_stage_enabled()) return;
    std::cerr << "[system_stage] run size=" << run.size() << std::endl;
    for (std::size_t i = 0; i < run.size(); ++i) {
        std::cerr << "  run[" << i << "] = " << run[i].to_string()
                  << std::endl;
    }
}

int arb_mid_cmp(arb_srcptr a, arb_srcptr b) {
    return arf_cmp(arb_midref(a), arb_midref(b));
}

}  // namespace

System::System(RationalMatrix de,
               std::vector<BoundarySpec> bc,
               bool point_at_infinity)
    : de_(std::move(de)),
      bc_singular_(std::move(bc)),
      bc_kind_(BCKind::Singular),
      point_at_infinity_(point_at_infinity) {
    if (!numeric::silent_mode()) {
        numeric::log_line("LoadSystem: integrals=" + std::to_string(de_.rows())
                 + ", point=" + (point_at_infinity_ ? "infinity" : "regular"));
    }
}

System::System(RationalMatrix de,
               std::vector<AcbValue> bc_regular,
               RationalComplex point)
    : de_(std::move(de)),
      bc_regular_(std::move(bc_regular)),
      bc_kind_(BCKind::Regular),
      point_at_infinity_(false),
      point_(std::move(point)) {
    if (!numeric::silent_mode()) {
        numeric::log_line("LoadSystem: integrals=" + std::to_string(de_.rows())
                 + ", point=regular");
    }
}

void System::inf_to_regular(const RationalComplex& x0, long prec) {
    if (!point_at_infinity_) {
        throw std::runtime_error("inf_to_regular: current boundary is not at infinity");
    }
    if (bc_kind_ != BCKind::Singular) {
        throw std::runtime_error("inf_to_regular: bc is not in singular form");
    }

    auto asy = calc_inf(de_, bc_singular_, prec);
    asy_exp_ = std::move(asy);

    AcbValue inv_x0_acb;
    AcbValue x0_acb = x0.to_acb(prec);
    acb_inv(inv_x0_acb.raw(), x0_acb.raw(), prec);

    bc_regular_.clear();
    bc_regular_.reserve(asy_exp_.size());
    for (const auto& term : asy_exp_) {
        bc_regular_.push_back(evaluate_asy_expansion(term, inv_x0_acb.raw(), prec));
    }

    bc_kind_ = BCKind::Regular;
    point_at_infinity_ = false;
    point_ = x0;

    dump_acb_vector("after_inf_to_regular", bc_regular_);

    if (!numeric::silent_mode()) numeric::log_line("InfToRegular: switched to regular point");
}

void System::regular_run(const std::vector<RationalComplex>& run, long prec) {
    if (point_at_infinity_) {
        throw std::runtime_error("regular_run: current boundary is at infinity");
    }
    if (run.empty()) {
        throw std::invalid_argument("regular_run: empty run list");
    }
    if (!(run.front() == point_)) {
        throw std::runtime_error("regular_run: first run point does not match current boundary");
    }

    std::vector<AcbValue> run_acb;
    run_acb.reserve(run.size());
    for (const auto& p : run) run_acb.push_back(p.to_acb(prec));

    dump_run(run);
    dump_acb_vector("before_regular_run", bc_regular_);
    auto new_bc = calc_run(de_, bc_regular_, run_acb, prec);
    bc_regular_ = std::move(new_bc);
    point_ = run.back();
    dump_acb_vector("after_regular_run", bc_regular_);

    if (!numeric::silent_mode()) numeric::log_line("RegularRun: arrived at last point");
}

std::vector<std::vector<AcbValue>>
System::regular_interpolation(const std::vector<RationalComplex>& samples, long prec) {
    if (point_at_infinity_) {
        throw std::runtime_error("regular_interpolation: current boundary is at infinity");
    }

    const std::vector<AcbValue> saved_bc = [&]() {
        std::vector<AcbValue> tmp;
        tmp.reserve(bc_regular_.size());
        for (const auto& v : bc_regular_) tmp.push_back(v.clone());
        return tmp;
    }();
    const RationalComplex saved_point = point_;

    auto poles = numeric::get_poles(de_);

    auto radius = [&]() -> AcbValue {
        AcbValue minv;
        bool have = false;
        AcbValue p_acb = point_.to_acb(prec);
        for (const auto& pole : poles) {
            AcbValue diff;
            AcbValue pole_acb = pole.to_acb(prec);
            acb_sub(diff.raw(), p_acb.raw(), pole_acb.raw(), prec);
            AcbValue mag;
            acb_abs(acb_realref(mag.raw()), diff.raw(), prec);
            arb_zero(acb_imagref(mag.raw()));
            if (!have
                    || arb_mid_cmp(acb_realref(mag.raw()),
                                   acb_realref(minv.raw())) < 0) {
                minv = mag.clone();
                have = true;
            }
        }
        if (!have) {
            AcbValue out;
            out.set_zero();
            return out;
        }
        AcbValue out;
        arb_div_ui(acb_realref(out.raw()), acb_realref(minv.raw()),
                   static_cast<unsigned long>(numeric::run_radius()), prec);
        arb_zero(acb_imagref(out.raw()));
        return out;
    };

    auto findnum = [&](const std::vector<RationalComplex>& sam) -> std::size_t {
        std::size_t k = 0;
        AcbValue r = radius();
        AcbValue p_acb = point_.to_acb(prec);
        while (k < sam.size()) {
            AcbValue diff;
            AcbValue s_acb = sam[k].to_acb(prec);
            acb_sub(diff.raw(), p_acb.raw(), s_acb.raw(), prec);
            AcbValue d;
            acb_abs(acb_realref(d.raw()), diff.raw(), prec);
            arb_zero(acb_imagref(d.raw()));
            if (arb_mid_cmp(acb_realref(d.raw()),
                            acb_realref(r.raw())) > 0) break;
            ++k;
        }
        return k;
    };

    std::vector<std::vector<AcbValue>> vsets;
    std::vector<RationalComplex> residues = samples;

    if (!residues.empty() && residues.front() == point_) {
        std::vector<AcbValue> bc_clone;
        bc_clone.reserve(bc_regular_.size());
        for (const auto& v : bc_regular_) bc_clone.push_back(v.clone());
        vsets.push_back(std::move(bc_clone));
        residues.erase(residues.begin());
    }

    while (!residues.empty()) {
        std::size_t number = findnum(residues);
        if (number == 0) {
            auto seg = run_segment(poles, point_, residues.front(), prec);
            if (seg.empty()) {
                throw std::runtime_error(
                    "regular_interpolation: cannot lay a contour to next sample");
            }
            regular_run(seg, prec);
            std::vector<AcbValue> bc_clone;
            bc_clone.reserve(bc_regular_.size());
            for (const auto& v : bc_regular_) bc_clone.push_back(v.clone());
            vsets.push_back(std::move(bc_clone));
            residues.erase(residues.begin());
        } else {
            auto nheq  = nh_equations(de_, EquationMode::Regular);
            auto nheqn = nh_equations_num(nheq, prec);
            AcbValue point_acb = point_.to_acb(prec);
            auto rule = calcx1x2(nheqn, bc_regular_, point_acb.raw(), prec);
            for (std::size_t i = 0; i < number; ++i) {
                AcbValue dh;
                AcbValue s_acb = residues[i].to_acb(prec);
                acb_sub(dh.raw(), s_acb.raw(), point_acb.raw(), prec);
                vsets.push_back(evaluate_taylor(rule, dh.raw(), prec));
            }
            bc_regular_.clear();
            for (const auto& v : vsets.back()) bc_regular_.push_back(v.clone());
            point_ = residues[number - 1];
            residues.erase(residues.begin(), residues.begin() + number);
        }
    }

    bc_regular_.clear();
    for (const auto& v : saved_bc) bc_regular_.push_back(v.clone());
    point_ = saved_point;

    return vsets;
}

void System::solve_asy_exp(long prec) {
    if (point_at_infinity_) {
        throw std::runtime_error("solve_asy_exp: current boundary is at infinity");
    }
    AcbValue x0_acb = point_.to_acb(prec);
    asy_exp_ = calc_zero(de_, bc_regular_, x0_acb.raw(), prec);
    if (debug_stage_enabled()) {
        std::cerr << "[system_stage] after_solve_asy_exp size="
                  << asy_exp_.size() << std::endl;
        for (std::size_t i = 0; i < asy_exp_.size(); ++i) {
            std::cerr << "  asy[" << i << "] = "
                      << pick_zero_rule_s(asy_exp_[i], prec).to_string(20)
                      << std::endl;
        }
    }
}

std::vector<AcbValue> System::pick_zero_solution(long prec) const {
    std::vector<AcbValue> out;
    out.reserve(asy_exp_.size());
    for (const auto& a : asy_exp_) out.push_back(pick_zero_rule_s(a, prec));
    // [bug U1] The calcx00 boundary-match solve goes through an exact-
    // rational Dixon solve (zero.cpp: rationalize_acb -> fmpq_mat_solve_dixon
    // -> arb_set_fmpq), which structurally drops the Arb radius.  bc_regular_
    // still carries the accumulated x_order truncation tail from
    // inf_to_regular + every calc_run step (see evaluate_taylor /
    // evaluate_expansion).  Re-inject that radius onto the picked-zero
    // values so the returned ball is honest about x_order.
    if (xorder_tail_bound_enabled()
            && bc_regular_.size() == out.size()) {
        mag_t r;
        mag_init(r);
        for (std::size_t i = 0; i < out.size(); ++i) {
            mag_max(r, arb_radref(acb_realref(bc_regular_[i].raw())),
                       arb_radref(acb_imagref(bc_regular_[i].raw())));
            acb_add_error_mag(out[i].raw(), r);
        }
        mag_clear(r);
    }
    return out;
}

std::vector<AcbValue> amflow(const RationalMatrix& de,
                               const std::vector<BoundarySpec>& bcs,
                               long prec) {
    const bool stage_time = std::getenv("AMFLOW_STAGE_TIMING") != nullptr;
    using clk = std::chrono::steady_clock;
    auto t0 = clk::now();
    auto mark = [&](const char* tag) {
        if (!stage_time) return;
        auto t1 = clk::now();
        std::cerr << "[amflow][stage] N=" << de.rows() << " " << tag << " "
                  << std::chrono::duration<double>(t1 - t0).count() << "s\n";
        t0 = t1;
    };
    if (!numeric::silent_mode()) numeric::log_line("AMFlow: start.");

    std::vector<BoundarySpec> bcs_clone;
    bcs_clone.reserve(bcs.size());
    for (const auto& b : bcs) bcs_clone.push_back(clone_boundary(b));

    System sys(de, std::move(bcs_clone), /*point_at_infinity=*/true);

    if (!numeric::silent_mode()) numeric::log_line("AMFlow: building integration contour");
    auto poles = numeric::get_poles(de);
    auto run = run_eta(poles, prec);
    if (run.empty()) {
        throw std::runtime_error("amflow: integration contour not generated");
    }
    if (!numeric::silent_mode()) {
        numeric::log_line("AMFlow: contour has " + std::to_string(run.size()) + " regular points");
    }
    mark("contour");

    if (!numeric::silent_mode()) numeric::log_line("AMFlow: solving near eta = Infinity");
    sys.inf_to_regular(run.front(), prec);
    mark("calc_inf");

    if (!numeric::silent_mode()) numeric::log_line("AMFlow: regular running");
    sys.regular_run(run, prec);
    mark("calc_run");

    if (!numeric::silent_mode()) numeric::log_line("AMFlow: solving near the last point");
    sys.solve_asy_exp(prec);
    mark("calc_zero");

    if (!numeric::silent_mode()) numeric::log_line("AMFlow: picking zero values");
    auto sol = sys.pick_zero_solution(prec);

    if (!numeric::silent_mode()) numeric::log_line("AMFlow: finished.");
    return sol;
}

// ===========================================================================
//  calcx00 audit: diagnostic helpers (env-gated callers).
// ===========================================================================

namespace {

// Scan a batched column set for non-finite (NaN/Inf mid or radius) acb
// entries; report the first few offenders and a per-column count.  Purely
// observational (AMF_DBG_X00).
void report_nonfinite_cols(const std::vector<std::vector<AcbValue>>& cols,
                           const char* tag) {
    long total = 0, shown = 0;
    for (std::size_t c = 0; c < cols.size(); ++c) {
        for (std::size_t i = 0; i < cols[c].size(); ++i) {
            if (!acb_is_finite(cols[c][i].raw())) {
                ++total;
                if (shown < 8) {
                    fprintf(stderr,
                            "[X00-scan] %s: NON-FINITE at col=%zu master=%zu\n",
                            tag, c, i);
                    ++shown;
                }
            }
        }
    }
    fprintf(stderr, "[X00-scan] %s: %ld non-finite entries total "
                    "(%zu cols x %zu masters)\n",
            tag, total, cols.size(),
            cols.empty() ? (std::size_t)0 : cols[0].size());
}

// Exact bc_last checkpoint loader (counterpart of the AMFLOW_X00_DUMP_BC
// writer in calc_zero_multi).  Format: line 1 "X00BC v1 <ncols> <N>", then
// per acb entry two lines (real arb, imag arb) in arb_dump_str format,
// column-major over entries within a column.
std::vector<std::vector<AcbValue>>
load_bc_checkpoint(const char* path, std::size_t ncols, std::size_t N) {
    std::FILE* f = std::fopen(path, "r");
    if (!f) {
        throw std::runtime_error(std::string("AMFLOW_X00_BC_FILE: cannot "
                                             "open ") + path);
    }
    char header[128];
    unsigned long fc = 0, fn = 0;
    if (std::fscanf(f, "%15s %15s %lu %lu ", header, header + 16, &fc, &fn)
            != 4 || fc != ncols || fn != N) {
        std::fclose(f);
        throw std::runtime_error("AMFLOW_X00_BC_FILE: header mismatch "
                                 "(wrong file / wrong system?)");
    }
    std::vector<std::vector<AcbValue>> out(ncols);
    for (auto& col : out) col.resize(N);
    char line[4096];
    for (std::size_t c = 0; c < ncols; ++c) {
        for (std::size_t i = 0; i < N; ++i) {
            for (int part = 0; part < 2; ++part) {
                if (!std::fgets(line, sizeof line, f)) {
                    std::fclose(f);
                    throw std::runtime_error(
                        "AMFLOW_X00_BC_FILE: truncated checkpoint");
                }
                std::size_t len = std::strlen(line);
                while (len && (line[len - 1] == '\n' || line[len - 1] == '\r'))
                    line[--len] = '\0';
                arb_ptr dst = part == 0 ? acb_realref(out[c][i].raw())
                                        : acb_imagref(out[c][i].raw());
                if (arb_load_str(dst, line) != 0) {
                    std::fclose(f);
                    throw std::runtime_error(
                        "AMFLOW_X00_BC_FILE: arb_load_str failed");
                }
            }
        }
    }
    std::fclose(f);
    return out;
}

}  // namespace

// ===========================================================================
//  amflow_multi  (fp-probe column batching)
// ===========================================================================

std::vector<std::vector<AcbValue>>
amflow_multi(const RationalMatrix& de,
             const std::vector<std::vector<BoundarySpec>>& bcs_cols,
             long prec) {
    const std::size_t ncols = bcs_cols.size();
    std::vector<std::vector<AcbValue>> out(ncols);
    if (ncols == 0) return out;

    const bool stage_time = std::getenv("AMFLOW_STAGE_TIMING") != nullptr;
    using clk = std::chrono::steady_clock;
    auto t0 = clk::now();
    auto mark = [&](const char* tag) {
        if (!stage_time) return;
        auto t1 = clk::now();
        std::cerr << "[amflow][stage] N=" << de.rows()
                  << " batched_cols=" << ncols << " " << tag << " "
                  << std::chrono::duration<double>(t1 - t0).count() << "s\n";
        t0 = t1;
    };

    if (!numeric::silent_mode()) {
        numeric::log_line("AMFlow: start (fp-probe batched, "
                          + std::to_string(ncols) + " cols).");
    }

    // calcx00 audit: scan the eta=inf boundary VALUES at
    // entry (pre-transport, CHEAP).  A NaN/Inf here means the divergence is
    // born in boundary assembly, before any transport.  AMF_DBG_X00 only.
    if (std::getenv("AMF_DBG_X00")) {
        long total = 0, shown = 0;
        for (std::size_t c = 0; c < ncols; ++c) {
            for (std::size_t i = 0; i < bcs_cols[c].size(); ++i) {
                for (std::size_t e = 0; e < bcs_cols[c][i].size(); ++e) {
                    const auto& be = bcs_cols[c][i][e];
                    if (!acb_is_finite(be.value.raw()) ||
                        !acb_is_finite(be.mu.raw())) {
                        ++total;
                        if (shown < 8) {
                            char* vr = arb_get_str(
                                acb_realref(be.value.raw()), 10, 0);
                            fprintf(stderr,
                                    "[X00-scan] bcs_cols entry NON-FINITE: "
                                    "col=%zu master=%zu entry=%zu value.re=%s "
                                    "mu_finite=%d\n",
                                    c, i, e, vr,
                                    acb_is_finite(be.mu.raw()) ? 1 : 0);
                            flint_free(vr);
                            ++shown;
                        }
                    }
                }
            }
        }
        fprintf(stderr, "[X00-scan] bcs_cols (eta=inf boundary specs): "
                        "%ld non-finite entries total (%zu cols)\n",
                total, ncols);
    }

    // ---- contour (shared) ----
    if (!numeric::silent_mode()) numeric::log_line("AMFlow: building integration contour");
    auto poles = numeric::get_poles(de);
    auto run = run_eta(poles, prec);
    if (run.empty()) {
        throw std::runtime_error("amflow_multi: integration contour not generated");
    }
    if (!numeric::silent_mode()) {
        numeric::log_line("AMFlow: contour has " + std::to_string(run.size())
                          + " regular points");
    }
    mark("contour");

    // ------------------------------------------------------------------
    //  calcx00 audit diagnostic re-entry hooks.  Both are
    //  env-gated and default-off; the env-OFF path below is byte-identical
    //  to the original sequence.
    //    AMFLOW_X00_SYNTH_BC=1   : skip the (expensive) eta=inf + regular
    //        transport entirely and hand CalcZero an all-ones bc_last.
    //        The CalcZero-side matrix (NormalizeMat / kernel rescue /
    //        asymptotic behavior / calcx00 kk-chains) is bc-independent,
    //        so this reproduces any STRUCTURAL Asys-side failure in
    //        seconds instead of ~1 h of transport.  DIAGNOSTIC ONLY —
    //        the returned values are meaningless.
    //    AMFLOW_X00_BC_FILE=path : load bc_last from a checkpoint written
    //        by AMFLOW_X00_DUMP_BC (see calc_zero_multi) — exact arb
    //        deserialization, i.e. a mid-pipeline re-entry at CalcZero
    //        without replaying transport.
    // ------------------------------------------------------------------
    const char* synth_env   = std::getenv("AMFLOW_X00_SYNTH_BC");
    const char* bcfile_env  = std::getenv("AMFLOW_X00_BC_FILE");
    const bool  synth_bc    = synth_env && synth_env[0] == '1';
    const std::size_t Nsys  = de.rows();

    // AMFLOW_X00_ONECOL=1 (DIAGNOSTIC ONLY): truncate
    // the batch to its first column so the eta=inf exec (where the
    // non-finite entries are born) reproduces in minutes instead of hours.
    // Output is a 1-column probe fragment — never ship it.
    std::vector<std::vector<BoundarySpec>> bcs_one;
    const std::vector<std::vector<BoundarySpec>>* bcs_use = &bcs_cols;
    std::size_t ncols_eff = ncols;
    {
        const char* oc = std::getenv("AMFLOW_X00_ONECOL");
        if (oc && *oc == '1' && ncols > 1) {
            numeric::log_line("AMFlow: [X00_ONECOL] DIAGNOSTIC: batch "
                              "truncated to column 0");
            bcs_one.push_back(clone_boundaries(bcs_cols[0]));
            bcs_use = &bcs_one;
            ncols_eff = 1;
        }
    }

    std::vector<std::vector<AcbValue>> bc_last;
    if (synth_bc || bcfile_env) {
        if (synth_bc) {
            numeric::log_line("AMFlow: [X00_SYNTH_BC] DIAGNOSTIC: transport "
                              "SKIPPED, bc_last = ones (values meaningless)");
            bc_last.resize(ncols);
            for (std::size_t c = 0; c < ncols; ++c) {
                bc_last[c].resize(Nsys);
                for (std::size_t i = 0; i < Nsys; ++i)
                    acb_one(bc_last[c][i].raw());
            }
        } else {
            numeric::log_line(std::string("AMFlow: [X00_BC_FILE] loading "
                              "bc_last checkpoint from ") + bcfile_env);
            bc_last = load_bc_checkpoint(bcfile_env, ncols, Nsys);
        }
        mark("calc_run");
    } else {

    // ---- CalcInf (shared region plans, per-column execute) ----
    if (!numeric::silent_mode()) numeric::log_line("AMFlow: solving near eta = Infinity");
    auto asys = calc_inf_multi(de, *bcs_use, prec);
    mark("calc_inf");

    // ---- inf -> regular per column (System::inf_to_regular tail) ----
    AcbValue x0_acb = run.front().to_acb(prec);
    AcbValue inv_x0_acb;
    acb_inv(inv_x0_acb.raw(), x0_acb.raw(), prec);

    std::vector<std::vector<AcbValue>> bc_reg(ncols_eff);
    for (std::size_t c = 0; c < ncols_eff; ++c) {
        bc_reg[c].reserve(asys[c].size());
        for (const auto& term_list : asys[c]) {
            bc_reg[c].push_back(
                evaluate_asy_expansion(term_list, inv_x0_acb.raw(), prec));
        }
    }
    if (!numeric::silent_mode()) {
        numeric::log_line("InfToRegular: switched to regular point [batched]");
    }
    if (std::getenv("AMF_DBG_X00")) {
        report_nonfinite_cols(bc_reg, "bc_reg (after InfToRegular)");
    }
    if (ncols_eff != ncols) {
        throw std::runtime_error(
            "AMFLOW_X00_ONECOL diagnostic complete (bc_reg scan above); "
            "this mode never produces shippable output");
    }

    // ---- CalcRun (shared nheqn + per-step re-expansion; per-column march) ----
    if (!numeric::silent_mode()) numeric::log_line("AMFlow: regular running");
    std::vector<AcbValue> run_acb;
    run_acb.reserve(run.size());
    for (const auto& p : run) run_acb.push_back(p.to_acb(prec));

    auto nheq  = nh_equations(de, EquationMode::Regular);
    auto nheqn = nh_equations_num(nheq, prec);
    bc_last = calc_run_multi(nheqn, bc_reg, run_acb, prec);
    if (!numeric::silent_mode()) numeric::log_line("RegularRun: arrived at last point");
    if (std::getenv("AMF_DBG_X00")) {
        report_nonfinite_cols(bc_last, "bc_last (after RegularRun)");
    }
    mark("calc_run");

    }  // end env-OFF transport path

    // ---- CalcZero (shared NormalizeMat; per-column calcx00) ----
    if (!numeric::silent_mode()) numeric::log_line("AMFlow: solving near the last point");
    AcbValue xb_acb = run.back().to_acb(prec);
    auto asy0 = calc_zero_multi(de, bc_last, xb_acb.raw(), prec);
    mark("calc_zero");

    // ---- pick zero values per column (System::pick_zero_solution) ----
    if (!numeric::silent_mode()) numeric::log_line("AMFlow: picking zero values");
    const bool tail = xorder_tail_bound_enabled();
    for (std::size_t c = 0; c < ncols; ++c) {
        out[c].reserve(asy0[c].size());
        for (const auto& a : asy0[c]) out[c].push_back(pick_zero_rule_s(a, prec));
        if (tail && bc_last[c].size() == out[c].size()) {
            mag_t r;
            mag_init(r);
            for (std::size_t i = 0; i < out[c].size(); ++i) {
                mag_max(r, arb_radref(acb_realref(bc_last[c][i].raw())),
                           arb_radref(acb_imagref(bc_last[c][i].raw())));
                acb_add_error_mag(out[c][i].raw(), r);
            }
            mag_clear(r);
        }
    }

    if (!numeric::silent_mode()) numeric::log_line("AMFlow: finished (batched).");
    return out;
}

}  // namespace amflow::ode
