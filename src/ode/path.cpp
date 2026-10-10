// SPDX-License-Identifier: MIT
// ode::path — implementation.
//

#include "amflow/ode/path.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>
#include <vector>

#include <flint/acb.h>
#include <flint/arb.h>
#include <flint/arf.h>
#include <flint/fmpq.h>
#include <flint/fmpz.h>

namespace amflow::ode {

using numeric::AcbValue;
using numeric::RationalComplex;

namespace {

int arb_mid_cmp(arb_srcptr a, arb_srcptr b) {
    return arf_cmp(arb_midref(a), arb_midref(b));
}

AcbValue rc_abs_value(const RationalComplex& z, long prec) {
    AcbValue zacb = z.to_acb(prec);
    AcbValue out;
    acb_abs(acb_realref(out.raw()), zacb.raw(), prec);
    arb_zero(acb_imagref(out.raw()));
    return out;
}

RationalComplex rc_add(const RationalComplex& a, const RationalComplex& b) {
    RationalComplex z;
    fmpq_add(z.re(), a.re(), b.re());
    fmpq_add(z.im(), a.im(), b.im());
    return z;
}

RationalComplex rc_sub(const RationalComplex& a, const RationalComplex& b) {
    RationalComplex z;
    fmpq_sub(z.re(), a.re(), b.re());
    fmpq_sub(z.im(), a.im(), b.im());
    return z;
}

RationalComplex rc_mul(const RationalComplex& a, const RationalComplex& b) {
    fmpq_t arbr, aibi, arbi, aibr;
    fmpq_init(arbr); fmpq_init(aibi); fmpq_init(arbi); fmpq_init(aibr);
    fmpq_mul(arbr, a.re(), b.re());
    fmpq_mul(aibi, a.im(), b.im());
    fmpq_mul(arbi, a.re(), b.im());
    fmpq_mul(aibr, a.im(), b.re());
    RationalComplex z;
    fmpq_sub(z.re(), arbr, aibi);
    fmpq_add(z.im(), arbi, aibr);
    fmpq_clear(arbr); fmpq_clear(aibi); fmpq_clear(arbi); fmpq_clear(aibr);
    return z;
}

RationalComplex rc_div(const RationalComplex& a, const RationalComplex& b) {
    fmpq_t br2, bi2, denom;
    fmpq_init(br2); fmpq_init(bi2); fmpq_init(denom);
    fmpq_mul(br2, b.re(), b.re());
    fmpq_mul(bi2, b.im(), b.im());
    fmpq_add(denom, br2, bi2);
    if (fmpq_is_zero(denom)) {
        fmpq_clear(br2); fmpq_clear(bi2); fmpq_clear(denom);
        throw std::domain_error("rc_div: division by zero");
    }
    fmpq_t arbr, aibi, aibr_neg, arbi_neg;
    fmpq_init(arbr); fmpq_init(aibi); fmpq_init(aibr_neg); fmpq_init(arbi_neg);
    fmpq_mul(arbr, a.re(), b.re());
    fmpq_mul(aibi, a.im(), b.im());
    fmpq_mul(arbi_neg, a.re(), b.im());
    fmpq_neg(arbi_neg, arbi_neg);
    fmpq_mul(aibr_neg, a.im(), b.re());
    RationalComplex z;
    fmpq_t numer_re, numer_im;
    fmpq_init(numer_re); fmpq_init(numer_im);
    fmpq_add(numer_re, arbr, aibi);
    fmpq_add(numer_im, arbi_neg, aibr_neg);
    fmpq_div(z.re(), numer_re, denom);
    fmpq_div(z.im(), numer_im, denom);
    fmpq_clear(numer_re); fmpq_clear(numer_im);
    fmpq_clear(arbr); fmpq_clear(aibi); fmpq_clear(aibr_neg); fmpq_clear(arbi_neg);
    fmpq_clear(br2); fmpq_clear(bi2); fmpq_clear(denom);
    return z;
}

void rc_abs2(fmpq_t out, const RationalComplex& z) {
    fmpq_t r2, i2;
    fmpq_init(r2); fmpq_init(i2);
    fmpq_mul(r2, z.re(), z.re());
    fmpq_mul(i2, z.im(), z.im());
    fmpq_add(out, r2, i2);
    fmpq_clear(r2); fmpq_clear(i2);
}

bool acb_real_to_fmpq(fmpq_t out, acb_srcptr x, int rationalize_digits, long prec) {
    // PRECISION-MISMATCH FIX (see src/ode/inf.cpp:acb_real_to_fmpq_local for
    // the rationale): cap rationalize_digits at what working_prec can actually
    // resolve.  Without this cap, rationalize_pre=100 with working_pre=120
    // bits captures binary representation noise as a "real" rational.
    long working_digits = static_cast<long>(prec * 0.301029995663981195L);
    if (working_digits >= 5) working_digits -= 5;
    if (working_digits < 1) working_digits = 1;
    if (rationalize_digits > working_digits) {
        rationalize_digits = static_cast<int>(working_digits);
    }
    // The imaginary-part test and the simplest-rational recovery live in
    // numeric::acb_real_to_fmpq_simplest (an arf threshold, so
    // rationalize_pre > ~308 does not underflow the way a double 10^-d did).
    return numeric::acb_real_to_fmpq_simplest(out, x, rationalize_digits);
}

AcbValue min_distance_to_poles(const RationalComplex& point,
                                 const std::vector<RationalComplex>& poles,
                                 long prec) {
    AcbValue best;
    bool have = false;
    for (const auto& p : poles) {
        RationalComplex diff = rc_sub(point, p);
        AcbValue d = rc_abs_value(diff, prec);
        if (!have
                || arb_mid_cmp(acb_realref(d.raw()),
                               acb_realref(best.raw())) < 0) {
            best = d.clone();
            have = true;
        }
    }
    return best;
}

}  // namespace

RationalComplex first_step(const std::vector<RationalComplex>& poles,
                            long prec) {
    if (poles.empty()
        || (poles.size() == 1 && poles[0].is_zero())) {
        return RationalComplex::from_fractions(1, 1, 0, 1);
    }
    AcbValue max_abs;
    bool have = false;
    for (const auto& p : poles) {
        AcbValue a = rc_abs_value(p, prec);
        if (!have
                || arb_mid_cmp(acb_realref(a.raw()),
                               acb_realref(max_abs.raw())) > 0) {
            max_abs = a.clone();
            have = true;
        }
    }
    if (!have || max_abs.is_zero()) {
        return RationalComplex::from_fractions(1, 1, 0, 1);
    }
    AcbValue val;
    arb_mul_ui(acb_realref(val.raw()), acb_realref(max_abs.raw()),
               static_cast<unsigned long>(numeric::run_radius()), prec);
    arb_zero(acb_imagref(val.raw()));
    fmpq_t q;
    fmpq_init(q);
    if (!acb_real_to_fmpq(q, val.raw(), numeric::rationalize_pre(), prec)) {
        fmpq_clear(q);
        throw std::runtime_error("first_step: rationalisation failed");
    }
    RationalComplex z;
    fmpq_set(z.re(), q);
    fmpq_clear(q);
    return z;
}

RationalComplex last_step(const std::vector<RationalComplex>& poles,
                           long prec) {
    if (poles.empty()
        || (poles.size() == 1 && poles[0].is_zero())) {
        return RationalComplex::from_fractions(1, 1, 0, 1);
    }
    AcbValue min_abs;
    bool have = false;
    for (const auto& p : poles) {
        if (p.is_zero()) continue;
        AcbValue a = rc_abs_value(p, prec);
        if (!have
                || arb_mid_cmp(acb_realref(a.raw()),
                               acb_realref(min_abs.raw())) < 0) {
            min_abs = a.clone();
            have = true;
        }
    }
    if (!have || min_abs.is_zero()) {
        return RationalComplex::from_fractions(1, 1, 0, 1);
    }
    AcbValue val;
    arb_div_ui(acb_realref(val.raw()), acb_realref(min_abs.raw()),
               static_cast<unsigned long>(numeric::run_radius()), prec);
    arb_zero(acb_imagref(val.raw()));
    fmpq_t q;
    fmpq_init(q);
    if (!acb_real_to_fmpq(q, val.raw(), numeric::rationalize_pre(), prec)) {
        fmpq_clear(q);
        throw std::runtime_error("last_step: rationalisation failed");
    }
    RationalComplex z;
    fmpq_set(z.re(), q);
    fmpq_clear(q);
    return z;
}

std::vector<RationalComplex>
run_unit(const std::vector<RationalComplex>& poles, long prec) {
    for (const auto& p : poles) {
        if (!fmpq_is_zero(p.im())) continue;
        int sgn = fmpq_sgn(p.re());
        if (sgn < 0) continue;
        fmpq_t one;
        fmpq_init(one);
        fmpq_set_si(one, 1, 1);
        int cmp = fmpq_cmp(p.re(), one);
        fmpq_clear(one);
        if (cmp <= 0) {
            return {};
        }
    }

    std::vector<RationalComplex> run;
    run.push_back(RationalComplex::from_fractions(0, 1, 0, 1));
    long max_steps = numeric::run_length();

    for (;;) {
        const RationalComplex& current = run.back();
        fmpq_t one;
        fmpq_init(one);
        fmpq_set_si(one, 1, 1);
        int cmp = fmpq_cmp(current.re(), one);
        fmpq_clear(one);
        if (cmp >= 0) break;

        AcbValue d = min_distance_to_poles(current, poles, prec);
        AcbValue step;
        arb_div_ui(acb_realref(step.raw()), acb_realref(d.raw()),
                   static_cast<unsigned long>(numeric::run_radius()), prec);
        arb_zero(acb_imagref(step.raw()));

        AcbValue current_real = current.to_acb(prec);
        AcbValue cand;
        acb_add(cand.raw(), current_real.raw(), step.raw(), prec);
        fmpq_t qcand;
        fmpq_init(qcand);
        if (!acb_real_to_fmpq(qcand, cand.raw(), numeric::rationalize_pre(), prec)) {
            fmpq_clear(qcand);
            return {};
        }
        fmpq_t one2;
        fmpq_init(one2);
        fmpq_set_si(one2, 1, 1);
        if (fmpq_cmp(qcand, one2) > 0) fmpq_set(qcand, one2);
        fmpq_clear(one2);

        RationalComplex next;
        fmpq_set(next.re(), qcand);
        fmpq_clear(qcand);

        run.push_back(std::move(next));
        if (static_cast<long>(run.size()) > max_steps) {
            return {};
        }
    }
    return run;
}

std::vector<RationalComplex>
run_segment(const std::vector<RationalComplex>& poles,
            const RationalComplex& ini,
            const RationalComplex& fin,
            long prec) {
    RationalComplex span = rc_sub(fin, ini);
    if (span.is_zero()) {
        std::vector<RationalComplex> out;
        out.push_back(ini);
        return out;
    }
    std::vector<RationalComplex> mapped;
    mapped.reserve(poles.size());
    for (const auto& p : poles) {
        RationalComplex d = rc_sub(p, ini);
        mapped.push_back(rc_div(d, span));
    }
    auto unit = run_unit(mapped, prec);
    if (unit.empty()) return {};
    std::vector<RationalComplex> out;
    out.reserve(unit.size());
    for (const auto& u : unit) {
        RationalComplex scaled = rc_mul(u, span);
        out.push_back(rc_add(scaled, ini));
    }
    return out;
}

std::vector<RationalComplex>
run_eta_direction_positive(const std::vector<RationalComplex>& poles,
                            long prec) {
    auto ini = first_step(poles, prec);
    auto fin = last_step(poles, prec);
    if (ini == fin) {
        std::vector<RationalComplex> out;
        out.push_back(ini);
        return out;
    }
    return run_segment(poles, ini, fin, prec);
}

std::vector<RationalComplex>
run_eta_direction(const std::vector<RationalComplex>& poles,
                   const RationalComplex& direction,
                   long prec) {
    if (direction.is_zero()) {
        throw std::invalid_argument("run_eta_direction: zero direction");
    }
    fmpq_t mag2;
    fmpq_init(mag2);
    rc_abs2(mag2, direction);
    AcbValue dir_acb = direction.to_acb(prec);
    AcbValue mag_acb;
    acb_abs(acb_realref(mag_acb.raw()), dir_acb.raw(), prec);
    arb_zero(acb_imagref(mag_acb.raw()));
    AcbValue dir_unit;
    acb_div(dir_unit.raw(), dir_acb.raw(), mag_acb.raw(), prec);
    fmpq_clear(mag2);

    fmpq_t qre, qim;
    fmpq_init(qre); fmpq_init(qim);
    AcbValue tmp_re_only, tmp_im_only;
    arb_set(acb_realref(tmp_re_only.raw()), acb_realref(dir_unit.raw()));
    arb_zero(acb_imagref(tmp_re_only.raw()));
    arb_set(acb_realref(tmp_im_only.raw()), acb_imagref(dir_unit.raw()));
    arb_zero(acb_imagref(tmp_im_only.raw()));
    if (!acb_real_to_fmpq(qre, tmp_re_only.raw(), numeric::rationalize_pre(), prec)
        || !acb_real_to_fmpq(qim, tmp_im_only.raw(), numeric::rationalize_pre(), prec)) {
        fmpq_clear(qre); fmpq_clear(qim);
        throw std::runtime_error("run_eta_direction: failed to rationalise unit direction");
    }
    RationalComplex dir;
    fmpq_set(dir.re(), qre);
    fmpq_set(dir.im(), qim);
    fmpq_clear(qre); fmpq_clear(qim);

    RationalComplex one = RationalComplex::from_fractions(1, 1, 0, 1);
    RationalComplex inv_dir = rc_div(one, dir);
    std::vector<RationalComplex> rotated;
    rotated.reserve(poles.size());
    for (const auto& p : poles) rotated.push_back(rc_mul(p, inv_dir));

    auto run = run_eta_direction_positive(rotated, prec);

    std::vector<RationalComplex> out;
    out.reserve(run.size());
    for (const auto& r : run) out.push_back(rc_mul(r, dir));
    return out;
}

std::vector<RationalComplex>
run_eta_direction(const std::vector<RationalComplex>& poles,
                   numeric::RunningOptions::Direction mode,
                   long prec) {
    long n_cand = numeric::run_candidate();
    std::vector<std::vector<RationalComplex>> candidates;
    candidates.reserve(2 * n_cand + 1);

    auto sample = [&](RationalComplex (*build)(long, long)) {
        for (long j = 0; j <= n_cand; ++j) {
            candidates.push_back(run_eta_direction(poles, build( j, n_cand), prec));
            if (j != 0) {
                candidates.push_back(run_eta_direction(poles, build(-j, n_cand), prec));
            }
        }
    };

    switch (mode) {
        case numeric::RunningOptions::Direction::Re:
            sample([](long j, long n) {
                return RationalComplex::from_fractions(1, 1, j, n);
            });
            break;
        case numeric::RunningOptions::Direction::Im:
            sample([](long j, long n) {
                return RationalComplex::from_fractions(-j, n, 1, 1);
            });
            break;
        case numeric::RunningOptions::Direction::NegRe:
            sample([](long j, long n) {
                return RationalComplex::from_fractions(-1, 1, -j, n);
            });
            break;
        case numeric::RunningOptions::Direction::NegIm:
            sample([](long j, long n) {
                return RationalComplex::from_fractions(j, n, -1, 1);
            });
            break;
        case numeric::RunningOptions::Direction::Custom:
            return {};
    }

    std::vector<RationalComplex>* best = nullptr;
    for (auto& c : candidates) {
        if (c.empty()) continue;
        if (best == nullptr || c.size() < best->size()) best = &c;
    }
    if (best == nullptr) return {};
    return std::move(*best);
}

std::vector<RationalComplex>
run_eta(const std::vector<RationalComplex>& poles, long prec) {
    return run_eta_direction(poles, numeric::run_direction(), prec);
}

}  // namespace amflow::ode
