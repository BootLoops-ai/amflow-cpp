#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz; code
# written by Claude (Anthropic) under his supervision.
"""Independent high-precision anchor for the 2-loop equal-mass sunrise.

J(p^2) = int d^dl1 d^dl2 / (i pi^{d/2})^2  1/((l1^2-m^2)(l2^2-m^2)((l1+l2-p)^2-m^2))
with d = 4-2eps, at Euclidean p^2 = -P^2 < 0 (so J is real).

Position-space (Bessel-moment) representation, exact in d:

    J(p^2=-P^2) = -2^d P^{-nu} m^{3nu} I(eps),    nu = d/2-1 = 1-eps,
    I(eps) = int_0^inf dx x^{3-d} J_nu(P x) K_nu(m x)^3 .

Derivation: Wick rotation (measure d^dl/(i pi^{d/2}) -> d^dl_E/pi^{d/2},
each of the 3 propagators contributes -1), Fourier representation of the
Euclidean propagator int d^dq e^{-iqx}/(q^2+m^2) = (2pi)^{d/2} x^{1-d/2}
m^{d/2-1} K_{d/2-1}(m x), convolution, d-dimensional angular integral.
References for the x-space / Bessel-moment route: S. Groote, J. G. Koerner,
A. A. Pivovarov, "On the evaluation of a certain class of Feynman diagrams in
x-space: sunrise-type topologies at any loop order", Annals Phys. 322 (2007)
2374 [hep-ph/0506286]; D. H. Bailey, J. M. Borwein, D. Broadhurst,
M. L. Glasser, "Elliptic integral evaluations of Bessel moments", J. Phys. A
41 (2008) 205203 [arXiv:0801.0891].

Numerics: split at x=1.
  [1,inf): direct quadrature (K_nu^3 ~ e^{-3mx} decay).
  [0,1]:   expand K_nu = pi/(2 sin(pi nu)) (I_{-nu} - I_nu) and J_nu in
           their everywhere-convergent ascending series; multiply series;
           integrate term-by-term:  int_0^1 x^beta dx = 1/(beta+1)
           (this IS the analytic continuation past the x->0 divergences
           of the d->4 limit; valid since the series converge on [0,1]).

Laurent coefficients in eps via Cauchy integrals on |eps| = 1/8
(nearest other poles of J(eps) at eps = +-1/2).

Self-checks: (a) digit stability under K/M/dps refinement is reported;
(b) compare against pySecDec / the engine externally.
"""
import argparse, json
from mpmath import mp, mpf, mpc, gamma, pi, sin, quad, exp, power, factorial, besselj, besselk, mpmathify


def series_I(alpha, z2_4_terms, K):
    """Ascending series coefficients of I_alpha(z) / (z/2)^alpha :
       sum_j (z^2/4)^j / (j! Gamma(alpha+j+1)); returns coeff list in u=(z/2)^2."""
    return [1 / (factorial(j) * gamma(alpha + j + 1)) for j in range(K + 1)]


def mul_u(a, b, K):
    return [sum(a[i] * b[k - i] for i in range(max(0, k - len(b) + 1), min(k, len(a) - 1) + 1))
            for k in range(K + 1)]


def I_small_x(eps, P, m, K):
    """int_0^1 x^{3-d} J_nu(Px) K_nu(mx)^3 dx by series, analytic continuation."""
    nu = 1 - eps
    d = 4 - 2 * eps
    c = pi / (2 * sin(pi * nu))
    # J_nu(Px) = (P/2)^nu x^nu * SJ(x),  SJ in u=(Px/2)^2 with alternating signs
    SJ = [(-1) ** j / (factorial(j) * gamma(nu + j + 1)) for j in range(K + 1)]
    # I_{+nu}(mx) = (m/2)^nu x^nu SP(x), I_{-nu} = (m/2)^{-nu} x^{-nu} SM(x)
    SP = series_I(nu, None, K)
    SM = series_I(-nu, None, K)
    total = mpc(0)
    # K_nu^3 = c^3 (I_- - I_+)^3 = c^3 sum_{r=0..3} C(3,r) (-1)^r I_-^{3-r} I_+^r
    from math import comb
    for r in range(4):
        # x-power: (3-d) + nu*( (r) - (3-r) ) + nu  [J factor]  = 3-d + nu(2r-2)
        base = 3 - d + nu * (2 * r - 2)
        # u-series: SM^{3-r} * SP^{r} * SJ, with u_m=(m x/2)^2 for I's, u_P=(Px/2)^2 for J
        # combine into series in x^2: coeff_k x^{2k}; prefactors (m/2)^{nu(r-(3-r))} (P/2)^nu
        ser = [mpc(1)] + [mpc(0)] * K
        for _ in range(3 - r):
            ser = mul_u(ser, [SM[j] * (m / 2) ** (2 * j) for j in range(K + 1)], K)
        for _ in range(r):
            ser = mul_u(ser, [SP[j] * (m / 2) ** (2 * j) for j in range(K + 1)], K)
        ser = mul_u(ser, [SJ[j] * (P / 2) ** (2 * j) for j in range(K + 1)], K)
        pref = c ** 3 * comb(3, r) * (-1) ** r \
            * power(m / 2, nu * (r - (3 - r))) * power(P / 2, nu)
        for k in range(K + 1):
            total += pref * ser[k] / (base + 2 * k + 1)   # int_0^1 x^{base+2k} dx
    return total


def I_large_x(eps, P, m):
    nu = 1 - eps
    d = 4 - 2 * eps

    def f(x):
        return power(x, 3 - d) * besselj(nu, P * x) * besselk(nu, m * x) ** 3
    # integrand decays e^{-3mx}; integrate [1, U] with U s.t. e^{-3mU} << 10^-dps
    U = mpf(mp.dps + 30) * mpf(2.303) / (3 * m) + 10
    return quad(f, [1, 4, 16, U])


def sunrise(eps, P, m):
    nu = 1 - eps
    d = 4 - 2 * eps
    I = I_small_x(eps, P, m, sunrise.K) + I_large_x(eps, P, m)
    return -power(2, d) * power(P, -nu) * power(m, 3 * nu) * I


sunrise.K = 60


def _circle_point(work):
    """Pool worker: one Cauchy-circle evaluation. work=(j,dps,P,m,K,rho,M)."""
    j, dps, P_s, m_s, K, rho_s, M = work
    mp.dps = dps
    sunrise.K = K
    th = 2 * pi * j / M
    e = mpmathify(rho_s) * exp(1j * th)
    v = sunrise(e, mpmathify(P_s), mpmathify(m_s))
    return j, str(v.real), str(v.imag)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dps", type=int, default=120)
    ap.add_argument("--orders", type=int, default=9, help="orders above eps^-2")
    ap.add_argument("--P2", type=str, default="1", help="Euclidean P^2 = -p^2 (>0)")
    ap.add_argument("--msq", type=str, default="1")
    ap.add_argument("--eps", type=str, default="1/1000")
    ap.add_argument("--K", type=int, default=0, help="series order (0=auto)")
    ap.add_argument("--M", type=int, default=0, help="Cauchy points (0=auto)")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    mp.dps = args.dps + 40
    P = mp.sqrt(mpmathify(args.P2.replace("/", "./") if "/" in args.P2 else args.P2))
    m = mp.sqrt(mpmathify(args.msq.replace("/", "./") if "/" in args.msq else args.msq))
    sunrise.K = args.K if args.K else max(60, int(1.1 * mp.dps / 2.4) + 40)

    n = args.orders
    rho = mpf(1) / 8
    M = args.M if args.M else max(64, 8 * (n + 6))
    # f(eps) values on the circle -- embarrassingly parallel
    import multiprocessing as mp_proc
    work = [(j, mp.dps, str(P), str(m), sunrise.K, str(rho), M) for j in range(M)]
    with mp_proc.Pool(min(M, max(4, mp_proc.cpu_count() - 6))) as pool:
        res = pool.map(_circle_point, work)
    fv = [None] * M
    for j, vr, vi in res:
        fv[j] = mpc(mpmathify(vr), mpmathify(vi))
    print(f"  {M} circle points done", flush=True)
    coeffs = {}
    for order in range(-2, n - 1):
        cval = mpc(0)
        for j in range(M):
            th = 2 * pi * j / M
            e = rho * exp(1j * th)
            cval += fv[j] * power(e, -order) / M
        coeffs[order] = cval

    # convergence diagnostic: recompute with M/2 points
    diag = {}
    for order in (-2, 0, n - 3):
        cval2 = mpc(0)
        for j in range(0, M, 2):
            th = 2 * pi * j / M
            e = rho * exp(1j * th)
            cval2 += fv[j] * power(e, -order) / (M // 2)
        diag[str(order)] = mp.nstr(abs(cval2 - coeffs[order]), 5)

    p, q = (args.eps.split("/") + ["1"])[:2]
    epsv = mpf(int(p)) / mpf(int(q))
    at_eps = sunrise(epsv, P, m)

    def fmt(x):
        return {"re": mp.nstr(mp.re(x), args.dps), "im": mp.nstr(mp.im(x), args.dps)}

    out = {
        "case": f"sunrise equal-mass p^2=-{args.P2}, msq={args.msq}",
        "convention": "d^dl/(i pi^{d/2}) per loop, no exp(eps gammaE)",
        "dps": args.dps, "series_K": sunrise.K, "cauchy_M": M,
        "orders": {str(k): fmt(v) for k, v in coeffs.items()},
        f"at_eps={args.eps}": fmt(at_eps),
        "halved_M_diff": diag,
    }
    with open(args.out, "w") as f:
        json.dump(out, f, indent=1)
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
