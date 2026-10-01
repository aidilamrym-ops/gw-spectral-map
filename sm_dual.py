#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Independent float64 engine for the Guinand-Weil block tau = W02 - WR - Wp.

Engine 2 (cheap tier) of the SPECTRAL_MAP_CERT two-tier design.

WHAT MAKES IT INDEPENDENT
  * Closed forms (digamma/trigamma at 1/4 + i*pi*n/L, geometric sums) are
    recomputed here with mpmath; matrix assembly and eigenvalue solving run
    in plain float64 + numpy/LAPACK.  Different arithmetic, different
    eigen-solver, different code path than the Arb implementation in
    GUINAND_WEIL/source_arb_ldlt_certify.py.  The shared element is the
    DEFINITION (the formulas), which is the spec both engines must satisfy.
  * Entry-level agreement Arb-vs-mpmath (1e-60) is separately established by
    the source module's own selftest; this module adds the eigen-level
    cross-check.

WHAT IT IS NOT
  * float64 results are CANDIDATES, never certificates.  They place Sturm
    shifts mu inside spectral gaps and cross-check certified counts.
    Every number that ends up in a claim must come from Arb (route S count
    or route E arb_mat.eig enclosure).

Scope: entry -> float64 rounding makes the matrix accurate to roughly
1e-13..1e-15 relative; eigvalsh adds LAPACK backward error.  Shift midpoints
derived from these eigenvalues must sit strictly inside a gap --
bracket_driver refuses mu within its ambiguity margin of any candidate
eigenvalue.

Formulas mirror the specification implemented in source_arb_ldlt_certify.py
(build_arb_tau, arb_closed_forms, _geom_sums, arb_J, arb_kappa).
"""
import math
import os

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

import mpmath as mp  # noqa: E402
import numpy as np  # noqa: E402


def prime_powers_up_to(c):
    """[(q, p)] for prime powers q = p^a <= c (same spec as source)."""
    primes = []
    x = 2
    while x <= c:
        if all(x % p for p in primes):
            primes.append(x)
        x += 1
    out = []
    for p in primes:
        q = p
        while q <= c:
            out.append((q, p))
            q *= p
    return out


def closed_forms(N, c, dps=60):
    """S, CC, XC (index 0..N) + L, kappa, J, pref02 in float64.

    mpmath computes the transcendental closed forms at `dps` decimal digits;
    results are cast to float64 afterwards (target accuracy ~1e-13).
    """
    with mp.workdps(dps):
        L = mp.log(c)
        PI = mp.pi
        quarter = mp.mpf(1) / 4
        psi_quarter = mp.digamma(quarter)
        thr = mp.mpf(10) ** (-(dps + 10))

        S = [0.0] * (N + 1)
        CC = [0.0] * (N + 1)
        XC = [0.0] * (N + 1)
        for n in range(N + 1):
            w = 2 * PI * n / L
            z = quarter + 1j * PI * n / L
            psi = mp.digamma(z)
            psi1 = mp.polygamma(1, z)
            # geometric sums (no ball widening needed: truncation << 1e-15)
            gS = gCC = gX1 = gX2 = mp.mpf(0)
            k = 0
            while True:
                c_k = mp.mpf(2 * k) + mp.mpf("0.5")
                e = mp.exp(-c_k * L)
                den = c_k * c_k + w * w
                gS += e / den
                if n != 0:
                    gCC += e * w * w / (c_k * den)
                gX1 += e * c_k / den
                gX2 += e * (c_k * c_k - w * w) / (den * den)
                if e < thr and k > 2:
                    break
                k += 1
            if n == 0:
                S[n] = 0.0
                CC[n] = 0.0
            else:
                S[n] = float(mp.im(psi) / 2 - w * gS)
                CC[n] = float(-(mp.re(psi) - psi_quarter) / 2 + gCC)
            XC[n] = float(mp.re(psi1) / 4 - L * gX1 - gX2)

        kappa = float(mp.log(4 * PI * (mp.e ** L - 1) / (mp.e ** L + 1))
                      + mp.euler)
        U = mp.exp(L / 2)
        J = float(-2 * mp.log(U + 1) + mp.log(U * U + 1)
                  + 2 * mp.atan(U) + mp.log(2) - PI / 2)
        pref02 = float(32 * L * mp.sinh(L / 4) ** 2)
        Lf = float(L)
    return dict(S=S, CC=CC, XC=XC, L=Lf, kappa=kappa, J=J, pref02=pref02,
                dps=dps)


def build_tau_f64(c, N, dps=60):
    """tau = W02 - WR - Wp as a float64 symmetric numpy matrix (DIM=2N+1)."""
    cf = closed_forms(N, c, dps)
    S, CC, XC = cf["S"], cf["CC"], cf["XC"]
    L, kappa, J, pref02 = cf["L"], cf["kappa"], cf["J"], cf["pref02"]

    pdata = prime_powers_up_to(c)
    weights = [math.log(p) * (q ** -0.5) for (q, p) in pdata]
    positions = [math.log(q) for (q, p) in pdata]

    PI = math.pi
    sp2 = 16 * PI * PI
    l2 = L * L
    DIM = 2 * N + 1
    M = np.empty((DIM, DIM), dtype=np.float64)

    def S_signed(nn):
        return S[nn] if nn >= 0 else -S[-nn]

    for i in range(DIM):
        n = i - N
        for j in range(i, DIM):
            m = j - N
            num = l2 - sp2 * m * n
            den = (l2 + sp2 * m * m) * (l2 + sp2 * n * n)
            W02 = pref02 * num / den
            if n == m:
                WR = kappa + 2 * CC[abs(n)] + J - (2.0 / L) * XC[abs(n)]
            else:
                WR = (S_signed(m) - S_signed(n)) / (PI * (n - m))
            Wp = 0.0
            for idx in range(len(weights)):
                y = positions[idx]
                if n == m:
                    q = 2.0 * (1.0 - y / L) * math.cos(2 * PI * n * y / L)
                else:
                    q = ((math.sin(2 * PI * m * y / L)
                          - math.sin(2 * PI * n * y / L)) / (PI * (n - m)))
                Wp += weights[idx] * q
            val = W02 - WR - Wp
            M[i, j] = val
            M[j, i] = val
    return M


def dual_eigenvalues(c, N, dps=60):
    """Sorted (ascending) float64 eigenvalues of tau(c, N)."""
    M = build_tau_f64(c, N, dps)
    return np.linalg.eigvalsh(M)


def w02_block(c, N, dps=60):
    """The W02 term alone as a matrix (for the rank-2 structural selftest)."""
    cf = closed_forms(N, c, dps)
    L, pref02 = cf["L"], cf["pref02"]
    PI = math.pi
    sp2 = 16 * PI * PI
    l2 = L * L
    DIM = 2 * N + 1
    W = np.empty((DIM, DIM), dtype=np.float64)
    for i in range(DIM):
        n = i - N
        for j in range(DIM):
            m = j - N
            W[i, j] = pref02 * (l2 - sp2 * m * n) / (
                (l2 + sp2 * m * m) * (l2 + sp2 * n * n))
    return W


def selftest():
    """Structural check: the W02 block is separable of rank 2.

    Algebra: (L^2 - 16pi^2 mn)/((L^2+16pi^2 m^2)(L^2+16pi^2 n^2))
           = L^2 * t1(m) t1(n) - 16pi^2 * t2(m) t2(n),
    t1(x) = 1/(L^2+16pi^2 x^2), t2(x) = x/(L^2+16pi^2 x^2).
    So every singular value beyond the 2nd must vanish to float64 noise.
    """
    N = 20
    W = w02_block(100, N)
    sv = np.linalg.svd(W, compute_uv=False)
    ratio = sv[2] / sv[0]
    ok = ratio < 1e-13
    print("[sm_dual] SELFTEST W02-rank2: sv[0]=%.6e sv[1]=%.6e sv[2]=%.6e "
          "ratio=%.3e -> %s" % (sv[0], sv[1], sv[2], ratio,
                                "PASS" if ok else "FAIL"), flush=True)
    return ok


if __name__ == "__main__":
    import sys
    sys.exit(0 if selftest() else 1)
