#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Route S: certified spectral counting  N(mu) = #{lambda_i < mu}.

MECHANISM (Sylvester inertia)
    A  real symmetric, mu exact (float64 -> dyadic, arb-exact):
    A' = A - mu*I  has eigenvalues (lambda_i - mu), hence
    N(mu) = n_neg(A') = number of strictly negative eigenvalues of A'.
    n_neg is PROVED by the source module's interval LDL^T
    (certified_inertia): every pivot ball strictly signed -> (n_pos, n_neg)
    proved, not observed.  An undetermined pivot aborts with
    undetermined=<index> -- fail-closed, the count is NEVER guessed.

WHY THIS IS THE SCALABLE RIGOROUS ROUTE
    route E (arb_mat.eig, full spectrum in one call) may be too slow or may
    fail to isolate at a given size/precision (measured precedent: 801x801
    @2851 bits ran 94 min without finishing; rump/vdhoeven_mourrain can
    raise NOT ISOLATED on radius-carrying matrices).  One LDL^T count is the
    same factorization OMEGA already runs, and bracket_driver turns counts
    into certified eigenvalue intervals.

PROVENANCE
    Each count carries the sha256 digest of the pivot transcript (the same
    digest pattern the source module writes into its certificates).

USAGE
    python sturm_count.py --selftest
    python sturm_count.py --bench --c 100 --N 200 --prec 2000 --lowprio
    python sturm_count.py --c 100 --N 200 --prec 2000 --mu 12.5 --json-out x.json
"""
import argparse
import datetime
import hashlib
import importlib.util
import json
import os
import platform
import sys
import time

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

HERE = os.path.dirname(os.path.abspath(__file__))
SOURCE_PATH = os.path.normpath(
    os.path.join(HERE, "..", "GUINAND_WEIL", "source_arb_ldlt_certify.py"))
RESULTS_DIR = os.path.join(HERE, "results")

_src_cache = None


def ts():
    return datetime.datetime.now().strftime("%H:%M:%S")


def lower_priority():
    """BelowNormal priority so an OMEGA sweep is never starved (log-only).

    ctypes x64 detail (measured): without explicit argtypes/restype the
    pseudo-handle of GetCurrentProcess is passed as 32-bit and
    SetPriorityClass silently returns 0.  With c_void_p types it succeeds
    (verified: GetPriorityClass -> 0x8000).
    """
    try:
        import ctypes
        k = ctypes.windll.kernel32
        k.GetCurrentProcess.restype = ctypes.c_void_p
        k.SetPriorityClass.argtypes = [ctypes.c_void_p, ctypes.c_uint]
        k.SetPriorityClass.restype = ctypes.c_int
        BELOW_NORMAL = 0x00008000
        handle = k.GetCurrentProcess()
        ok = k.SetPriorityClass(handle, BELOW_NORMAL)
        return bool(ok)
    except Exception as exc:  # pragma: no cover - platform guard
        print("[warn] could not lower priority: %s" % exc, flush=True)
        return False


def load_source():
    """Import GUINAND_WEIL/source_arb_ldlt_certify.py by path (read-only use)."""
    global _src_cache
    if _src_cache is None:
        if not os.path.exists(SOURCE_PATH):
            raise SystemExit("source module not found: %s" % SOURCE_PATH)
        spec = importlib.util.spec_from_file_location(
            "source_arb_ldlt_certify", SOURCE_PATH)
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        _src_cache = mod
    return _src_cache


def shift_matrix(A, DIM, mu):
    """B = A - mu*I as a fresh arb_mat (mu: python float -> arb, exact dyadic)."""
    from flint import arb, arb_mat
    m = arb(mu)          # python float is dyadic -> exact ball at ctx.prec
    B = arb_mat(DIM, DIM)
    for i in range(DIM):
        for j in range(i, DIM):
            v = A[i, j]
            if i == j:
                v = v - m
            B[i, j] = v
            B[j, i] = v
    return B


def count_shifted(src, A, DIM, mu, heartbeat=0):
    """Certified N(mu) on an already-built matrix A.  Returns a record dict."""
    B = shift_matrix(A, DIM, mu)
    t0 = time.time()
    n_pos, n_neg, undet, transcript = src.certified_inertia(B, DIM, heartbeat)
    t_ldlt = time.time() - t0
    digest = hashlib.sha256("\n".join(transcript).encode()).hexdigest()
    return dict(mu=float(mu),
                n_pos=int(n_pos), n_neg=int(n_neg),
                undetermined=(None if undet is None else int(undet)),
                ldlt_seconds=round(t_ldlt, 2),
                pivot_digest_sha256=digest,
                pivot_count=len(transcript))


def selftest():
    """Two parts, both fail-closed:
    (a) the source module's own entry-agreement selftest (Arb vs mpmath 1e-60);
    (b) certified counts at chosen shifts vs the independent float64 engine.
    """
    c, N, prec = 13, 8, 300
    print("[%s] STURM selftest c=%d N=%d prec=%d" % (ts(), c, N, prec),
          flush=True)
    ok = True

    src = load_source()
    if not src.selftest(c=c, N=N, prec=prec):
        print("  FAIL part (a): source entry selftest", flush=True)
        ok = False
    else:
        print("  PASS part (a): source entry selftest (Arb vs mpmath 1e-60)",
              flush=True)

    import sm_dual
    A, DIM = src.build_arb_tau(c, N, prec)
    evals = sm_dual.dual_eigenvalues(c, N)
    if len(evals) != DIM:
        print("  FAIL part (b): dual eigenvalue count %d != DIM %d"
              % (len(evals), DIM), flush=True)
        return False

    # shift set: below-spectrum, above-spectrum, mu=0, every gap midpoint.
    # Vetting rule (same spirit as bracket_driver): if any dual eigenvalue
    # sits within `ambig` of mu, the float64 placement cannot decide that
    # mu's side of the eigenvalue -> the pair (arb count vs dual count) is
    # NOT checkable there and the point is skipped as AMBIGUOUS, never
    # counted as engine disagreement.  (Small-c matrices carry eigenvalues
    # near 1e-16, far below float64 resolution of this spectrum -- Arb at
    # 300+ bits is the only engine that can decide those signs.)
    scale = max(1.0, abs(float(evals[0])), abs(float(evals[-1])))
    ambig = 1e-6 * scale
    mus = [float(evals[0] - 10 * scale), 0.0]
    for i in range(len(evals) - 1):
        mus.append(0.5 * float(evals[i] + evals[i + 1]))
    mus.append(float(evals[-1] + 10 * scale))

    ok_b = True
    n_match = checked = skipped = 0
    for mu in mus:
        margin = float(abs(evals - mu).min())
        if margin < ambig:
            skipped += 1
            continue
        checked += 1
        expected = int((evals < mu).sum())
        rec = count_shifted(src, A, DIM, mu, heartbeat=0)
        if rec["undetermined"] is not None:
            print("  FAIL part (b): undetermined pivot %s at mu=%r"
                  % (rec["undetermined"], mu), flush=True)
            ok_b = False
            continue
        if rec["n_neg"] != expected:
            print("  FAIL part (b): mu=%r certified n_neg=%d != dual %d"
                  % (mu, rec["n_neg"], expected), flush=True)
            ok_b = False
        else:
            n_match += 1
    if checked < 3:
        print("  FAIL part (b): only %d checkable shifts (non-vacuous "
              "minimum is 3)" % checked, flush=True)
        ok_b = False
    print("  %s part (b): %d/%d checkable certified counts match the float64 "
          "engine (%d skipped as ambiguous < ambig %.2e)"
          % ("PASS" if ok_b else "FAIL", n_match, checked, skipped, ambig),
          flush=True)
    ok = ok and ok_b

    # structural: sm_dual W02-rank2 check rides along (cheap, validates Notes)
    if not sm_dual.selftest():
        ok = False

    print("[%s] STURM SELFTEST %s" % (ts(), "PASS" if ok else "FAIL"),
          flush=True)
    return ok


def np_count_below(evals, mu):
    return int((evals < mu).sum())


def bench(args):
    src = load_source()
    print("[%s] G2 BENCH build c=%d N=%d prec=%d ..." %
          (ts(), args.c, args.N, args.prec), flush=True)
    t0 = time.time()
    A, DIM = src.build_arb_tau(args.c, args.N, args.prec)
    t_build = time.time() - t0
    print("[%s] built DIM=%d in %.1f s; count mu=%r ..." %
          (ts(), DIM, t_build, args.mu), flush=True)
    rec = count_shifted(src, A, DIM, args.mu, heartbeat=0)
    rec.update(stage="g2_baseline", c=args.c, N=args.N, dimension=DIM,
               prec_bits=args.prec, build_seconds=round(t_build, 1),
               date=datetime.datetime.now().isoformat(timespec="seconds"),
               python_version=platform.python_version(),
               script=os.path.basename(__file__))
    undet = rec["undetermined"]
    if undet is None:
        print("[%s] COUNT DONE: n_pos=%d n_neg=%d (ldlt %.1f s, build %.1f s)"
              % (ts(), rec["n_pos"], rec["n_neg"], rec["ldlt_seconds"],
                 t_build), flush=True)
        if args.mu == 0.0 and args.c == 100 and args.N == 200:
            match = (rec["n_pos"] == 401 and rec["n_neg"] == 0)
            print("  vs manuscript certificate (c=100,N=200): expected "
                  "n_pos=401 n_neg=0 -> %s"
                  % ("MATCH" if match else "MISMATCH (fail-closed)"),
                  flush=True)
            rec["manuscript_expected"] = dict(n_pos=401, n_neg=0)
            rec["manuscript_match"] = bool(match)
    else:
        print("[%s] UNDETERMINED at pivot %s -- raise --prec (fail-closed)"
              % (ts(), undet), flush=True)

    os.makedirs(RESULTS_DIR, exist_ok=True)
    out = args.json_out or os.path.join(
        RESULTS_DIR, "g2_sturm_c%d_N%d_p%d_mu%g.json"
        % (args.c, args.N, args.prec, args.mu))
    with open(out, "w") as f:
        json.dump(rec, f, indent=2)
    print("wrote %s" % out, flush=True)
    if undet is not None:
        return 2
    if rec.get("manuscript_match") is False:
        return 1
    return 0


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--c", type=int, default=100)
    p.add_argument("--N", type=int, default=200)
    p.add_argument("--prec", type=int, default=2000)
    p.add_argument("--mu", type=float, default=0.0)
    p.add_argument("--selftest", action="store_true")
    p.add_argument("--bench", action="store_true")
    p.add_argument("--json-out", type=str, default="")
    p.add_argument("--lowprio", action="store_true",
                   help="run BelowNormal (required while OMEGA is alive)")
    args = p.parse_args()

    if args.lowprio:
        lowered = lower_priority()
        print("[%s] priority BelowNormal: %s" % (ts(), lowered), flush=True)

    if args.selftest:
        sys.exit(0 if selftest() else 1)
    sys.exit(bench(args))


if __name__ == "__main__":
    main()
