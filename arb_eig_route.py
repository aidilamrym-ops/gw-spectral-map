#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Route E: full rigorous spectrum via arb_mat.eig on a NATIVE-Arb matrix.

WHY ROUTE E CAN SUCCEED WHERE OMEGA v1 FAILED
  v1 built its matrix in mpmath and handed the digits to arb as DECIMAL
  TEXT: (1) slow -- arb_mat.eig on 801x801 @2851 bits ran 94 minutes
  without finishing and was killed; (2) the mpmath origin of the digits was
  NOT enclosed, a real rigour caveat.  Here build_arb_tau computes every
  entry as a genuine Arb ball end-to-end, so eig enclosures are rigorous
  FOR THIS MATRIX.

WHEN ROUTE E IS THE WRONG TOOL
  On radius-carrying matrices the isolation algorithms can raise: precedent
  recorded in GUINAND_WEIL/gw_arb_measured.py -- rump and vdhoeven_mourrain
  both NOT ISOLATED at prec 1024 and 4096 on the (100,40) matrix.
  Exit code 2 = route E unavailable at this size/precision -> fall back to
  route S (sturm_count), which never separates eigenvalues, only counts
  them.  That is an honest tool limitation, never massaged away.

EXIT CODES
  0 = all DIM eigenvalues isolated, strictly positive (consistent pos-def)
  1 = ran but isolation incomplete / enclosures contain zero (fail-closed)
  2 = eig raised for every algorithm (route E unavailable here)
"""
import argparse
import datetime
import json
import os
import platform
import sys
import time

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from sturm_count import lower_priority, load_source, ts, RESULTS_DIR  # noqa: E402


def run_route_e(c, N, prec, algs):
    src = load_source()
    print("[%s] ROUTE E build c=%d N=%d prec=%d ..." % (ts(), c, N, prec),
          flush=True)
    t0 = time.time()
    A, DIM = src.build_arb_tau(c, N, prec)
    t_build = time.time() - t0
    print("[%s] built DIM=%d in %.1f s" % (ts(), DIM, t_build), flush=True)

    results = []
    exit_code = 2
    for alg in algs:
        rec = dict(algorithm=alg)
        t1 = time.time()
        try:
            ev = A.eig(algorithm=alg)
        except Exception as exc:
            rec.update(raised=True,
                       error="%s: %s" % (type(exc).__name__, str(exc)[:200]),
                       seconds=round(time.time() - t1, 1))
            print("  %-22s RAISED after %.1f s: %s: %s"
                  % (alg, rec["seconds"], type(exc).__name__,
                     str(exc)[:120]), flush=True)
            results.append(rec)
            continue

        n = len(ev)
        n_contains_zero = 0
        n_strict_pos = 0
        lam_min_s = lam_max_s = ""
        try:
            reals = [w.real for w in ev]
            for r in reals:
                if bool(r.contains(0)):
                    n_contains_zero += 1
                if bool(r > 0):
                    n_strict_pos += 1
            lo = min(reals, key=lambda r: float(r.lower()))
            hi = max(reals, key=lambda r: float(r.upper()))
            lam_min_s = str(lo)[:120]
            lam_max_s = str(hi)[:120]
        except Exception as exc:
            rec["stats_error"] = "%s: %s" % (type(exc).__name__, str(exc)[:200])

        rec.update(raised=False, seconds=round(time.time() - t1, 1),
                   n_returned=int(n), dimension=int(DIM),
                   n_contains_zero=int(n_contains_zero),
                   n_strictly_positive=int(n_strict_pos),
                   lambda_min=lam_min_s, lambda_max=lam_max_s)
        fully = (n == DIM and n_contains_zero == 0 and n_strict_pos == DIM)
        rec["fully_isolated_positive_definite"] = bool(fully)
        print("  %-22s %.1f s: %d/%d isolated, contains0=%d, strictPos=%d"
              % (alg, rec["seconds"], n, DIM, n_contains_zero, n_strict_pos),
              flush=True)
        results.append(rec)
        if fully:
            exit_code = 0
        elif exit_code == 2:
            exit_code = 1

    out_rec = dict(stage="g2_route_e", c=c, N=N, prec_bits=prec,
                   dimension=DIM, build_seconds=round(t_build, 1),
                   date=datetime.datetime.now().isoformat(timespec="seconds"),
                   python_version=platform.python_version(),
                   script=os.path.basename(__file__),
                   algorithms=results)
    return exit_code, out_rec


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--c", type=int, default=100)
    p.add_argument("--N", type=int, default=200)
    p.add_argument("--prec", type=int, default=2000)
    p.add_argument("--algs", type=str, default="rump,vdhoeven_mourrain")
    p.add_argument("--json-out", type=str, default="")
    p.add_argument("--lowprio", action="store_true")
    args = p.parse_args()

    if args.lowprio:
        print("[%s] priority BelowNormal: %s" % (ts(), lower_priority()),
              flush=True)

    code, rec = run_route_e(args.c, args.N, args.prec,
                            [a.strip() for a in args.algs.split(",") if a.strip()])
    os.makedirs(RESULTS_DIR, exist_ok=True)
    out = args.json_out or os.path.join(
        RESULTS_DIR, "g2_route_e_c%d_N%d_p%d.json"
        % (args.c, args.N, args.prec))
    with open(out, "w") as f:
        json.dump(rec, f, indent=2)
    print("wrote %s" % out, flush=True)
    print("[%s] ROUTE E EXIT %d" % (ts(), code), flush=True)
    sys.exit(code)


if __name__ == "__main__":
    main()
