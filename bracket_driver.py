#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Level-0 certified spectral brackets for tau(c, N) -- vetted Sturm counts.

TWO TIERS AT WORK
  dual engine (sm_dual, float64): proposes shift positions and provides the
  cross-check counts.  Arb (source interval LDL^T): the ONLY authority --
  every recorded count is proved, never sampled.

PIPELINE
  1. dual eigenvalues -> planned shifts: gap midpoints + two outer bounds.
  2. VETTING: a gap midpoint whose margin to the nearest dual eigenvalue is
     below ambig (= ambig_rel * spectral scale) is flagged TOO_TIGHT and no
     count is attempted there.  Float64 cannot place a shift inside a gap it
     cannot resolve; the flag is a statement about the CHEAP ENGINE ONLY --
     Arb is not implicated.  (Small-c matrices carry eigenvalues around
     1e-16 and this family's lambda_min can reach ~1e-102, so near-zero
     clusters are expected, not exceptional.)
  3. One matrix build; one certified count per planned shift; JSONL
     checkpoint per count -> power-loss safe resume.
  4. FAIL-CLOSED cross-check on every counted shift: certified n_neg must
     equal the dual engine's expected count, else CONFLICT (recorded, no
     bracket certified through it, exit 1).
  5. SEGMENTATION into a contiguous partition of (outer_lo, outer_hi]:
       * two adjacent clean counts differing by exactly 1 -> INDIVIDUAL
         bracket (lambda_i in (lo, hi]), dual lambda cross-checked inside;
       * a run of TOO_TIGHT gaps between two clean counts -> CLUSTER
         bracket: Arb proves exactly k = count(hi) - count(lo) eigenvalues
         live in (lo, hi]; the dual engine must agree on k (boundary margins
         are clean there, so the dual count is checkable even inside).
  6. Exit 0 iff the segments cover ALL DIM eigenvalues with zero conflicts
     and zero undetermined pivots.  Cluster brackets are honest level-0
     certificates ("k eigenvalues, width W"); refining them by Sturm
     bisection inside the zone is LEVEL-1 work, queued until OMEGA is
     released.

EXIT CODES
  0 = full coverage certified, zero conflicts, zero undetermined
  1 = conflict, coverage gap, or cluster/dual count mismatch (recorded)
  2 = undetermined pivot (fail-closed; raise --prec)

The full campaign at c=100, N=200 is QUEUED until OMEGA is released
(single-thread rule).  --selftest runs the whole pipeline at c=13, N=8.
"""
import argparse
import datetime
import json
import os
import sys
import time

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from sturm_count import (RESULTS_DIR, count_shifted, load_source,  # noqa: E402
                         lower_priority, ts)
import sm_dual  # noqa: E402


def mu_key(mu):
    return "%.17g" % mu


def plan_shifts(evals, ambig_rel=1e-6):
    """Planned shifts with expected dual counts (0-based ordering)."""
    DIM = len(evals)
    scale = max(1.0, abs(float(evals[0])), abs(float(evals[-1])))
    ambig = ambig_rel * scale
    shifts = []
    gap0 = float(evals[1] - evals[0]) if DIM > 1 else scale
    mu_lo = float(evals[0] - max(abs(gap0), 10 * ambig))
    shifts.append(dict(kind="outer_lo", mu=mu_lo, expected=0))
    for i in range(DIM - 1):
        mu = 0.5 * float(evals[i] + evals[i + 1])
        margin = min(mu - float(evals[i]), float(evals[i + 1]) - mu)
        sh = dict(kind="gap", index=i, mu=mu, expected=i + 1, margin=margin)
        if margin < ambig:
            sh["flag"] = "TOO_TIGHT"
        shifts.append(sh)
    gap_last = float(evals[-1] - evals[-2]) if DIM > 1 else scale
    mu_hi = float(evals[-1] + max(abs(gap_last), 10 * ambig))
    shifts.append(dict(kind="outer_hi", mu=mu_hi, expected=DIM))
    return shifts, scale, ambig


def load_done(out_path):
    done = {}
    if os.path.exists(out_path):
        with open(out_path, "r") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                try:
                    rec = json.loads(line)
                    done[rec["key"]] = rec
                except Exception:
                    continue  # torn line after power loss: redo that shift
    return done


def segment(shifts, done, evals, tol):
    """Contiguous partition of (outer_lo, outer_hi] into certified segments.

    Returns (segments, problems).  Each segment: lo, hi, k, kind
    ("individual"|"cluster"), optionally dual_lambda (individual).
    """
    segs, problems = [], []
    prev = None            # last shift with a certified count
    zone_lo = None         # boundary shift of a pending TOO_TIGHT run
    cum = 0                # eigenvalues already covered by segs
    for sh in shifts:
        if sh.get("flag") == "TOO_TIGHT":
            if zone_lo is None:
                zone_lo = prev
                if zone_lo is None:
                    problems.append("TOO_TIGHT run starts before any "
                                    "certified count")
                    return segs, problems
            continue
        rec = done.get(mu_key(sh["mu"]))
        if rec is None or rec.get("status") != "OK":
            break          # conflict or uncounted beyond this point
        if zone_lo is not None:
            k = rec["n_neg"] - zone_lo["n_neg"]
            lo, hi = zone_lo["mu"], rec["mu"]
            n_dual = int(((evals >= lo) & (evals < hi)).sum())
            if k < 1:
                problems.append("cluster zone (mu=%r..%r) certified k=%d"
                                % (lo, hi, k))
            elif n_dual != k:
                problems.append("cluster zone mu=%r..%r: Arb k=%d != dual "
                                "count %d" % (lo, hi, k, n_dual))
            segs.append(dict(kind="cluster", lo=lo, hi=hi, k=int(k),
                             lo_count=int(zone_lo["n_neg"]),
                             hi_count=int(rec["n_neg"])))
            cum += int(k)
            zone_lo = None
        elif prev is not None:
            k = rec["n_neg"] - prev["n_neg"]
            lo, hi = prev["mu"], rec["mu"]
            if k == 1:
                segs.append(dict(kind="individual", lo=lo, hi=hi, k=1,
                                 lo_count=int(prev["n_neg"]),
                                 hi_count=int(rec["n_neg"])))
                cum += 1
            else:
                problems.append("adjacent clean counts differ by %d "
                                "(mu=%r..%r), expected 1 without a zone"
                                % (k, lo, hi))
                break
        prev = dict(mu=sh["mu"], n_neg=rec["n_neg"])
    # cross-checks: dual eigenvalues inside their segments
    for s in segs:
        if s["kind"] == "individual" and "lo_count" in s:
            idx = s["lo_count"]                      # 0-based eigenvalue
            lam = float(evals[idx])
            if not (s["lo"] - tol <= lam <= s["hi"] + tol):
                problems.append("dual lambda[%d]=%r outside certified "
                                "bracket (%r, %r]" % (idx, lam, s["lo"],
                                                      s["hi"]))
            s["dual_lambda"] = lam
    return segs, problems


def run(c, N, prec, out_path, ambig_rel, lowprio=False, heartbeat_every=10):
    if lowprio:
        print("[%s] priority BelowNormal: %s" % (ts(), lower_priority()),
              flush=True)

    print("[%s] PLAN dual engine c=%d N=%d ..." % (ts(), c, N), flush=True)
    evals = sm_dual.dual_eigenvalues(c, N)
    DIM = len(evals)
    shifts, scale, ambig = plan_shifts(evals, ambig_rel)
    tol = 1e-9 * scale
    n_too = sum(1 for s in shifts if s.get("flag") == "TOO_TIGHT")
    print("[%s] DIM=%d, %d planned shifts (%d TOO_TIGHT -> cluster zones)"
          % (ts(), DIM, len(shifts), n_too), flush=True)

    done = load_done(out_path)
    if done:
        print("[%s] resume: %d counts already certified in %s"
              % (ts(), len(done), out_path), flush=True)
    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)

    src = load_source()
    print("[%s] building matrix once (c=%d, N=%d, prec=%d) ..."
          % (ts(), c, N, prec), flush=True)
    t0 = time.time()
    A, DIM_A = src.build_arb_tau(c, N, prec)
    t_build = time.time() - t0
    if DIM_A != DIM:
        print("[%s] FATAL dimension mismatch dual=%d arb=%d"
              % (ts(), DIM, DIM_A), flush=True)
        return 2
    print("[%s] built DIM=%d in %.1f s" % (ts(), DIM, t_build), flush=True)

    conflicts, undetermined = [], None
    t_run0 = time.time()
    n_done_before = len(done)
    with open(out_path, "a") as out_f:
        for k, sh in enumerate(shifts):
            if sh.get("flag") == "TOO_TIGHT":
                continue
            key = mu_key(sh["mu"])
            rec = done.get(key)
            if rec is None or rec.get("status") not in ("OK", "CONFLICT"):
                cnt = count_shifted(src, A, DIM, sh["mu"], heartbeat=0)
                status = "OK"
                if cnt["undetermined"] is not None:
                    status = "UNDETERMINED"
                elif cnt["n_neg"] != sh["expected"]:
                    status = "CONFLICT"
                rec = dict(key=key, mu=sh["mu"], kind=sh["kind"],
                           expected=sh["expected"], status=status,
                           n_neg=cnt["n_neg"], n_pos=cnt["n_pos"],
                           undetermined=cnt["undetermined"],
                           ldlt_seconds=cnt["ldlt_seconds"],
                           pivot_digest_sha256=cnt["pivot_digest_sha256"],
                           date=datetime.datetime.now().isoformat(
                               timespec="seconds"))
                out_f.write(json.dumps(rec) + "\n")
                out_f.flush()
                done[key] = rec
            if rec["status"] == "UNDETERMINED":
                undetermined = rec
                break
            if rec["status"] == "CONFLICT":
                conflicts.append(rec)
            if (k + 1) % heartbeat_every == 0 or k == len(shifts) - 1:
                el = time.time() - t_run0
                print("  [%s] shift %d/%d  elapsed=%.0fs  conflicts=%d"
                      % (ts(), k + 1, len(shifts), el, len(conflicts)),
                      flush=True)

    if undetermined is not None:
        print("[%s] FAIL-CLOSED: undetermined pivot at mu=%r (raise --prec)"
              % (ts(), undetermined["mu"]), flush=True)
        return 2

    segs, problems = segment(shifts, done, evals, tol)
    for p in problems:
        print("  CONFLICT: %s" % p, flush=True)

    covered = sum(s["k"] for s in segs)
    individuals = [s for s in segs if s["kind"] == "individual"]
    clusters = [s for s in segs if s["kind"] == "cluster"]
    for s in individuals:
        loc = float(evals[min(s["lo_count"] + 1, DIM - 1)]
                    - evals[s["lo_count"]]) or 1.0
        s["width"] = s["hi"] - s["lo"]
        s["width_over_gap"] = s["width"] / loc
    for s in clusters:
        s["width"] = s["hi"] - s["lo"]

    full = (covered == DIM and not problems and not conflicts)
    out_json = out_path.replace(".jsonl", "_brackets.json")
    summary = dict(stage="level0_brackets", c=c, N=N, prec_bits=prec,
                   dimension=DIM, build_seconds=round(t_build, 1),
                   planned_shifts=len(shifts),
                   counts_certified=len([r for r in done.values()
                                         if r.get("status") == "OK"]),
                   counts_from_previous_runs=n_done_before,
                   conflicts=len(conflicts) + len(problems),
                   too_tight_zones=n_too,
                   individual_brackets=len(individuals),
                   cluster_brackets=len(clusters),
                   eigenvalues_covered=covered,
                   coverage_of_dimension=covered / float(DIM),
                   all_certified=bool(full),
                   max_width_over_gap=max(
                       [s["width_over_gap"] for s in individuals] or [None]),
                   clusters=[dict(lo=s["lo"], hi=s["hi"], k=s["k"],
                                  width=s["width"]) for s in clusters],
                   date=datetime.datetime.now().isoformat(timespec="seconds"),
                   dual_eigenvalues=[float(x) for x in evals],
                   brackets=segs)
    with open(out_json, "w") as f:
        json.dump(summary, f, indent=2)
    print("[%s] coverage %d/%d  individual=%d  clusters=%d  conflicts=%d "
          "-> %s" % (ts(), covered, DIM, len(individuals), len(clusters),
                     len(conflicts) + len(problems), out_json), flush=True)
    print("[%s] BRACKETS %s" % (ts(), "CERTIFIED" if full else "INCOMPLETE"),
          flush=True)
    return 0 if full else 1


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--c", type=int, default=100)
    p.add_argument("--N", type=int, default=200)
    p.add_argument("--prec", type=int, default=2000)
    p.add_argument("--out", type=str, default="",
                   help="JSONL checkpoint (default results/brackets_...jsonl)")
    p.add_argument("--ambig-rel", type=float, default=1e-6,
                   help="gap vetting margin, relative to spectral scale")
    p.add_argument("--selftest", action="store_true",
                   help="full pipeline at c=13 N=8 prec=300")
    p.add_argument("--lowprio", action="store_true")
    args = p.parse_args()

    c, N, prec = args.c, args.N, args.prec
    if args.selftest:
        c, N, prec = 13, 8, 300
    out = args.out or os.path.join(
        RESULTS_DIR, "brackets_c%d_N%d_p%d.jsonl" % (c, N, prec))
    sys.exit(run(c, N, prec, out, args.ambig_rel, lowprio=args.lowprio))


if __name__ == "__main__":
    main()
