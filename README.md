# SPECTRAL_MAP_CERT — certified eigenvalue spectrum of the Guinand–Weil block Q_Nc

Tier-3 "clear map" project: turn the certified LDL^T machinery already used by
OMEGA-CORE into a **certified spectral staircase** — rigorous eigenvalue
brackets for the Guinand–Weil matrix `tau = W02 - WR - Wp` — and feed them to
the S1/S2 Montgomery pair-correlation test with honest error bars.

Origin: 2026-10-01 discussion with the Architect ("extreme-precision data ->
a clear map -> hidden patterns / anomalies / spectral regularity along the
critical line").  Background and construction findings are persisted in
[NOTES.md](NOTES.md).

## Status (2026-10-01)

| Gate | Result |
|------|--------|
| `py_compile` (4 files) | OK |
| `sm_dual.py` selftest (W02 rank-2) | PASS — sv2/sv0 = 3.7e-17 |
| `sturm_count.py --selftest` | PASS — part (a) source Arb-vs-mpmath 1e-60; part (b) 10/10 checkable certified counts match float64 engine (9 skipped as ambiguous, vetted) |
| `bracket_driver.py --selftest` | PASS — coverage 17/17, 8 individual + 1 cluster (k=9), 0 conflicts, JSONL resume verified |
| G2 baseline (DIM=401) | **MEASURED — see table below**; precision ladder found the minimum certifiable precision, manuscript certificate (401, 0) reproduced at TWO precisions |
| Full campaign (c=100, N=200) | **QUEUED** — runs only after OMEGA-CORE is released + `gw_verify_results.py` (brain #1001702), single-thread, BelowNormal |

### G2 measured baseline (2026-10-01, single thread, BelowNormal, during OMEGA sweep)

Route S precision ladder at `c=100, N=200, mu=0` (one fresh build per rung; artifacts in `results/g2_sturm_*.json`):

| prec (bits) | build (s) | certified LDL^T count | verdict |
|---|---|---|---|
| 2000 | 339.5 | undetermined at pivot 186 | fail-closed, too tight |
| 4000 | 1074.0 | undetermined at pivot 303 | fail-closed, too tight |
| **6000** | **1764.9** | **184.4 s -> n_pos=401, n_neg=0** | **CERTIFIED, manuscript MATCH** (minimum) |
| 9000 | 3246.2 | 305.2 s -> n_pos=401, n_neg=0 | CERTIFIED, manuscript MATCH (cross-check) |

- **Independent reproduction**: our toolchain re-derives the published finite
  certificate `c=100, N=200 -> (n_pos=401, n_neg=0)` at both 6000 and 9000 bits.
- **Route E** (`arb_mat.eig`, same matrix, prec 2000): build OK, eig did **not
  return after 7040 s** -> killed under the OMEGA resource rule; honest record
  in `results/g2_route_e_c100_N200_p2000_KILLED.json`. Retest on an idle
  machine is queued; route S is the selected path meanwhile.
- **Campaign cost model (measured, not estimated)**: one build (1764.9 s,
  reused across shifts) + ~184 s per certified count at 6000 bits
  -> 401 shifts ≈ **20-21 h single-thread**, checkpointed JSONL, resumable.
  Escalation to 9000 bits for any shift whose pivot comes back undetermined
  (measured ~305 s/count if needed).
- **Contention note (measured)**: overlapping this baseline with OMEGA-CORE
  cost OMEGA ~15% cumulative (per-pivot rate 7.9e-5 vs 5.1e-5 s/unit during
  overlap; two of the benchmark processes were observed at AboveNormal
  despite an internal BelowNormal request and were force-corrected). The
  campaign itself must run on an idle machine — which is exactly why it is
  queued behind OMEGA.

## Two routes, two engines

```
                        cheap tier (candidates)           rigorous tier (certificates)
                     sm_dual: float64 rebuild +          Arb (python-flint), read-only
                     numpy/LAPACK eigvalsh               import of source_arb_ldlt_certify
                              |                                        |
   route E --------------------------------------------------------------------> arb_mat.eig
                              |              full spectrum, one call                    |
   route S: planned shifts --> vetted Sturm counts: certified_inertia(A - mu*I) --------
                              |              N(mu) = n_neg, proved (Sylvester)          |
                              v
             segments = individual brackets (lo,hi]  +  cluster brackets (exactly k in (lo,hi])
                              |
                              v
             S1/S2 statistics with bracket error bars at sigma_u = 1.0 (queued)
```

- **Route E** runs `arb_mat.eig` on a *natively built* Arb matrix — OMEGA v1's
  caveats (mpmath digits handed over as decimal text; 94-minute timeout at
  801x801) do not apply here.  If every algorithm raises, exit code 2 records
  the tool limitation honestly and route S is the answer.
- **Route S** never separates eigenvalues, only counts them — the same
  interval LDL^T OMEGA runs.  `mu` is float64 -> dyadic -> exact in Arb.
- **Vetting**: a shift within `ambig` of a dual eigenvalue is not attempted
  (float64 cannot place it); the zone gets a *cluster certificate* from two
  endpoint counts instead.  This is not an edge case: this family's
  `lambda_min` can reach ~1e-102 (measured precedent at (100,40)), far below
  float64 resolution.

## Usage

```bash
# selftests (seconds, safe while OMEGA runs)
python sm_dual.py
python sturm_count.py --selftest --lowprio
python bracket_driver.py --selftest --lowprio

# G2 baseline: one build + one certified count (DIM=401)
python sturm_count.py --bench --c 100 --N 200 --prec 2000 --mu 0 --lowprio

# route E benchmark (same size)
python arb_eig_route.py --c 100 --N 200 --prec 2000 --lowprio

# level-0 campaign: one count per vetted shift, JSONL checkpoint, resumable
# (QUEUED until OMEGA is released)
python bracket_driver.py --c 100 --N 200 --prec 2000 --lowprio
```

Exit codes: `0` pass / fully certified, `1` fail-closed finding recorded
(conflict, coverage gap, mismatch), `2` undetermined pivot or route-E tool
limitation (raise `--prec`, never nudge).

## Honest bounds (read before quoting anything)

1. **This is a finite matrix.**  Prime cutoff `c = 100`, basis `|n| <= N`.
   Its spectrum is not "the spectrum of zeta"; truncation error is a model
   accuracy issue that more bits do NOT cure.
2. **Eigenvalue <-> zero identification is a hypothesis, not a theorem.**
   What the code proves is inertia and brackets of `tau`.  The statistical
   bridge to zero spacings is exactly what the S1/S2 instrument tests — and
   its current honest verdict is **S1 = INCONCLUSIVE, S2 = INCONSISTENT WITH
   MONTGOMERY at this resolution**, with the binding constraint measured to
   be the unfolding bandwidth (power >= 90% only at `sigma_u = 1.0`).
3. **No RH / Weil-positivity / prime-counting claims.**  Each run certifies
   one finite instance; a negative pivot would be a falsification *signal*
   requiring analysis (off-line zero vs construction artefact), never an
   automatic refutation.
4. **Cluster brackets are level-0**: "exactly k eigenvalues in (lo, hi]",
   width W.  Refining inside a zone by Sturm bisection (level-1) is queued
   and needs an explicit budget — naive refinement to `lambda_min ~ 1e-102`
   costs ~300 counts per chain (see NOTES.md C.4).
5. **Float64 tier is never evidence.**  It proposes and cross-checks;
   certificates come only from Arb (or arb_mat.eig enclosures).
6. **Resource rules**: everything here is single-thread + BelowNormal; the
   full campaign waits until OMEGA-CORE is released (brain #1001702).  The
   source module is imported **read-only**; `D:\gw_ckpt` is never touched.

## Files

| File | Role |
|------|------|
| `sm_dual.py` | independent float64 engine (mpmath closed forms + LAPACK) + W02 rank-2 selftest |
| `sturm_count.py` | route S: certified `N(mu)`, selftest, G2 bench runner |
| `arb_eig_route.py` | route E: `arb_mat.eig` full-spectrum benchmark with honest NOT-ISOLATED handling |
| `bracket_driver.py` | level-0 campaign: vetted shifts -> counts -> individual/cluster brackets, JSONL resume |
| `NOTES.md` | Montgomery–Odlyzko background + Q_Nc construction findings + design log |
| `results/` | selftest evidence, G2 JSONs, bracket JSONLs/checkpoints |

Provenance: matrix construction and interval LDL^T come from
`../GUINAND_WEIL/source_arb_ldlt_certify.py` (imported read-only), the same
module OMEGA-CORE certifies with.
