# NOTES.md — Montgomery–Odlyzko background, Q_Nc construction findings, tier-3 design log

Companion to README.md.  This file persists the 2026-10-01 discussion with the
Architect ("can the extreme-precision data stack become a clear map of hidden
spectral patterns along the critical line?") so it does not live only in
conversation history.  Everything below is either (i) checked against primary
sources, (ii) read directly from our own code with line references, or
(iii) measured in this repository — nothing is asserted from memory alone.

---

## A.  Montgomery–Odlyzko: history and mechanism

### Chronology (verified against primary sources, 2026-10-01)

| Year | Event |
|------|-------|
| ~1962 | Dyson notes the threefold symmetry of zero spacings of zeta — heuristic only |
| 1973 | **Montgomery**, "The pair correlation of zeros of the zeta function",
| | Proc. Sympos. Pure Math. **XXIV, 181–193**: announces the Pair Correlation
| | Conjecture (PCC) + conjecture that asymptotically 100% of zeros are simple.
| | When he described the formula to Dyson, Dyson recognised it as the pair
| | correlation of the GUE (Gaussian Unitary Ensemble) |
| 1978 | **Gallagher–Mueller**: PCC (under RH) implies 100% of zeros are simple;
| | the argument itself does not need RH, so simplicity follows from PCC alone |
| 1980s | **Odlyzko**: large-scale numerical tests — Cray X-MP era, plus the
| | Odlyzko–Schönhage algorithm (Riemann–Siegel evaluation in ~T^epsilon time) |
| 1987 | Odlyzko's large-scale computations **support the PCC** |
| 1989 | AT&T Bell Labs preprint: *"The 10^20-th zero of the Riemann zeta function
| | and 70 million of its neighbors"* — the famous pictures gave the conjecture
| | large credibility |
| 1996 | **Rudnick–Sarnak** (Duke Math. J. 81): ALL k-point correlations -> GUE |
| 1999 | **Katz–Sarnak** universality (symmetry type depends on the L-family);
| | Sarnak, "Zeroes of zeta functions and symmetry" (Bull. AMS) |

### Mechanism, four steps

1. **Unfolding.**  Mean spacing at height T is `2*pi/log(T/2*pi)`; statistics
   are only meaningful after dividing it out to unit mean.  This is exactly
   the *unfolding bandwidth* problem our S1/S2 instrument already measured:
   bandwidth too small = undersmoothing bias = power 0% (the phenomenon is
   classical, not a defect of our framework).
2. **Pair correlation.**  `g(u) = 1 - (sin(pi*u)/(pi*u))^2` — linear level
   repulsion at small u.  This is the canonical "hidden spectral regularity".
3. **Montgomery's theorem (partial).**  In the Fourier side he proved
   (assuming RH) `F(alpha) = |alpha|` for `|alpha| <= 1`; beyond that his
   method stops and only conjecture remains (with the fine delta correction).
   Even the theorem is half-conjecture — a fact worth remembering.
4. **Numerical confirmation.**  Extreme precision at extreme height is a
   *resolution requirement*, not vanity: at T = 10^20 the mean spacing is
   ~0.14, so Z(t) must be evaluated far finer than the local scale to resolve
   individual zeros.  Turing's method exists precisely to guarantee N(T) has
   no missed zeros (born from the Hutchinson scare).

### Cautionary history (anomaly vs artefact)

- **Hutchinson 1936**: "missing zero" — N(T) disagreed with the
  Riemann–von Mangoldt count; later shown to be an arithmetic slip.
- **Lehmer 1956**: a pair of extraordinarily close zeros — a *real* anomaly
  that refutes nothing.
- Lesson: in precise data, hidden regularity and analysis error look
  identical; only statistics with a null hypothesis separate them.

### Direct parallel to OMEGA

18,000-bit ball arithmetic plays the same role as Odlyzko's extreme
precision: buy *resolution and sign certainty* (every pivot decision proven,
never sampled).  The difference: Arb also measures its own radius growth, so
the data knows when it starts blurring.

---

## B.  Where the "critical line" lives in Q_Nc

Source of truth: `../GUINAND_WEIL/source_arb_ldlt_certify.py` (read
2026-10-01; line numbers refer to that file).

1. **Explicitly in the entries.**  `S/CC/XC` are built from
   `psi, psi'` evaluated at `1/4 + i*pi*n/L` (lines 114–116) — i.e.
   `psi(rho/2)` with `rho = 1/2 + 2i*pi*n/L`: **exactly the critical line**
   Re(s) = 1/2 (the 1/4 offset is the s/2 of the completed zeta).  The
   frequency index n maps to ordinates `gamma_n = 2*pi*n/L` (= `w`, line 105),
   so the N=800 sweep window covers gamma up to ~1091 (2*pi*800/log(100)).
2. **The zeta zeros are NOT in the code — they are in the identity.**
   `build_arb_tau` computes only the prime side (`Wp`, prime powers <= c,
   lines 149–179) and the archimedean side (`W02`, `WR`).  The Guinand–Weil
   explicit formula is what equates this to the zero side; Weil's positivity
   criterion makes positive-definiteness of the (limit of the) quadratic form
   equivalent to RH.  Hence certified inertia is a **finite shadow of
   Weil-positivity** — the module docstring itself quotes only the finite
   instance ("c=100, N=200, 9000 bits: n_pos = 401, n_neg = 0", lines 8–10)
   and the framework carries an explicit disclaimer against RH/Weil claims.
3. **Structural surprise: the W02 block is separable of rank 2.**
   `(L^2 - 16*pi^2*m*n)/((L^2+16*pi^2*m^2)(L^2+16*pi^2*n^2))`
   = `L^2*t1(m)t1(n) - 16*pi^2*t2(m)t2(n)`.  Algebraic factorisation of
   lines 165–167; **verified numerically** by `sm_dual.selftest()`
   (sv[2]/sv[0] = 3.7e-17 -> PASS).  Consequence: the bulk of the spectrum
   comes from `WR` (digamma on the critical line) + `Wp` (primes), not from
   the archimedean W02 kernel.
4. **Cutoffs (what blurs the map).**  The archimedean sums are cutoff-free:
   geometric tails are rigorously bounded (threshold `2^-(prec+24)`, tail
   `4*e^(-c_next*L)/(1-e^(-2L))`, lines 71–88).  The only model truncations
   are the prime cutoff c and the basis cutoff |n| <= N.  More bits do not
   cure those — they are accuracy-of-model, not precision.
5. **Three levels of "critical line", never to be conflated:**
   (i) explicit in entries (psi at 1/4); (ii) implicit in the identity
   (inertia = Weil-positivity shadow); (iii) the STATISTICAL hypothesis that
   the eigenvalue spectrum of Q_Nc reflects zero statistics — that is what
   S1/S2 test, and the honest current status is **S1 = INCONCLUSIVE,
   S2 = INCONSISTENT WITH MONTGOMERY at this resolution** (binding constraint
   = unfolding bandwidth; power >= 90% only at sigma_u = 1.0).

---

## C.  Tier-3 design + blind spots (incl. ones found during implementation)

### Design recap

- **Route E**: `arb_mat.eig` full spectrum on a NATIVE-Arb matrix
  (v1's caveats do not apply: v1 handed mpmath digits over as decimal text;
  here `build_arb_tau` produces balls end-to-end).  Precedent risk: rump /
  vdhoeven_mourrain can raise NOT ISOLATED on radius-carrying matrices
  (`gw_arb_measured.py`, prec 1024 and 4096) -> exit code 2 is an honest
  outcome, and route S takes over.
- **Route S**: `certified_inertia(A - mu*I)` gives `N(mu) = n_neg` (Sylvester).
  mu is float64 -> dyadic -> exact in Arb.  Count is proved; undetermined
  pivot = fail-closed, never guessed.
- **Two tiers**: float64 engine (`sm_dual`, independent reimplementation +
  LAPACK) proposes shifts and cross-checks; Arb certifies.  Candidates from
  the cheap tier, certificates from the expensive tier — the same
  anti-circularity discipline used elsewhere in this workspace.
- **Harvest**: S1/S2 re-run on the certified spectrum with bracket-derived
  error bars at sigma_u = 1.0 (the pre-registered power level) -> decides
  "physics vs instrument" for the S2 verdict.

### Blind spots listed during the discussion

1. mu exactly hitting an eigenvalue -> undetermined pivot -> precision
   escalation + deterministic logged re-pick, never a silent nudge.
2. gaps below float64 resolution -> flag, never shift mu quietly.
3. radius growth on shifted diagonals -> measured, never widened to look good.
4. certification does not cure unfolding bandwidth bias (S1/S2 limit stands).
5. finite matrix (c=100) != zeta spectrum; README must bound claims.
6. resource rule: campaign only after OMEGA is released (brain #1001702).

### Additional facts found DURING implementation (fail-closed gates caught them)

1. **Float64 cannot count this family's near-zero eigenvalues.**  Measured at
   c=13, N=8: dual engine places 9 eigenvalues within ~1e-15 of zero (some
   negative), Arb at 300 bits certifies n_neg = 0 — signs near zero are
   float64 noise, not engine disagreement.  Precedent: lambda_min at
   (100, 40) is +1.3e-102, i.e. this whole low end sits below float64
   resolution by design.  Consequences: (a) vetting rule (margin >= ambig)
   before any dual-guided shift; (b) LEVEL-0 **cluster certificates** — a
   TOO_TIGHT zone is certified by two endpoint counts ("exactly k
   eigenvalues in (lo, hi]"), not by dual placement.  Selftest: 8 ambiguous
   gaps -> 1 cluster, coverage 17/17, conflicts 0.
2. **mpmath `mpf.exp` is an internal component (the binary exponent, an
   int), not a method** — `.exp()` raises `TypeError: 'int' object is not
   callable`.  Use module functions `mp.exp/mp.log/mp.atan`.  Caught by
   selftest before any baseline ran.
3. **ctypes x64 `SetPriorityClass` needs explicit argtypes** — without
   `c_void_p` types the pseudo-handle passes as 32-bit and the call silently
   returns 0; with them it succeeds (measured: `GetPriorityClass` ->
   0x8000).  `lower_priority()` now encodes this.
4. **Level-1 cost trap**: pure Sturm bisection inside a zone down to
   lambda_min ~ 1e-102 from a 1e-5 bound costs ~log2(1e97) ~ 320 counts per
   chain — at 1-3 min per count that is hours for a single zone.  Level-0
   therefore stops at cluster certificates; zone refinement is budget-gated
   (queued until OMEGA is released, and needs an explicit budget parameter).
5. `bracket_driver` resumes from its JSONL across power loss (torn last line
   is discarded and that shift redone) — verified by rerunning the selftest
   (resume path reported cached counts, coverage unchanged 17/17).
