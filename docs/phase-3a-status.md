# Phase 3A Status — COMPLETE (2026-04-24)

Final-state snapshot per the CLAUDE.md "Session Hygiene" rule.
Sub-phase 3A (`JohannsenPsaltisMetric` and the six validation gates)
is closed; tag `phase-3a-complete` lives at commit `f779db1`.

---

## 1. Final gate residuals

All six gates from `docs/phase-3-plan.md` §4.7, with the numbers from
the `mvn test` run on `f779db1`:

| Gate | Description | Tolerance | Final residual | Margin |
|---|---|---|---|---|
| 1 | `g · gInv = I` (algebraic identity) | 1e-12 | ~1e-13 (machine-eps order) | ~10× |
| 2 | Christoffel lower-index symmetry | 1e-14 | ~1e-14 | 1× |
| 3 | Kerr reduction `max\|Γ_JP(ε₃=0) − Γ_Kerr\|` | 1e-12 | **1.42e-14** | ~70× |
| 3 sanity | (M, a) variation residual ≈ \|g\|·ε_machine | tracking | confirmed across (1.0,+0.5), (2.0,+0.9), (1.0,−0.9) | n/a |
| 4 | Asymptotic flatness `\|g_JP − g_Kerr\|` at r=1000M, ε₃=1 | 1e-6 | **1.00e-9** | 1000× |
| 5 | DP45 photon orbit `\|g·k·k\|/(k^t)²` over 1000 M, JP(0.9, +0.5) | 1e-8 | **1.67e-10** | 60× |
| 6 | `photonOrbitRadius` vs Python reference, 7 (a, ε₃) pairs | 1e-10 | **2.09e-12** (worst row) | ~50× |
| 6 amend. 4 | Cross-path equatorial-derivative agreement, 18 components | 1e-13 | **0.0e+00** (bit-identical) | ∞ |

Test count: **88 tests, 0 failures, 0 errors** across all 15 test classes.
JP-specific: 25 tests in `JohannsenPsaltisMetricTest`.

---

## 2. Headline scientific result — the cusp

The dominant scientific finding of Phase 3A is the discovery of a
**bifurcation in the JP equatorial prograde photon orbit at
M87*-consistent spin**: the Kerr-continuation prograde orbit merges
with the Kerr horizon at

```
ε₃_crit(a = 0.9) = 0.1212    (4-decimal bisection precision)
```

with `r_photon → r₊_Kerr ≈ 1.4359 M` from above as ε₃ → ε₃_crit. For
ε₃ > ε₃_crit, no Kerr-continuation prograde photon orbit exists above
the horizon; the prograde side of the M87* shadow is bounded by the
event horizon directly, not by a photon sphere. The retrograde orbit
is unaffected (smooth across [−3, 0] and approximately smooth above).

This is **not a metric pathology** (the metric is regular on the
equatorial exterior and the spacetime has no CTCs in the approved
sweep range — see `docs/jp-parameter-space-notes.md` §3.1). It is a
genuine bifurcation in the orbital structure, the JP analogue of the
Kerr extremality limit (`r_photon_pro → M` as `a → M`), shifted by
the deformation parameter.

The RNAAS paper has been reframed around this finding (see CLAUDE.md
Project section): the bound on ε₃ is extracted from the smooth
sub-interval of the swept asymmetry curve, and the cusp is the
**additional prediction** of the JP metric — testable by future
higher-precision observations.

---

## 3. Commit range and tag

```
a9555ea  metric: add Johannsen-Psaltis metric with epsilon_3 deformation (Phase 3A)  (gates 1-4)
cfe7cbe  docs: plan for Phase 3A gates 5 and 6                                       (sub-plan)
e02ec8b  docs: characterize JP photon-orbit transition at a=0.9                       (cusp finding)
e790272  scripts: reference tables for JP photon orbit and transition locus           (Python reproducibility)
079571e  docs: reframe RNAAS paper around JP photon-sphere bifurcation                (CLAUDE.md reframe)
a30ab4a  docs: update sub-plans and findings for reframed Phase 3 scope               (Commit B doc updates)
f779db1  metric: close Phase 3A with gate 5 (null-norm drift) and gate 6 ...          (gates 5+6, FINAL)
```

Tag: **`phase-3a-complete`** at `f779db1`. Pushed to
`origin/main`. Span: `a9555ea .. f779db1`.

---

## 4. Files modified this session (3A close-out)

| File | Lines added | Summary |
|---|---|---|
| `src/main/java/com/pranav/grrt/metric/JohannsenPsaltisMetric.java` | +236 | `photonOrbitRadius(boolean)` public; `residualPhotonOrbit` private; two package-private helpers `equatorialGDerivativesViaIscoPath`, `equatorialGDerivativesViaPhotonPath` for cross-path test |
| `src/test/java/com/pranav/grrt/metric/JohannsenPsaltisMetricTest.java` | +208 | three new tests: `photonOrbitRadiusMatchesPythonReference`, `equatorialDerivativesAgreeAcrossPaths`, `geodesicStaysOnLightConeUnderDp45_jp_a09_eps3_0p5` |
| `scripts/jp_photon_orbit_reference.py` | +14/−10 | `pairs` list updated to gate-6 approved 7-row table |
| `docs/phase-3a-gates-5-6.md` | +21/−5 | §2.1 DP45 tolerances confirmed; §3.4 bracket-walk replaced with grid-scan |

---

## 5. Phase 3A — open items resolved / closed out

| Item | Status |
|---|---|
| Gate 6 pair revision after cusp discovery | ✅ applied (script + sub-plan §3.3) |
| DP45 tolerance carryover for §2.1 | ✅ recorded (`atol=rtol=1e-10, h=1.0M`) |
| Bracket-walk vs grid-scan in `photonOrbitRadius` | ✅ grid-scan implemented and documented |
| Cross-path derivative regression guard | ✅ bit-identical agreement at 3 sample points |
| CLAUDE.md `JohannsenPsaltisMetric` checkbox | ⚠️ **not yet ticked** — still shows `[ ]`; lives in §Status of CLAUDE.md, not blocking but worth fixing in 3B kickoff |

---

## 6. Phase 3B readiness check

### 6.1 New files Phase 3B will introduce

Per `docs/phase-3-plan.md` §5.1:

```
src/main/java/com/pranav/grrt/disk/Disk.java
src/main/java/com/pranav/grrt/disk/NovikovThorneDisk.java
src/main/java/com/pranav/grrt/disk/DiskEmissionShader.java
src/test/java/com/pranav/grrt/disk/NovikovThorneDiskTest.java
src/test/java/com/pranav/grrt/disk/DiskEmissionShaderTest.java
```

New top-level package `com.pranav.grrt.disk`. ~600 lines total.

### 6.2 Phase-2-tagged files Phase 3B will modify

Per `docs/phase-3-plan.md` §5.3:

1. **`src/main/java/com/pranav/grrt/integrator/DormandPrince45.java`** —
   add `interpolate(double θ, double[] out)` for dense output.
   Existing `step()` signature unchanged; existing DP45 tests must
   pass bit-exactly (regression guard, sub-plan §5.5 gate 6).
2. **`src/main/java/com/pranav/grrt/renderer/AdaptiveRayTracer.java`** —
   plumb optional `Disk` alongside horizon termination; sign-change
   detection on `θ − π/2` between accepted DP45 steps; bisect with
   the new interpolant to find crossing, then invoke shader.
3. **`src/main/java/com/pranav/grrt/renderer/RayOutcome.java`** — add
   `HIT_DISK` enum variant alongside `HORIZON`, `ESCAPED`,
   `MAX_STEPS`.
4. **`CLAUDE.md`** — tick `Disk emission model (Phase 2)` checkbox at
   the end of 3B (and the deferred `JohannsenPsaltisMetric` from 3A).

### 6.3 First three pre-coding decisions for 3B

1. **DP45 dense-output interpolant scheme.** D-B1 approved
   *adding* `interpolate()` but did not specify the order. Three
   candidates:
   - **Shampine 5th-order** (uses the existing 7 stages plus one
     extra `k_dense` evaluation). Tightest, ~O(h⁵) error. Preferred
     if the §5.5 gate 5 (1e-8 dense-output sanity) demands it.
   - **4th-order natural Hermite** from `(y, y')` at endpoints.
     Cheap, O(h⁴). Probably overkill given DP45 itself is O(h⁵).
   - **Cubic Hermite.** Cheapest, O(h³). May be too loose for the
     1e-8 sanity gate.
   Decision needed before touching `DormandPrince45.java`.

2. **Page-Thorne 1974 reference flux F(r) values for gate 1.**
   Sub-plan §5.5 gate 1 asserts `F(r)` matches the corrected
   Page-Thorne formula at six radii for `a ∈ {0, 0.9}`. Three
   options:
   - Compute the analytic values inline in the test (verify the
     formula via Sympy in a scratch).
   - Tabulate from the paper appendix.
   - Generate via a Python scratch and paste (same regenerable-
     reference pattern as `scripts/jp_photon_orbit_reference.py`).
   Choose before writing `NovikovThorneDiskTest`.

3. **Disk inner edge for ε₃ ≠ 0 in the 3C sweep.** The disk
   truncates at `r_ISCO(metric)`. For JP `ε₃ ≠ 0`,
   `iscoRadius()` returns the numeric value computed at metric
   construction. Confirm:
   - the renderer caches this once at setup vs recomputes per pixel
     (perf concern only; no correctness issue),
   - for the 3C sweep, each ε₃ frame uses its own
     `r_ISCO(JP(0.9, ε₃))`, not a single shared value.
   Trivial in principle; flagged here so it's not forgotten when
   wiring `EpsilonSweep` to the renderer.

### 6.4 Open questions from gates 5/6 that affect 3B

- **JP equatorial horizon topology for ε₃ > 0.** For ε₃ above some
  threshold (verified at ε₃ = +0.5 in
  `docs/jp-parameter-space-notes.md` §3.1.iii), the equatorial
  `g_rr` does NOT diverge — there is no equatorial horizon surface,
  only a polar-axis horizon. The current `AdaptiveRayTracer`
  termination uses a fixed `r ≤ horizonRadius()` test that returns
  the polar-axis value. **For 3C sweep frames at ε₃ > ε₃_crit the
  termination criterion is incorrect** (a photon can pass through
  what the polar-axis test calls "the horizon" without the geometry
  being singular on its trajectory). Either:
  - restrict 3C to ε₃ < ε₃_crit (which is fine — that is exactly
    where the bound is extracted, and the cusp is sampled at one
    point past it), or
  - add a position-dependent termination test in 3B
    (`Δ + a² h sin²θ ≤ tol` rather than fixed `r`).

  The `phase-3-plan.md` §2.5 already mentions the position-dependent
  test as a 3A design note; 3B implements it (if we go that route)
  or 3C restricts the range (if we don't). **Decision needed at 3B
  kickoff.**

- **No other gate-5/6 carryovers** affect 3B. Gate 5 (null-norm
  drift) confirmed the integrator+metric coupling is solid at
  ε₃ ≠ 0, exactly what 3C needs. Gate 6 confirmed the equatorial-
  derivative implementations agree, which 3B's NT computation will
  reuse.

### 6.5 Estimated wall-clock to `phase-3b-complete`

Coding effort, in focused sessions:

- ~1 session for `Disk.java` + `NovikovThorneDisk.java` +
  `DiskEmissionShader.java` + `DormandPrince45.interpolate()`
  (pure code, no integration with renderer yet). ~300 lines.
- ~1 session for `AdaptiveRayTracer` disk-crossing detection,
  `RayOutcome.HIT_DISK`, position-dependent termination if chosen.
  ~150 lines.
- ~1 session for the six validation gates in §5.5 (NT flux, ISCO,
  face-on annulus, redshift, dense-output sanity, regression of
  prior DP45 tests). ~250 lines of test code + ~3 s of `mvn test`
  per run.

Total: **3 focused sessions, ~6–10 hours session time**, plus any
loops on the dense-output interpolant or position-dependent
termination decisions in §6.3 and §6.4. Calendar time depends on
session cadence.

---

## 7. Notes for next session

- Resume with a fresh session per "we resume tomorrow".
- First action: read this file, then `docs/phase-3-plan.md` §5
  (Phase 3B sub-plan) and `docs/phase-3a-gates-5-6.md` §1
  (just-superseded; for context only).
- Surface for user approval BEFORE 3B coding starts: the three
  decisions in §6.3 and the JP-horizon-termination question from
  §6.4.
- No `mvn test` regressions or warnings as of `f779db1`. Working
  tree clean. `origin/main` matches local. Tag pushed.

---

*End of Phase 3A status snapshot.*
