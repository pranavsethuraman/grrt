# Phase 3B Status — COMPLETE (2026-04-27)

Final-state snapshot per the CLAUDE.md "Session Hygiene" rule. Sub-phase
3B (`Disk` package, Novikov-Thorne disk emission, DP45 dense output,
position-dependent JP horizon termination) is closed; tag
`phase-3b-nt-disk` lives at commit `23526b3`.

---

## 1. Final gate residuals

All six gates from `docs/phase-3-plan.md` §5.5 (gate scoreboard format
matching `docs/phase-3a-status.md` §1):

| Gate | Description | Tolerance | Final residual | Margin |
|---|---|---|---|---|
| 1 | NT radial flux F(r) vs Page-Thorne 1974 ref, 12 rows × 2 spins | 1e-6 rel | **1.20e-14 rel** | ~8 orders |
| 2 | Kerr a=0.9 prograde ISCO ≈ 2.32088 (Bardeen 1972) | 1e-4 | exact via `iscoRadius` override | n/a |
| 3 | Face-on Schwarzschild disk centroid at 256² | 0.5 px | **0.0 px** (exact) | ∞ |
| 4 | Redshift formula consistency (Schwarzschild + Kerr ISCO face-on) | 1e-12 | machine precision | n/a |
| 5 | DP45 dense output at θ=0.5 vs half-step direct, Kerr photon orbit | 1e-8 | **1.00e-11** | 3 orders |
| 6 | Phase 2 + 3A regression (88 prior tests) | bit-exact | all pass | n/a |
| 3B.3 | JP `isInsideHorizon` override at both ε₃ signs (equator + axis) | algebraic | F-sign verified at hand-derived points | n/a |

Test count: **108 tests, 0 failures, 0 errors** across all 17 test
classes. New tests added in 3B (relative to `phase-3a-complete`):

- `NovikovThorneDiskTest` — 9 tests
- `DiskEmissionShaderTest` — 5 tests
- `DormandPrince45KerrTest` — 3 added (gate 5 + interpolate negative cases)
- `RendererTest.faceOnSchwarzschildDiskCentroidIsNearOrigin` — 1 (gate 3)
- `JohannsenPsaltisMetricTest.isInsideHorizonOverridesDefault*` — 2

Net 3B test delta: **+20** (88 → 108).

---

## 2. Commit chain and tag

```
0c906c7  disk: add Novikov-Thorne disk model and DP45 dense output (Phase 3B.1)
c0c856d  renderer: plumb disk crossing into AdaptiveRayTracer with HIT_DISK (Phase 3B.2)
23526b3  metric: JP isInsideHorizon override unit test, close Phase 3B (Phase 3B.3)
```

Tag: **`phase-3b-nt-disk`** at `23526b3`. Pushed to `origin/main`. Span:
`f779db1..23526b3`. Per-commit file counts: 11 / 7 / 2 = 20 file
modifications across 3B.

Pre-3B.1 supporting commits (already pushed before 3B.1 landed):

```
915cd66  docs: tick JohannsenPsaltisMetric in CLAUDE.md status
45abbd5  docs: relabel disk emission model from Phase 2 to Phase 3
4030069  docs: record approved Phase 3B sub-plan
```

---

## 3. Files created and modified across 3B

### 3.1 Created (8)

- `scripts/page_thorne_reference.py` — Sympy + scipy.integrate.quad
  reference for the gate-1 NT-flux table.
- `src/main/java/com/pranav/grrt/disk/Disk.java` — interface.
- `src/main/java/com/pranav/grrt/disk/NovikovThorneDisk.java` —
  Page-Thorne 1974 eq. (15n) implementation.
- `src/main/java/com/pranav/grrt/disk/DiskEmissionShader.java` —
  bolometric Planck g⁴σT⁴/π shader.
- `src/test/java/com/pranav/grrt/disk/NovikovThorneDiskTest.java` —
  gates 1, 2, JP smoke, validation tests.
- `src/test/java/com/pranav/grrt/disk/DiskEmissionShaderTest.java` —
  gate-4 redshift formula tests.
- `docs/phase-3b-plan.md` — committed as `4030069` before 3B.1.
- `docs/phase-3b-status.md` — this file.

### 3.2 Modified (8)

- `src/main/java/com/pranav/grrt/integrator/DormandPrince45.java` —
  Hairer contd5 5th-order dense output, lazy yPrev/yNext/hAccepted/
  denseValid state, new `interpolate(θ, out)`.
- `src/main/java/com/pranav/grrt/metric/Metric.java` — two new
  default methods: `iscoRadius`, `isInsideHorizon`.
- `src/main/java/com/pranav/grrt/metric/KerrMetric.java` —
  `@Override` on existing `iscoRadius`.
- `src/main/java/com/pranav/grrt/metric/JohannsenPsaltisMetric.java` —
  `@Override` on existing `iscoRadius`; new override of
  `isInsideHorizon` with position-dependent test.
- `src/main/java/com/pranav/grrt/renderer/RayOutcome.java` — added
  `HIT_DISK`.
- `src/main/java/com/pranav/grrt/renderer/AdaptiveRayTracer.java` —
  Disk-aware constructor overload, sign-change disk-crossing
  detection with `interpolate`-driven bisection, switched termination
  to `metric.isInsideHorizon`.
- `src/main/java/com/pranav/grrt/disk/DiskEmissionShader.java` —
  `shade()` dispatches on `HIT_DISK` (returns `intensity(state)`).
- `src/test/java/com/pranav/grrt/integrator/DormandPrince45KerrTest.java` —
  gate-5 dense-output test + 2 interpolate negative-case tests.
- `src/test/java/com/pranav/grrt/renderer/RendererTest.java` —
  gate-3 face-on Schwarzschild disk render test.
- `src/test/java/com/pranav/grrt/validation/KerrShadowTest.java` —
  added `case HIT_DISK -> Float.NaN` to keep exhaustive switch valid.
- `src/test/java/com/pranav/grrt/metric/JohannsenPsaltisMetricTest.java` —
  two new `isInsideHorizon` override tests.
- `CLAUDE.md` — ticked `Disk emission model (Phase 3)` checkbox;
  added interface-extension paragraph documenting the two new
  `Metric` default methods (`iscoRadius`, `isInsideHorizon`).

---

## 4. `Metric` interface extensions added in Phase 3B

Documented in CLAUDE.md §Coding rules. Two new `default` methods:

| Method | Default body | Overridden by | Sub-phase |
|---|---|---|---|
| `double iscoRadius(boolean prograde)` | throws `UnsupportedOperationException` | `KerrMetric`, `JohannsenPsaltisMetric` | 3B.1 |
| `boolean isInsideHorizon(double[] x, double tol)` | `x[1] - horizonRadius() < tol` | `JohannsenPsaltisMetric` (position-dependent) | 3B.2 |

The `isInsideHorizon` default is bit-exact with the prior
`AdaptiveRayTracer` `rNow < rCapture` check, so all Phase 2 + 3A
regression tests pass unchanged.

---

## 5. Open items resolved / closed out

| Item | Status |
|---|---|
| Page-Thorne 1974 reference table generation | ✅ committed in 3B.1 (`scripts/page_thorne_reference.py`) |
| DP45 dense-output 5th-order interpolant | ✅ Hairer contd5 in 3B.1, gate 5 at 1e-11 |
| Disk inner edge per-frame caching | ✅ NT constructor caches `r_ISCO` once |
| `RayOutcome.HIT_DISK` consumer audit | ✅ all consumers reviewed; KerrShadowTest switch arm extended |
| Position-dependent JP horizon termination | ✅ `isInsideHorizon` override + 2 unit tests in 3B.3 |
| JP(0.9, +0.20) full-render smoke | 🟡 deferred to 3C — falls out naturally as one of 30 sweep frames |
| `Metric` interface extension docs | ✅ CLAUDE.md §Coding rules updated 3B.3 |

---

## 6. Phase 3C readiness check

### 6.1 New files Phase 3C will introduce

Per `docs/phase-3-plan.md` §6.1:

```
src/main/java/com/pranav/grrt/analysis/RingExtractor.java
src/main/java/com/pranav/grrt/analysis/EpsilonSweep.java
src/test/java/com/pranav/grrt/analysis/RingExtractorTest.java
src/test/java/com/pranav/grrt/analysis/EpsilonSweepIT.java
```

New top-level package `com.pranav.grrt.analysis`. ~600 lines total.

### 6.2 Phase-3B-tagged files Phase 3C will modify

Per `docs/phase-3-plan.md` §6.2:

1. **`pom.xml`** — surefire `excludedGroups=slow` so `mvn test` stays
   under ~30 s and the full sweep runs only under
   `mvn verify -PrunSlow`.
2. **`CLAUDE.md`** — tick `EHT comparison pipeline` checkbox at end of 3C.
3. **`.gitignore`** — confirm `output/` is ignored (already is; spot-check).

### 6.3 First three pre-coding decisions for 3C

1. **NovikovThorneDisk surfaceFlux caching strategy.** 3B.1 left
   `surfaceFlux(r)` as adaptive Simpson per call (~10 µs each). The
   gate-3 face-on Schwarzschild render at 256² took 1.13 s with
   ~32k disk hits — so per-pixel cost is ~30 µs (mostly disk shader,
   not flux integration). At 512² with ~130k disk hits per frame and
   30 sweep frames, total disk-shader time scales to ~2.4 minutes
   wall-clock — comparable to the ray-tracing budget but not
   dominating. **Decision: profile after the first sweep frame; if
   surfaceFlux dominates, add a 1000-point log-spaced spline cache in
   the constructor (no public API change).** Otherwise keep on-the-fly.

2. **Slow-test categorization.** The full sweep takes ~12 minutes, far
   too long for `mvn test`. The plan says "surefire excludedGroups=slow
   so mvn test stays fast and the full sweep runs under mvn verify
   -PrunSlow". **Decision: tag `EpsilonSweepIT` with JUnit 5 `@Tag("slow")`;
   configure surefire to exclude `slow` by default; add a Maven profile
   `runSlow` that drops the exclusion.** Standard pattern.

3. **Asymmetry metric: peak-intensity radius per bin tie-breaking.**
   Per `phase-3-plan.md` §6.9, "Use radially-outer peak per bin; flag
   bins where multiple peaks exist within ~2 px." For face-on
   Schwarzschild (no asymmetry) every bin has a single peak; for
   inclined renders the n=1 lensed ring may cross the direct image at
   some azimuths, producing two peaks. **Decision: outermost peak
   wins; emit a per-frame WARN (`System.out.println` with
   `[ring-warn]` prefix) when any bin has multiple peaks within 2 px,
   and store the affected bin count in a CSV column for diagnostics.**

### 6.4 Open questions from 3B that affect 3C

- **JP equatorial horizon at ε₃ > ε₃_crit:** the `isInsideHorizon`
  override is unit-tested but not yet exercised by an end-to-end
  render. The 3C sweep includes ε₃ ∈ {+0.13, +0.20, +1.00} which sit
  past the cusp; rendering these frames will be the de-facto smoke
  test. If a frame produces unexpected NaNs or "max-steps" outcomes
  beyond the expected post-cusp count, investigate `isInsideHorizon`
  before continuing the sweep.
- **Camera at face-on inclination:** for Phase 3C primary inclination
  i = 17°, the existing Camera at non-axis θ is fine. The face-on
  smoke (i = 1e-4 rad) was a 3B-only test convenience; not a 3C
  concern.
- **Resumability semantics:** CSV schema in `phase-3-plan.md` §6.6
  uses `git_sha`; confirm before sweep starts whether this is the
  abbreviated 7-char or full 40-char SHA. **Default: 7-char abbrev**
  (matches `git log --oneline` and is enough to disambiguate within
  the project's commit volume).

### 6.5 Estimated wall-clock to `phase-3c-sweep`

Coding effort (focused sessions):

- ~1 session for `RingExtractor.java` + `RingExtractorTest.java`
  (synthetic-circle + ellipse gates 1, 2 + bin-convergence gate 6).
  Pure analysis code; no integration.
- ~1 session for `EpsilonSweep.java` + `EpsilonSweepIT.java`
  (CSV resumability gate 7 + Kerr/JP consistency gate 3).
- ~1 session for the actual sweep run (`mvn verify -PrunSlow`):
  - 30 frames × ~24 s/frame at 512² = **~12 min** wall-clock.
  - Plus diagnostics: ring-extraction post-process, CSV emit.
- ~1 session for figure generation in the 3D paper sub-phase
  (separate, follows 3C tag).

Total to `phase-3c-sweep`: **3 focused sessions, ~6–10 hours session
time, plus ~12 min sweep wall-clock**. Calendar time depends on
session cadence.

### 6.6 Performance budget at 3C

| Step | Budget per frame | 30-frame total |
|---|---|---|
| AdaptiveRayTracer + Disk render at 512² | ~24 s | ~12 min |
| RingExtractor (180 bins) | < 1 s | ~30 s |
| FITS write | ~50 ms | ~1.5 s |
| CSV emit | < 10 ms | < 1 s |
| **Total per sweep run** | **~25 s** | **~13 min** |

`mvn test` (no slow tests) budget unchanged: < 30 s, currently ~28.6 s
at `phase-3b-nt-disk`.

---

## 7. Notes for next session

- Resume with a fresh session per "we resume tomorrow".
- First action: read this file, then `docs/phase-3-plan.md` §6
  (Phase 3C sub-plan) and `docs/jp-parameter-space-notes.md` (cusp
  context for the post-+0.12 ε₃ rendering).
- Surface for user approval BEFORE 3C coding starts: the three
  decisions in §6.3 (caching, slow-test tagging, asymmetry tie-break)
  and the open question on `git_sha` length (§6.4).
- Working tree clean. `origin/main` matches local. Tag pushed.

---

*End of Phase 3B status snapshot.*
