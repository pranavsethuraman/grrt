# Phase 3C Sub-plan — Ring asymmetry extraction + ε₃ sweep

Companion to `docs/phase-3-plan.md` §6 and `docs/jp-parameter-space-notes.md`
§5.2. Builds the analysis pipeline that turns rendered FITS frames into
the `δ_r/⟨r⟩(ε₃)` curve underlying the RNAAS bound, and runs the
approved 16-point single-curve sweep at `a = 0.9` for both inclinations.

Parent commit: `23526b3` (`phase-3b-nt-disk`). Sub-phase 3C lands as
three commits leading to tag `phase-3c-sweep`. Drafted 2026-04-28
following the 3C-readiness review in `docs/phase-3b-status.md` §6.

This plan supersedes `docs/phase-3-plan.md` §6 where the two overlap.

---

## 1. Scope

Sub-phase 3C delivers four additions and zero `Metric`/`Integrator`
interface extensions:

1. New `com.pranav.grrt.analysis` package: `RingExtractor` (image →
   per-azimuthal-bin peak-intensity radius), `CircularityMetric`
   (radius array → mean / RMS / peak-to-peak / Fourier m=1), and
   `EpsilonSweep` (orchestrator: per-frame metric + disk + render +
   extract + CSV row).
2. Surefire group split: `@Tag("slow")` integration tests excluded
   from default `mvn test`, runnable via `mvn verify -PrunSlow`.
3. Resumable CSV-driven sweep, idempotent under interrupt/restart.
4. Final figure-quality data product `output/sweep.csv` populated by
   the production sweep run (32 frames at 512², ~13 min wall-clock).

No edits to `Metric`, `Integrator`, `Camera`, or any existing
`renderer/` class. The `Disk` / `DiskEmissionShader` plumbing landed
in 3B is exercised end-to-end here for the first time across a
parameter sweep.

---

## 2. Decisions to confirm at 3C kickoff

The five carry-overs from `docs/phase-3b-status.md` §§6.3–6.4. Each
needs explicit user approval before 3C.1 coding starts.

### 2.1 Class split: `RingExtractor` vs `CircularityMetric`

Parent plan §6.1 lists two analysis classes (`RingExtractor`,
`EpsilonSweep`). User-requested 3C scope adds a third:
`CircularityMetric`, isolating the pure-math reduction from the
image-space extraction.

**Proposed contracts:**

```java
public final class RingExtractor {
    public RingExtractor(int numBins, double pixelToM, double cxPx, double cyPx);
    public double[] extract(float[][] image, int width, int height);   // length numBins, NaN allowed for empty bins
    public int multipeakBinCount(float[][] image, int width, int height, double pxThreshold);
}

public final class CircularityMetric {
    public record Result(double meanR, double deltaRrms, double deltaRp2p,
                         double fourierM1Amplitude, double fourierM1Phase,
                         int validBins) {}
    public static Result compute(double[] perBinRadii);
}

public final class EpsilonSweep {
    public record Config(double spin, double[] epsilonGrid, double[] inclinationDegrees,
                         int resolution, double rOuter, Path csvPath, Path fitsDir);
    public void runSweep(Config cfg);
    public static void main(String[] args);   // production sweep, hard-coded config
}
```

`RingExtractor` holds image-space geometry (centroid, pixel scale,
bin count) so it can be reused across frames. `CircularityMetric`
is stateless. `EpsilonSweep` owns CSV r/w + per-frame orchestration.

**Why split?** Three independently-testable concerns: image-space
peak finding, vector-space statistics, and sweep orchestration.
Splitting also matches the gate decomposition (gates 1–2 exercise
the round-trip; gates with fixed synthetic radius arrays exercise
`CircularityMetric` in isolation; gates 3 and 7 exercise the
orchestrator without re-testing math).

**Decision needed:** approve the three-class split, or collapse
`CircularityMetric` back into `RingExtractor` (matching parent plan
§6.1). The split is the recommendation here.

### 2.2 NT `surfaceFlux` caching: profile-driven, default off

Per `phase-3b-status.md` §6.3 decision 1: 3B left
`NovikovThorneDisk.surfaceFlux(r)` as adaptive Simpson per call
(~10 µs each). At 512² with ~130k disk hits per frame and 32 sweep
frames, total disk-shader time projects to ~2.4 min wall-clock —
non-dominant.

**Proposal:** profile during the first three sweep frames in 3C.3.
If `surfaceFlux` exceeds 30% of frame wall-clock, add a 1000-point
log-spaced spline cache in `NovikovThorneDisk`'s constructor (no
public API change, no test break). Otherwise keep the per-call
adaptive Simpson.

**Decision needed:** profile-then-decide is a reasonable default. If
user wants to add the spline cache pre-emptively to keep the sweep
strictly under 12 min, do it as a fourth commit (3C.0 before 3C.1).

### 2.3 Slow-test tagging: `@Tag("slow")` + Maven `runSlow` profile

Per `phase-3b-status.md` §6.3 decision 2: standard pattern.

`pom.xml` modification (added in commit 3C.1):

```xml
<plugin>
  <artifactId>maven-surefire-plugin</artifactId>
  <version>3.5.2</version>
  <configuration>
    <excludedGroups>slow</excludedGroups>
  </configuration>
</plugin>

<profiles>
  <profile>
    <id>runSlow</id>
    <build>
      <plugins>
        <plugin>
          <artifactId>maven-surefire-plugin</artifactId>
          <configuration>
            <excludedGroups combine.self="override"></excludedGroups>
          </configuration>
        </plugin>
      </plugins>
    </build>
  </profile>
</profiles>
```

`EpsilonSweepIT` carries `@Tag("slow")`. Default `mvn test` budget
remains < 30 s. `mvn verify -PrunSlow` runs the full IT suite
(~75 s, see §8.1).

**Decision needed:** confirm the surefire group key is `slow` (vs.
`integration`, `e2e`, etc.). `slow` matches the `phase-3b-status.md`
language and is the recommendation.

### 2.4 Multi-peak tie-break: outermost peak, WARN print, CSV column

Per `phase-3b-status.md` §6.3 decision 3: outermost peak per bin
wins. `RingExtractor` emits a `[ring-warn]`-prefixed
`System.out.println` when any bin has multiple peaks within 2 px,
and returns the affected bin count via `multipeakBinCount(...)`.
`EpsilonSweep` writes that count to the `multipeak_bins` CSV column
for diagnostics.

This matters most past the cusp (§9.2): bins straddling the
prograde shadow edge will likely show a disk-edge peak alongside a
faint photon-ring residual. Counting them lets us flag frames where
RingExtractor's interpretation of "the ring" has shifted from
"photon ring" to "disk inner edge".

**Decision needed:** the 2 px threshold is a guess (≈ 0.04 M at
512² with ±15 M field of view). Confirm or override.

### 2.5 `git_sha` length and dirty-tree behavior

Per `phase-3b-status.md` §6.4: 7-char abbreviation matches
`git log --oneline`.

**Proposal:** `EpsilonSweep` runs `git rev-parse --short=7 HEAD` at
process startup. If the working tree is dirty
(`git diff-index --quiet HEAD --` non-zero), append `-dirty`
(e.g. `23526b3-dirty`). Each CSV row records the `git_sha` valid at
the time of that frame's render — if the user does
modify-rebuild-resume mid-sweep, different rows will have different
SHAs, which is auditably correct.

**Decision needed:** confirm 7-char + `-dirty` suffix.

---

## 3. Sweep grid

Per `docs/jp-parameter-space-notes.md` §5.2 (approved 2026-04-24).
Single-curve sweep at `a = 0.9` for both inclinations:

```
Negative side (smooth, 8 points):
  ε₃ ∈ {−2.5, −1.5, −1.0, −0.5, −0.2, −0.1, −0.05, 0.0}

Cusp approach (smooth, 4 points):
  ε₃ ∈ {+0.05, +0.08, +0.10, +0.11}

Cusp itself (1 point):
  ε₃ = +0.12                     (just below ε₃_crit ≈ 0.1212)

Transition regime (3 points):
  ε₃ ∈ {+0.13, +0.20, +1.00}     (above ε₃_crit; prograde shadow at horizon)

Inclinations: i ∈ {17°, 60°}     (primary + reference)
Resolution: 512²
```

Total: **16 ε₃ × 2 inclinations = 32 frames at 512²**.

### 3.1 Off-by-one with `phase-3-plan.md` §6.5 / `jp-parameter-space-notes.md` §5.2

Both prior docs claim "15 points × 2 inclinations = 30 frames" but
list 16 values. The 16-point list is authoritative: dropping any
single ε₃ would either remove a smooth-regime data point that the
RNAAS bound interpolation depends on, or thin the cusp-resolution
neighborhood `{+0.10, +0.11, +0.12}` that the headline figure
needs. **Recommend rendering all 16 points (32 frames).**

Wall-clock impact: ~13 min vs ~12 min. Negligible.

**Decision needed:** approve 32-frame sweep at 512², or prune to
exactly 15 points (and which one to drop).

### 3.2 Vertical-line annotation on the final figure

`paper/scripts/make_figures.py` (Phase 3D) draws a dashed vertical
line at `ε₃_crit = 0.1212` so the cusp is unambiguous in print. The
sweep itself does not annotate; it just produces the data.

---

## 4. CSV schema

Path: `output/sweep.csv` (gitignored — see §5).

Header (line 1, fixed):

```
epsilon_3,inclination_deg,resolution,mean_r,delta_r_rms,delta_r_p2p,fourier_m1_amp,fourier_m1_phase,multipeak_bins,render_wallclock_s,git_sha,timestamp_iso
```

Field reference:

| Column | Type | Source | Format / units |
|---|---|---|---|
| `epsilon_3` | double | grid input | `%+.4f` (e.g. `+0.1200`, `−2.5000`) |
| `inclination_deg` | double | grid input | `%.1f` (e.g. `17.0`) |
| `resolution` | int | grid input | `512` (sweep) or `64`–`256` (IT) |
| `mean_r` | double | `CircularityMetric.Result.meanR` | `%.6e`, units M |
| `delta_r_rms` | double | `Result.deltaRrms` | `%.6e`, units M |
| `delta_r_p2p` | double | `Result.deltaRp2p` | `%.6e`, units M |
| `fourier_m1_amp` | double | `Result.fourierM1Amplitude` | `%.6e`, units M |
| `fourier_m1_phase` | double | `Result.fourierM1Phase` | `%.6e`, radians |
| `multipeak_bins` | int | `RingExtractor.multipeakBinCount` | integer count, `[0, numBins]` |
| `render_wallclock_s` | double | `System.nanoTime` delta around the renderer call | `%.3f` |
| `git_sha` | string | `git rev-parse --short=7 HEAD` (+ `-dirty` if applicable) | 7 chars or `7chars-dirty` |
| `timestamp_iso` | string | `Instant.now().toString()` UTC | ISO 8601 (e.g. `2026-04-28T15:30:00Z`) |

Total: **12 columns**. Parent plan §6.6 listed 10; the two additions
(`fourier_m1_phase`, `multipeak_bins`) are required by §2.4
(diagnostic) and §6.7 gate 1 isotropy check (Fourier phase must be
defined for the synthetic-circle test even if amplitude is zero).

`fourier_m1_amp` is the absolute value of the discrete-Fourier m=1
coefficient computed over `numBins = 180` bins:

```
c_k = (1/N) Σᵢ rᵢ exp(−2π i k iₛ / N)         iₛ = bin index
fourierM1Amp = |c₁|
fourierM1Phase = arg(c₁)
```

`fourier_m1_amp / mean_r` is the EHT-comparable Fourier m=1
fractional amplitude; the ratio is computed in 3D
`ConsistencyBound`, not stored.

Empty bins (NaN in `RingExtractor.extract` output) are **excluded**
from `mean_r`, `delta_r_rms`, `delta_r_p2p`, and Fourier sums.
`validBins` in the `Result` record reports the count, but is not
persisted to CSV (derivable from `multipeak_bins` and `numBins`).

---

## 5. Resumability design

### 5.1 File policy

`.gitignore` currently has `output/*.fits` and `output/*.png`. CSV
needs to be added. Modification in 3C.1:

```diff
 output/*.fits
 output/*.png
+output/*.csv
 !output/.gitkeep
 !output/first_shadow.png
```

`output/sweep.csv` is regenerable from any `phase-3c-sweep`-tagged
checkout; not committed.

### 5.2 Restart algorithm

```
on EpsilonSweep.runSweep(cfg):
  1. ensure cfg.csvPath parent exists
  2. if cfg.csvPath does not exist:
       write header line; treat completedKeys as empty
     else:
       read CSV; verify header matches §4 exactly (fail fast on mismatch)
       parse all data rows into completedKeys = set of (eps3, incl, res) triples,
         with eps3 keyed by `String.format("%+.4f", v)` for exact match
  3. for (eps3, incl) in cartesian(cfg.epsilonGrid, cfg.inclinationDegrees):
       key = (formattedEps3, formattedIncl, cfg.resolution)
       if key in completedKeys: continue
       build JohannsenPsaltisMetric(spin, eps3)
       build NovikovThorneDisk(metric, cfg.rOuter, 0.001)
       build Camera(metric, incl)
       build AdaptiveRayTracer(metric, disk)
       wallStart = System.nanoTime()
       render image (FITS) → path_for_key
       wallEnd = System.nanoTime()
       radii = ringExtractor.extract(image)
       result = CircularityMetric.compute(radii)
       multipeak = ringExtractor.multipeakBinCount(image, ...)
       gitSha = currentGitSha()    // includes -dirty if applicable
       append CSV row, fsync, close
       completedKeys.add(key)
  4. log "sweep complete: N rendered, M skipped"
```

### 5.3 Crash invariants

- **CSV header mismatch on resume → fail fast.** Schema migration is
  out of scope for 3C.
- **FITS file present but no CSV row → re-render.** A crash between
  FITS write and CSV append loses ~24 s of work; acceptable. The
  re-rendered FITS overwrites bit-for-bit (same metric, same
  numerics, deterministic).
- **CSV row present but no FITS file → keep CSV row, skip render.**
  The CSV is the figure-of-merit; FITS files are reproducible
  artifacts. If the user wants the FITS regenerated, they delete the
  CSV row.
- **Concurrent runs against the same CSV are unsupported.** No file
  locking. Documented in `EpsilonSweep` Javadoc.

### 5.4 Per-frame FITS path

```
output/sweep_a0.9_eps{epsilon_3:+.4f}_i{inclination_deg:.1f}_res{resolution}.fits
```

Example: `output/sweep_a0.9_eps+0.1200_i17.0_res512.fits`. Stable
across re-runs at the same key.

---

## 6. Sub-phase decomposition (3 commits, 1 tag)

### 6.1 Commit 3C.1 — Analysis core + Maven slow-test infrastructure

Self-contained: no renderer or sweep code yet. Pure analysis +
build config.

**Files created:**

| File | Est. lines | Purpose |
|---|---|---|
| `src/main/java/com/pranav/grrt/analysis/RingExtractor.java` | ~180 | Image-space geometry, per-bin peak finder, multi-peak counter |
| `src/main/java/com/pranav/grrt/analysis/CircularityMetric.java` | ~90 | Pure stats: mean, RMS dispersion, peak-to-peak, Fourier m=1 |
| `src/test/java/com/pranav/grrt/analysis/RingExtractorTest.java` | ~200 | Gates 1, 2, 6 |
| `src/test/java/com/pranav/grrt/analysis/CircularityMetricTest.java` | ~140 | Synthetic radius arrays: pure-m1, pure-m2, constant, two-mode mix |

**Files modified:**

| File | Modification |
|---|---|
| `pom.xml` | Add surefire `excludedGroups=slow` and `runSlow` profile per §2.3 |
| `.gitignore` | Add `output/*.csv` per §5.1 |

**Test budget:** `mvn test` < 30 s (existing 28.6 s + ~1 s for new
fast tests). No slow tests yet.

### 6.2 Commit 3C.2 — Sweep driver + integration tests

**Files created:**

| File | Est. lines | Purpose |
|---|---|---|
| `src/main/java/com/pranav/grrt/analysis/EpsilonSweep.java` | ~320 | Orchestrator: CSV r/w, key matching, per-frame loop, `main(String[])` |
| `src/test/java/com/pranav/grrt/analysis/EpsilonSweepIT.java` | ~280 | Gates 3, 4 (non-blocking), 7. `@Tag("slow")` |

**Files modified:** none.

**Test budget:**
- `mvn test`: unchanged (no new fast tests).
- `mvn verify -PrunSlow`: ~75 s additional (gate 3 ~50 s, gate 4
  ~15 s, gate 7 ~3 s, fixed overhead ~7 s).

If 3C.2 ends up under ~150 lines added (unlikely), fold into 3C.1.

### 6.3 Commit 3C.3 — Run production sweep, post-validate, tag

**Files modified:**

| File | Modification |
|---|---|
| `CLAUDE.md` | Tick `EHT comparison pipeline` checkbox per `phase-3-plan.md` §6.2 + `phase-3b-status.md` §6.2 |

**Manual / non-source steps:**

1. `mvn package` (build the executable).
2. `mvn exec:java -Dexec.mainClass=com.pranav.grrt.analysis.EpsilonSweep`
   (or equivalent entry point — exact invocation locked at coding
   time). Runtime ≈ 13 min.
3. Visual inspection of the 32 emitted FITS frames — particularly
   `eps+0.1200`, `eps+0.1300`, `eps+0.2000`, `eps+1.0000` at
   `i17.0` to confirm prograde-shadow structure matches §9.2
   expectations.
4. **Gate 5 monotonicity inspection (manual):** plot
   `delta_r_rms / mean_r` vs `epsilon_3` on the smooth segment
   `[−2.5, +0.10]` at `i17.0`. Verify monotonicity. If
   non-monotonic, halt and investigate before tagging. Optional
   automation: `scripts/check_monotonicity.py` — defer if hand
   inspection is enough.
5. `git tag phase-3c-sweep && git push origin main --tags`.

`output/sweep.csv` is **not** committed (§5.1).

---

## 7. Validation gates

Seven gates per `phase-3-plan.md` §6.7, refined for the three-class
split and adjusted tolerances.

| # | Test | Tolerance | Lives in | Slow? |
|---|---|---|---|---|
| 1 | RingExtractor + CircularityMetric round-trip on synthetic perfect circle (rasterized at `r = 5.5 M` on 512²) | `delta_r_rms / mean_r < 1e-10` | `RingExtractorTest` | no |
| 2 | RingExtractor + CircularityMetric round-trip on synthetic ellipse (`a = 5.5 M`, `b = 5.2 M`, 512²) | `delta_r_rms / mean_r` matches analytic to 1e-4 relative | `RingExtractorTest` | no |
| 3 | **Kerr/JP consistency** at `(a=0.9, ε₃=0, i=17°, 512²)`: JP sweep row vs direct Kerr render of same frame | `\|Δ(δ_r/⟨r⟩)\| < 0.1 pp` (parent plan correction 4) | `EpsilonSweepIT` | yes |
| 4 | External shadow-diameter cross-check: `KerrMetric(a=0.9)` shadow at `i=17°` vs published GRay value (Chan et al. 2013 Tab. 1) | within 2% (secondary, **non-blocking** — flag failure, do not fail the build) | `EpsilonSweepIT` | yes |
| 5 | Smooth-regime monotonicity: `δ_r/⟨r⟩` monotonic in `ε₃` on `[−2.5, +0.10]` at `i = 17°` from production CSV | manual inspection (§6.3 step 4); knee at `ε₃ ≈ +0.12` is EXPECTED, not a failure | n/a (post-sweep) | manual |
| 6 | Bin-count convergence: 90/180/360 bins on a single saved 256² Kerr+NT frame | 5% relative agreement on `δ_r/⟨r⟩` | `RingExtractorTest` | no |
| 7 | Resumability: 3-frame synthetic sweep at 64², kill after frame 1, restart, complete; final CSV identical to single-shot run | exact on numeric columns; `timestamp_iso` and `render_wallclock_s` may differ | `EpsilonSweepIT` | yes |

Plus six `CircularityMetric`-isolated unit tests (not numbered as
gates; in `CircularityMetricTest`):

| Synthetic input | Expected |
|---|---|
| Constant `rᵢ = 5.5` for all bins | `delta_r_rms = delta_r_p2p = fourier_m1_amp = 0` (all 1e-15) |
| Pure m=1: `rᵢ = R₀(1 + A cos(2π i/N))` with `R₀ = 5.5`, `A = 0.05` | `delta_r_rms = R₀ A / √2`, `fourier_m1_amp = R₀ A / 2` (1e-12) |
| Pure m=2: `rᵢ = R₀(1 + A cos(4π i/N))` | `delta_r_rms = R₀ A / √2`, `fourier_m1_amp ≈ 0` (1e-12) |
| Two-mode mix: m=1 at A₁=0.05, m=2 at A₂=0.03 | RMS = R₀√(A₁²+A₂²)/√2; m=1 amplitude isolates A₁ (1e-12) |
| With NaN-filled bins (10 of 180) | `validBins = 170`; statistics computed on 170 only |
| All-NaN array | `validBins = 0`; all metrics NaN; no exception |

### 7.1 Gate 3 numerical detail (most important)

At `ε₃ = 0`, `JohannsenPsaltisMetric(0.9, 0)` should reduce
component-wise to `KerrMetric(0.9)` at machine precision (verified
in 3A gate 3). Therefore any 512² render through both metrics with
the same camera, integrator tolerances, and disk should produce
**bit-identical** images, modulo floating-point reassociation in
the JP closed-form contractions.

The 0.1 pp tolerance budgets up to ~5e-4 relative drift in
`δ_r/⟨r⟩` (assuming `δ_r/⟨r⟩ ~ 0.05` for face-on Kerr at i=17°).
That margin absorbs any reassociation drift. If gate 3 fails, the
likely cause is a bug in `JohannsenPsaltisMetric.geodesicAcceleration`
that 3A grid-tests didn't catch — investigate before continuing
the sweep.

### 7.2 Gate 4 reference value

Chan, Psaltis, Özel 2013 Tab. 1: Kerr `a = 0.9`, `i = 17°` shadow
half-diameter ≈ 5.07 M (cross-checked against `KerrShadowTest` gate
in Phase 2). Gate 4 compares the mean `RingExtractor` radius on a
disk-emission-only Kerr frame (no JP) to this published value. The
2% tolerance accommodates differences between "shadow boundary"
(GRay, photon-orbit-grazing) and "ring peak" (RingExtractor,
intensity-weighted). Non-blocking: if gate fails, log a WARN and
proceed; the bound extraction does not depend on this comparison.

---

## 8. Wall-clock estimates (M3, 10 threads)

### 8.1 Test suite

| Stage | Wall-clock |
|---|---|
| `mvn test` after 3C.1 (108 + ~30 fast tests) | ~30 s |
| `mvn test` after 3C.2 (slow IT excluded) | ~30 s |
| `mvn verify -PrunSlow` after 3C.2 | ~30 s + 75 s = ~105 s |

### 8.2 Production sweep (3C.3)

| Per-frame breakdown | Time |
|---|---|
| `AdaptiveRayTracer` + `NovikovThorneDisk` render at 512² | ~24 s |
| `RingExtractor.extract` (180 bins, image scan) | ~0.5 s |
| `RingExtractor.multipeakBinCount` | ~0.3 s |
| `CircularityMetric.compute` (180-element FFT-free m=1) | < 1 ms |
| FITS write | ~50 ms |
| CSV append + fsync | ~5 ms |
| **Per-frame total** | **~25 s** |

| Sweep total | Frames | Wall-clock |
|---|---|---|
| 16 ε₃ × 2 inclinations | 32 | **~13 min** |

### 8.3 Coding effort

Per `phase-3b-status.md` §6.5 estimate: ~3 focused sessions, ~6–10
hours session time, plus ~13 min sweep wall-clock. Unchanged.

This plan does not change the 3C calendar projection.

---

## 9. Failure modes

### 9.1 General (parent plan §6.9)

| Failure | Detection | Mitigation |
|---|---|---|
| "Photon ring" ambiguity (critical curve vs resolved ring) | n/a | D-C3 locked: dominant bright ring (n=0 + n=1 unresolved); documented in `RingExtractor` Javadoc |
| Centroid drift with large `\|ε₃\|` | Visual inspection of post-cusp frames | `RingExtractor` accepts a fixed centroid; pass intensity-weighted centroid from a one-pass image scan with 3σ clip — falls back to image origin if clip diverges |
| Azimuthal binning bias | Gate 6 (bin convergence at 90/180/360) | 180 bins default; gate fails the build if convergence > 5% relative |
| Multi-valued azimuthal profile (n=1 ring crosses direct image) | `[ring-warn]` print + `multipeak_bins` CSV column | Outermost peak per bin (§2.4); audit affected frames manually |
| Sweep crash mid-run | Gate 7 (resumability IT) | Idempotent CSV; FITS re-rendered if missing, CSV row preserved if present (§5.3) |
| `output/` not gitignored for CSV | `.gitignore` modification in 3C.1 | Spot-check `git status` after first sweep run |
| RingExtractor returns NaN/zero radii on a high-`\|ε₃\|` frame | Per-frame `validBins` check; `EpsilonSweep` aborts if `validBins < 0.5 · numBins` | Investigate `RingExtractor` peak threshold or interpret as post-cusp shadow-edge regime (§9.2) |

### 9.2 Cusp-specific failure modes (`ε₃ ≈ +0.12` and beyond)

This is the headline risk for 3C. Per `docs/jp-parameter-space-notes.md`,
the prograde Kerr-continuation photon orbit at `a = 0.9` merges with
the Kerr horizon at `ε₃_crit ≈ 0.1212`. The sweep grid samples this
transition densely from below (`+0.10, +0.11, +0.12`) and then
crosses it (`+0.13, +0.20, +1.00`).

| ε₃ region | Expected physical state | Expected numerical state | Watch for |
|---|---|---|---|
| `[+0.05, +0.11]` (smooth, near-cusp) | Prograde photon orbit shrinks: `r_ph ∈ [1.515, 1.451]` (parent §2.2 of `jp-parameter-space-notes.md`) | All gates pass; `δ_r/⟨r⟩` monotone in ε₃ | Sharp gradient in `δ_r/⟨r⟩` near ε₃ = +0.11; smooth, not discontinuous |
| `+0.12` (just below `ε₃_crit`) | `r_ph ≈ 1.4377 M`, margin ~0.0018 M to `r₊_K = 1.4359 M` | Photons grazing this orbit have very narrow capture conditions; integrator must resolve the ~0.13% margin | Watch for `RayOutcome.MAX_STEPS` count > a few % of pixels; if so, tighten `atol`/`rtol` per-pixel for this frame only (does not break determinism) |
| `+0.13` | `r_ph` does not exist; prograde shadow bounded by JP equatorial horizon | First end-to-end exercise of the 3B `isInsideHorizon` override (per `phase-3b-status.md` §6.4 — "de-facto smoke test") | NaN propagation, unexpected `MAX_STEPS` rates, prograde-side per-bin radii landing on the disk inner edge (~`r_ISCO`) instead of a photon ring |
| `+0.20`, `+1.00` | Same as `+0.13`, deeper into transition regime | Asymmetric shadow with retrograde photon ring at `r ≈ 3.88 M` and prograde horizon edge | Multi-peak count likely jumps on prograde bins (disk inner edge + horizon edge); `δ_r/⟨r⟩` jumps to 0.5+ regime |

**Specific recovery paths:**

1. **NaN in any output column at any post-cusp ε₃.** Investigate
   `JohannsenPsaltisMetric.isInsideHorizon` first — verify
   `Δ(r) + a² h(r,θ) sin²θ` evaluation matches the closed form in
   `phase-3-plan.md` §2.5. The 3B unit tests covered the override
   at the equator and axis only; this is the first full-coverage
   stress. Do not relax tolerances; root-cause.

2. **`MAX_STEPS` > 10% of pixels at `ε₃ = +0.12`.** Likely the
   integrator is taking many tiny steps near the photon-orbit /
   horizon coincidence point. Two responses:
   - Increase `maxSteps` (currently set in `AdaptiveRayTracer`
     constructor) for that frame only.
   - Tighten `atol`/`rtol` from default `1e-10` to `1e-12` for that
     frame only.
   Do not relax. Either response is per-frame; the rest of the
   sweep keeps the default.

3. **`RingExtractor` returns NaN/zero on prograde bins for
   `ε₃ ≥ +0.13`.** This is a semantic mismatch, not a bug. The
   "ring" no longer exists on the prograde side; what's there is
   the disk inner edge. Two responses, in order:
   - Confirm visually (open the FITS file) that the prograde shadow
     boundary is at the JP equatorial horizon and not pathological.
   - Decide whether to extract "the brightest pixel per bin" (which
     becomes the disk inner edge for these bins) or to leave the
     bin NaN. Recommendation: **extract the brightest pixel per
     bin always**; the resulting per-bin profile is exactly what
     the headline figure needs, and its `δ_r/⟨r⟩` correctly
     captures the qualitative shadow-structure change. The
     `multipeak_bins` count will reflect the regime change.

4. **Gate 5 (manual monotonicity) fails on the smooth segment
   `[−2.5, +0.10]`.** Halt before tagging `phase-3c-sweep`. Likely
   causes (in priority order):
   - JP `geodesicAcceleration` bug in a regime not covered by 3A
     gates. Re-run 3A `Kerr-reduction` test at `ε₃ = 0` and the
     `closed-form prediction` test at the failing ε₃.
   - RingExtractor centroid-clip bug. Run gate 6 on the failing
     frame at multiple bin counts; if convergence breaks, root-cause.
   - DP45 step-size adapter pathology at a specific `(ε₃, r)`.
     Re-run with tighter `atol`.

5. **Visual inspection (3C.3 step 3) reveals the `+1.00` frame is
   visually pathological** (ringing artifacts, shadow fragmentation,
   etc.). Per `jp-parameter-space-notes.md` §3.1, the metric is
   regular and CTC-free at `(0.9, +1.00)`. Pathology in the image
   would point to a numerical issue — most likely insufficient
   integrator robustness in the high-`h` regime. Triage path:
   review `1 + h` minimum across the frame (should be ~0.34 per
   §3.1 analysis); confirm `Δ + a² h sin²θ` does not change sign
   except at the horizon. Do not silently accept the frame.

### 9.3 New failure modes added in 3C

| Failure | Mitigation |
|---|---|
| `output/sweep.csv` from a prior schema is stale | Header validation on resume (§5.3) — fail fast |
| Two `EpsilonSweep` processes against the same CSV | Document as unsupported; no locking |
| `git rev-parse` not available (e.g. tarball checkout) | Fall back to `"unknown"` for `git_sha` |
| Resume after a `git rebase` changes the recorded SHA mid-sweep | Acceptable: each row records the SHA valid at its render time. Audit trail intact. |

---

## 10. Commit + tag structure

```
[3C.1] analysis: add RingExtractor and CircularityMetric with slow-test infrastructure
[3C.2] analysis: add EpsilonSweep driver with resumable CSV and integration tests
[3C.3] sweep: render 32-frame ε₃ × inclination sweep; close Phase 3C
       (followed by `git tag phase-3c-sweep`)
```

Author: configured git identity, no attribution trailers per
CLAUDE.md "Agent behavior expectations".

If 3C.2 ends up under ~150 lines, fold into 3C.1.

---

## 11. Pre-coding checklist

Each row needs explicit user confirmation before 3C.1 lands.

| # | Item | §  | Default proposal |
|---|---|---|---|
| 1 | Three-class split: `RingExtractor`, `CircularityMetric`, `EpsilonSweep` | 2.1 | Approve split |
| 2 | NT `surfaceFlux` caching: profile-then-decide | 2.2 | Defer to mid-3C.3 unless user wants pre-emptive cache |
| 3 | Slow-test surefire group key | 2.3 | `slow` |
| 4 | Multi-peak threshold | 2.4 | 2 px |
| 5 | `git_sha` length and dirty handling | 2.5 | 7 chars + `-dirty` suffix |
| 6 | 16-point sweep grid (vs 15-point literal) | 3.1 | Render all 16 (32 frames, ~13 min) |
| 7 | CSV column count (12, including `fourier_m1_phase` and `multipeak_bins`) | 4 | Approve 12-column schema |
| 8 | `output/*.csv` gitignore addition | 5.1 | Approve |
| 9 | Tick `EHT comparison pipeline` checkbox at end of 3C | 6.3 | Approve (per `phase-3b-status.md` §6.2) |
| 10 | Gate 5 (smooth-regime monotonicity) as manual inspection only | 6.3, 7 | Approve; defer `scripts/check_monotonicity.py` |

---

## 12. Pre-3C reading order (for fresh sessions)

1. `docs/phase-3b-status.md` — Phase 3B close-out and 3C readiness
   check (§6 in particular).
2. `docs/phase-3c-plan.md` — this file (the authoritative 3C plan).
3. `docs/phase-3-plan.md` §6 — parent plan, superseded by this
   sub-plan where they overlap.
4. `docs/jp-parameter-space-notes.md` §§4–5 — cusp physics + sweep
   grid rationale, especially for understanding the post-`ε₃_crit`
   image regime.
5. `CLAUDE.md` — repo rules. No `Metric` interface change in 3C, so
   the "ask when touching Metric" caveat does not bind.

---

## 13. Deferred / out of scope for 3C

- **`ConsistencyBound` class and EHT-Paper-VI ingestion.** Phase 3D.
  3C produces the data; 3D extracts the bound.
- **Frequency-channel emission** beyond bolometric. Phase 3D / out
  of scope for the RNAAS.
- **Per-pixel adaptive tolerance.** Default `atol`/`rtol = 1e-10`
  uniform across pixels for the production sweep. Per-frame override
  is allowed (§9.2 recovery path 2) but not pre-emptive.
- **NT spline-cache constructor argument.** Default off; added in a
  4th commit only if §2.2 profile result demands it.
- **Multi-process sweep parallelism.** Single-process,
  `parallelStream` over pixels only. The bulk is renderer-bound; CSV
  serialization is fast enough.
- **`paper/` directory.** Phase 3D.
- **External shadow-diameter cross-check expansion** beyond the
  single Chan-Psaltis-Özel reference value. Gate 4 is a sanity
  check, not a literature survey.

---

*End of sub-plan.*
