# Phase 3 Plan (approved 2026-04-24)

Authoritative reference for Phase 3 of grrt. This document supersedes any
conflicting summary recovered from `/compact` output or prior session
transcripts. The plan is approved; execution order is fixed.

End goal: an RNAAS Research Note measuring the M87* photon-ring asymmetry
constraint on a single-parameter Johannsen-Psaltis deviation from Kerr.

---

## 1. Fixed scientific decisions

Locked by user at Phase 3 kickoff. Do not revisit without explicit
re-approval.

| Decision | Value | Source / rationale |
|---|---|---|
| Deformation parameter | `ε₃` | JP 2011 notation retained |
| Deformation form | `h(r, θ) = ε₃ M³ r / Σ²` | Johannsen 2013 / Psaltis et al. 2020; θ-dependent (D-A1) |
| Spin | `a = 0.9 M` | High-spin regime, M87*-consistent |
| Primary inclination | `i = 17°` | M87* jet axis (Walker et al. 2018) |
| Reference inclination | `i = 60°` | Methodological comparator on final figure |
| Scientific bound | From 17° curve only | 60° curve is illustrative only |
| Disk | Novikov-Thorne thin, optically thick | `r ∈ [r_ISCO(a, ε₃), 20 M]` |
| Asymmetry metric (primary) | `δ_r/⟨r⟩ ≡ √mean((rᵢ − ⟨r⟩)²) / ⟨r⟩` | RMS; correction 3 |
| Asymmetry metric (secondary) | peak-to-peak, Fourier m=1 | Stored in CSV for robustness |
| EHT benchmark | Fourier m=1 fractional amplitude, Paper VI Table 7 | D-D1; Psaltis 2020 comparator |

---

## 2. Metric form (D-A1, Johannsen 2013 ansatz)

### 2.1 Structural ansatz

Johannsen & Psaltis 2011, Phys. Rev. D **83**, 124015, equation (1),
with the deformation scalar taken from Johannsen 2013
(Phys. Rev. D **88**, 044002):

```
h(r, θ) = ε₃ M³ r / Σ²                  (dimensionless)
Σ(r, θ) = r² + a² cos²θ
Δ(r)    = r² − 2Mr + a²
```

### 2.2 Metric components (Boyer-Lindquist-like, signature −+++)

```
g_tt  = −(1 − 2Mr/Σ)(1 + h)
g_tφ  = −(2 M a r sin²θ / Σ)(1 + h)
g_rr  =  Σ (1 + h) / (Δ + a² h sin²θ)
g_θθ  =  Σ
g_φφ  =  sin²θ · [ (r² + a²) + (2 M a² r sin²θ / Σ)(1 + h) ]
```

All other components are zero by stationarity (`∂_t`) and axisymmetry
(`∂_φ`).

### 2.3 Inverse metric (closed form)

The `(r, θ)` block is diagonal; the `(t, φ)` block is 2×2 and inverted
analytically. Let

```
D = g_tt · g_φφ − g_tφ²
```

Then

```
g^{tt}  =  g_φφ / D
g^{φφ}  =  g_tt / D
g^{tφ}  = −g_tφ / D
g^{rr}  =  (Δ + a² h sin²θ) / [Σ (1 + h)]
g^{θθ}  =  1 / Σ
```

No runtime 4×4 matrix inversion. Closed-form expressions will be emitted
by `scripts/derive_jp_christoffels.py`.

### 2.4 Kerr reduction

Setting `ε₃ = 0` implies `h = 0` and every component above collapses
to the Kerr form implemented in `KerrMetric.java`. This is tested at
machine precision in 3A gate 3.

### 2.5 Horizon location (θ-dependent)

The event horizon is the outermost root of `g_rr → ∞`, which is

```
Δ(r) + a² h(r, θ) sin²θ = 0
```

On the polar axis (`sin θ = 0`) this reduces to `Δ = 0` → Kerr horizon
`r₊ = M + √(M² − a²) ≈ 1.4359 M` for `a = 0.9 M`. Off-axis, the
horizon radius depends on `ε₃`:

- `ε₃ > 0`: equatorial horizon moves inward (`a² h > 0`)
- `ε₃ < 0`: equatorial horizon moves outward

`Metric.horizonRadius()` returns a single scalar. For JP we will return
the θ-maximum radius (equatorial for `ε₃ < 0`, polar for `ε₃ ≥ 0`) as
the safe outer bound. Actual ray termination in `AdaptiveRayTracer`
will use a position-dependent test `Δ(r) + a² h(r, θ) sin²θ ≤ tol`
rather than a fixed `r`. This is a 3A design note; see section 4.1.

---

## 3. Coordinate pathology verification (gate 2)

**Claim:** `1 + h(r, θ) > 0` for all `ε₃ ∈ [−2, +2]`, `a = 0.9 M`,
`r ≥ r₊ = 1.4359 M`, `θ ∈ [0, π]`.

### 3.1 Analytic argument

Write `h = ε₃ · r / Σ²` (in `M = 1`). Take `η ≡ cos²θ ∈ [0, 1]`,
so `Σ = r² + a² η`. Then

```
∂h/∂η = −2 ε₃ · r · a² / (r² + a² η)³
```

which is sign-definite in `η` at fixed `(r, ε₃, a)`. So `h` is monotonic
in `η`, and `|h|` achieves its extrema at `η = 0` (equator) or `η = 1`
(pole). Checking both:

**Equator** (`η = 0`, `Σ = r²`):
```
h_eq(r) = ε₃ · r / r⁴ = ε₃ / r³
|h_eq| is maximized at minimum r = r₊ = 1.4359 M:
|h_eq|_max = |ε₃| / r₊³ = 2 / 2.962 = 0.6753
```
Therefore `1 + h_eq ∈ [0.3247, 1.6753]`. ✓

**Pole** (`η = 1`, `Σ = r² + a²`):
```
h_pol(r) = ε₃ · r / (r² + a²)²
d(r / (r² + a²)²)/dr = [(r² + a²) − 4r²] / (r² + a²)³ = (a² − 3r²)/(r² + a²)³
```
Zero at `r = a/√3 = 0.520 M`, which lies inside `r₊`. So for
`r ∈ [r₊, ∞)` the function `h_pol(r)` is monotonically decreasing in `r`,
and its maximum magnitude is at `r = r₊`:
```
|h_pol|_max = 2 · 1.4359 / (1.4359² + 0.81)² = 2 · 1.4359 / 2.872² = 0.3479
```
Therefore `1 + h_pol ∈ [0.6521, 1.3479]`. ✓

**Conclusion.** Across the entire sweep range (ε₃ ∈ [−2, +2],
a = 0.9 M, r ≥ r₊, θ ∈ [0, π]), the factor `1 + h` stays in
`[0.3247, 1.6753]`. It never reaches zero; there is no `(1 + h)`-driven
coordinate pathology. ✓

### 3.2 Reproducibility

A numerical cross-check over a grid is embedded in
`scripts/derive_jp_christoffels.py` (gate 3 artifact). Running the
script prints `coordinate pathology grid check: PASSED (min 1+h = ...)`
before the symbolic derivation begins. The analytic bound above is
the authoritative record; the grid check is belt-and-suspenders.

### 3.3 Separate concern: `Δ + a² h sin²θ` sign

`1 + h > 0` does **not** rule out pathology in `g_rr`, which also
depends on `Δ + a² h sin²θ`. For `ε₃ < 0` this expression vanishes
at `r > r₊_Kerr` on the equator (the JP horizon sits outside the Kerr
horizon). This is **not** a coordinate pathology in the problematic
sense — it is the true JP event horizon, and ray integration should
terminate there. See section 4.1 for the termination strategy. The
present gate is specifically about `1 + h > 0`, which has been
verified.

---

## 4. Sub-phase 3A: JohannsenPsaltisMetric

### 4.1 Scope

Single-file `Metric` implementation for Johannsen-Psaltis deformed
Kerr with the `h(r, θ) = ε₃ M³ r / Σ²` scalar. The `Metric` interface
is unchanged (honors CLAUDE.md "ask when touching Metric" rule).

### 4.2 New files

```
src/main/java/com/pranav/grrt/metric/JohannsenPsaltisMetric.java
src/test/java/com/pranav/grrt/metric/JohannsenPsaltisMetricTest.java
scripts/derive_jp_christoffels.py           # committed in gate 3, pre-3A
```

### 4.3 New interfaces

None. JP is a drop-in `Metric`.

### 4.4 Files modified

- `CLAUDE.md` at 3A completion: tick `JohannsenPsaltisMetric` checkbox;
  add a JP landmarks subsection (horizon θ-dependence, equatorial photon
  orbit closed form).

Nothing else — the architecture guarantees JP is a single-file
extension.

### 4.5 Christoffels (correction 1)

Derived symbolically via Sympy, not by hand.
`scripts/derive_jp_christoffels.py`:

1. Defines symbolic `t, r, θ, φ, M, a, ε₃`.
2. Constructs `h`, `Σ`, `Δ`, then the five non-zero `g_μν` components
   from section 2.2.
3. Computes the closed-form inverse (section 2.3) and verifies
   `g · g⁻¹ = I` symbolically.
4. Computes all 40 independent Christoffels
   `Γ^α_{μν} = ½ g^{ασ} (∂_μ g_{σν} + ∂_ν g_{σμ} − ∂_σ g_{μν})`
   (stationarity + axisymmetry set `∂_t = ∂_φ = 0`; only `∂_r` and
   `∂_θ` survive).
5. Applies `sympy.cse` for common-subexpression elimination and emits
   Java-paste-ready expressions keyed to a pre-allocated scratch array.
6. At `ε₃ = 0`, prints a numerical spot-check against hand-derived
   Kerr Christoffels at `(r=10M, θ=π/3)` and asserts agreement to
   1e-12. This catches transcription errors in the Sympy input.

The CLAUDE.md "derivation first" rule is satisfied by the committed
script plus a pointer from `JohannsenPsaltisMetric`'s Javadoc. The
manual hand-derivation that KerrMetric embeds is replaced by the
algorithmically-generated, CSE-optimized Sympy output for JP.

**Python environment:** the script requires `python3` and `sympy`.
Not a runtime dependency of the Java project; Maven build is
unaffected. Script runtime: ~15–30 s on M3.

### 4.6 `geodesicAcceleration` optimization

After Sympy emits the Christoffels, hand-fuse the non-zero contractions
into the four output components (same approach as
`KerrMetric.geodesicAcceleration`). Regression-test the optimized
version against the default contraction-through-`christoffel()` at
1e-12 on a 16-point `(r, θ, k)` grid.

### 4.7 Validation gates (JUnit 5, AAA pattern)

| # | Test | Tolerance |
|---|---|---|
| 1 | `g · g⁻¹ = I` on 64-point `(r, θ, ε₃)` grid | 1e-12 |
| 2 | Christoffel symmetry `Γ^α_{μν} = Γ^α_{νμ}` | 1e-12 |
| 3 | **Kerr reduction**: `JP(a=0.9, ε₃=0)` vs `KerrMetric(a=0.9)` — `g`, `g⁻¹`, `Γ` componentwise on shared grid | 1e-12 |
| 4 | Asymptotic flatness: `g_μν(r=1000 M)` vs Minkowski | 1e-6 |
| 5 | Null-norm drift: DP45 photon orbit from `(r₀=50 M, b=5.5 M)` in `JP(0.9, +0.5)`, `\|g_μν k^μ k^ν\| / (k^t)²` over 1000 M affine length | < 1e-6 |
| 6 | **Closed-form prediction**: equatorial prograde circular photon orbit `r_ph(a, ε₃)` at six `(a, ε₃)` pairs incl. `(0.9, 0)`, `(0.9, ±0.5)`, `(0.9, ±1)` | 1e-10 |
| 7 | `geodesicAcceleration` optimized vs default contraction | 1e-12 |

Test 6 derives `r_ph` from the equatorial quartic
`(d/dr)[g_tt + 2 L g_tφ + L² g_φφ]_eq = 0` with the null condition,
solved numerically from a closed-form expression derived in the Sympy
script and pasted into the test. This satisfies CLAUDE.md test rule 4
("at least one closed-form physical prediction").

### 4.8 Wall-clock (M3)

`mvn test -Dtest=JohannsenPsaltisMetricTest` — budget < 5 s. No
rendering in 3A.

### 4.9 Primary failure modes

| Failure | Detection | Mitigation |
|---|---|---|
| Sign error in `h'(r)` or `∂h/∂θ` propagated into Γ | Kerr-reduction test (gate 3) catches any term that does not vanish at `ε₃ = 0` | Regenerate via Sympy; check Sympy inputs against section 2.2 |
| Coordinate-chart drift near JP horizon | Null-norm gate 5; ray termination guard | Use `Δ + a² h sin²θ ≤ tol` termination, not fixed `r` |
| NaN propagation from `1 + h` near zero | Analytic verification (section 3) rules this out in sweep range; also `IllegalArgumentException` guard | Fail loudly if caller passes `ε₃` outside approved range |
| `geodesicAcceleration` fused form wrong | Gate 7 against default | Regenerate from Sympy |

### 4.10 Tag at completion

`phase-3a-jp-metric`

---

## 5. Sub-phase 3B: Novikov-Thorne disk emission

Absorbs the Phase-2 leftover "disk emission model" checkbox — this is
the precondition for 3C.

### 5.1 New files

```
src/main/java/com/pranav/grrt/disk/Disk.java
src/main/java/com/pranav/grrt/disk/NovikovThorneDisk.java
src/main/java/com/pranav/grrt/disk/DiskEmissionShader.java
src/test/java/com/pranav/grrt/disk/NovikovThorneDiskTest.java
src/test/java/com/pranav/grrt/disk/DiskEmissionShaderTest.java
```

### 5.2 New interfaces

`Disk` (small, local to `disk/` package):

```java
public interface Disk {
    boolean crossedEquator(double[] xPrev, double[] xCurr);
    double[] keplerianFourVelocity(double[] xOnDisk, Metric m);
    double temperature(double r, Metric m);
    double rIsco();
    double rOuter();
}
```

No change to `Metric`, `Integrator`, or `Camera`.

### 5.3 Files modified

1. `src/main/java/com/pranav/grrt/integrator/DormandPrince45.java` —
   add `interpolate(double θ, double[] out)` for dense output
   (D-B1 approved). `step()` signature unchanged. Existing DP45
   tests must pass bit-exactly.
2. `src/main/java/com/pranav/grrt/renderer/AdaptiveRayTracer.java` —
   optional `Disk` injection; between accepted DP45 steps, test
   sign change on `θ − π/2`; on sign change, bisect with the new
   interpolant to find crossing, check `rIsco ≤ r(λ*) ≤ rOuter`,
   invoke shader.
3. `src/main/java/com/pranav/grrt/renderer/RayOutcome.java` — add
   `HIT_DISK` variant (enum additions only; existing callers
   handled via exhaustive switch + default).
4. `CLAUDE.md` — tick disk-model checkbox.

### 5.4 Emission model

Novikov-Thorne 1973 thin optically-thick disk in the equatorial plane
of the current `Metric`. Surface radial flux from Page-Thorne 1974
eq. (15n), using the metric's own `Metric.g(x)` to build the
angular-velocity and specific-energy profiles — no Kerr-specific code
paths. Local effective temperature:

```
T(r) = [ F(r) / σ_SB ]^(1/4)
```

Shader emission:

```
I_obs = g⁴ · B_bol(T(r_emit))      (bolometric)
g     = −k_μ u_obs^μ / (−k_ν u_emit^ν)
```

`B_bol(T) = σ_SB T⁴ / π`. Frequency-channel extension is deferred.

### 5.5 Validation gates

| # | Test | Tolerance |
|---|---|---|
| 1 | NT radial flux `F(r)` vs Page-Thorne 1974 eq. (15n) at 6 radii, `a ∈ {0, 0.9}` | 1e-6 |
| 2 | ISCO: `a = 0.9` prograde → 2.3209 M (Bardeen 1972 Table 1) | 1e-6 |
| 3 | Face-on Schwarzschild disk at 256²: symmetric annulus, centroid at origin | 0.5 px |
| 4 | Inner-ring redshift, `a=0.9`, `i=85°` edge-on vs analytic `g = 1/√(−g_tt_eq(r_ISCO))` | 1e-4 |
| 5 | DP45 dense-output interpolant at `θ_mid` vs DP45 half-step direct integration on Kerr photon orbit | 1e-8 |
| 6 | All prior DP45 tests (Phase 2) pass bit-exactly | exact |

### 5.6 Wall-clock (M3, 10 threads)

Phase-2 Kerr + BinaryShader, 256²: ~3 s. NT shader adds per-step
equatorial check + Planck evaluation at hit (~2×). Estimates:

- 256², `i = 17°`, Kerr + NT: ~6 s/frame
- 512², `i = 17°`, Kerr + NT: ~24 s/frame

These are per-frame; 3C multiplies by 26 sweep points.

### 5.7 Primary failure modes

| Failure | Mitigation |
|---|---|
| Missed equator crossing when DP45 step straddles disk | Cap step when `\|θ − π/2\| < 0.1 rad`; sign-change on both endpoints |
| Plunging-region emission | NT truncates at `r_ISCO`; sub-ISCO rays treated as reaching horizon (no emission) |
| Divergent Page-Thorne correction near ISCO | Terminate at `r_ISCO`, no smoothing |
| `AdaptiveRayTracer` regression breaking Phase-2 binary shader | `HIT_DISK` is a new outcome; existing shader paths unchanged; DP45 tests run bit-exactly (gate 6) |

### 5.8 Tag at completion

`phase-3b-nt-disk`

---

## 6. Sub-phase 3C: Ring asymmetry extraction + ε₃ sweep

### 6.1 New files

```
src/main/java/com/pranav/grrt/analysis/RingExtractor.java
src/main/java/com/pranav/grrt/analysis/EpsilonSweep.java
src/test/java/com/pranav/grrt/analysis/RingExtractorTest.java
src/test/java/com/pranav/grrt/analysis/EpsilonSweepIT.java
```

### 6.2 Files modified

- `pom.xml` — surefire `excludedGroups=slow` so `mvn test` stays fast
  and the full sweep runs under `mvn verify -PrunSlow`.
- `CLAUDE.md` — tick EHT comparison pipeline.
- `.gitignore` — ensure `output/` stays ignored (already is, confirm).

### 6.3 Asymmetry metric (correction 3)

**Primary:** RMS, defined as

```
δ_r/⟨r⟩ ≡ √( mean_i [ (rᵢ − ⟨r⟩)² ] ) / ⟨r⟩
```

where `rᵢ` is the peak-intensity radius in azimuthal bin `i` (180 bins
default) and `⟨r⟩ = mean_i rᵢ`.

**Secondary** (for robustness, stored in CSV):
- peak-to-peak: `(r_max − r_min) / ⟨r⟩`
- Fourier m=1 fractional amplitude: `|c₁| / c₀` where
  `cₖ = (1/N) Σᵢ rᵢ exp(−2π i k i / N)`

### 6.4 Ring definition (D-C3)

Dominant bright ring = direct image + unresolved n=1 lensed ring,
sampled at 180 azimuthal bins with per-bin peak-intensity radius. This
is explicitly **not** the asymptotic critical curve (n → ∞), which
would require sub-pixel-resolved 2048² + adaptive sampling and is out
of scope for this RNAAS. The distinction is recorded in both the
`RingExtractor` Javadoc and the paper's Methods section.

### 6.5 Sweep grid (correction 2)

```
ε₃ ∈ {−2, −1, −0.5, −0.2, −0.1, −0.05, 0, +0.05, +0.1, +0.2, +0.5, +1, +2}
     (13 points, denser near zero where the EHT bound is extracted)

i  ∈ {17°, 60°}                (both at 512²; correction to D-C2)
```

Total: 26 frames at 512².

### 6.6 Sweep driver (resumability requirement)

`EpsilonSweep` reads `output/sweep.csv` on startup. Each completed row
is keyed by `(ε₃, inclination_deg, resolution)`. Rows present are
skipped; rows absent are computed. Progress is flushed after every
frame. Crash/interrupt can be restarted without re-running completed
renders. CSV schema:

```
epsilon_3, inclination_deg, resolution,
mean_r, delta_r_rms, delta_r_p2p, fourier_m1,
render_wallclock_s, git_sha, timestamp_iso
```

`timestamp_iso` is ISO 8601 UTC (e.g. `2026-04-24T10:15:00Z`).

### 6.7 Validation gates

| # | Test | Tolerance |
|---|---|---|
| 1 | RingExtractor on synthetic perfect circle (`r = 5.5` on 512² raster) | `δ_r/⟨r⟩ < 1e-10` |
| 2 | RingExtractor on synthetic ellipse `(a=5.5, b=5.2)` vs analytic RMS | 1e-4 relative |
| 3 | **Kerr/JP consistency (correction 4)**: at `(a=0.9, ε₃=0, i=17°, 512²)`, JP sweep row vs direct Kerr render (same a, i, resolution) | `\|Δ(δ_r/⟨r⟩)\| < 0.1 pp` |
| 4 | External shadow-diameter cross-check vs GRay/Kerr published values | within 2% (secondary sanity, non-blocking) |
| 5 | Monotonicity: `δ_r/⟨r⟩` monotonic in `\|ε₃\|` at i=17° | manual inspection; non-monotonic halts sweep |
| 6 | Bin-count convergence: 90/180/360 bins on same frame | 5% relative agreement |
| 7 | Resumability: interrupt after 5 frames, restart, complete — final CSV identical to uninterrupted run (modulo `timestamp_iso`) | exact on numeric columns |

Gate 3 is the strong form requested in correction 4. Gate 4 is a
secondary sanity check only.

### 6.8 Wall-clock (M3, 10 threads)

- 26 frames × ~24 s/frame at 512² ≈ **~10.5 min** wall-clock for the
  full sweep.
- Ring extraction: < 1 s/frame. Negligible.
- I/O: 26 × ~2 MB FITS + ~5 KB CSV, under gitignored `output/`.

Resumable, so cost of interruption is bounded.

### 6.9 Primary failure modes

| Failure | Mitigation |
|---|---|
| "Photon ring" ambiguity (critical curve vs resolved ring) | D-C3 locked; explicitly documented |
| Centroid drift with large `\|ε₃\|` | Per-frame intensity-weighted centroid with 3σ clip, not origin |
| Azimuthal binning bias | 180 bins with convergence check gate 6 |
| Multi-valued azimuthal profile (n=1 ring crosses direct image) | Use radially-outer peak per bin; flag bins where multiple peaks exist within ~2 px |
| Sweep crash during 10-min run | Resumable CSV (gate 7) |

### 6.10 Tag at completion

`phase-3c-sweep`

---

## 7. Sub-phase 3D: EHT M87* consistency bound + RNAAS manuscript

### 7.1 New files

```
paper/manuscript.tex
paper/refs.bib
paper/Makefile
paper/figures/asymmetry_vs_epsilon3.pdf          # regenerable
paper/figures/image_gallery.pdf                   # regenerable
paper/scripts/make_figures.py                     # reads output/sweep.csv
src/main/java/com/pranav/grrt/analysis/ConsistencyBound.java
src/test/java/com/pranav/grrt/analysis/ConsistencyBoundTest.java
```

### 7.2 Files modified

- `README.md` — "Reproducing the paper" section with
  `mvn verify -PrunSlow && cd paper && make`.
- `.gitignore` — add `paper/figures/*.pdf`, `paper/*.aux`,
  `paper/*.log`, `paper/*.out`, `paper/*.pdf`. Source `.tex`, `.bib`,
  `Makefile`, and figure-generation scripts tracked (D-D2).
- `CLAUDE.md` — tick EHT comparison and RNAAS paper.

### 7.3 EHT benchmark (D-D1)

Primary benchmark: **Fourier m=1 fractional amplitude** from EHT 2019
Paper VI Table 7. Rationale: this is the quantity Psaltis et al. 2020
uses to derive their bounds, and the RNAAS paper's methodological
comparator is Psaltis 2020. Using the same convention eliminates a
conversion step in the comparison section.

Conversion requirement: our sweep primary metric is RMS
`δ_r/⟨r⟩`. For the bound inversion we convert RMS to Fourier m=1
fractional amplitude using the relation for a perturbed circle
`r(φ) = ⟨r⟩ (1 + Σ_k A_k cos(k φ + φ_k))`:

```
(δ_r/⟨r⟩)²_RMS = ½ Σ_k A_k²
```

If higher harmonics are negligible (gate in 7.4), `A_1 ≈ √2 · δ_r_RMS / ⟨r⟩`.
`ConsistencyBound` computes both and uses the measured Fourier m=1
directly from the CSV (stored in 6.3) rather than the RMS-derived
approximation. Conversion factor documented and unit-tested.

### 7.4 Validation gates

| # | Test | Tolerance |
|---|---|---|
| 1 | `ConsistencyBound` on synthetic monotonic curve — inversion vs analytic | 1e-10 |
| 2 | RMS-to-Fourier-m1 conversion on synthetic pure-m1 profile | 1e-10 |
| 3 | Higher-harmonic content on 17° Kerr frame: `\|A_k\|/A_1 < 0.2` for `k ≥ 2` | non-blocking flag |
| 4 | `make paper` produces a valid PDF | build succeeds |
| 5 | All figures regenerable from `output/sweep.csv` alone | zero hand-artifacts |
| 6 | RNAAS limits: body ≤ 1000 words, figures ≤ 2, refs ≤ 15 | strict |

### 7.5 Wall-clock (M3)

- `ConsistencyBound` computation + figure generation: < 1 min.
- Manuscript writing: human-in-loop, not time-budgeted.

### 7.6 Primary failure modes

| Failure | Mitigation |
|---|---|
| Our `δ_r/⟨r⟩` convention ≠ EHT's | Direct Fourier m=1 comparison via D-D1 conversion; both stored in CSV |
| Systematic disk-model dependence | Acknowledge as caveat in paper; do not over-claim |
| LaTeX toolchain missing | Detect `pdflatex`/`tectonic`/`latexmk` in `paper/Makefile`; flag at build time |
| RNAAS length overrun | Compress iteratively; drop the reference-60° curve to its own inset if needed |

### 7.7 Tag at completion

`phase-3-complete` when the compiled PDF exists and gates pass.
Submission to RNAAS is a separate human action.

---

## 8. Pre-3A gating sequence

In order. 3A coding begins only after all three:

1. **`docs/phase-3-plan.md` committed** (this file).
   Commit message: `docs: record approved Phase 3 plan`.
2. **`1 + h(r, θ) > 0` verified** for `ε₃ ∈ [−2, +2]`, `a = 0.9 M`,
   `r ≥ r₊`, `θ ∈ {0, π/2}` and by extension for all `θ`.
   Verified analytically in section 3 of this document; numerical
   belt-and-suspenders check embedded in the Sympy script.
3. **`scripts/derive_jp_christoffels.py` committed.**

---

## 9. Decisions log (final)

| ID | Decision | Status | Value |
|---|---|---|---|
| D-A1 | JP deformation form | Locked | Johannsen 2013 `h(r,θ) = ε₃ M³ r / Σ²` |
| D-B1 | DP45 dense-output edit | Approved | Add `interpolate()`; `step()` untouched; Phase-2 tests bit-exact |
| D-C1 | Sweep grid | Superseded by correction 2 | 13 denser-near-zero points |
| D-C2 | Resolution | Revised | 512² at both inclinations |
| D-C3 | Ring definition | Approved | Dominant bright ring, 180 bins, per-bin peak |
| D-C4 | Asymmetry metric | Superseded by correction 3 | RMS primary; p2p + m1 secondary |
| D-D1 | EHT benchmark | Locked | Fourier m=1, Paper VI Table 7 |
| D-D2 | `paper/` in git | Locked | Yes for source; PDF/build-artifacts excluded |

---

## 10. Corrections log

| # | Correction | Applied to |
|---|---|---|
| 1 | Christoffels via Sympy, not hand | §4.5 |
| 2 | 13-point denser-near-zero grid | §6.5 |
| 3 | RMS primary asymmetry metric | §6.3, §1 |
| 4 | Kerr/JP consistency < 0.1 pp primary; 2% external secondary | §6.7 gate 3 |

---

*End of plan.*
