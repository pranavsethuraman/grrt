# Phase 3A Gates 5 & 6 — Sub-plan

Companion to `docs/phase-3-plan.md` §4.7. Closes the remaining two
validation gates for `JohannsenPsaltisMetric` before moving to 3B
(Novikov-Thorne disk emission).

Parent commit: `a9555ea` (Phase 3A gates 1–4). Gates 5 and 6 land as
one follow-up commit, after which sub-phase 3A is tagged
`phase-3a-complete`.

---

## 1. Scope

Two validation gates; one commit. Both exercise physics that gates
1–4 cannot see:

- **Gate 5** catches ε₃-specific sign errors in the Christoffels that
  cancel at ε₃ = 0 (so gate 3 Kerr-reduction misses them) but that
  cause drift in long-integration orbits.
- **Gate 6** cross-checks a closed-form physical prediction
  (equatorial photon orbit radius) at non-zero ε₃ against an
  independent Python derivation.

Both live in the existing `JohannsenPsaltisMetricTest` class plus one
new standalone Python script.

---

## 2. Gate 5 — Null-norm drift on a DP45 orbit

### 2.1 Setup

- Metric: `JohannsenPsaltisMetric(M=1, a=0.9, ε₃=+0.5)`.
- Integrator: `DormandPrince45` with default tolerances (reltol 1e-8,
  abstol 1e-10 per Phase 2 conventions — to confirm from existing
  DP45 test defaults before coding).
- Initial position: `r₀ = 50 M`, `θ = π/2`, `φ = 0`, `t = 0`.
- Initial momentum: equatorial photon with impact parameter
  `b = L/E = 5.5 M`, inward radial component.
- Integration length: `λ ∈ [0, 1000 M]`.

### 2.2 Initial-momentum construction

For a photon with conserved `E = -k_t`, `L = k_φ`, and `b = L/E`:

```
k_t = -E        (lower index)
k_φ = +L = b E  (lower index)
k^t = -g^{tt} E + g^{tφ} L
k^φ = -g^{tφ} E + g^{φφ} L
k^θ = 0         (equatorial)
k^r = -√[(-k^t g_tt k^t − 2 k^t g_tφ k^φ − g_φφ (k^φ)²) / g_rr]
                (negative root = inward)
```

Choose `E = 1` (affine normalization is arbitrary for null geodesics).
Initial null-norm is zero exactly by construction (to machine
precision at the starting point).

**Radicand safety.** The helper computes the radicand
`R = (−k^t g_tt k^t − 2 k^t g_tφ k^φ − g_φφ (k^φ)²) / g_rr` **before**
taking the square root, and asserts `R ≥ 0`. If `R < 0`, the chosen
`(r₀, b)` does not admit a real inward photon direction (typical
cause: `|b|` below the critical-curve window at that radius). The
helper throws `IllegalStateException` with a message naming `r₀`,
`b`, and the numerical value of `R`. Silent NaN propagation into the
integrator is worse than a loud throw — the failure mode would be a
late-trajectory blow-up with no clear attribution.

### 2.3 Success criterion

```
max over accepted DP45 steps of
    |g_μν k^μ k^ν| / (k^t)²   <   1e-8
```

Normalizing by `(k^t)²` absorbs the overall affine-parameter
rescaling. The 1e-8 bound is a ~10× margin over the Phase 2 Kerr
DP45 baseline of ~1e-9 on comparable orbits (see
`DormandPrince45KerrTest`). A real ε₃-odd Christoffel sign bug
produces drift orders of magnitude larger than baseline and is
caught cleanly at 1e-8; the looser 1000× margin (1e-6) from the
Phase 3 plan §4.7 could hide a small-but-real bug and is superseded
by this tighter tolerance.

### 2.4 New code

- New test method in `JohannsenPsaltisMetricTest`:
  `geodesicStaysOnLightConeUnderDp45_jp_a09_eps3_0p5`.
- Helper for initial-momentum construction: local static method in
  the test class (not added to the main source — test-only concern).
- No change to `JohannsenPsaltisMetric` itself.
- Uses existing `DormandPrince45`, `Metric.nullNorm`.

### 2.5 Wall clock

Single orbit at ~5000 DP45 steps for 1000 M of affine length
(estimated from existing `DormandPrince45KerrTest` step counts at
similar accuracy). ~100 ms at JIT steady state. Test budget: < 2 s.

### 2.6 Primary failure modes

- Wrong sign on an ε₃-odd Christoffel term (would cancel at ε₃ = 0,
  invisible to gate 3). Manifests as monotonic drift, typically
  O(1e-3) to O(1) by λ = 1000 M. Gate 5 catches this.
- Transposed index in `fillChristoffelScratch` scratch layout. Would
  show as early-step failure or rapid drift in one coordinate.
- Coupling error between `geodesicAcceleration` and the DP45 `rhs`
  contraction (already tested at Kerr but not at ε₃ ≠ 0). Any such
  error manifests as drift proportional to |ε₃|.

---

## 3. Gate 6 — Closed-form photon orbit cross-check

### 3.1 Approach

Three independent computations of the equatorial photon-orbit radius
`r_photon(a, ε₃)` must agree. Two separate Java code paths compute
the equatorial metric and its r-derivatives; gate 6 validates both
against the Python reference. Each path has a distinct role:

- **Python (Sympy + scipy):** builds the JP metric symbolically,
  forms the equatorial Christoffels via the standard
  `½ g^{ασ} (…)` formula at θ = π/2, and solves the "null +
  circular" system numerically. Authoritative reference.
- **Java path A (ISCO-derived):** uses the equatorial metric
  components and r-derivatives that `eCircular` already implements.
  Validated at the Kerr limit via the ISCO test
  (`zeroEpsilonIscoMatchesKerrBardeenPressTeukolsky`).
- **Java path B (photon-orbit-derived):** independently re-derives
  and implements the same equatorial metric components and
  r-derivatives in `photonOrbitRadius`. Same underlying math,
  separate expressions, kept deliberately duplicate so a regression
  in one path cannot silently propagate to the other.

If Python, path A, and path B all agree at 1e-10 for
`r_photon(pro)` and `r_photon(retro)`, a bug must exist simultaneously
in all three to escape detection.

An optional helper test `equatorialDerivativesAgreeAcrossPaths`
directly asserts that paths A and B compute identical values for
`∂_r g_tt, ∂_r g_tφ, ∂_r g_φφ` at three `(r, a, ε₃)` sample points
within 1e-13 — turning the independence property from a structural
guarantee into a runtime check. See §3.4.

### 3.2 Physics

Equatorial (θ = π/2, cos²θ = 0) with ω ≡ k^φ / k^t:

```
(A) null:     g_tt + 2 g_tφ ω + g_φφ ω² = 0
                  → ω_± = (−g_tφ ± √(g_tφ² − g_tt g_φφ)) / g_φφ

(B) circular: ∂_r g_tt + 2 ∂_r g_tφ ω + ∂_r g_φφ ω² = 0
                  (from Γ^r_{αβ} k^α k^β = 0, factoring −½ g^{rr})
```

`r_photon_prograde` and `r_photon_retrograde` are the roots of
`F(r) ≡ ∂_r g_tt + 2 ∂_r g_tφ ω_±(r) + ∂_r g_φφ ω_±(r)² = 0` for the
two branches of ω.

### 3.3 Seven (a, ε₃) pairs

```
  a     ε₃      notes
-----------------------------------------------------------------
  0.0   0.0    Schwarzschild,     r_photon = 3 M (exact)
  0.9   0.0    Kerr a=0.9,        r_pro = 1.5578, r_retro = 3.9103
  0.9  +0.5    JP positive small  (primary sweep point)
  0.9  −0.5    JP negative small  (primary sweep point)
  0.9  +1.0    JP positive large  (boundary of sweep density)
  0.9  −1.0    JP negative large  (boundary of sweep density)
  0.5  +0.5    intermediate spin  (coverage for intermediate-a code paths)
```

The pair `(0.0, 0.0)` anchors to the Schwarzschild exact value 3 M,
providing an independent sanity check before the JP-specific rows.
The `(0.9, 0)` row cross-checks against Bardeen 1972 closed form
already printed by `BardeenShadowTest` (`r_ph_pro = 1.5578546274`,
`r_ph_retro = 3.9102679391`) to 1e-10. The four JP-ε₃ rows at
`a = 0.9` validate the deformation at the primary sweep spin. The
`(0.5, +0.5)` row closes a coverage gap: every other nonzero-ε₃ row
is at `a = 0.9`, so a bug specific to intermediate-a code paths
would go undetected without this row.

### 3.4 New code

- `scripts/jp_photon_orbit_reference.py`: Sympy + scipy, ~70 lines.
  - Builds the JP metric symbolically (reuses metric setup from
    `derive_jp_christoffels.py` by copy; no module dependency).
  - For each (a, ε₃), constructs `F(r)` at θ = π/2 symbolically,
    lambdifies, and root-finds with `scipy.optimize.brentq` on
    `[1.01 · r₊, 10 M]` (prograde) and `[r₊, 10 M]` (retrograde).
  - Emits a whitespace-delimited table to stdout with header
    `# a eps3 r_photon_pro r_photon_retro`.
  - Committed alongside `derive_jp_christoffels.py`.
- `JohannsenPsaltisMetric.photonOrbitRadius(boolean prograde)`:
  new public method, ~60 lines.
  - Analytic equatorial `g_tt, g_tφ, g_φφ, ∂_r g_{…}` —
    **independently implemented from `eCircular`, not a reuse.**
    The two Java code paths (ISCO via `eCircular`, photon-orbit via
    `photonOrbitRadius`) must compute the equatorial metric and its
    r-derivatives in separate expressions. Same underlying math,
    deliberately duplicate code, so a regression in one path is
    caught by `equatorialDerivativesAgreeAcrossPaths` rather than
    silently propagating. Javadoc cross-references `eCircular` and
    this sub-plan.
  - Bracket walk identical in structure to `iscoRadius` (start at
    `1.1 r₊`, ×1.05 until null discriminant is real), then bisect
    `F(r) = 0` on `[rLo, 10 M]`.
- Two test methods in `JohannsenPsaltisMetricTest`:
  - `photonOrbitRadiusMatchesPythonReference`. Hardcodes the
    seven-row reference table as `double[][]` (the Python script's
    output, spot-checked by a human before the commit; the script
    itself is NOT run at `mvn test` time). For each row, asserts
    Java `photonOrbitRadius(pro)` and `photonOrbitRadius(retro)`
    agree with the tabulated value to 1e-10.
  - `equatorialDerivativesAgreeAcrossPaths`. Calls two package-
    private helpers on `JohannsenPsaltisMetric` that expose the
    equatorial derivative computations of path A (ISCO-derived) and
    path B (photon-orbit-derived), and asserts element-wise equality
    of `{g_tt, g_tφ, g_φφ, ∂_r g_tt, ∂_r g_tφ, ∂_r g_φφ}` to 1e-13
    at three `(r, a, ε₃)` sample points.

To support the helper test, two new package-private methods on
`JohannsenPsaltisMetric`:

```
double[] equatorialGDerivativesViaIscoPath(double r);
double[] equatorialGDerivativesViaPhotonPath(double r);
```

Each returns `{g_tt, g_tφ, g_φφ, ∂_r g_tt, ∂_r g_tφ, ∂_r g_φφ}`
computed by its respective path's expressions. Production code does
not depend on either helper; they exist solely as test-visibility
hooks. `eCircular` and `photonOrbitRadius` call their own internal
inline expressions (which the helpers mirror line-for-line), not
these package-private methods, so the runtime hot paths remain
unchanged.

### 3.5 Why hardcode the reference table, not load it at test time?

- Avoids Python as a test-time dependency of the Maven build.
- Keeps the test hermetic (no external file I/O).
- The Python script is committed and rerunnable if the metric form
  ever changes; its output becomes a new reference table that the
  developer pastes into the test. Regeneration policy mirrors the
  CSE-block pattern in `fillChristoffelScratch`.

### 3.6 Success criterion

```
for each of 6 rows:
    |photonOrbitRadius_JP(pro)   − ref_pro|   < 1e-10
    |photonOrbitRadius_JP(retro) − ref_retro| < 1e-10
```

### 3.7 Wall clock

- Python script: ~12 s total (Sympy symbolic setup + 14 root-finds
  — 7 rows × 2 branches). Run once, output pasted by hand.
- Java test: 14 bracket walks + bisections at < 100 iterations each.
  < 60 ms; test budget < 1 s for both gate-6 tests combined
  (`photonOrbitRadiusMatchesPythonReference` +
  `equatorialDerivativesAgreeAcrossPaths`).

### 3.8 Primary failure modes

- Python/Java disagreement beyond 1e-10. Root cause likely in the
  analytic equatorial r-derivatives in `eCircular` (the only
  JP-specific code not exercised by gate 5). Debug by comparing
  `dgtt, dgtp, dgpp` between the two implementations at a shared
  point.
- Wrong ω-branch selection (`prograde ↔ +`, `retrograde ↔ −` at
  a > 0). Catches via the Kerr `(0.9, 0)` row: must match
  `r_photon_pro = 1.5578` and `r_photon_retro = 3.9103` — mixing
  branches would swap these.
- Bracket failure on large `|ε₃|` if `r_photon` leaves
  `[1.01 r₊, 10 M]`. Mitigate by widening `rHi` to `20 M` if needed
  (ISCO bracketing uses `20 M`; the photon orbit is always inside
  the ISCO).

---

## 4. Commit structure

Single commit, proposed message:

```
metric: close Phase 3A with gate 5 (null-norm drift) and gate 6 (photon orbit)

Gate 5: null-norm drift of a DP45 photon orbit in JP(a=0.9, eps3=+0.5)
stays below 1e-6 over 1000 M of affine length. Catches eps3-odd
Christoffel sign errors that Kerr reduction at eps3=0 misses.

Gate 6: equatorial photon orbit radius cross-checked at six (a, eps3)
pairs against an independent Sympy+scipy derivation in
scripts/jp_photon_orbit_reference.py. Agreement to 1e-10.

New: JohannsenPsaltisMetric.photonOrbitRadius(boolean).
New: scripts/jp_photon_orbit_reference.py (regenerable reference table).
New: two test methods in JohannsenPsaltisMetricTest.

Refs: docs/phase-3a-gates-5-6.md, docs/phase-3-plan.md §4.7.
```

Files modified / added:

- `scripts/jp_photon_orbit_reference.py` (new, ~60 lines).
- `src/main/java/com/pranav/grrt/metric/JohannsenPsaltisMetric.java`
  (new `photonOrbitRadius` method, ~50 lines).
- `src/test/java/com/pranav/grrt/metric/JohannsenPsaltisMetricTest.java`
  (two new test methods, ~80 lines total).
- `CLAUDE.md` (tick the `JohannsenPsaltisMetric` checkbox at the end
  of this sub-phase).

After commit + `mvn test` green, tag `phase-3a-complete`.

---

## 5. Pre-coding checklist

Before writing any Java for gates 5/6:

1. Confirm DP45 default tolerances from the existing
   `DormandPrince45Test` / `DormandPrince45KerrTest` so gate 5 uses
   the same settings and the 1e-6 drift budget is apples-to-apples
   with Phase 2.
2. Verify the Kerr photon-orbit closed-form values cited in §3.3
   match the `BardeenShadowTest` printed values at a = 0.9
   (1.5578546274 and 3.9102679391 — already in its console log).
3. Run `scripts/jp_photon_orbit_reference.py` once, eyeball the
   table for plausibility (prograde monotonic in ε₃, retrograde
   ditto, `(0.0, 0.0)` → exactly 3 M to double precision).
4. Only then write `photonOrbitRadius` and the tests.

---

## 6. Deferred

Nothing. Gates 5 and 6 complete sub-phase 3A as specified by the
parent plan. After this, sub-phase 3B (Novikov-Thorne disk emission)
begins per `docs/phase-3-plan.md` §5.

*End of sub-plan.*
