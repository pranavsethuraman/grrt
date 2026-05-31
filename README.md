# grrt

**General Relativistic Ray Tracer** — a from-scratch Java 21 implementation for simulating null geodesics in curved spacetime, producing black hole shadow and accretion disk images under parameterised deviations from the Kerr metric.

Built to support a Research Notes of the AAS (RNAAS) paper on photon-ring structure in the Johannsen-Psaltis spacetime at M87*-consistent spin.

## Highlights

- Four spacetime metrics: Minkowski, Schwarzschild, Kerr (Boyer-Lindquist), and Johannsen-Psaltis (JP) with arbitrary &epsilon;&#8323; deformation
- Adaptive Dormand-Prince 4(5) integrator with Shampine 5th-order dense-output interpolation for precise disk-crossing detection
- Novikov-Thorne thin-disk emission model with Bardeen-Press-Teukolsky circular-orbit formulas and bolometric Planck intensity
- Pluggable metric architecture: each spacetime is a drop-in replacement with hand-derived, Sympy-verified Christoffel symbols
- 161 tests validating metric identities, Christoffel symmetries, null-norm conservation, shadow geometry, photon-orbit radii, disk flux, and ring extraction to machine precision
- Reproducible 28-frame &epsilon;&#8323; sweep at 512&times;512 with resumable CSV output

## Key results

**Photon-orbit bifurcation.** At spin a = 0.9 M, the JP prograde equatorial photon orbit merges with the Kerr horizon at &epsilon;&#8323;_crit = 0.1212. Beyond this threshold the prograde shadow edge is set by the event horizon, not a photon sphere.

**ISCO disappearance.** The prograde innermost stable circular orbit ceases to exist for &epsilon;&#8323; &isin; (0.13, 0.20) at a = 0.9 M. Above this range, the standard Novikov-Thorne thin-disk model has no prograde branch.

**Monotonic circularity trend.** The ring-circularity metric &delta;r/&#10216;r&#10217; increases monotonically from 0.121 to 0.198 across &epsilon;&#8323; &isin; [-2.5, -0.2] at inclination i = 17&deg;, confirming that increasingly negative JP deformations produce measurably rounder shadows.

**Extraction precision limit.** Integer-pixel peak finding floors at ~12% fractional radius dispersion, well above the EHT M87* circularity of ~5.5% (Paper VI &sect;7.4). Sub-pixel extraction methods are identified as the required next step for EHT-level constraints on &epsilon;&#8323;.

## Architecture

```
src/main/java/com/pranav/grrt/
  metric/         Metric interface + 4 implementations (Schwarzschild, Kerr, Minkowski, JP)
  integrator/     RK4 + Dormand-Prince 4(5) with adaptive step control and dense output
  camera/         Observer tetrad and pixel-to-ray mapping
  renderer/       Parallel ray tracer with configurable shaders and disk-crossing detection
  disk/           Novikov-Thorne disk model + bolometric emission shader
  analysis/       Ring extractor, circularity metric, epsilon sweep driver, consistency bound
  io/             FITS image writer

paper/            RNAAS manuscript source (AASTeX 6.3.1), figures, and build system
scripts/          Python reference-value generators (Sympy + scipy)
docs/             Phase plans, status snapshots, and parameter-space analysis
```

Each metric supplies its own closed-form Christoffel symbols. The integrator, camera, renderer, and disk never branch on metric type. Swapping Schwarzschild for Kerr for JP is a one-line constructor change.

## Build and test

Requires Java 21+ and Maven 3.9+.

```bash
git clone https://github.com/pranavsethuraman/grrt.git
cd grrt

mvn test                    # 134 fast tests (~25 s)
mvn verify -PrunSlow        # full suite including integration tests (161 tests, ~55 s)
```

## Reproducing the paper

The RNAAS note and its figures are fully reproducible from the committed sweep data.

```bash
# 1. Regenerate output/sweep.csv + FITS frames (slow, ~2 h):
mvn verify -PrunSlow

# 2. Build figures (no LaTeX needed):
cd paper
make figures

# 3. Compile the manuscript (requires a LaTeX engine):
make             # autodetects latexmk / pdflatex / tectonic
make wordcount   # RNAAS body word-count gate (<= 1000 words)
```

`make figures` provisions a local virtual environment and writes both PDFs deterministically (byte-identical across runs via SOURCE_DATE_EPOCH).

## Validation

Every phase ships with quantitative validation gates. Selected highlights:

| Gate | Tolerance | Achieved |
|---|---|---|
| Metric inverse identity g &middot; g&#8315;&#185; = I | 10&#8315;&#185;&#178; | ~10&#8315;&#185;&#179; |
| Kerr reduction of JP at &epsilon;&#8323; = 0 | 10&#8315;&#185;&#178; | 1.42 &times; 10&#8315;&#185;&#8308; |
| Null-norm drift over 1000 M (JP, a = 0.9, &epsilon;&#8323; = 0.5) | 10&#8315;&#8312; | 1.67 &times; 10&#8315;&#185;&#8304; |
| NT flux vs Page-Thorne 1974 reference (12 radii) | 10&#8315;&#8310; | 1.20 &times; 10&#8315;&#185;&#8308; |
| DP45 dense output at midstep vs half-step direct | 10&#8315;&#8312; | 1.00 &times; 10&#8315;&#185;&#185; |
| Face-on disk centroid symmetry at 256&times;256 | 0.5 px | 0.0 px |
| JP photon-orbit radius vs independent Sympy derivation (7 pairs) | 10&#8315;&#185;&#8304; | 2.09 &times; 10&#8315;&#185;&#178; |

Full gate details are recorded in the phase status documents under `docs/`.

## Phase progression

| Phase | Scope | Tag |
|---|---|---|
| 1 | Schwarzschild shadow, RK4 integrator, FITS output | `phase-1-complete` |
| 2 | Kerr metric, adaptive DP45, polar-axis handling | `phase-2-complete` |
| 3A | Johannsen-Psaltis metric, 6 validation gates, cusp discovery | `phase-3a-complete` |
| 3B | Novikov-Thorne disk, dense-output interpolation, disk-crossing detection | `phase-3b-nt-disk` |
| 3C | Ring extraction, circularity analysis, 28-frame production sweep | `phase-3c-sweep` |
| 3D | Consistency bound, figures, RNAAS manuscript | `phase-3-complete` |

## References

- Event Horizon Telescope Collaboration (2019). First M87 Event Horizon Telescope Results. I-VI. *ApJL*, 875.
- Johannsen, T. and Psaltis, D. (2011). Metric for rapidly spinning black holes suitable for strong-field tests of the no-hair theorem. *Phys. Rev. D*, 83, 124015.
- Novikov, I. D. and Thorne, K. S. (1973). Astrophysics of Black Holes. In *Black Holes (Les Astres Occlus)*.
- Page, D. N. and Thorne, K. S. (1974). Disk-Accretion onto a Black Hole. I. *ApJ*, 191, 499.

## Licence

MIT. See `LICENSE`.
