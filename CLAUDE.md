# CLAUDE.md

Context for Claude Code sessions in this repository. Read this before writing code.

---

## Project

**grrt** is a general relativistic ray tracer in Java 21, built with Maven. It
simulates null geodesics in curved spacetime to produce black hole images. The
end goal is a Research Notes of the AAS (RNAAS) paper measuring photon ring
asymmetry vs. a parameterized deviation from Kerr (Johannsen-Psaltis formalism),
compared against Event Horizon Telescope observations of M87*.

**Phases**
1. Schwarzschild ray tracer with fixed-step RK4. (current)
2. Kerr metric, adaptive Dormand-Prince, disk emission.
3. Johannsen-Psaltis metric, EHT comparison, paper.

**Hardware target:** MacBook Pro M3, 16 GB. macOS aarch64.

---

## Commands

Run from the repo root (`~/grrt`):

| Command | Purpose |
|---|---|
| `mvn compile` | Compile main sources |
| `mvn test` | Compile + run all tests |
| `mvn test -Dtest=SchwarzschildMetricTest` | Run one test class |
| `mvn test -Dtest=RK4Test#radialPhotonInMinkowskiIsStraightLine` | Run one method |
| `mvn package` | Produce jar in `target/` |
| `mvn clean` | Wipe `target/` |

Build must be green before any commit. Never commit with failing tests.

---

## Repo layout
grrt/
├── pom.xml                 Maven config. Java 21, junit-jupiter 5.11.3, nom-tam-fits 1.21.2.
├── CLAUDE.md               This file.
├── .cursorrules            Short pointer for Cursor. CLAUDE.md is the source of truth.
├── src/
│   ├── main/java/com/pranav/grrt/
│   │   ├── metric/         Metric interface + implementations
│   │   ├── integrator/     ODE solvers
│   │   ├── camera/         Observer tetrad, pixel → 4-momentum
│   │   ├── renderer/       Pixel loop, FITS orchestration
│   │   └── io/             FITS I/O via nom-tam-fits
│   └── test/java/com/pranav/grrt/
│       └── <mirrors main>  JUnit 5 tests
└── output/                 Rendered FITS images (gitignored)
The architecture is four classes talking through small interfaces. Each
`Metric` is a drop-in replacement; the integrator, camera, and renderer
never branch on metric type. Preserving this is how Kerr and JP get added
as single-file extensions later.

---

## Physics conventions (do not drift from these)

- **Signature:** (-, +, +, +)
- **Units:** geometrized, G = c = 1. Mass M has dimensions of length; default M = 1.
- **Coordinates:** (t, r, θ, φ) with indices 0, 1, 2, 3.
- **State vector (length 8):** `[t, r, θ, φ, k^t, k^r, k^θ, k^φ]`
- **Affine parameter:** λ. Trace null geodesics forward in λ (h > 0) from
  the observer; by time-symmetry of the geodesic equation, the traced
  curve is the reverse of the physical photon's path. Intersections with
  horizon, disk, or celestial sphere identify the source.
- **Null condition:** `g_{μν} k^μ k^ν = 0`. Monitored for numerical drift.

**Schwarzschild landmarks (used in tests):**
- Horizon: r = 2M
- Photon sphere: r = 3M
- Critical impact parameter: b_crit = 3√3 M ≈ 5.196 M

**Kerr landmarks (Phase 2):**
- Outer horizon: r₊ = M + √(M² - a²)
- Ergosphere (equator): r_ergo = 2M
- Prograde ISCO varies with spin

If you are uncertain about a sign, index ordering, or unit, stop and ask.
Silent guessing in GR will introduce bugs that take hours to find.

---

## Coding rules

- Java 21. Use records where they help, not cosmetically.
- State vectors are `double[]`. No `Vector4` class. No boxing.
- No DI frameworks. No Spring. Plain constructors.
- No new dependencies without approval. `nom-tam-fits` and JUnit are it.
- All public APIs have Javadoc stating physical meaning, units, and valid
  domain (e.g. "requires r > 2M").
- Christoffels are computed in closed form per metric, derived by hand, with
  the derivation in Javadoc. Never numerical differentiation.
- `Metric.geodesicAcceleration` has a correctness default but every concrete
  metric overrides it with a sparsity-exploiting version.
- Fail loudly on invalid input (`IllegalArgumentException`). Do not let NaN
  propagate silently.

**Hot-path rules (integrator inner loop):**
- No allocation inside `Integrator.step`. Use pre-allocated scratch arrays.
- No lambdas, streams, or boxing.
- Integrators hold scratch buffers and are **not** thread-safe. One instance
  per worker thread.

**Parallelism:** pixel loop only. Parallel streams or `ForkJoinPool`.

---

## Testing rules

- JUnit 5. Test files mirror source: `src/test/java/.../FooTest.java` for
  `src/main/java/.../Foo.java`.
- Every new `Metric` requires tests for:
  1. Metric × inverse = identity
  2. Christoffel symmetry in lower indices
  3. Asymptotic behavior (flat limit, horizon limit)
  4. At least one closed-form physical prediction (e.g. circular photon orbit)
- Every new `Integrator` requires:
  1. Straight-line test in Minkowski
  2. Closed-orbit preservation in Schwarzschild
  3. Null-condition drift bound
  4. Input validation
- Tolerances:
  - Algebraic identities: `1e-12`
  - Short integration (< 1 orbital period): `1e-6`
  - Long integration: state explicit tolerance in the test name

**When a test fails, debug the cause. Do not relax the tolerance.**

---

## Status (update as phases complete)

- [x] Maven scaffold, Java 21, dependencies
- [x] `Metric` interface
- [x] `SchwarzschildMetric` with closed-form Christoffels + optimized acceleration
- [x] `MinkowskiMetric` for integrator testing
- [x] `Integrator` interface + fixed-step `RK4`
- [x] `Camera` (static observer tetrad in Schwarzschild, pixel → null 4-momentum)
- [ ] `Renderer` (pixel loop, FITS output, parallel)
- [ ] Validation: shadow radius converges to 3√3 M
- [ ] `KerrMetric` (Phase 2)
- [ ] Adaptive `DormandPrince45` integrator (Phase 2)
- [ ] Disk emission model (Phase 2)
- [ ] `JohannsenPsaltisMetric` (Phase 3)
- [ ] EHT comparison pipeline (Phase 3)
- [ ] RNAAS paper

---

## Agent behavior expectations

- **Produce compiling, tested code.** No placeholders, no pseudocode, no
  "TODO: implement this". If a dependency is missing, say so instead of
  inventing it.
- **Never add a dependency silently.** Surface it for approval first.
- **Never relax a test tolerance to make a test pass.** Find the cause.
- **Run `mvn test` before claiming a change is done.** If you cannot run it
  in this environment, say so explicitly.
- **State complexity.** Numerical routines get a time/space complexity note
  in Javadoc.
- **Derive, don't guess.** For new metrics, show the Christoffel derivation
  in the Javadoc before coding the tensor.
- **Prefer fewer files.** This codebase is intentionally small.
- **Ask when a change would touch the `Metric` interface.** It is the
  load-bearing abstraction for Phases 2-3.

---

## Key references

Not stored in repo; cite as needed.
1. Hartle — *Gravity*, Ch. 9 (Schwarzschild geodesics)
2. Chan, Psaltis, Özel (2013) — *GRay* ray tracer methodology
3. Johannsen & Psaltis (2011) — bumpy metric formalism
4. EHT Collaboration (2019, Paper V) — M87* imaging results
5. Bronzwaer et al. (2018) — *RAPTOR*, a modern GRRT for comparison
