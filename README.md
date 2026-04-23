# grrt

General relativistic ray tracer in Java for null geodesic integration in stationary, axisymmetric spacetimes. Built from scratch to render black hole shadow images and, ultimately, quantify how photon ring observables constrain deviations from the Kerr geometry against Event Horizon Telescope data.

## Scientific motivation

The shadow of a black hole (the region of a sky image where null geodesics terminate on the event horizon) is sensitive to the metric outside the horizon. For a Kerr black hole of spin `a`, the shadow is non-circular: frame dragging flattens the prograde edge and produces a characteristic D-shape at high inclination (Bardeen 1973). Parametrized deviations from Kerr (e.g., Johannsen-Psaltis 2011) distort the shadow in predictable, testable ways. Quantitative comparison of rendered shadows against the EHT images of M87* (EHT Collaboration 2019) places model-independent bounds on such deviations, complementing analyses already reported by the EHT consortium.

This repository provides an independent, open, from-first-principles implementation of the geometric ray-tracing component of that comparison. It is deliberately simple enough to audit, fast enough to sweep a parameter space on a workstation, and validated against closed-form analytic benchmarks at every stage.

## Current status

The code is developed in phases, each gated by numerical validation against an analytic reference. Milestones are tagged.

| Phase | Scope | Gate | Status |
|-------|-------|------|--------|
| 1 | Schwarzschild metric, RK4 integrator, pinhole camera, binary shader | Shadow radius matches `3√3 M` to sub-pixel precision | Complete (tag `phase-1-complete`) |
| 2 | Kerr metric (Boyer-Lindquist), Dormand-Prince 5(4) adaptive integrator, polar-axis handling, RayTracer strategy split | Kerr D-shape horizontal diameter at `a = 0.9 M, i = π/2` matches Bardeen (1973) to the discretization floor at 256² and 1024²; E, L drift at machine precision | Complete (tag `phase-2-complete`) |
| 3 | Johannsen-Psaltis parameterized metric, accretion disk emission, EHT M87* comparison | Photon ring observables as a function of JP deviation parameters; consistency bounds from M87* | Planned |

Details of the validation suite for each phase, including the specific tolerances, initial conditions, and measured drift/error values, are recorded in the commit history and in the test sources under `src/test/java/com/pranav/grrt/`.

## Architecture

Four core abstractions, each with a minimal interface:

- `Metric` provides `g_{μν}(x)`, `g^{μν}(x)`, Christoffel symbols, horizon radius, and an optimized `geodesicAcceleration` override where exploitable sparsity warrants it. Implementations: `SchwarzschildMetric`, `KerrMetric`.
- `Integrator` / `AdaptiveIntegrator` provide fixed-step (`RK4`) and adaptive (`DormandPrince45`) solvers for the null geodesic ODE. The adaptive interface returns a `StepStatus` record carrying acceptance, proposed next step, and error norm; tolerances are passed per call, not stored on the instance.
- `RayTracer` is a strategy interface over the per-ray integration loop. `FixedStepRayTracer` preserves Phase 1 semantics; `AdaptiveRayTracer` adds FSAL-aware polar-axis reflection.
- `Camera` and `Renderer` implement a pinhole observer in an orthonormal tetrad at large `r`, with a parallelized pixel loop via `Supplier<RayTracer>` and `ThreadLocal` state.

Further implementation notes (coordinate conventions, sign choices, and the derivation paper trail for each metric) are in `CLAUDE.md`.

## Build and run

Requirements: Java 21, Maven 3.9+.

```bash
mvn test        # run the validation suite
mvn package     # build the jar
```

The test suite is the primary entry point. Shadow renders are produced by the validation tests under `src/test/java/com/pranav/grrt/validation/` (Phase 1 Schwarzschild, Phase 2 Kerr D-shape at two resolutions) and written to `output/` as PNGs. A standalone rendering CLI is not part of Phase 1 or 2 scope.

## Reproducibility

Each tagged phase is a self-contained validation checkpoint. To reproduce a phase's published results:

```bash
git checkout phase-2-complete
mvn test
```

Rendered images and the full numeric output of the validation gates are written to stdout and to `output/`.

## License

MIT. See `LICENSE`.

## Citation

If this code supports work leading to a publication, please cite the repository by its tag. A DOI and a formal citation block will be added when the first associated manuscript is submitted.

```
Sethuraman, P. (2026). grrt: a general relativistic ray tracer in Java.
https://github.com/pranavsethuraman/grrt  (tag: phase-2-complete)
```

## References

- Bardeen, J. M. (1973). *Timelike and null geodesics in the Kerr metric.* In *Black Holes*, ed. C. DeWitt and B. S. DeWitt, Gordon and Breach.
- Chandrasekhar, S. (1983). *The Mathematical Theory of Black Holes.* Oxford University Press.
- Dormand, J. R. and Prince, P. J. (1980). *A family of embedded Runge-Kutta formulae.* J. Comput. Appl. Math. 6, 19.
- Event Horizon Telescope Collaboration (2019). *First M87 Event Horizon Telescope Results.* ApJL 875, L1 through L6.
- Hairer, E., Nørsett, S. P., and Wanner, G. (1993). *Solving Ordinary Differential Equations I: Nonstiff Problems* (2nd ed.), Springer.
- Johannsen, T. and Psaltis, D. (2011). *A metric for rapidly spinning black holes suitable for strong-field tests of the no-hair theorem.* Phys. Rev. D 83, 124015.

## See also

- `CLAUDE.md` for coordinate conventions, derivation notes, and contributor guidance.
