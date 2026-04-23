# General Relativistic Ray Tracing for Black Hole Imagery: Photon Ring Asymmetry as a Probe of Non-Kerr Spacetime

## Abstract

This project implements a general relativistic (GR) backward ray tracer from scratch in Java to simulate the appearance of accreting supermassive black holes. Null geodesics are integrated in the Kerr spacetime and a parameterized non-Kerr "bumpy" generalization using the classical 4th-order Runge-Kutta (RK4) method on the first-order form of the geodesic equations. For each pixel on the observer screen, a photon is traced backward in affine parameter until it either crosses the event horizon, intersects the accretion disk, or escapes to infinity. The resulting image is emitted in FITS format for quantitative comparison against Event Horizon Telescope (EHT) observations of M87*.

The scientific contribution centers on measuring the **photon ring asymmetry** as a function of the Johannsen-Psaltis (JP) deviation parameter `epsilon_3`. The hypothesis under test is that the observed asymmetry of the M87* photon ring places a measurable bound on `epsilon_3`, providing a null test of the Kerr hypothesis. A short-form submission to Research Notes of the AAS (RNAAS) is the intended publication target.

Project goals:

- Full Schwarzschild and Kerr null-geodesic integration, validated against closed-form results (photon sphere at `r = 3M`, ISCO at `r = 6M` for Schwarzschild; prograde/retrograde ISCO formulas for Kerr)
- Parameterized non-Kerr metric via the Johannsen-Psaltis formalism
- Geometrically thin, optically thick Novikov-Thorne accretion disk model for initial imaging
- FITS-format image output compatible with DS9 and standard EHT analysis tooling
- Photon ring asymmetry extraction, with comparison against published EHT M87* constraints

> **Status:** in development. Tables marked **[TBD]** are placeholders that will be populated as the corresponding experiments complete.

---

## Step-by-Step Reproduction Guide

This section is written so that any reader (including RNAAS reviewers) can replicate every result in this project from a clean machine. Follow the steps in order.

### Prerequisites

| Requirement    | Version           | Notes                                          |
|----------------|-------------------|------------------------------------------------|
| Java (JDK)     | 21 LTS or newer   | Records, pattern matching, Vector API          |
| Maven          | 3.9 or newer      | Build system                                   |
| Python         | 3.10 or newer     | Plotting, photon ring analysis                 |
| DS9            | 8.5 or newer      | FITS image inspection (optional)               |
| IntelliJ IDEA  | 2024.1 or newer   | Recommended IDE (optional)                     |
| GMP / JBLAS    | not required      | Native linear algebra is not used              |

**Operating system:** Linux, macOS, or Windows 10/11. No platform-specific dependencies.
**Hardware:** A multi-core CPU is required for reasonable render times. A modern discrete GPU is optional and can accelerate batch image production by one to two orders of magnitude through the OpenCL / CUDA backend (Step 5b).

### Step 0: Clone and Build

```bash
git clone https://github.com/<user>/gr-raytracer.git
cd gr-raytracer
mvn clean package
```

Produces `target/raytracer-<version>.jar` plus a `bin/raytracer` launch script.
**Runtime:** < 60 s first build (dependency resolution), < 15 s thereafter.

### Step 1: Metric Validation (Pure Geometry)

Verifies that each implemented `Metric` implementation satisfies analytic consistency checks before any integration is attempted. No images are produced in this step.

```bash
mvn test -Dtest=MetricValidationTest
```

**Checks:**

- Schwarzschild metric components `g_tt`, `g_rr`, `g_theta_theta`, `g_phi_phi` match closed form to machine precision
- Kerr metric reduces to Schwarzschild as `a -> 0`
- Johannsen-Psaltis metric reduces to Kerr as `epsilon_3 -> 0`
- Christoffel symbols satisfy the torsion-free identity `Gamma^a_bc = Gamma^a_cb`
- All metrics are asymptotically flat: `g_{mu nu} -> eta_{mu nu}` as `r -> infinity`

**Runtime:** < 5 s.

### Step 2: Integrator Validation (Null Geodesic Conservation)

Traces photons in Schwarzschild and Kerr spacetimes and verifies that the conserved quantities drift by less than a prescribed tolerance over a long integration arc.

```bash
bin/raytracer validate-integrator --metric schwarzschild --steps 200000
bin/raytracer validate-integrator --metric kerr --spin 0.9 --steps 200000
```

**Checks:**

- Null condition: `g_{mu nu} k^mu k^nu` stays within `1e-10` of zero
- Energy `E = -k_t` conserved to `1e-9` relative
- Axial angular momentum `L_z = k_phi` conserved to `1e-9` relative
- Carter constant `Q` (Kerr only) conserved to `1e-8` relative

**Outputs:** `data/validation/conservation_<metric>.csv`
**Runtime:** ~30 s per metric on a single core.

### Step 3: First Schwarzschild Image

Produces the canonical black hole shadow image for a non-rotating black hole viewed from the equatorial plane.

```bash
bin/raytracer render \
  --metric schwarzschild \
  --inclination 85 \
  --resolution 1024x1024 \
  --output images/schwarzschild_i85.fits
```

**Outputs:**

- `images/schwarzschild_i85.fits` (FITS image, intensity in cgs)
- `images/schwarzschild_i85_meta.json` (ray statistics, metric parameters)

**Validation:** shadow radius must match `r_shadow = 3 sqrt(3) M` to sub-pixel precision.
**Runtime:** ~2 min on 16 cores at 1024x1024.

### Step 4: Kerr Shadow Image Sweep

Generates shadow images over a grid of spin values and observer inclinations, reproducing the classical Bardeen (1973) shadow shapes for validation against published contours.

```bash
bin/raytracer batch \
  --metric kerr \
  --spin 0.0,0.5,0.9,0.99 \
  --inclination 17,45,75,85 \
  --resolution 2048x2048 \
  --output images/kerr/
```

**Outputs:** `images/kerr/a<spin>_i<incl>.fits` for each (spin, inclination) pair.
**Validation:** extracted shadow outlines overlaid on Bardeen analytic contours in `figures/kerr_shadow_validation.png`.
**Runtime:** ~40 min total on 16 cores for the 16-image grid.

### Step 5: Johannsen-Psaltis Non-Kerr Sweep

The core scientific experiment. Sweeps the JP deviation parameter `epsilon_3` over a physically motivated range at M87*-like inclination, holding spin fixed at the EHT-inferred value.

```bash
bin/raytracer batch \
  --metric jp \
  --spin 0.9 \
  --epsilon3 -5,-2,-1,0,1,2,5 \
  --inclination 17 \
  --resolution 4096x4096 \
  --output images/jp/
```

**Outputs:** `images/jp/eps3_<value>_i17.fits` for each `epsilon_3` value.
**Runtime:** ~2 h total on 16 cores at 4096x4096. GPU backend reduces this to ~10 min (see Step 5b).

### Step 5b: GPU Acceleration (Optional)

For batch work, the OpenCL integrator kernel ports the RK4 step to the GPU.

```bash
bin/raytracer batch \
  --backend opencl \
  --device 0 \
  [...same args as Step 5...]
```

**Runtime:** ~12 min for the full `epsilon_3` sweep on a modern discrete GPU.

### Step 6: EHT M87* Comparison (Photon Ring Asymmetry)

Extracts the photon ring asymmetry `A = (r_max - r_min) / (r_max + r_min)` from each JP image and computes a chi-squared against the EHT 2019 and 2024 M87* constraints.

```bash
python analysis/photon_ring_asymmetry.py \
  --input images/jp/ \
  --eht-constraints data/eht/m87_2019_ring_constraints.json \
  --output figures/asymmetry_vs_eps3.png
```

**Outputs:**

- `figures/asymmetry_vs_eps3.png` (main paper figure)
- `data/results/jp_asymmetry.csv` (per-image asymmetry with uncertainty)
- `data/results/chi2_vs_eps3.json` (chi-squared curve, 1-sigma and 2-sigma bounds on `epsilon_3`)

**Runtime:** < 30 s.

### Quick Reference: Complete Reproduction in One Block

```bash
git clone https://github.com/<user>/gr-raytracer.git
cd gr-raytracer
mvn clean package

# Validation
mvn test -Dtest=MetricValidationTest
bin/raytracer validate-integrator --metric schwarzschild --steps 200000
bin/raytracer validate-integrator --metric kerr --spin 0.9 --steps 200000

# Imaging
bin/raytracer render  --metric schwarzschild --inclination 85  --resolution 1024x1024 --output images/schwarzschild_i85.fits
bin/raytracer batch   --metric kerr --spin 0.0,0.5,0.9,0.99 --inclination 17,45,75,85 --resolution 2048x2048 --output images/kerr/
bin/raytracer batch   --metric jp   --spin 0.9 --epsilon3 -5,-2,-1,0,1,2,5 --inclination 17 --resolution 4096x4096 --output images/jp/

# Analysis
python analysis/photon_ring_asymmetry.py --input images/jp/ --eht-constraints data/eht/m87_2019_ring_constraints.json --output figures/asymmetry_vs_eps3.png
```

### Expected Runtimes

| Step | Command                                   | Approximate Time (16 cores) | Approximate Time (GPU) |
|------|-------------------------------------------|-----------------------------|------------------------|
| 0    | `mvn clean package`                       | < 60 s (cold), < 15 s (warm)| N/A                    |
| 1    | `mvn test -Dtest=MetricValidationTest`    | < 5 s                       | N/A                    |
| 2    | Conservation tests (2 metrics)            | ~1 min                      | N/A                    |
| 3    | Schwarzschild 1024^2 image                | ~2 min                      | ~15 s                  |
| 4    | Kerr shadow sweep (16 images, 2048^2)     | ~40 min                     | ~4 min                 |
| 5    | JP sweep (7 images, 4096^2)               | ~2 h                        | ~12 min                |
| 6    | Photon ring asymmetry analysis            | < 30 s                      | N/A                    |

Total wall time for full reproduction (CPU only): approximately 3 hours. GPU: approximately 25 minutes.

---

## Key Results

### Metric and Integrator Validation

| Test                                    | Tolerance | Measured  |
|-----------------------------------------|-----------|-----------|
| Null condition drift                    | 1e-10     | **[TBD]** |
| Energy conservation (Schwarzschild)     | 1e-9      | **[TBD]** |
| Energy conservation (Kerr, a=0.9M)      | 1e-9      | **[TBD]** |
| Angular momentum conservation           | 1e-9      | **[TBD]** |
| Carter constant conservation (Kerr)     | 1e-8      | **[TBD]** |
| Schwarzschild shadow radius vs 3*sqrt(3)*M | 1 pixel | **[TBD]** |

### Kerr Shadow Reproduction

| Spin `a/M` | Inclination | Shadow deviation from Bardeen contour |
|------------|-------------|----------------------------------------|
| 0.0        | 85 deg      | **[TBD]**                              |
| 0.5        | 85 deg      | **[TBD]**                              |
| 0.9        | 85 deg      | **[TBD]**                              |
| 0.99       | 85 deg      | **[TBD]**                              |

### Photon Ring Asymmetry vs JP Deviation Parameter (Main Result)

| `epsilon_3` | Ring asymmetry `A` | Delta chi^2 vs EHT M87* |
|-------------|--------------------|--------------------------|
| -5.0        | **[TBD]**          | **[TBD]**                |
| -2.0        | **[TBD]**          | **[TBD]**                |
| -1.0        | **[TBD]**          | **[TBD]**                |
|  0.0 (Kerr) | **[TBD]**          | **[TBD]**                |
| +1.0        | **[TBD]**          | **[TBD]**                |
| +2.0        | **[TBD]**          | **[TBD]**                |
| +5.0        | **[TBD]**          | **[TBD]**                |

**1-sigma bound on `epsilon_3`:** **[TBD]**
**2-sigma bound on `epsilon_3`:** **[TBD]**

---

## Methodology

### The Null Geodesic Equation

Photons follow null geodesics in curved spacetime. In affine-parameter form:

```
d^2 x^a / d lambda^2 + Gamma^a_{bc} (dx^b / d lambda) (dx^c / d lambda) = 0
```

subject to the null condition `g_{ab} (dx^a / d lambda)(dx^b / d lambda) = 0`. This is recast as a 1st-order system in 8 variables `(x^a, k^a) = (t, r, theta, phi, k^t, k^r, k^theta, k^phi)` and integrated with classical RK4.

### Backward Ray Tracing

Rather than tracing photons from emitter to camera, each camera pixel emits a ray traveling *backward* in affine parameter. This is efficient because only rays that reach the camera contribute to the image. Termination conditions:

1. `r < r_horizon + epsilon` -> photon captured by hole (pixel set to 0)
2. Ray crosses the accretion disk plane within `r_isco <= r <= r_out` -> intensity sampled from the disk emission model
3. `r > r_far` with outgoing radial velocity -> photon escapes to infinity (CMB background or 0)
4. Step count exceeds `MAX_STEPS` -> marked as stuck (rare, logged for diagnostics)

### Adaptive Step Size

Fixed-step RK4 wastes work near the horizon and under-resolves tight orbits. A simple error controller based on a step-doubling estimate (Richardson extrapolation between one full step and two half steps) adjusts the affine-parameter step by:

```
dl_new = safety * dl_old * (tol / err)^(1/5)
```

with `safety = 0.9`, `tol = 1e-9`, and hard bounds `[1e-4, 1e-1] M`.

### Kerr Metric (Boyer-Lindquist)

```
ds^2 = -(1 - 2Mr/Sigma) dt^2
       - (4 M a r sin^2(theta) / Sigma) dt dphi
       + (Sigma / Delta) dr^2
       + Sigma d theta^2
       + (r^2 + a^2 + 2 M a^2 r sin^2(theta) / Sigma) sin^2(theta) dphi^2
```

with `Sigma = r^2 + a^2 cos^2(theta)` and `Delta = r^2 - 2 M r + a^2`. Christoffel symbols are computed analytically, not numerically, to avoid finite-difference error near the horizon where `Delta -> 0`.

### Johannsen-Psaltis "Bumpy" Metric

The JP metric introduces four deviation functions `h_t`, `h_r`, `h_theta`, `h_phi` parameterized by `epsilon_n`. For the single-parameter subfamily used here (`epsilon_3` only, all others set to zero), the metric reduces to Kerr when `epsilon_3 = 0` and retains the Carter constant structure, allowing an identical integration pipeline. The full metric expression is implemented in `src/main/java/metric/JohannsenPsaltisMetric.java` and verified against the original 2011 paper.

### Photon Ring Asymmetry Metric

Following the EHT 2019 Paper VI convention:

```
A = (d_max - d_min) / (d_max + d_min)
```

where `d_max` and `d_min` are the maximum and minimum angular diameters of the photon ring, extracted by:

1. Locating the brightness peak along radial cuts at 1-degree azimuthal spacing
2. Fitting a smoothed ellipse to the resulting (x, y) peak locus
3. Reading off the major and minor axes

Uncertainties are estimated by bootstrap resampling over azimuthal bins (1000 resamples).

### Why Implement from Scratch (vs `gyoto`, `GRay2`, `bhlight`)

Existing codes are excellent but two project goals motivate a clean-room implementation:

1. **Educational:** full understanding of every numerical choice, which is the stated point of the RNAAS submission
2. **Parameterized metric freedom:** easier to plug in bespoke non-Kerr metrics without fighting a large foreign codebase

Validation against `gyoto` output for the Kerr case is a planned sanity check (see Future Work).

---

## Architecture

```
  +-------------+        +---------------+         +----------------+
  |   Camera    |  rays  |   Integrator  |  query  |     Metric     |
  |             +------->+               +-------->+                |
  | - FOV       |        | - RK4 step    |         | - Schwarzschild|
  | - Resolution|        | - Adaptive dl |         | - Kerr         |
  | - Position  |        | - Termination |         | - JP non-Kerr  |
  | - Tetrad    |        |   checks      |         |                |
  +------+------+        +-------+-------+         +-------+--------+
         |                       |                         |
         | initial (x^a, k^a)    | g_{ab}, Gamma^a_{bc}    |
         |                       |                         |
         v                       v                         |
  +------+-----------------------+-----+                   |
  |               Renderer             |<------------------+
  |                                    |
  | - Accumulate flux per pixel        |
  | - Query disk emission model        |
  | - Apply redshift factor g^4        |
  | - Write FITS image                 |
  +------------------+-----------------+
                     |
                     v
               +-----+-----+
               |   FITS    |
               |   image   |
               +-----+-----+
                     |
                     v
            +--------+---------+           +----------------+
            | Python analysis  |<----------+ EHT M87* data  |
            | - Ring extraction|           | (constraints)  |
            | - chi^2 vs EHT   |           +----------------+
            +--------+---------+
                     |
                     v
            +--------+---------+
            | figures/         |
            | paper tables     |
            +------------------+
```

**Core classes:**

- `Metric` (interface): `double[] metric(double[] x)`, `double[][] christoffel(double[] x)`
- `Integrator`: owns the RK4 step, error controller, and termination logic; depends only on `Metric`
- `Camera`: constructs an orthonormal tetrad at the observer, emits initial conditions for each pixel
- `Renderer`: drives the pixel loop, accumulates intensity, writes FITS via `nom.tam.fits`

**Concurrency:** the outer pixel loop is embarrassingly parallel. CPU backend uses a `ForkJoinPool` with per-thread integrator instances to avoid contention. GPU backend ports the inner RK4 step to OpenCL with one work item per pixel.

**Precision:** all geometry is `double`. No 32-bit float is used anywhere in the integrator. Vector API (`jdk.incubator.vector`) is used selectively for Christoffel contractions.

---

## Repository Structure

```
gr-raytracer/
|
|-- README.md                     # This file
|-- pom.xml                       # Maven build
|-- LICENSE                       # MIT
|
|-- src/main/java/
|   |-- metric/
|   |   |-- Metric.java           # Interface: metric tensor + Christoffel symbols
|   |   |-- SchwarzschildMetric.java
|   |   |-- KerrMetric.java
|   |   `-- JohannsenPsaltisMetric.java
|   |
|   |-- integrator/
|   |   |-- Integrator.java       # RK4 step + adaptive error controller
|   |   |-- Geodesic.java         # Typed record for (x^a, k^a)
|   |   `-- Terminator.java       # Horizon / disk / infinity checks
|   |
|   |-- camera/
|   |   |-- Camera.java           # Observer position, tetrad, pixel -> ray
|   |   `-- Tetrad.java           # Orthonormal frame construction
|   |
|   |-- disk/
|   |   |-- AccretionDisk.java    # Geometrically thin, optically thick
|   |   `-- NovikovThorne.java    # Emission profile
|   |
|   |-- renderer/
|   |   |-- Renderer.java         # Pixel loop, FITS output
|   |   `-- RedshiftShader.java   # g^4 intensity factor
|   |
|   |-- backend/
|   |   |-- CpuBackend.java       # ForkJoinPool
|   |   `-- OpenClBackend.java    # GPU kernel dispatch
|   |
|   `-- cli/
|       `-- Main.java             # Argument parsing, subcommands
|
|-- src/test/java/
|   |-- MetricValidationTest.java
|   |-- ConservationTest.java
|   `-- ShadowContourTest.java
|
|-- analysis/
|   |-- photon_ring_asymmetry.py  # Ring extraction + chi^2 vs EHT
|   |-- plot_shadows.py           # Kerr shadow validation figure
|   `-- fits_io.py                # Thin astropy wrapper
|
|-- data/
|   |-- eht/
|   |   `-- m87_2019_ring_constraints.json
|   |-- validation/               # Conservation-test CSVs
|   `-- results/                  # Final asymmetry and chi^2 tables
|
|-- images/
|   |-- schwarzschild_i85.fits
|   |-- kerr/                     # (spin, inclination) grid
|   `-- jp/                       # epsilon_3 sweep
|
|-- figures/                      # Publication PNGs
|
`-- paper/
    |-- rnaas.tex                 # RNAAS short-form manuscript
    `-- refs.bib
```

---

## Output Figures

| File                                             | Contents                                                         |
|--------------------------------------------------|------------------------------------------------------------------|
| `figures/schwarzschild_i85.png`                  | Canonical Schwarzschild shadow at 85 deg inclination             |
| `figures/kerr_shadow_validation.png`             | Simulated Kerr shadows overlaid on Bardeen analytic contours     |
| `figures/kerr_shadow_grid.png`                   | (spin, inclination) image grid                                   |
| `figures/jp_image_grid.png`                      | Image grid across `epsilon_3` values                             |
| `figures/asymmetry_vs_eps3.png`                  | **Main result:** ring asymmetry vs `epsilon_3` with EHT band     |
| `figures/chi2_vs_eps3.png`                       | chi-squared curve with 1- and 2-sigma bounds                     |
| `figures/conservation_schwarzschild.png`         | Null / E / L drift vs affine parameter (Schwarzschild)           |
| `figures/conservation_kerr.png`                  | Null / E / L / Q drift (Kerr, a=0.9M)                            |

---

## Troubleshooting

**"Integrator produced NaN near horizon"**
The adaptive step size has a hard floor of `dl_min = 1e-4 M`. Near extremal Kerr (`a -> 1`), `Delta` can drop below working precision. Either (a) cap `a` at 0.998, as is standard in the astrophysics literature, or (b) switch to Kerr-Schild coordinates via `--coords kerr-schild`, which are horizon-penetrating.

**"FITS image all zeros"**
Most commonly a tetrad orientation bug: the camera was pointing away from the hole. Confirm the observer radius is large (default 1000 M) and that `--inclination` is in degrees, not radians.

**"Kerr shadow does not match Bardeen contour"**
Double-check that the spin sign convention matches Bardeen's. This implementation uses `a > 0` for prograde-spinning black holes and places the spin axis along `+z`. Inclination is measured from the spin axis.

**OpenCL kernel fails to compile**
The GPU backend requires OpenCL 1.2 or newer with `cl_khr_fp64`. On consumer NVIDIA cards, enable double precision explicitly. Falls back to CPU cleanly if `--backend opencl` is requested but not available.

**Maven build fails with "class file version 65.0"**
You are on an older JDK. This project requires Java 21 (class version 65). Run `java -version` and upgrade.

---

## Assumptions and Limitations

1. **Geometrically thin disk.** The Novikov-Thorne model ignores disk thickness and self-obscuration. This is appropriate for low-accretion-rate sources like M87* but would fail for, e.g., Sgr A* ADAFs.
2. **Single deviation parameter.** Only `epsilon_3` is varied; all other JP deviation functions are held at zero. A fuller analysis would marginalize over the full deviation vector.
3. **Fixed spin and inclination.** Spin is taken from the EHT-preferred value and inclination is fixed at the jet-inferred 17 deg. The constraint on `epsilon_3` is conditional on these priors.
4. **Optically thick disk.** Photons are terminated on first disk crossing. No radiative transfer is performed; intensity comes from the local disk emission model times the redshift factor `g^4`.
5. **No polarimetry.** EHT polarization results are not used.

---

## Future Work

1. **Cross-validation against `gyoto`** for the Kerr case, to rule out systematic bugs before trusting JP results
2. **Full JP deviation sweep** over the four-parameter family, with MCMC to produce posterior bounds
3. **Novikov-Thorne to RIAF** disk model for Sgr A* application
4. **Polarimetric ray tracing** with parallel transport of the polarization four-vector
5. **GPU kernel autotuning** across register budgets and local memory layouts

---

## Citation

Intended publication venue: Research Notes of the AAS.

```bibtex
@article{<lastname>2026raytracer,
  author  = {<Author>},
  title   = {Photon Ring Asymmetry as a Probe of Non-Kerr Spacetime with a
             From-Scratch General Relativistic Ray Tracer},
  journal = {Research Notes of the AAS},
  year    = {2026},
  note    = {In preparation}
}
```

**License:** MIT.
