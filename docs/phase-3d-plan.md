# Phase 3D Sub-plan — EHT consistency bound + RNAAS manuscript

Drafted 2026-05-26 as the successor to `docs/phase-3c-status.md`.
Parent commit: `9def498` (`phase-3c-sweep`). Sub-phase 3D lands as
three commits leading to tag `phase-3-complete`. Supersedes
`docs/phase-3-plan.md` §7 where the two overlap.

This plan covers only what was deferred from 3C: the algorithmic
inversion (`ConsistencyBound`), the figure-generation script
(`make_figures.py`), and the AASTeX manuscript. The scientific
inputs (`output/sweep.csv`, 28 rendered frames) already exist at
the parent commit; 3D does not re-render.

---

## 1. Scope

Sub-phase 3D delivers three additions and zero
`Metric`/`Integrator`/`Renderer` changes:

1. New `com.pranav.grrt.analysis.ConsistencyBound` (pure-math
   inversion of `δ_r/⟨r⟩(ε₃)` against an EHT amplitude band) and
   its unit test.
2. New `paper/` directory: `Makefile`, `scripts/make_figures.py`,
   `refs.bib`, `manuscript.tex`. Reads only `output/sweep.csv`
   (and the 28 FITS frames for the gallery); no hand-edited
   plot data.
3. Repo housekeeping: `.gitignore` (PDF build artifacts), `README.md`
   ("Reproducing the paper" section), `CLAUDE.md` (tick the final
   checkbox).

No edits to any class under `src/main/java/com/pranav/grrt/{metric,
integrator,camera,renderer,io}`. The `Metric` interface is left
exactly as it stands at `9def498`.

---

## 2. Pre-approved decisions (lock; do not revisit)

These five decisions are the user-approved deltas relative to
`docs/phase-3-plan.md` §7 and `docs/phase-3c-status.md` §9.4.

### 2.1 Bound interpolation domain: `[−2.5, −0.2]` only

`ConsistencyBound` inverts `δ_r/⟨r⟩ → ε₃` on the strictly-monotone
sub-interval `[−2.5, −0.2]` at `i = 17°` (five sweep points). The
small-|ε₃| plateau `[−0.1, +0.08]` (six points clustered at
δ_r/⟨r⟩ ≈ 0.195 ± 0.001) is collapsed into a single "near-zero"
data cluster and used only to flag the upper-bound case "EHT
amplitude lies within the plateau band → upper bound is
unconstrained from above by this sweep" (§5.2 below).

Rationale: at i = 17° the extractor noise floor on the small-|ε₃|
side is wider than the ε₃-induced variation. Inverting through the
plateau would produce numerically valid but physically meaningless
ε₃ estimates. The caveat is documented in the manuscript.

### 2.2 Cusp framing: two-sided bound + ISCO disappearance as separate prediction

The original plan §7.3 envisioned a third result — a qualitative
shadow-structure change at ε₃_crit detectable in the azimuthal-
average asymmetry. The 3C cusp-neighborhood data (status doc §4)
shows no such jump in the 180-bin azimuthal-average δ_r/⟨r⟩
metric. The paper therefore reframes:

| Result | Source |
|---|---|
| Two-sided bound on ε₃ at `a = 0.9, i = 17°` | Smooth segment, EHT m=1 inversion |
| ISCO disappears at `ε₃_ISCO ∈ (+0.13, +0.20)` | Future-observable, separate prediction; cited from 3C sweep skips (status doc §6) |

**No sharp-cusp claim** is made in the manuscript. The ε₃_crit
landmark appears in Figure 1 as a vertical dashed line for
context, with text noting that detecting the cusp directly
requires a cusp-sensitive diagnostic beyond the scope of this
note.

### 2.3 EHT Paper VI Table 7 m=1 amplitude hardcoded

`ConsistencyBound` reads no external EHT file. The Paper VI Table 7
fractional m=1 amplitude (mid value + 1σ half-width) is committed
as two `private static final double` constants at the top of
`ConsistencyBound.java`, each with an inline `@see` comment citing
the table number, page, and value. Updating EHT inputs is a
one-line source edit + recompile.

Rationale: a single integer-and-two-floats lookup does not justify
a CSV parser or a new dependency, and pins the bound to a
specific, version-controlled EHT data version.

### 2.4 Two figures, single column

| Figure | Content |
|---|---|
| `asymmetry_vs_epsilon3.pdf` | δ_r/⟨r⟩ vs ε₃ for `i = {17°, 60°}`; horizontal shaded band for EHT m=1 ± 1σ; vertical dashed line at ε₃_crit ≈ +0.12; vertical shaded band over `(+0.13, +0.20)` for ε₃_ISCO; vertical solid lines at the resulting two-sided bound on the i = 17° curve |
| `image_gallery.pdf` | 1×4 grid of FITS frames at ε₃ ∈ {−2.5, −1.0, 0.0, +0.13}, `i = 17°` only, shared log color scale |

No third figure (no Schwarzschild reference inset, no
i = 60°-only inset). The reference-60° curve lives on Figure 1
as a second line.

### 2.5 Manuscript ≤ 1000 words, AASTeX 6.3.1

Single-column RNAAS template (`\documentclass[RNAAS]{aastex631}`).
Body word budget:

| Section | Budget (words) |
|---|---|
| Abstract | ≤ 80 |
| §1 Introduction | ≤ 160 |
| §2 Methods | ≤ 240 |
| §3 Results | ≤ 320 |
| §4 Discussion & caveats | ≤ 160 |
| Acknowledgments | ≤ 40 |
| **Total** | **≤ 1000** |

References ≤ 15. Figures ≤ 2. `paper/Makefile` enforces the word
count as a build-time gate (see §5 gate 6).

---

## 3. Files

### 3.1 Created (8)

```
src/main/java/com/pranav/grrt/analysis/ConsistencyBound.java
src/test/java/com/pranav/grrt/analysis/ConsistencyBoundTest.java
paper/manuscript.tex
paper/refs.bib
paper/Makefile
paper/scripts/make_figures.py
paper/scripts/requirements.txt
paper/figures/.gitkeep
```

### 3.2 Modified (3)

- `README.md` — add a "Reproducing the paper" section invoking
  `mvn verify -PrunSlow` (regenerates `output/sweep.csv` if absent)
  and `cd paper && make`.
- `.gitignore` — add `paper/figures/*.pdf`, `paper/*.aux`,
  `paper/*.log`, `paper/*.out`, `paper/*.pdf`,
  `paper/scripts/__pycache__/`.
- `CLAUDE.md` — tick the final `RNAAS paper` checkbox in the
  Status section.

### 3.3 Read-only inputs

- `output/sweep.csv` — 28 rows, 12 columns (schema in
  `docs/phase-3c-plan.md` §4). Both `delta_r_rms` and
  `fourier_m1_amp` columns present.
- `output/sweep_a0.9_eps*.fits` — 4 frames for the gallery
  (gallery selection chosen from 16 available at i = 17°).

---

## 4. Contracts for the three named files

### 4.1 `ConsistencyBound.java`

```java
public final class ConsistencyBound {

    /** EHT 2019 Paper VI Table 7, M87* m=1 fractional amplitude.
     *  Mid value to be filled at implementation time from the
     *  published table; this plan does not pre-commit a number. */
    private static final double EHT_M1_AMPLITUDE_MID = /* TBD */ 0.0;
    private static final double EHT_M1_AMPLITUDE_SIGMA = /* TBD */ 0.0;

    /** Domain on which the i=17° sweep is strictly monotone in
     *  δ_r/⟨r⟩. Hardcoded from sweep inspection at commit 9def498. */
    private static final double EPS3_DOMAIN_LO = -2.5;
    private static final double EPS3_DOMAIN_HI = -0.2;

    public record SweepPoint(double epsilon3, double meanR,
                             double deltaRrms, double fourierM1Amp) {}

    public enum BoundStatus { TWO_SIDED, UPPER_OPEN_PLATEAU,
                              LOWER_OPEN_BELOW_PATHOLOGY,
                              PLATEAU_CONSISTENT }

    public record Bound(double epsilonLower, double epsilonUpper,
                        BoundStatus status, String caveat) {}

    /** Reads i=17° rows from output/sweep.csv, filters to the
     *  monotone domain, inverts EHT band → (ε₃_low, ε₃_high)
     *  via piecewise-linear interpolation in Fourier m=1 space. */
    public static Bound invertFromCsv(Path csv);

    /** Converts CSV RMS column to Fourier-m=1 fractional amplitude
     *  using A_1 = √2 · (δ_r/⟨r⟩)_RMS when higher harmonics are
     *  negligible (cf. parent §7.3). Returns the measured Fourier
     *  m=1 directly when the column is present; falls back to the
     *  RMS conversion otherwise. */
    public static double[] toFourierM1(double[] rmsArray,
                                       double[] fourierColumn);
}
```

Interpolation is piecewise-linear in `(ε₃, A_1)` space on the five
domain points. Cubic spline is rejected: five-point splines can
ring near the right endpoint, and the 3C status doc shows the
extractor noise floor is ~1 %, well above any benefit from a
smoother interpolant.

ISCO-disappearance prediction is **not** an output of
`ConsistencyBound` — it is a constant quoted in the manuscript
from the 3C sweep-skip evidence (status doc §6.1). Keeping it
outside the algorithmic class avoids encoding a numerical value
that was determined by inspection rather than by a method this
class implements.

### 4.2 `make_figures.py`

Python 3.11+, dependencies pinned in `paper/scripts/requirements.txt`:

```
numpy>=2.0
matplotlib>=3.9
pandas>=2.2
astropy>=6.1
```

Single CLI entry point: `python make_figures.py` (no flags) reads
`../../output/sweep.csv` and `../../output/sweep_*.fits` relative
to the script location, writes `../figures/asymmetry_vs_epsilon3.pdf`
and `../figures/image_gallery.pdf`. Deterministic output (fixed
random seed for any jitter, fixed font size, frozen `matplotlib`
style). Re-running on identical inputs reproduces byte-identical
PDFs (gate 5).

The EHT band on Figure 1 is drawn from the same two constants the
Java class uses, redeclared at the top of the Python file with a
matching comment. A docstring-level note flags that these two
sources must be edited together.

### 4.3 `manuscript.tex`

AASTeX 6.3.1 RNAAS class. Section skeleton:

```latex
\documentclass[RNAAS]{aastex631}
\begin{document}
\title{An EHT M87* Photon-Ring Asymmetry Bound on a Single-Parameter
       Johannsen--Psaltis Deviation from Kerr}
\author{...}
\begin{abstract} ... \end{abstract}
\section{Introduction}        % §1, ≤160 words
\section{Methods}             % §2, ≤240 words; cites grrt repo commit
\section{Results}             % §3, ≤320 words; embeds both figures
\section{Discussion}          % §4, ≤160 words; caveats per §5.1-5.3 of 3C status
\acknowledgments              % ≤40 words
\bibliography{refs}
\end{document}
```

The Methods section names this repo's commit hash and the sweep CSV
git_sha (`d2e1da2`-equivalent per 3C status §2.3) so the work is
exactly reproducible.

---

## 5. Validation gates

Six gates. 1–3 are JUnit (run by `mvn verify -PrunSlow`); 4–6 are
build-time `make` targets in `paper/Makefile`.

| # | Gate | Tolerance / criterion |
|---|---|---|
| 1 | `ConsistencyBound.invertFromCsv` on a synthetic linear monotone CSV (no plateau) returns analytic `(ε₃_low, ε₃_high)` for a chosen synthetic EHT band | 1e-10 (analytic linear inversion is exact under linear interpolation) |
| 2 | `ConsistencyBound.toFourierM1` on a synthetic pure-cos(φ) ring profile recovers `A_1 = √2 · (δ_r/⟨r⟩)_RMS` | 1e-10 |
| 3 | `ConsistencyBound.invertFromCsv` raises `IllegalStateException` when EHT band lies entirely above the plateau (`PLATEAU_CONSISTENT`) and when it lies entirely below `δ_r/⟨r⟩(−2.5)` (`LOWER_OPEN_BELOW_PATHOLOGY`); returns `UPPER_OPEN_PLATEAU` when the upper-band crossing would require the plateau region | status-flag inspection + caveat-string non-empty |
| 4 | `make figures` produces both PDFs with byte-identical SHA-256 across two consecutive runs on the same `output/sweep.csv` | exact |
| 5 | `make paper` produces a valid PDF via `pdflatex` (or `tectonic` / `latexmk` if `pdflatex` missing); `paper/Makefile` autodetects | exit code 0; PDF page count ≥ 2 |
| 6 | `make wordcount` parses `manuscript.tex` (strips `\section*`, `\title`, `\acknowledgments`, math, comments) and reports body word count ≤ 1000; counts figures ≤ 2 and `\bibitem` refs ≤ 15 | strict |

Gates 1–3 use a stable synthetic CSV fixture committed under
`src/test/resources/analysis/synthetic_sweep.csv` so the unit
tests do not depend on `output/sweep.csv`. Gate 4's byte-stability
is the regenerability criterion from parent §7.4 gate 5.

Non-blocking flag (mirrors parent §7.4 gate 3): on the i = 17°
ε₃ = 0 sweep row, `|A_k / A_1| < 0.2` for `k ≥ 2`. Reported in
the manuscript Discussion section, not gated.

---

## 6. Commit chain

Three commits on top of `9def498`:

```
9def498  docs: phase-3c-complete status snapshot and 3D readiness  (parent)
<3D.1>   analysis: add ConsistencyBound with EHT m1 inversion (Phase 3D.1)
<3D.2>   paper: add make_figures.py and Makefile producing both PDFs (Phase 3D.2)
<3D.3>   paper: AASTeX manuscript and repo housekeeping (Phase 3D.3)
```

Tag at the third commit: **`phase-3-complete`**.

### 6.1 Commit 3D.1 — ConsistencyBound + tests

Files: `ConsistencyBound.java`, `ConsistencyBoundTest.java`,
`src/test/resources/analysis/synthetic_sweep.csv`.

Acceptance: gates 1–3 green via `mvn verify -PrunSlow`. Existing
141 tests still pass. No changes outside the `analysis` package.

### 6.2 Commit 3D.2 — paper scaffold + figures

Files: `paper/Makefile`, `paper/scripts/make_figures.py`,
`paper/scripts/requirements.txt`, `paper/figures/.gitkeep`,
`paper/refs.bib` (skeleton — citation entries added incrementally
in 3D.3).

Acceptance: `cd paper && make figures` produces both PDFs from
`output/sweep.csv` alone (gate 4). Visual sanity by user before
3D.3.

### 6.3 Commit 3D.3 — manuscript + housekeeping + tag

Files: `paper/manuscript.tex`, `paper/refs.bib` (filled),
modifications to `README.md`, `.gitignore`, `CLAUDE.md`.

Acceptance: gates 5–6 green via `cd paper && make`. The compiled
PDF is byte-identical across two consecutive `make clean && make`
invocations. Tag `phase-3-complete` placed on this commit only
after the user reviews the rendered PDF.

---

## 7. Wall-clock estimate (M3, 16 GB)

| Item | Estimate |
|---|---|
| 3D.1 `ConsistencyBound` + tests (~200 LOC + ~120 LOC tests) | one session, ~2 h |
| 3D.2 `make_figures.py` + `Makefile` + figure tuning | one session, ~2–3 h (figure aesthetics dominate) |
| 3D.3 manuscript first draft | one session, ~3 h (drafting + LaTeX wrangling) |
| 3D.3 manuscript revision passes (human-in-loop) | 2–3 short sessions, ~1 h each, not time-budgeted |
| Build / gate verification at each commit | ~15 min per commit |
| **Critical path to `phase-3-complete`** | **3 focused sessions + revision time** |

Compute time for the entire 3D phase is < 1 min wall-clock — no
new rendering, only CSV reads and PDF generation.

---

## 8. Primary failure modes

| Failure | Likelihood | Mitigation |
|---|---|---|
| EHT m=1 band entirely above plateau → no constraint on upper bound | likely (the plateau sits at A_1 ≈ √2 · 0.196 ≈ 0.277, well above any plausible EHT value) | This is the expected outcome under decision 2.1; `BoundStatus.UPPER_OPEN_PLATEAU` is the designed return; manuscript text already framed for it |
| LaTeX toolchain (`pdflatex` / `tectonic` / `latexmk`) absent on the user's machine | moderate | `paper/Makefile` autodetects all three and emits a clear `make: *** no LaTeX engine found` if none present; not a code bug |
| Python deps not installed | moderate | `paper/scripts/requirements.txt` plus a `make deps` target invoking `pip install -r scripts/requirements.txt`; README documents this in the "Reproducing the paper" section |
| Word count overrun | moderate | `make wordcount` is gate 6; iterative compression in 3D.3 revision passes; drop adjective-heavy phrasing before dropping content |
| Plateau-cluster interpretation challenged in peer review | possible | Discussion §4 explicitly states the plateau is the extractor noise floor (per status doc §5.1) and frames the upper bound as conservative; future work bullet mentions sub-pixel parabolic peak finding (per `docs/phase-3c-plan.md` §13) as a path to tighter upper bounds |
| EHT amplitude lies below `δ_r/⟨r⟩(−2.5)` → bound runs off the left edge | unlikely at a = 0.9 (the curve at ε₃ = −2.5 is A_1 ≈ 0.17, comparable to plausible EHT values) | `BoundStatus.LOWER_OPEN_BELOW_PATHOLOGY` returned with a caveat; would prompt extending the sweep toward ε₃_pathology ≈ −2.97 as 3E work, not blocking the note |

---

## 9. Tag at completion

`phase-3-complete` when:
- All 6 validation gates green at HEAD.
- `output/sweep.csv` unchanged from parent (no re-render).
- `paper/manuscript.pdf` exists locally and the user has reviewed it.
- `CLAUDE.md` Status section has the final checkbox ticked.

RNAAS submission is a separate human action and is **not** part of
the tag criteria.

---

## 10. Open questions to surface during 3D.1

These do not block planning, but the implementer (likely Claude
under user direction) should raise them before writing code:

1. **EHT numerical values.** What exact `(A_1_mid, A_1_sigma)`
   should be hardcoded in `ConsistencyBound.java`? Paper VI Table 7
   has multiple amplitude estimates across the imaging methods;
   the manuscript text should match the constants. User to supply
   or approve a specific row reference.
2. **Authorship list and ORCID.** Required by AASTeX. User to
   provide.
3. **Acknowledgment of computing resources.** Required by AAS
   journals; user to supply or approve "MacBook Pro M3 with 16 GB"
   as the literal text.
4. **`refs.bib` citation keys.** Convention preferred? `lastname:year`
   (e.g. `eht2019:paper-vi`, `johannsen:2011`) is the
   recommendation; user to approve or override.

Each is a ~30-second response from the user and prevents a
mid-implementation pause.

---

*End of Phase 3D plan.*
