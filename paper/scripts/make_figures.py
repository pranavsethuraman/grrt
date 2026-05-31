#!/usr/bin/env python3
"""Generate the two RNAAS figures for the grrt Phase 3D note.

Reads ``../../output/sweep.csv`` and the per-frame FITS images
(``../../output/sweep_a0.9_eps*_i17.0_res512.fits``) relative to this
script's location, and writes

  * ``../figures/asymmetry_vs_epsilon3.pdf`` -- ring circularity
    ``delta_r/<r>`` vs the Johannsen-Psaltis deviation ``eps3`` for
    i = 17 deg and 60 deg, with the EHT 2019 Paper VI Sec. 7.4 / Fig. 18
    circularity band shown as a horizontal reference. The band falls
    *below* both curves: the i = 17 deg extractor noise floor (~0.12)
    exceeds the EHT circularity (~0.055), so no two-sided eps3 bound is
    placed (the honest-null result of Phase 3D; see
    ``ConsistencyBound`` Javadoc).
  * ``../figures/image_gallery.pdf`` -- a 1x4 gallery of the i = 17 deg
    disk images at eps3 in {-2.5, -1.0, 0.0, +0.13} on a shared log scale.

No CLI flags: ``python make_figures.py`` does everything. Output is
deterministic -- re-running on identical inputs yields byte-identical
PDFs (docs/phase-3d-plan.md Sec. 5, gate 4). Determinism relies on
``SOURCE_DATE_EPOCH`` (fixed PDF timestamps), the bundled DejaVu/serif
fonts, and a frozen rcParams block.

EHT constants below are duplicated from
``ConsistencyBound.EHT_CIRCULARITY_MID`` / ``_SIGMA``. **Edit both
together** -- the Java bound and this figure must cite the same number.
"""

import os

# Must be set before matplotlib is imported so the PDF backend stamps a
# fixed CreationDate/ModDate (reproducibility, gate 4). 2010-01-01 UTC.
os.environ.setdefault("SOURCE_DATE_EPOCH", "1262304000")

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("pdf")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import LogNorm  # noqa: E402
from astropy.io import fits  # noqa: E402

# --- EHT 2019 Paper VI, Sec. 7.4 "The Circular Shapes of the M87 Images"
#     and Figure 18: fractional diameter spread peaks at ~0.05-0.06.
#     Keep in sync with ConsistencyBound.EHT_CIRCULARITY_MID / _SIGMA. ---
EHT_CIRC_MID = 0.055
EHT_CIRC_SIGMA = 0.005

# Johannsen-Psaltis landmarks at a = 0.9 (docs/jp-parameter-space-notes.md).
EPS3_CRIT = 0.1212            # prograde photon-orbit bifurcation (cusp)
EPS3_ISCO_LO, EPS3_ISCO_HI = 0.13, 0.20   # ISCO-disappearance band
EPS3_DOMAIN_LO, EPS3_DOMAIN_HI = -2.5, -0.2  # monotone inversion sub-interval

SPIN = 0.9
RESOLUTION = 512
GALLERY_EPS = (-2.5, -1.0, 0.0, 0.13)
GALLERY_INCL = 17.0
GALLERY_DECADES = 3.0         # log color scale spans this many decades below max

SCRIPT_DIR = Path(__file__).resolve().parent
OUTPUT_DIR = SCRIPT_DIR.parent.parent / "output"
FIG_DIR = SCRIPT_DIR.parent / "figures"
CSV_PATH = OUTPUT_DIR / "sweep.csv"


def configure_matplotlib() -> None:
    """Freeze a deterministic, LaTeX-friendly style."""
    plt.rcParams.update({
        "font.family": "serif",
        "font.size": 9,
        "axes.titlesize": 9,
        "axes.labelsize": 9,
        "legend.fontsize": 7.5,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "axes.linewidth": 0.8,
        "lines.linewidth": 1.1,
        "pdf.fonttype": 42,
        "figure.dpi": 200,
        "savefig.dpi": 200,
    })


def fits_path(eps: float, incl: float) -> Path:
    """Frame path in the EpsilonSweep naming convention."""
    return OUTPUT_DIR / (
        f"sweep_a{SPIN:.1f}_eps{eps:+.4f}_i{incl:.1f}_res{RESOLUTION}.fits"
    )


def load_curve(df: pd.DataFrame, incl: float):
    """Return (eps3, delta_r/<r>) sorted by eps3 for one inclination."""
    sub = df[np.isclose(df["inclination_deg"], incl)].sort_values("epsilon_3")
    eps = sub["epsilon_3"].to_numpy()
    circ = (sub["delta_r_rms"] / sub["mean_r"]).to_numpy()
    return eps, circ


def make_asymmetry_figure(df: pd.DataFrame) -> Path:
    fig, ax = plt.subplots(figsize=(3.4, 2.9))

    # ISCO-disappearance band and cusp line (context for results 1-2).
    ax.axvspan(EPS3_ISCO_LO, EPS3_ISCO_HI, color="#cfe3f3", alpha=0.85,
               zorder=0, label=r"$\epsilon_3^{\rm ISCO}\in(0.13,0.20)$")
    ax.axvline(EPS3_CRIT, ls="--", lw=0.9, color="0.35", zorder=1)

    # EHT 2019 circularity band -- horizontal reference; falls below curves.
    # § = section sign; matplotlib's default text renderer has no \S.
    ax.axhspan(EHT_CIRC_MID - EHT_CIRC_SIGMA, EHT_CIRC_MID + EHT_CIRC_SIGMA,
               color="0.78", alpha=0.9, zorder=0,
               label="EHT 2019 §7.4 circularity")

    series = [
        (17.0, "o", "#1b4965", r"$i=17^\circ$"),
        (60.0, "s", "#bc4749", r"$i=60^\circ$"),
    ]
    for incl, marker, color, label in series:
        eps, circ = load_curve(df, incl)
        ax.plot(eps, circ, marker=marker, ms=3.2, color=color, label=label)

    ax.set_xlim(-2.65, 0.26)
    ax.set_ylim(0.0, 0.52)
    ax.set_xlabel(r"$\epsilon_3$")
    ax.set_ylabel(r"$\delta_r / \langle r \rangle$")

    # Cusp label: vertical, set between the two curves so it clears the
    # top legend strip and the markers.
    ax.text(EPS3_CRIT - 0.05, 0.295, r"$\epsilon_3^{\rm crit}\approx0.12$",
            rotation=90, ha="center", va="center", fontsize=6.5, color="0.35")

    ax.legend(loc="upper center", frameon=False, ncol=2, fontsize=7,
              columnspacing=1.2, handletextpad=0.4, borderaxespad=0.3)
    fig.tight_layout(pad=0.4)
    out = FIG_DIR / "asymmetry_vs_epsilon3.pdf"
    fig.savefig(out)
    plt.close(fig)
    return out


def make_gallery_figure() -> Path:
    frames = []
    for eps in GALLERY_EPS:
        path = fits_path(eps, GALLERY_INCL)
        if not path.exists():
            raise FileNotFoundError(f"gallery frame missing: {path}")
        frames.append((eps, np.asarray(fits.getdata(path), dtype=float)))

    positive = np.concatenate([d[d > 0.0].ravel() for _, d in frames])
    if positive.size == 0:
        raise ValueError("all gallery frames are non-positive; cannot log-scale")
    vmax = float(positive.max())
    vmin = vmax / (10.0 ** GALLERY_DECADES)
    norm = LogNorm(vmin=vmin, vmax=vmax)

    fig, axes = plt.subplots(1, len(frames), figsize=(6.8, 1.95))
    im = None
    for ax, (eps, data) in zip(axes, frames):
        im = ax.imshow(np.clip(data, vmin, None), origin="lower",
                       norm=norm, cmap="inferno", interpolation="nearest")
        ax.set_title(rf"$\epsilon_3={eps:+.2f}$")
        ax.set_xticks([])
        ax.set_yticks([])
    cbar = fig.colorbar(im, ax=axes, fraction=0.046, pad=0.02)
    cbar.set_label("disk intensity (arb., log)", fontsize=7.5)
    cbar.ax.tick_params(labelsize=7)
    out = FIG_DIR / "image_gallery.pdf"
    fig.savefig(out, bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)
    return out


def main() -> int:
    if not CSV_PATH.exists():
        raise FileNotFoundError(
            f"sweep CSV not found: {CSV_PATH}\n"
            "Run `mvn verify -PrunSlow` (or the EpsilonSweep main) first."
        )
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    configure_matplotlib()

    df = pd.read_csv(CSV_PATH)
    fig1 = make_asymmetry_figure(df)
    fig2 = make_gallery_figure()
    print(f"wrote {fig1}")
    print(f"wrote {fig2}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
