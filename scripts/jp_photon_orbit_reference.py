#!/usr/bin/env python3
"""
scripts/jp_photon_orbit_reference.py
====================================

Reference values for the equatorial photon-orbit radius in the
Johannsen-Psaltis metric at seven (a, eps3) pairs. Output is pasted by
hand into JohannsenPsaltisMetricTest.java as the reference table for
gate 6 (see docs/phase-3a-gates-5-6.md §3.3 and §3.4).

Method
------

1. Build the JP metric symbolically in Sympy, substitute theta = pi/2
   to get the equatorial metric components g_tt(r), g_tphi(r), g_phiphi(r)
   as functions of (M, a, eps3, r).
2. Differentiate each symbolically w.r.t. r to obtain ∂_r g_{αβ}. This
   path is deliberately INDEPENDENT of the Java `eCircular` closed-form
   r-derivatives (sub-plan §3.1 three-layer validation).
3. For each (a, eps3) pair, find r such that the combined
   null + circular condition at theta = pi/2 is satisfied:
     null     : g_tt + 2 g_tφ ω + g_φφ ω² = 0
                 → ω_± = (−g_tφ ± √(g_tφ² − g_tt g_φφ)) / g_φφ
     circular : Γ^r_{tt} + 2 Γ^r_{tφ} ω + Γ^r_{φφ} ω² = 0
                 → at equator with diagonal (r, θ) block of g^{-1}:
                   ∂_r g_tt + 2 ∂_r g_tφ ω + ∂_r g_φφ ω² = 0
   r_photon is the root of F(r) ≡ LHS of (circular) with ω substituted
   from (null). Prograde uses the + root, retrograde the − root.
4. Emit a whitespace-delimited table to stdout. Columns:
     a   eps3   r_photon_pro   r_photon_retro

Requires: python3, sympy, numpy, scipy. Same dev-only toolchain as
scripts/derive_jp_christoffels.py.

Runtime: ~15 s on Apple M3.

Usage
-----

    python3 scripts/jp_photon_orbit_reference.py

Regeneration policy
-------------------

If the JP metric form ever changes, rerun this script and paste the
new table into JohannsenPsaltisMetricTest.photonOrbitRadiusMatchesPythonReference.
Same pattern as scripts/derive_jp_christoffels.py.
"""

import sys

try:
    import numpy as np
    import sympy as sp
    from scipy.optimize import brentq
except ImportError as exc:
    sys.stderr.write(
        f"ERROR: required package missing ({exc.name}).\n"
        "Install with: pip3 install sympy numpy scipy\n"
    )
    sys.exit(1)


# ---------------------------------------------------------------------------
# Symbolic JP metric at theta = pi/2
# ---------------------------------------------------------------------------

r, theta = sp.symbols('r theta', real=True, positive=True)
M, a, eps3 = sp.symbols('M a eps3', real=True)

cos_t = sp.cos(theta)
sin_t = sp.sin(theta)
sin2 = sin_t ** 2

Sigma = r ** 2 + a ** 2 * cos_t ** 2
Delta = r ** 2 - 2 * M * r + a ** 2
h_scalar = eps3 * M ** 3 * r / Sigma ** 2
one_plus_h = 1 + h_scalar

g_tt_sym = -(1 - 2 * M * r / Sigma) * one_plus_h
g_tp_sym = -(2 * M * a * r * sin2 / Sigma) * one_plus_h
g_pp_sym = sin2 * ((r ** 2 + a ** 2)
                    + (2 * M * a ** 2 * r * sin2 / Sigma) * one_plus_h)

# Substitute equator
eq_subs = {theta: sp.pi / 2}
g_tt_eq = g_tt_sym.subs(eq_subs)
g_tp_eq = g_tp_sym.subs(eq_subs)
g_pp_eq = g_pp_sym.subs(eq_subs)

dg_tt_eq = sp.diff(g_tt_eq, r)
dg_tp_eq = sp.diff(g_tp_eq, r)
dg_pp_eq = sp.diff(g_pp_eq, r)

syms = (M, a, eps3, r)
gtt_fn  = sp.lambdify(syms, g_tt_eq,  'numpy')
gtp_fn  = sp.lambdify(syms, g_tp_eq,  'numpy')
gpp_fn  = sp.lambdify(syms, g_pp_eq,  'numpy')
dgtt_fn = sp.lambdify(syms, dg_tt_eq, 'numpy')
dgtp_fn = sp.lambdify(syms, dg_tp_eq, 'numpy')
dgpp_fn = sp.lambdify(syms, dg_pp_eq, 'numpy')


# ---------------------------------------------------------------------------
# Null + circular residual F(r) and root finder
# ---------------------------------------------------------------------------

def residual(r_val: float, M_val: float, a_val: float, eps3_val: float,
             prograde: bool) -> float:
    """F(r) at the chosen prograde/retrograde branch of omega.

    Raises ValueError if the null-condition discriminant is negative
    at this r (no real omega; physically inside the photon trapping
    window for that branch).
    """
    gtt  = gtt_fn(M_val,  a_val, eps3_val, r_val)
    gtp  = gtp_fn(M_val,  a_val, eps3_val, r_val)
    gpp  = gpp_fn(M_val,  a_val, eps3_val, r_val)
    dgtt = dgtt_fn(M_val, a_val, eps3_val, r_val)
    dgtp = dgtp_fn(M_val, a_val, eps3_val, r_val)
    dgpp = dgpp_fn(M_val, a_val, eps3_val, r_val)

    disc = gtp * gtp - gtt * gpp
    if disc < 0:
        raise ValueError(
            f"null discriminant negative: r={r_val}, a={a_val}, "
            f"eps3={eps3_val}, disc={disc}"
        )
    sign = +1.0 if prograde else -1.0
    omega = (-gtp + sign * np.sqrt(disc)) / gpp
    return dgtt + 2.0 * dgtp * omega + dgpp * omega * omega


def photon_orbit(M_val: float, a_val: float, eps3_val: float,
                 prograde: bool) -> float:
    """Root-find F(r) = 0 by scanning a log-spaced grid for sign changes.

    Returns NaN if no sign change of F is found on the branch — this
    indicates the photon orbit on that branch (prograde/retrograde)
    does not exist for this (a, eps3), which is a real physical feature
    of the JP metric at large |eps3| (see Johannsen 2013 §IV).
    """
    if a_val * a_val >= M_val * M_val:
        raise ValueError(f"non-subextremal spin: a={a_val}, M={M_val}")
    r_plus = M_val + np.sqrt(M_val * M_val - a_val * a_val)

    # Scan log-spaced grid from well below the Kerr horizon (the JP
    # equatorial horizon can sit inside r_plus_Kerr for eps3 > 0) up to
    # 15 M. 500 points resolves roots to ~0.6% before brentq refines.
    r_min = 0.3 * M_val
    r_max = 15.0 * M_val
    rs = np.logspace(np.log10(r_min), np.log10(r_max), 500)

    valid_rs = []
    valid_fs = []
    for r_val in rs:
        try:
            f_val = residual(r_val, M_val, a_val, eps3_val, prograde)
        except ValueError:
            continue
        if np.isfinite(f_val):
            valid_rs.append(r_val)
            valid_fs.append(f_val)

    if len(valid_rs) < 2:
        return float('nan')

    # Find sign changes. Take the outermost — it corresponds to the
    # photon orbit proper; inner sign changes (if any) are inside the
    # trapping region and not relevant here.
    brackets = []
    for i in range(len(valid_fs) - 1):
        if valid_fs[i] * valid_fs[i + 1] < 0:
            brackets.append((valid_rs[i], valid_rs[i + 1]))

    if not brackets:
        return float('nan')

    r_lo, r_hi = brackets[-1]
    return brentq(
        lambda rr: residual(rr, M_val, a_val, eps3_val, prograde),
        r_lo, r_hi, xtol=1e-14, rtol=1e-14, maxiter=200,
    )


# ---------------------------------------------------------------------------
# Reference-table pairs and output
# ---------------------------------------------------------------------------

pairs = [
    (0.0,  0.0),    # Schwarzschild, r_ph = 3 M (exact)
    (0.9,  0.0),    # Kerr a=0.9, BardeenShadowTest cross-check
    (0.9, +0.5),    # primary sweep point
    (0.9, -0.5),    # primary sweep point
    (0.9, +1.0),    # sweep-density boundary
    (0.9, -1.0),    # sweep-density boundary
    (0.5, +0.5),    # intermediate-a coverage (amendment 3)
]

print("# Johannsen-Psaltis equatorial photon orbit reference (M = 1)")
print("# Generated by scripts/jp_photon_orbit_reference.py")
print("# Metric form: Johannsen 2013 h(r, theta) = eps3 * M^3 * r / Sigma^2")
print("# Columns: a  eps3  r_photon_pro  r_photon_retro")
print("#")

for (a_val, eps3_val) in pairs:
    r_pro   = photon_orbit(1.0, a_val, eps3_val, prograde=True)
    r_retro = photon_orbit(1.0, a_val, eps3_val, prograde=False)
    print(f"  {a_val:+.4f}  {eps3_val:+.4f}  "
          f"{r_pro:.12f}  {r_retro:.12f}")

print()
print("Done.")
