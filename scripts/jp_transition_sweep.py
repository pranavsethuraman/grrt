#!/usr/bin/env python3
"""
scripts/jp_transition_sweep.py — Phase 3A gate-6 side-investigation.

At what eps_3 does the Kerr-continuation prograde equatorial photon
orbit cease to exist above the Kerr outer horizon r_plus_Kerr?

Sweeps a in {0.5, 0.7, 0.8, 0.9} at M=1, eps_3 from 0 in 0.01 steps,
reports the largest eps_3 for which r_photon_pro > r_plus_Kerr, then
bisects the transition to 1e-4 precision.

Scratch script: not part of the committed reference toolchain, not
needed at runtime. Kept alongside derive_jp_christoffels.py and
jp_photon_orbit_reference.py for provenance of Phase 3A parameter-
space analysis.
"""

import sys

try:
    import numpy as np
    import sympy as sp
    from scipy.optimize import brentq
except ImportError as exc:
    sys.stderr.write(f"ERROR: missing package {exc.name}\n")
    sys.exit(1)


# ---------------------------------------------------------------------------
# Symbolic setup (copy from jp_photon_orbit_reference.py — deliberately
# duplicate so this scratch script has no cross-file dependency)
# ---------------------------------------------------------------------------

r_sym, theta_sym = sp.symbols('r theta', real=True, positive=True)
M_sym, a_sym, eps3_sym = sp.symbols('M a eps3', real=True)
cos_t = sp.cos(theta_sym)
sin_t = sp.sin(theta_sym)
sin2 = sin_t ** 2
Sigma = r_sym ** 2 + a_sym ** 2 * cos_t ** 2
Delta = r_sym ** 2 - 2 * M_sym * r_sym + a_sym ** 2
h_scalar = eps3_sym * M_sym ** 3 * r_sym / Sigma ** 2
one_plus_h = 1 + h_scalar

g_tt_sym = -(1 - 2 * M_sym * r_sym / Sigma) * one_plus_h
g_tp_sym = -(2 * M_sym * a_sym * r_sym * sin2 / Sigma) * one_plus_h
g_pp_sym = sin2 * ((r_sym ** 2 + a_sym ** 2)
                    + (2 * M_sym * a_sym ** 2 * r_sym * sin2 / Sigma) * one_plus_h)

eq_subs = {theta_sym: sp.pi / 2}
g_tt_eq = g_tt_sym.subs(eq_subs)
g_tp_eq = g_tp_sym.subs(eq_subs)
g_pp_eq = g_pp_sym.subs(eq_subs)
dg_tt_eq = sp.diff(g_tt_eq, r_sym)
dg_tp_eq = sp.diff(g_tp_eq, r_sym)
dg_pp_eq = sp.diff(g_pp_eq, r_sym)

syms = (M_sym, a_sym, eps3_sym, r_sym)
gtt_fn  = sp.lambdify(syms, g_tt_eq,  'numpy')
gtp_fn  = sp.lambdify(syms, g_tp_eq,  'numpy')
gpp_fn  = sp.lambdify(syms, g_pp_eq,  'numpy')
dgtt_fn = sp.lambdify(syms, dg_tt_eq, 'numpy')
dgtp_fn = sp.lambdify(syms, dg_tp_eq, 'numpy')
dgpp_fn = sp.lambdify(syms, dg_pp_eq, 'numpy')


def residual(r_val, M_val, a_val, eps3_val, prograde):
    gtt  = gtt_fn(M_val,  a_val, eps3_val, r_val)
    gtp  = gtp_fn(M_val,  a_val, eps3_val, r_val)
    gpp  = gpp_fn(M_val,  a_val, eps3_val, r_val)
    dgtt = dgtt_fn(M_val, a_val, eps3_val, r_val)
    dgtp = dgtp_fn(M_val, a_val, eps3_val, r_val)
    dgpp = dgpp_fn(M_val, a_val, eps3_val, r_val)
    disc = gtp * gtp - gtt * gpp
    if disc < 0:
        raise ValueError
    sign = +1.0 if prograde else -1.0
    omega = (-gtp + sign * np.sqrt(disc)) / gpp
    return dgtt + 2.0 * dgtp * omega + dgpp * omega * omega


def photon_orbit_above_kerr_horizon(M_val, a_val, eps3_val, prograde=True):
    """r_photon in (r_plus_Kerr, 15M), NaN if no root there."""
    r_plus = M_val + np.sqrt(M_val * M_val - a_val * a_val)
    rs = np.linspace(r_plus * 1.0001, 15.0 * M_val, 1000)
    valid_rs, valid_fs = [], []
    for r_val in rs:
        try:
            f_val = residual(r_val, M_val, a_val, eps3_val, prograde)
            if np.isfinite(f_val):
                valid_rs.append(r_val)
                valid_fs.append(f_val)
        except ValueError:
            continue
    if len(valid_rs) < 2:
        return float('nan')
    brackets = [(valid_rs[i], valid_rs[i + 1])
                for i in range(len(valid_fs) - 1)
                if valid_fs[i] * valid_fs[i + 1] < 0.0]
    if not brackets:
        return float('nan')
    r_lo_b, r_hi_b = brackets[-1]
    try:
        return brentq(
            lambda rr: residual(rr, M_val, a_val, eps3_val, prograde),
            r_lo_b, r_hi_b, xtol=1e-14, rtol=1e-14, maxiter=200)
    except Exception:
        return float('nan')


def find_eps3_crit(a_val, eps3_cap=3.0, coarse_step=0.01, tol=1e-4):
    """Largest eps_3 at which the Kerr-continuation prograde photon orbit
    exists above r_plus_Kerr. Coarse sweep locates the transition, then
    bisect to tol."""
    eps3 = 0.0
    last_valid = 0.0
    last_valid_r = None
    first_invalid = None
    while eps3 <= eps3_cap:
        r_pro = photon_orbit_above_kerr_horizon(1.0, a_val, eps3, True)
        if np.isnan(r_pro):
            first_invalid = eps3
            break
        last_valid = eps3
        last_valid_r = r_pro
        eps3 += coarse_step
    if first_invalid is None:
        return None, None
    lo = last_valid
    hi = first_invalid
    while hi - lo > tol:
        mid = 0.5 * (lo + hi)
        r_mid = photon_orbit_above_kerr_horizon(1.0, a_val, mid, True)
        if np.isnan(r_mid):
            hi = mid
        else:
            lo = mid
            last_valid_r = r_mid
    return lo, last_valid_r


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

print("# JP Kerr-continuation prograde photon-orbit transition")
print("# Generated by scripts/jp_transition_sweep.py")
print("# M = 1. eps_3_crit = largest eps_3 with r_photon_pro > r_plus_Kerr.")
print("# Transition bisected to tol = 1e-4.")
print("#")
print("#   a      eps3_crit    r_photon_at_crit    r_plus_Kerr")
print("# --------------------------------------------------------")

for a in [0.5, 0.7, 0.8, 0.9]:
    r_plus = 1.0 + np.sqrt(1.0 - a * a)
    eps3_crit, r_at_crit = find_eps3_crit(a)
    if eps3_crit is None:
        print(f"  {a:.3f}    (no transition in [0, 3])    —              {r_plus:.6f}")
    else:
        print(f"  {a:.3f}    {eps3_crit:.4f}       {r_at_crit:.6f}        {r_plus:.6f}")

print()

# Full sweep at a=0.9 for documentation (the user's primary spin)
print("# Full sweep at a = 0.9, eps_3 ∈ [0.00, 0.50], step 0.01:")
print("# eps_3   r_photon_pro   (NaN = no Kerr-continuation orbit above r_plus_Kerr)")
for i in range(51):
    eps3 = round(0.01 * i, 4)
    r_pro = photon_orbit_above_kerr_horizon(1.0, 0.9, eps3, True)
    r_str = f"{r_pro:.10f}" if np.isfinite(r_pro) else "NaN"
    print(f"  {eps3:+.4f}   {r_str}")
