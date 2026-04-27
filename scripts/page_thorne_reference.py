#!/usr/bin/env python3
"""
scripts/page_thorne_reference.py
================================

Reference values for the Page-Thorne 1974 radial energy flux F(r, a)
from a Novikov-Thorne 1973 thin accretion disk in Kerr equatorial.
Output is pasted by hand into NovikovThorneDiskTest.java as the
reference table for Phase 3B gate 1 (see docs/phase-3b-plan.md §2.2
and §4 gate 1).

Flux formula (Page-Thorne 1974, eq. 15n; cf. Krolik 1999 §7.5)
--------------------------------------------------------------

    F(r, a) = -1 / (4 pi r sqrt(Delta))
              * (dOmega/dr) / (E - Omega L)^2
              * integral_{r_ISCO}^{r} (E - Omega L) dL/dr' dr'

with M = 1, accretion rate Mdot normalized to 1. F is a pure-numerical
function of (r, a) in arbitrary geometrized units that
NovikovThorneDisk.surfaceFlux must replicate.

Bardeen-Press-Teukolsky 1972 closed forms for Kerr equatorial prograde
circular orbits (M = 1):

    Omega(r, a) = 1 / (r^{3/2} + a)
    E(r, a)     = (r^2 - 2 r + a sqrt(r))
                  / (r sqrt(r^2 - 3 r + 2 a sqrt(r)))
    L(r, a)     = sqrt(r) (r^2 - 2 a sqrt(r) + a^2)
                  / (r sqrt(r^2 - 3 r + 2 a sqrt(r)))

ISCO from Bardeen 1972 closed form:

    Z1   = 1 + (1 - a^2)^{1/3} ((1+a)^{1/3} + (1-a)^{1/3})
    Z2   = sqrt(3 a^2 + Z1^2)
    rms  = 3 + Z2 - sqrt((3 - Z1)(3 + Z1 + 2 Z2))   (prograde)

Method
------

1. Build Omega, E, L, Delta symbolically in Sympy as functions of (r, a).
2. Differentiate symbolically to obtain dOmega/dr, dL/dr.
3. Lambdify all six functions for fast numerical evaluation.
4. For each (a, r) test point:
   a. Compute r_ISCO from Bardeen 1972.
   b. scipy.integrate.quad the integrand (E - Omega L) dL/dr from
      r_ISCO + delta to r, where delta = 1e-12 keeps quad's adaptive
      sampler away from any sqrt-tail at the endpoint denominator.
   c. Multiply by the prefactor.
5. Emit a whitespace-delimited table to stdout. Columns:
     a   r_over_M   F_pageThorne

Requires: python3, sympy, numpy, scipy. Same dev-only toolchain as
scripts/jp_photon_orbit_reference.py.

Runtime: ~10 s on Apple M3.

Usage
-----

    python3 scripts/page_thorne_reference.py

Regeneration policy
-------------------

If NovikovThorneDisk.surfaceFlux ever changes form, rerun this script
and paste the new table into NovikovThorneDiskTest.java. Same pattern
as scripts/jp_photon_orbit_reference.py.
"""

import sys

try:
    import numpy as np
    import sympy as sp
    from scipy.integrate import quad
except ImportError as exc:
    sys.stderr.write(
        f"ERROR: required package missing ({exc.name}).\n"
        "Install with: pip3 install sympy numpy scipy\n"
    )
    sys.exit(1)


# ---------------------------------------------------------------------------
# Symbolic Kerr equatorial circular-orbit quantities (M = 1)
# ---------------------------------------------------------------------------

r, a = sp.symbols('r a', real=True, positive=True)

# Bardeen-Press-Teukolsky 1972 prograde equatorial circular orbits.
sqrt_r = sp.sqrt(r)
common_inner = r ** 2 - 3 * r + 2 * a * sqrt_r
common_den = r * sp.sqrt(common_inner)

Omega_sym = 1 / (r ** sp.Rational(3, 2) + a)
E_sym = (r ** 2 - 2 * r + a * sqrt_r) / common_den
L_sym = sqrt_r * (r ** 2 - 2 * a * sqrt_r + a ** 2) / common_den

Delta_sym = r ** 2 - 2 * r + a ** 2

dOmega_sym = sp.diff(Omega_sym, r)
dL_sym = sp.diff(L_sym, r)

syms = (r, a)
omega_fn   = sp.lambdify(syms, Omega_sym,  'numpy')
e_fn       = sp.lambdify(syms, E_sym,      'numpy')
l_fn       = sp.lambdify(syms, L_sym,      'numpy')
delta_fn   = sp.lambdify(syms, Delta_sym,  'numpy')
domega_fn  = sp.lambdify(syms, dOmega_sym, 'numpy')
dl_fn      = sp.lambdify(syms, dL_sym,     'numpy')


# ---------------------------------------------------------------------------
# ISCO closed form (Bardeen 1972)
# ---------------------------------------------------------------------------

def r_isco_prograde(a_val: float) -> float:
    """Prograde Kerr ISCO radius in units of M = 1 (Bardeen 1972)."""
    if not (0.0 <= a_val < 1.0):
        raise ValueError(f"require 0 <= a < 1; got a = {a_val}")
    Z1 = 1.0 + (1.0 - a_val ** 2) ** (1.0 / 3.0) * (
        (1.0 + a_val) ** (1.0 / 3.0) + (1.0 - a_val) ** (1.0 / 3.0)
    )
    Z2 = np.sqrt(3.0 * a_val ** 2 + Z1 ** 2)
    return 3.0 + Z2 - np.sqrt((3.0 - Z1) * (3.0 + Z1 + 2.0 * Z2))


# ---------------------------------------------------------------------------
# Page-Thorne 1974 eq. (15n) flux
# ---------------------------------------------------------------------------

def page_thorne_integrand(rp: float, a_val: float) -> float:
    """Inner-integral integrand: (E - Omega L) dL/dr at r' = rp."""
    return (e_fn(rp, a_val) - omega_fn(rp, a_val) * l_fn(rp, a_val)) * dl_fn(rp, a_val)


def page_thorne_flux(r_val: float, a_val: float) -> float:
    """F(r, a) per the formula in the module docstring.

    The integrand is regular on (r_ISCO, infty); the denominator
    sqrt(r^2 - 3r + 2a sqrt(r)) in E and L vanishes at the photon
    orbit, well inside ISCO, so the lower limit is safely above any
    singularity. We still bump by 1e-12 to keep quad's adaptive
    sampler away from any endpoint sqrt-tail.
    """
    r_ms = r_isco_prograde(a_val)
    if r_val < r_ms:
        raise ValueError(
            f"r = {r_val} below ISCO (r_ms = {r_ms:.6f}) for a = {a_val}"
        )

    # Inner integral from ISCO to r.
    integral, abserr = quad(
        page_thorne_integrand,
        r_ms + 1e-12, r_val,
        args=(a_val,),
        epsabs=1e-14, epsrel=1e-12, limit=200,
    )

    if abserr > 1e-10:
        sys.stderr.write(
            f"WARNING: quad abserr = {abserr:.3e} at r = {r_val}, a = {a_val}\n"
        )

    # Prefactor: -1 / (4 pi r sqrt(Delta)) * dOmega/dr / (E - Omega L)^2
    delta_r = delta_fn(r_val, a_val)
    if delta_r <= 0.0:
        raise ValueError(
            f"Delta <= 0 at r = {r_val}, a = {a_val}: {delta_r}"
        )
    sqrt_delta = np.sqrt(delta_r)
    e_val = e_fn(r_val, a_val)
    l_val = l_fn(r_val, a_val)
    om_val = omega_fn(r_val, a_val)
    dom_val = domega_fn(r_val, a_val)

    prefactor = -1.0 / (4.0 * np.pi * r_val * sqrt_delta) \
        * dom_val / (e_val - om_val * l_val) ** 2

    return prefactor * integral


# ---------------------------------------------------------------------------
# Reference grid: 6 radii for each of two spins
# ---------------------------------------------------------------------------

# Schwarzschild (a = 0): r_ISCO = 6 M; sample from just outside ISCO out
# to 20 M.
schwarz_radii = [6.5, 7.0, 8.0, 10.0, 15.0, 20.0]

# Kerr a = 0.9: r_ISCO ≈ 2.32088 M; sample from just outside out to 20 M.
kerr_radii = [2.5, 3.0, 4.0, 6.0, 10.0, 20.0]

print("# Page-Thorne 1974 eq. (15n) radial flux from Novikov-Thorne disk")
print("# Generated by scripts/page_thorne_reference.py")
print("# Convention: Krolik 1999 sec. 7.5; M = 1; Mdot = 1; geometrized units")
print("# F(r, a) = -1/(4 pi r sqrt(Delta)) * dOmega/dr / (E - Omega L)^2")
print("#          * integral_{r_ISCO}^{r} (E - Omega L) dL/dr' dr'")
print("# Columns: a  r_over_M  F_pageThorne")
print("#")

print(f"# r_ISCO(a=0.0) = {r_isco_prograde(0.0):.10f}")
print(f"# r_ISCO(a=0.9) = {r_isco_prograde(0.9):.10f}")
print("#")

for a_val, radii in [(0.0, schwarz_radii), (0.9, kerr_radii)]:
    for r_val in radii:
        f_val = page_thorne_flux(r_val, a_val)
        print(f"  {a_val:+.4f}  {r_val:7.3f}   {f_val:.15e}")

print()
print("Done.")
