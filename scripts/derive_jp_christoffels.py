#!/usr/bin/env python3
"""
scripts/derive_jp_christoffels.py
=================================

Symbolic derivation of the Johannsen-Psaltis metric Christoffel symbols
for grrt Phase 3A.  The deformation form is Johannsen 2013 / Psaltis 2020:

    h(r, theta) = eps3 * M^3 * r / Sigma^2
    Sigma       = r^2 + a^2 cos^2(theta)
    Delta       = r^2 - 2 M r + a^2

The metric ansatz is Johannsen-Psaltis 2011 eq. (1) with this single-parameter
h(r, theta).  See docs/phase-3-plan.md sections 2 and 4 for the full plan.

WHAT THIS SCRIPT DOES
---------------------

Five phases, all printed to stdout:

  1. Coordinate pathology grid check -- numerically scan 1 + h across
     eps3 in [-2, +2], a = 0.9, r in [r+, 100 M], theta in [0, pi].
     Prints PASSED/FAILED.  Aborts if 1 + h <= 0 anywhere in the range.
     This is the belt-and-suspenders counterpart to the analytic bound
     recorded in docs/phase-3-plan.md section 3.1.

  2. Symbolic construction of g and the closed-form inverse g^{-1};
     verification of g * g^{-1} = I at three numerical test points.

  3. Derivation of all non-zero Christoffel symbols
       Gamma^alpha_{mu nu} = (1/2) g^{alpha sigma}
             * (d_mu g_{sigma nu} + d_nu g_{sigma mu} - d_sigma g_{mu nu}).

  4. Kerr-reduction cross-check at eps3 = 0 against an independently built
     pure-Kerr Christoffel pipeline at (r = 10 M, theta = pi/3) to 1e-12.

  5. Common-subexpression elimination via sympy.cse and emission of the
     resulting expressions as Java-syntax comments ready for manual paste
     into JohannsenPsaltisMetric.java.

REQUIREMENTS
------------

  python3, sympy, numpy.  Not runtime dependencies of the Maven build.
  Installation: `pip3 install sympy numpy` (or use a venv).

RUNTIME
-------

  ~15-45 s on Apple M3.  The expensive step is the Christoffel loop
  plus sympy.cse.  Intermediate simplification is skipped to stay fast.

USAGE
-----

  python3 scripts/derive_jp_christoffels.py
  python3 scripts/derive_jp_christoffels.py > jp_christoffels.txt

REGENERATION POLICY
-------------------

The JohannsenPsaltisMetric Java source is a handwritten translation of
this script's output.  When the metric form changes, re-run this script
and rewrite the Java by hand using the updated expressions as the
authoritative reference.  Do not auto-generate Java from this script
-- the translation is trivial once CSE has been applied, and a
human-reviewed Java implementation keeps the hot-path code readable.
"""

import sys
import time

try:
    import numpy as np
    import sympy as sp
except ImportError as exc:
    sys.stderr.write(
        f"ERROR: required package missing ({exc.name}).\n"
        "Install with: pip3 install sympy numpy\n"
    )
    sys.exit(1)


# ---------------------------------------------------------------------------
# Phase 1 -- coordinate pathology grid check
# ---------------------------------------------------------------------------

def pathology_check() -> None:
    """Scan 1 + h over a (eps3, r, theta) grid; abort if it reaches 0.

    This is the numerical counterpart of the analytic bound in
    docs/phase-3-plan.md section 3.1.  The analytic bound is
    authoritative; this check catches arithmetic transcription errors.
    """
    print("=" * 72)
    print("Phase 1 -- coordinate pathology grid check")
    print("=" * 72)

    a_val = 0.9
    m_val = 1.0
    r_plus = m_val + np.sqrt(m_val ** 2 - a_val ** 2)

    eps_values = np.array([-2.0, -1.0, -0.5, -0.2, -0.1, -0.05,
                           0.0,
                           0.05, 0.1, 0.2, 0.5, 1.0, 2.0])
    r_values = np.linspace(r_plus, 100.0, 400)
    theta_values = np.linspace(0.0, np.pi, 101)

    h_min = np.inf
    h_max = -np.inf
    arg_min = None

    for eps in eps_values:
        for rv in r_values:
            for th in theta_values:
                sigma = rv * rv + a_val * a_val * np.cos(th) ** 2
                h = eps * (m_val ** 3) * rv / (sigma * sigma)
                val = 1.0 + h
                if val < h_min:
                    h_min = val
                    arg_min = (float(eps), float(rv), float(th))
                if val > h_max:
                    h_max = val

    status = "PASSED" if h_min > 0.0 else "FAILED"
    print(f"  grid: eps3 in [-2, 2] (13 pts), r in [{r_plus:.4f}, 100] (400 pts), theta in [0, pi] (101 pts)")
    print(f"  1 + h range: [{h_min:.6f}, {h_max:.6f}]")
    print(f"  min at (eps3, r, theta) = {arg_min}")
    print(f"  coordinate pathology grid check: {status}")
    print()

    if h_min <= 0.0:
        sys.exit(f"ABORT: 1 + h reached {h_min:.6e}; gate 2 numerical check failed.")


pathology_check()


# ---------------------------------------------------------------------------
# Phase 2 -- symbolic metric + closed-form inverse
# ---------------------------------------------------------------------------

print("=" * 72)
print("Phase 2 -- symbolic metric and inverse")
print("=" * 72)

# symbolic coordinates and parameters
t, r, th, phi = sp.symbols('t r theta phi', real=True)
M, a, eps3 = sp.symbols('M a eps3', real=True)
coords = [t, r, th, phi]

cos_th = sp.cos(th)
sin_th = sp.sin(th)
sin2_th = sin_th ** 2

Sigma = r ** 2 + a ** 2 * cos_th ** 2
Delta = r ** 2 - 2 * M * r + a ** 2
h_scalar = eps3 * M ** 3 * r / Sigma ** 2
one_plus_h = 1 + h_scalar

# JP 2011 eq. (1) components with Johannsen 2013 h(r, theta)
g_tt = -(1 - 2 * M * r / Sigma) * one_plus_h
g_tp = -(2 * M * a * r * sin2_th / Sigma) * one_plus_h
g_rr = Sigma * one_plus_h / (Delta + a ** 2 * h_scalar * sin2_th)
g_thth = Sigma
g_pp = sin2_th * ((r ** 2 + a ** 2)
                  + (2 * M * a ** 2 * r * sin2_th / Sigma) * one_plus_h)

g = sp.Matrix([
    [g_tt, 0,      0,       g_tp],
    [0,    g_rr,   0,       0],
    [0,    0,      g_thth,  0],
    [g_tp, 0,      0,       g_pp],
])

# closed-form inverse: 2x2 (t, phi) block + diagonal (r, theta)
D_tp = g_tt * g_pp - g_tp ** 2
g_inv = sp.Matrix([
    [g_pp / D_tp,   0,         0,           -g_tp / D_tp],
    [0,             1 / g_rr,  0,           0],
    [0,             0,         1 / g_thth,  0],
    [-g_tp / D_tp,  0,         0,           g_tt / D_tp],
])

# numerical spot-checks of g * g_inv = I
print("  Numerical g * g^-1 = I check at three points:")
test_points = [
    dict(M_v=1.0, a_v=0.9, eps_v=0.5,  r_v=10.0, th_v=sp.pi / 4),
    dict(M_v=1.0, a_v=0.5, eps_v=-1.0, r_v=5.0,  th_v=sp.pi / 2),
    dict(M_v=1.0, a_v=0.9, eps_v=0.0,  r_v=20.0, th_v=sp.pi / 3),
]
for pt in test_points:
    subs = {M: pt['M_v'], a: pt['a_v'], eps3: pt['eps_v'],
            r: pt['r_v'], th: pt['th_v']}
    prod = (g.subs(subs) * g_inv.subs(subs)).evalf()
    err = max(abs(float(prod[i, j] - (1 if i == j else 0)))
              for i in range(4) for j in range(4))
    print(f"    a={pt['a_v']}, eps3={pt['eps_v']:+.2f}, r={pt['r_v']}, "
          f"theta={pt['th_v']}: max |err| = {err:.2e}")
    if err > 1e-10:
        sys.exit(f"ABORT: g * g^-1 != I at {subs}; max err = {err:.2e}")
print("  OK")
print()


# ---------------------------------------------------------------------------
# Phase 3 -- Christoffel symbols
# ---------------------------------------------------------------------------

def build_christoffels(g_matrix: sp.Matrix, g_matrix_inv: sp.Matrix,
                       label: str) -> dict:
    """Return {(alpha, mu, nu): expression} with mu <= nu, non-zero only.

    Uses stationarity and axisymmetry: only d_r (index 1) and d_theta
    (index 2) contribute.
    """
    print(f"  [{label}] differentiating metric components...")
    t0 = time.time()
    dg = {}
    for mu in (1, 2):
        for sigma in range(4):
            for nu in range(4):
                dg[(mu, sigma, nu)] = sp.diff(g_matrix[sigma, nu], coords[mu])
    print(f"    partials ready ({time.time() - t0:.1f} s)")

    print(f"  [{label}] assembling Gamma^alpha_{{mu nu}}...")
    t1 = time.time()
    gamma = {}
    for alpha in range(4):
        for mu in range(4):
            for nu in range(mu, 4):
                acc = sp.S.Zero
                for sigma in range(4):
                    gi = g_matrix_inv[alpha, sigma]
                    if gi == 0:
                        continue
                    term = sp.S.Zero
                    if (mu, sigma, nu) in dg:
                        term += dg[(mu, sigma, nu)]
                    if (nu, sigma, mu) in dg:
                        term += dg[(nu, sigma, mu)]
                    if (sigma, mu, nu) in dg:
                        term -= dg[(sigma, mu, nu)]
                    acc += gi * term
                expr = acc / 2
                if expr != 0:
                    gamma[(alpha, mu, nu)] = expr
    print(f"    {len(gamma)} non-zero Christoffels "
          f"({time.time() - t1:.1f} s)")
    return gamma


print("=" * 72)
print("Phase 3 -- JP Christoffel symbols")
print("=" * 72)
Gamma = build_christoffels(g, g_inv, "JP")
print()


# ---------------------------------------------------------------------------
# Phase 4 -- Kerr reduction cross-check
# ---------------------------------------------------------------------------

print("=" * 72)
print("Phase 4 -- Kerr reduction at eps3 = 0")
print("=" * 72)

# Independent Kerr pipeline.
g_K = sp.Matrix([
    [-(1 - 2 * M * r / Sigma), 0, 0, -(2 * M * a * r * sin2_th / Sigma)],
    [0, Sigma / Delta,         0, 0],
    [0, 0,                     Sigma, 0],
    [-(2 * M * a * r * sin2_th / Sigma), 0, 0,
     sin2_th * ((r ** 2 + a ** 2) + 2 * M * a ** 2 * r * sin2_th / Sigma)],
])
D_K = g_K[0, 0] * g_K[3, 3] - g_K[0, 3] ** 2
g_K_inv = sp.Matrix([
    [g_K[3, 3] / D_K, 0, 0, -g_K[0, 3] / D_K],
    [0, 1 / g_K[1, 1], 0, 0],
    [0, 0, 1 / g_K[2, 2], 0],
    [-g_K[0, 3] / D_K, 0, 0, g_K[0, 0] / D_K],
])

Gamma_K = build_christoffels(g_K, g_K_inv, "Kerr")

# Numerical comparison at eps3 = 0, (r = 10 M, theta = pi/3, a = 0.9, M = 1)
subs_num = {M: 1.0, a: 0.9, r: 10.0, th: sp.pi / 3, eps3: 0.0}
subs_num_K = {k: v for k, v in subs_num.items() if k is not eps3}

all_keys = set(Gamma.keys()) | set(Gamma_K.keys())
max_err = 0.0
max_err_key = None
for key in all_keys:
    jp_val = float(Gamma.get(key, sp.S.Zero).subs(subs_num).evalf())
    k_val = float(Gamma_K.get(key, sp.S.Zero).subs(subs_num_K).evalf())
    diff = abs(jp_val - k_val)
    if diff > max_err:
        max_err = diff
        max_err_key = key

print()
print(f"  Compared {len(all_keys)} Christoffels at eps3=0, (a=0.9, r=10, theta=pi/3)")
print(f"  max |Gamma_JP(0) - Gamma_Kerr| = {max_err:.2e} at {max_err_key}")
if max_err > 1e-12:
    sys.exit(f"ABORT: Kerr reduction failed, max_err = {max_err:.2e}")
print("  Kerr reduction: PASSED")
print()


# ---------------------------------------------------------------------------
# Phase 5 -- CSE and emit Java-syntax expressions
# ---------------------------------------------------------------------------

print("=" * 72)
print("Phase 5 -- CSE and Java-syntax emission")
print("=" * 72)

exprs = list(Gamma.values())
keys = list(Gamma.keys())

print("  Applying sympy.cse...")
t2 = time.time()
cse_subs, cse_reduced = sp.cse(exprs, optimizations='basic')
print(f"    {len(cse_subs)} common subexpressions, "
      f"{len(cse_reduced)} Christoffels ({time.time() - t2:.1f} s)")
print()


def to_java(expr: sp.Expr) -> str:
    """Translate a Sympy expression into Java-compatible source.

    Uses Sympy's C99 printer as a starting point and rewrites the
    function calls and integer powers to idiomatic Java.
    """
    s = sp.ccode(expr, standard='C99')
    replacements = [
        ('sin(', 'Math.sin('),
        ('cos(', 'Math.cos('),
        ('tan(', 'Math.tan('),
        ('sqrt(', 'Math.sqrt('),
        ('fabs(', 'Math.abs('),
        ('pow(', 'Math.pow('),
    ]
    for old, new in replacements:
        s = s.replace(old, new)
    return s


NAMES = {0: 't', 1: 'r', 2: 'th', 3: 'phi'}

print("// " + "=" * 70)
print("// JP Christoffels -- generated by scripts/derive_jp_christoffels.py")
print("// Metric form: h(r, theta) = eps3 * M^3 * r / Sigma^2 (Johannsen 2013)")
print("// Signature: (-, +, +, +), coordinates (t, r, theta, phi).")
print("// Regenerate when the metric form changes; translation to Java is")
print("// manual (see script docstring REGENERATION POLICY).")
print("// " + "=" * 70)
print()
print("// Common subexpressions:")
for sym, expr in cse_subs:
    print(f"final double {sym} = {to_java(expr)};")
print()
print("// Non-zero Christoffels (symmetric in lower indices):")
for key, expr in zip(keys, cse_reduced):
    alpha, mu, nu = key
    label = f"G_{NAMES[alpha]}__{NAMES[mu]}_{NAMES[nu]}"
    print(f"final double {label} = {to_java(expr)};")

print()
print("// End generated Christoffel block.")
print()
print("Done.")
