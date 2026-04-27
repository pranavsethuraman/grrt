package com.pranav.grrt.metric;

/**
 * Johannsen-Psaltis metric: a single-parameter deviation from Kerr
 * parameterized by ε₃, in Boyer-Lindquist-like coordinates (t, r, θ, φ).
 *
 * <h2>Metric form (Johannsen 2013 / Psaltis et al. 2020)</h2>
 *
 * Deformation scalar:
 * <pre>
 *   h(r, θ) = ε₃ M³ r / Σ²,        Σ = r² + a² cos²θ
 *   Δ       = r² − 2 M r + a²
 * </pre>
 *
 * Line element with this h (Johannsen-Psaltis 2011, Phys. Rev. D 83, 124015
 * eq. 1):
 * <pre>
 *   ds² = -(1 − 2Mr/Σ)(1 + h) dt²
 *         − (4 M a r sin²θ / Σ)(1 + h) dt dφ
 *         + Σ (1 + h) / (Δ + a² h sin²θ) dr²
 *         + Σ dθ²
 *         + sin²θ [(r² + a²) + (2 M a² r sin²θ / Σ)(1 + h)] dφ²
 * </pre>
 *
 * Setting ε₃ = 0 ⇒ h = 0 collapses every component to the Kerr line element
 * implemented by {@link KerrMetric}.
 *
 * <h2>Inverse metric (closed form)</h2>
 *
 * The (r, θ) block is diagonal. The (t, φ) 2×2 block is inverted directly
 * using the determinant D = g_{tt} g_{φφ} − g_{tφ}². For the JP form above
 * this gives (where 1+h ≠ 0 holds across the approved sweep range):
 * <pre>
 *   g^{rr}   = (Δ + a² h sin²θ) / [Σ (1 + h)]
 *   g^{θθ}   = 1/Σ
 *   g^{tt}   =  g_{φφ} / D
 *   g^{φφ}   =  g_{tt} / D
 *   g^{tφ}   = -g_{tφ} / D
 * </pre>
 * No 4×4 generic inversion is performed, matching the pattern used by
 * {@link KerrMetric}.
 *
 * <h2>Component derivation sketches</h2>
 *
 * <ul>
 *   <li>{@code g_tt = -(1 − 2Mr/Σ)(1+h)} — JP eq. 1 multiplies the Kerr
 *       {@code g_tt^K = -(1 − 2Mr/Σ)} by the conformal-like factor (1+h).</li>
 *   <li>{@code g_tφ = -(2Mar sin²θ/Σ)(1+h)} — same (1+h) multiplication on
 *       the Kerr frame-dragging term.</li>
 *   <li>{@code g_rr = Σ(1+h)/(Δ + a² h sin²θ)} — the (1+h) in the numerator
 *       is inherited from the ansatz; the {@code a² h sin²θ} in the
 *       denominator is the JP correction whose vanishing (at ε₃ &lt; 0)
 *       locates the θ-dependent JP horizon outside the Kerr horizon.</li>
 *   <li>{@code g_θθ = Σ} — unmodified; the deformation leaves the θ direction
 *       intact.</li>
 *   <li>{@code g_φφ = sin²θ[(r²+a²) + (2Ma²r sin²θ/Σ)(1+h)]} — the (1+h)
 *       multiplies only the frame-dragging piece of {@code g_φφ^K}; the
 *       radial piece {@code (r²+a²) sin²θ} stays unmodified.</li>
 * </ul>
 *
 * <h2>Horizon</h2>
 *
 * The true event horizon of the JP-deformed geometry is the outermost root of
 * {@code Δ + a² h sin²θ = 0}, which is θ-dependent. On the polar axis
 * {@code (sin θ = 0)} this reduces to {@code Δ = 0} → Kerr
 * {@code r₊ = M + √(M² − a²)}. For {@code ε₃ < 0} the equatorial horizon
 * moves outward of the Kerr horizon; for {@code ε₃ > 0} it moves inward.
 * {@link #horizonRadius()} returns the polar-axis value as the interface
 * requires a single scalar; ray termination in the tracer uses a
 * position-dependent test (see Phase 3 plan §2.5).
 *
 * <h2>Coordinate-pathology range</h2>
 *
 * For the approved sweep range {@code ε₃ ∈ [−2, +2]}, {@code a = 0.9 M},
 * {@code r ≥ r₊_Kerr}, {@code θ ∈ [0, π]}, the factor {@code (1 + h)} stays
 * in {@code [0.324, 1.676]}: it never reaches zero and no (1+h)-driven
 * coordinate pathology arises. This is verified analytically in
 * {@code docs/phase-3-plan.md} §3.1 and numerically by the grid scan in
 * {@code scripts/derive_jp_christoffels.py}. The constructor enforces the
 * equivalent pointwise bound {@code 1 + ε₃ M³ / r₊³ > 0} given the
 * current {@code (M, a)} pair.
 *
 * <h2>Christoffel symbols</h2>
 *
 * The 20 non-zero Christoffels were derived symbolically in
 * {@code scripts/derive_jp_christoffels.py} with Sympy and cross-checked
 * against an independent Kerr pipeline at {@code ε₃ = 0},
 * {@code (r = 10 M, θ = π/3)} to machine precision. Common-subexpression
 * elimination was applied and the resulting block is pasted verbatim into
 * {@link #fillChristoffelScratch}. Regenerate by re-running the Sympy
 * script; do not hand-edit.
 *
 * <h2>Hot-path discipline</h2>
 *
 * {@link #geodesicAcceleration} is the integrator inner-loop hot path and
 * is fully non-allocating: it computes the 20 Christoffels via
 * {@link #fillChristoffelScratch} into a 20-double array local to the call
 * (escape-analysis eligible) and contracts them directly into the provided
 * {@code out}. This matches the pattern used by
 * {@link KerrMetric#geodesicAcceleration}. {@link #g}, {@link #gInv}, and
 * {@link #christoffel} allocate their return values as the {@link Metric}
 * interface contract mandates (each returns a freshly allocated array).
 *
 * <h2>References</h2>
 *
 * <ul>
 *   <li>Johannsen &amp; Psaltis (2011), Phys. Rev. D 83, 124015 — metric ansatz.</li>
 *   <li>Johannsen (2013), Phys. Rev. D 88, 044002 — θ-dependent h form.</li>
 *   <li>Psaltis et al. (2020), Phys. Rev. Lett. 125, 141104 — EHT bound
 *       comparator.</li>
 *   <li>{@code scripts/derive_jp_christoffels.py} — authoritative
 *       derivation, rerun to regenerate the Christoffel block.</li>
 *   <li>{@code docs/phase-3-plan.md} §2, §3, §4 — full plan and validation
 *       gates.</li>
 * </ul>
 */
public final class JohannsenPsaltisMetric implements Metric {

    // Christoffel scratch layout (filled by fillChristoffelScratch).
    // Matches the emission order of scripts/derive_jp_christoffels.py.
    private static final int G_T_TR   = 0;
    private static final int G_T_TTH  = 1;
    private static final int G_T_RPH  = 2;
    private static final int G_T_THPH = 3;
    private static final int G_R_TT   = 4;
    private static final int G_R_TPH  = 5;
    private static final int G_R_RR   = 6;
    private static final int G_R_RTH  = 7;
    private static final int G_R_THTH = 8;
    private static final int G_R_PHPH = 9;
    private static final int G_TH_TT   = 10;
    private static final int G_TH_TPH  = 11;
    private static final int G_TH_RR   = 12;
    private static final int G_TH_RTH  = 13;
    private static final int G_TH_THTH = 14;
    private static final int G_TH_PHPH = 15;
    private static final int G_PH_TR   = 16;
    private static final int G_PH_TTH  = 17;
    private static final int G_PH_RPH  = 18;
    private static final int G_PH_THPH = 19;

    private static final int CHRISTOFFEL_SCRATCH_SIZE = 20;

    private final double M;
    private final double a;
    private final double eps3;

    /**
     * @param mass     gravitational mass M in geometrized units (G = c = 1);
     *                 must be positive and finite
     * @param spin     Kerr spin parameter a in geometrized units; must be
     *                 finite with |a| &lt; M (subextremal)
     * @param epsilon3 JP deformation parameter ε₃; must be finite and satisfy
     *                 {@code 1 + ε₃ M³ / r₊³ > 0} (the equatorial-horizon
     *                 pathology bound — see class Javadoc)
     * @throws IllegalArgumentException on any boundary violation
     */
    public JohannsenPsaltisMetric(double mass, double spin, double epsilon3) {
        if (mass <= 0.0 || !Double.isFinite(mass)) {
            throw new IllegalArgumentException(
                    "Mass must be positive and finite: " + mass);
        }
        if (!Double.isFinite(spin)) {
            throw new IllegalArgumentException("Spin must be finite: " + spin);
        }
        if (Math.abs(spin) >= mass) {
            throw new IllegalArgumentException(
                    "Spin must be subextremal |a| < M: a = " + spin + ", M = " + mass);
        }
        if (!Double.isFinite(epsilon3)) {
            throw new IllegalArgumentException(
                    "epsilon3 must be finite: " + epsilon3);
        }
        double rPlus = mass + Math.sqrt(mass * mass - spin * spin);
        double hMin = epsilon3 * mass * mass * mass / (rPlus * rPlus * rPlus);
        if (1.0 + hMin <= 0.0) {
            throw new IllegalArgumentException(
                    "Coordinate pathology: 1 + h reaches " + (1.0 + hMin)
                    + " at (r = r+, equator) for epsilon3 = " + epsilon3
                    + ", a = " + spin + ", M = " + mass);
        }
        this.M = mass;
        this.a = spin;
        this.eps3 = epsilon3;
    }

    /** Unit-mass JP hole (M = 1) with the given spin and ε₃. */
    public JohannsenPsaltisMetric(double spin, double epsilon3) {
        this(1.0, spin, epsilon3);
    }

    @Override public double mass() { return M; }

    /** @return Kerr spin parameter a in geometrized units */
    public double spin() { return a; }

    /** @return JP deformation parameter ε₃ */
    public double epsilon3() { return eps3; }

    /**
     * Outer event horizon on the polar axis, {@code r₊ = M + √(M² − a²)}.
     *
     * <p>The true JP horizon is θ-dependent (the outermost root of
     * {@code Δ + a² h sin²θ = 0}). Off-axis shifts are second-order in ε₃
     * and are deliberately ignored by this single-scalar accessor per
     * {@code docs/phase-3-plan.md} §2.5. Ray integration should use a
     * position-dependent termination test when ε₃ ≠ 0.
     */
    @Override
    public double horizonRadius() {
        return M + Math.sqrt(M * M - a * a);
    }

    /**
     * Innermost stable circular equatorial (ISCO) timelike-orbit radius.
     *
     * <p>Root-find on {@code dE_circ/dr = 0} where {@code E_circ(r)} is the
     * specific energy of an equatorial circular timelike orbit at radius r.
     * For a stationary axisymmetric metric with equatorial components
     * {@code g_tt(r), g_tφ(r), g_φφ(r)} the circular-orbit condition
     * {@code Γ^r_{μν} u^μ u^ν = 0} gives
     * <pre>
     *   Ω_± = (−g_tφ' ± √(g_tφ'² − g_tt' g_φφ')) / g_φφ'
     *   E   = -(g_tt + g_tφ Ω) / √(-g_tt - 2 g_tφ Ω - g_φφ Ω²)
     * </pre>
     * with the + root for prograde (co-rotating) orbits in the {@code a > 0}
     * convention. The ISCO is the r at which {@code dE/dr = 0}.
     *
     * <p>Solver: bracket + bisection. The starting lower bound
     * {@code 1.1 · r₊} may sit inside the photon sphere (for Schwarzschild
     * prograde, or for retrograde at any spin where the retrograde photon
     * sphere lies well outside the horizon). In that region no timelike
     * circular orbit exists and {@code eCircular} throws. We walk {@code rLo}
     * outward by a factor of 1.05 until {@code dEdrCircular} returns
     * without error, then bisect {@code [rLo, 20 M]} to relative tolerance
     * 1e-12 or 100 iterations. After the walk, {@code dE/dr} is negative
     * in the (unstable) circular-orbit region {@code r_photon < r < r_ISCO}
     * and positive outside {@code r > r_ISCO}, so the bisection sign change
     * locates the ISCO. Kerr (ε₃ = 0) recovery is verified against the
     * Bardeen-Press-Teukolsky closed form in the test suite.
     *
     * <p>{@code dE/dr} is computed by 4-point centered finite differences
     * with step {@code h = 1e-3 M}. Truncation error is {@code O(h⁴) ~ 1e-12};
     * round-off contribution {@code O(ε/h) ~ 2e-13}. The 2-point scheme
     * previously used here gave only {@code ~2e-10} precision, which after
     * dividing by the (small) retrograde {@code d²E/dr²} pushed the
     * retrograde r_ISCO error up to {@code ~1.5e-8}. The 4-point scheme
     * brings both prograde and retrograde r_ISCO below {@code 1e-10}.
     * This is not a hot path; called once during renderer setup.
     *
     * @param prograde true for co-rotating ISCO, false for counter-rotating
     * @return ISCO radius in geometrized units
     * @throws IllegalStateException if the bracket fails
     */
    @Override
    public double iscoRadius(boolean prograde) {
        double rHi = 20.0 * M;
        double rLo = 1.1 * horizonRadius();

        // Walk rLo outward until a valid timelike circular orbit exists.
        // Required for retrograde at any spin (retrograde photon sphere
        // sits well outside the horizon) and for Schwarzschild prograde
        // (photon sphere at 3M > 1.1 r+ = 2.2M).
        double fLo;
        while (true) {
            try {
                fLo = dEdrCircular(rLo, prograde);
                break;
            } catch (IllegalStateException ex) {
                rLo *= 1.05;
                if (rLo >= rHi) {
                    throw new IllegalStateException(
                            "ISCO bracket failed: no timelike circular orbit below r = "
                            + rHi + " (a = " + a + ", eps3 = " + eps3
                            + ", prograde = " + prograde + ")");
                }
            }
        }
        double fHi = dEdrCircular(rHi, prograde);
        if (Math.signum(fLo) == Math.signum(fHi) && fLo != 0.0) {
            throw new IllegalStateException(
                    "ISCO bracket failed: dE/dr has same sign at r = "
                    + rLo + " (" + fLo + ") and r = " + rHi + " (" + fHi + ")");
        }
        for (int iter = 0; iter < 100; iter++) {
            double rMid = 0.5 * (rLo + rHi);
            double fMid = dEdrCircular(rMid, prograde);
            if (Math.signum(fMid) == Math.signum(fLo)) {
                rLo = rMid; fLo = fMid;
            } else {
                rHi = rMid;
            }
            if (rHi - rLo < 1e-12 * Math.abs(rMid)) return rMid;
        }
        return 0.5 * (rLo + rHi);
    }

    /**
     * 4-point centered finite-difference approximation of {@code dE/dr}
     * along the equatorial circular-orbit sequence:
     * <pre>
     *   f'(r) ≈ [-f(r+2h) + 8 f(r+h) - 8 f(r-h) + f(r-2h)] / (12 h) + O(h⁴)
     * </pre>
     * Step {@code h = 1e-3 M} is near-optimal. Total error is
     * {@code |f⁽⁵⁾| h⁴/30 + 1.5 ε_m |f|/h}; minimizing gives
     * {@code h_opt = (11.25 ε_m |f|/|f⁽⁵⁾|)^(1/5) ≈ 1.2e-3 M} for
     * {@code ε_m ≈ 2.22e-16} and well-scaled {@code E(r)} near ISCO radii.
     * Residual error at this step is {@code ~3e-13}, which after division
     * by {@code |d²E/dr²(r_ISCO)|} gives r_ISCO precision well below 1e-10
     * for both prograde and retrograde Kerr at a = 0.9.
     *
     * <p>The naive 2-point stencil at {@code h = 1e-6 M} had roundoff-limited
     * precision {@code ~2e-10} in dE/dr and produced {@code ~1.5e-8} error
     * in retrograde r_ISCO (shallow {@code d²E/dr²}); that scheme is
     * superseded by this one. General rule: optimal h for a central
     * difference of order p scales as {@code ε_m^(1/(p+1))}; 2-point gives
     * {@code ~6e-6}, 4-point gives {@code ~1e-3}, 6-point gives {@code ~6e-3}.
     */
    private double dEdrCircular(double r, boolean prograde) {
        double h = 1e-3 * M;
        double fm2 = eCircular(r - 2.0 * h, prograde);
        double fm1 = eCircular(r - h,       prograde);
        double fp1 = eCircular(r + h,       prograde);
        double fp2 = eCircular(r + 2.0 * h, prograde);
        return (-fp2 + 8.0 * fp1 - 8.0 * fm1 + fm2) / (12.0 * h);
    }

    /**
     * Diagnostic hook: returns the post-walk lower bracket {@code rLo} that
     * {@link #iscoRadius(boolean)} would hand to its bisection loop. Package
     * private; used by the test suite to confirm that the walking-bracket
     * lands at distinct positions for prograde vs retrograde without
     * exposing the private {@link #dEdrCircular}.
     */
    double iscoBracketLowerBound(boolean prograde) {
        double rHi = 20.0 * M;
        double rLo = 1.1 * horizonRadius();
        while (true) {
            try {
                dEdrCircular(rLo, prograde);
                return rLo;
            } catch (IllegalStateException ex) {
                rLo *= 1.05;
                if (rLo >= rHi) {
                    throw new IllegalStateException(
                            "iscoBracketLowerBound: walk exceeded rHi = " + rHi);
                }
            }
        }
    }

    private double eCircular(double r, boolean prograde) {
        double h = eps3 * M * M * M / (r * r * r);
        double hPrime = -3.0 * eps3 * M * M * M / (r * r * r * r);
        double onePlusH = 1.0 + h;
        double twoMOverR = 2.0 * M / r;
        double twoMaOverR = 2.0 * M * a / r;
        double gtt = -(1.0 - twoMOverR) * onePlusH;
        double gtp = -twoMaOverR * onePlusH;
        double gpp = (r * r + a * a) + (2.0 * M * a * a / r) * onePlusH;
        double dgtt = -(2.0 * M / (r * r)) * onePlusH - (1.0 - twoMOverR) * hPrime;
        double dgtp = (2.0 * M * a / (r * r)) * onePlusH - twoMaOverR * hPrime;
        double dgpp = 2.0 * r + 2.0 * M * a * a * (hPrime / r - onePlusH / (r * r));
        double disc = dgtp * dgtp - dgtt * dgpp;
        if (disc < 0.0) {
            throw new IllegalStateException(
                    "No circular orbit at r = " + r + "; discriminant = " + disc);
        }
        double sqrtDisc = Math.sqrt(disc);
        double sgn = prograde ? +1.0 : -1.0;
        double omega = (-dgtp + sgn * sqrtDisc) / dgpp;
        double denomSq = -gtt - 2.0 * gtp * omega - gpp * omega * omega;
        if (denomSq <= 0.0) {
            throw new IllegalStateException(
                    "No timelike circular orbit at r = " + r + "; denom² = " + denomSq);
        }
        return -(gtt + gtp * omega) / Math.sqrt(denomSq);
    }

    // ------------------------------------------------------------------
    // Equatorial photon-orbit radius (independent code path from eCircular)
    // ------------------------------------------------------------------

    /**
     * Residual {@code F(r) = ∂_r g_tt + 2 ∂_r g_tφ ω + ∂_r g_φφ ω²} for
     * the equatorial circular-null condition, evaluated on the prograde
     * or retrograde branch of {@code ω(r)} from the null condition
     * {@code g_tt + 2 g_tφ ω + g_φφ ω² = 0}. F vanishes at the photon-
     * orbit radius.
     *
     * <p>Equatorial g and ∂_r g are written inline, INDEPENDENTLY of the
     * expressions used in {@link #eCircular} (the ISCO helper). Same
     * math, separate code paths — see
     * {@code docs/phase-3a-gates-5-6.md} §3.4 amendment 4. The
     * package-private {@link #equatorialGDerivativesViaPhotonPath}
     * mirrors this method's expressions line-for-line for the
     * cross-path agreement test.
     *
     * @throws IllegalStateException if the null discriminant is
     *         negative (no real ω at this r) or {@code g_φφ ≤ 0}
     */
    private double residualPhotonOrbit(double r, boolean prograde) {
        // Equatorial: Σ = r²; h = ε₃ M³/r³; h' = -3 ε₃ M³/r⁴.
        final double r2 = r * r;
        final double r3 = r2 * r;
        final double r4 = r3 * r;
        final double M3 = M * M * M;
        final double Ma = M * a;
        final double Ma2 = Ma * a;
        final double e = eps3;

        final double oneMinusTwoMOverR = 1.0 - 2.0 * M / r;
        final double twoMOverR2 = 2.0 * M / r2;
        final double minus3eM3OverR4 = -3.0 * e * M3 / r4;
        final double eM3OverR3 = e * M3 / r3;
        final double onePlusH_ = 1.0 + eM3OverR3;

        final double gtt_ = -oneMinusTwoMOverR * onePlusH_;
        final double gtp_ = -(2.0 * Ma / r) * onePlusH_;
        final double gpp_ = (r2 + a * a) + (2.0 * Ma2 / r) * onePlusH_;

        final double dgtt_ = -twoMOverR2 * onePlusH_ - oneMinusTwoMOverR * minus3eM3OverR4;
        final double dgtp_ = (2.0 * Ma / r2) * onePlusH_ - (2.0 * Ma / r) * minus3eM3OverR4;
        final double dgpp_ = 2.0 * r + 2.0 * Ma2 * (minus3eM3OverR4 / r - onePlusH_ / r2);

        final double disc = gtp_ * gtp_ - gtt_ * gpp_;
        if (disc < 0.0) {
            throw new IllegalStateException(
                    "Null discriminant negative at r = " + r
                    + " (a = " + a + ", eps3 = " + eps3 + "): disc = " + disc);
        }
        if (gpp_ <= 0.0) {
            throw new IllegalStateException(
                    "g_φφ non-positive at r = " + r + ": " + gpp_);
        }
        final double sgn = prograde ? +1.0 : -1.0;
        final double omega = (-gtp_ + sgn * Math.sqrt(disc)) / gpp_;
        return dgtt_ + 2.0 * dgtp_ * omega + dgpp_ * omega * omega;
    }

    /**
     * Equatorial photon-orbit radius for the prograde or retrograde
     * branch.
     *
     * <p>Solves {@code F(r) = 0} where F is the combined null + circular
     * null-orbit residual evaluated by {@link #residualPhotonOrbit}.
     * The combined system does not reduce to a low-degree polynomial
     * for the JP metric (the {@code h(r,θ) = ε₃ M³ r/Σ²} deformation
     * introduces nontrivial r-dependence beyond Kerr's quartic), so we
     * locate the root numerically:
     *
     * <ol>
     *   <li>Scan a 500-point linear grid on
     *       {@code [1.001 r₊, 10 M]} for points where the null
     *       discriminant is non-negative (i.e., the ω-branch is real).</li>
     *   <li>Identify the OUTERMOST sign change of F(r) on the valid grid.</li>
     *   <li>Bisect the bracketing interval to relative tolerance
     *       {@code 1e-12} or 100 iterations.</li>
     * </ol>
     *
     * <p><b>Why a grid scan, not the iscoRadius walk-from-1.1·r₊
     * pattern.</b> For high-spin prograde photon orbits the simple-walk
     * start point {@code 1.1 r₊} can land just OUTSIDE
     * {@code r_photon} (e.g., for Kerr {@code a = 0.9}: {@code r_photon
     * ≈ 1.558 M}, {@code 1.1 r₊ ≈ 1.580 M}), so the bisection bracket
     * {@code [1.1 r₊, 10 M]} contains no F sign change and the bisection
     * fails. The grid scan is robust for both branches across the
     * approved sweep range. This matches the pattern used in
     * {@code scripts/jp_photon_orbit_reference.py}.
     *
     * <p><b>Implementation independence.</b> The equatorial g and ∂_r g
     * expressions inside {@link #residualPhotonOrbit} are written
     * independently of those used inside {@link #eCircular}. The
     * package-private helpers
     * {@link #equatorialGDerivativesViaIscoPath} and
     * {@link #equatorialGDerivativesViaPhotonPath} expose the two
     * paths' expressions for the cross-path agreement test (see
     * {@code docs/phase-3a-gates-5-6.md} §3.4 amendment 4).
     *
     * @param prograde true for the + branch of ω (co-rotating with spin
     *                 at large r); false for retrograde
     * @return r_photon in geometrized units
     * @throws IllegalStateException if no F sign change is found in
     *         the search interval (the Kerr-continuation orbit does not
     *         exist above the Kerr horizon for this (a, ε₃))
     */
    public double photonOrbitRadius(boolean prograde) {
        final int N = 500;
        final double rLoStart = 1.001 * horizonRadius();
        final double rHi = 10.0 * M;

        double[] rs = new double[N];
        double[] fs = new double[N];
        int valid = 0;
        for (int i = 0; i < N; i++) {
            double r = rLoStart + (rHi - rLoStart) * i / (N - 1);
            double f;
            try {
                f = residualPhotonOrbit(r, prograde);
            } catch (IllegalStateException ex) {
                continue;
            }
            if (Double.isFinite(f)) {
                rs[valid] = r;
                fs[valid] = f;
                valid++;
            }
        }
        if (valid < 2) {
            throw new IllegalStateException(
                    "photonOrbitRadius: insufficient valid grid points (a = "
                    + a + ", eps3 = " + eps3 + ", prograde = " + prograde + ")");
        }

        int bestI = -1;
        for (int i = 0; i < valid - 1; i++) {
            if (fs[i] * fs[i + 1] < 0.0) {
                bestI = i;
            }
        }
        if (bestI < 0) {
            throw new IllegalStateException(
                    "photonOrbitRadius: no F sign change in [" + rLoStart + ", "
                    + rHi + "] (a = " + a + ", eps3 = " + eps3
                    + ", prograde = " + prograde
                    + ") — Kerr-continuation orbit may not exist; see "
                    + "docs/jp-parameter-space-notes.md");
        }

        double rLo = rs[bestI];
        double rHi2 = rs[bestI + 1];
        double fLo = fs[bestI];
        for (int iter = 0; iter < 100; iter++) {
            double rMid = 0.5 * (rLo + rHi2);
            double fMid = residualPhotonOrbit(rMid, prograde);
            if (Math.signum(fMid) == Math.signum(fLo)) {
                rLo = rMid; fLo = fMid;
            } else {
                rHi2 = rMid;
            }
            if (rHi2 - rLo < 1e-12 * Math.abs(rMid)) return rMid;
        }
        return 0.5 * (rLo + rHi2);
    }

    // ------------------------------------------------------------------
    // Test-visibility hooks: package-private equatorial-derivative paths
    // ------------------------------------------------------------------

    /**
     * Returns {@code {g_tt, g_tφ, g_φφ, ∂_r g_tt, ∂_r g_tφ, ∂_r g_φφ}}
     * on the equator (θ = π/2) computed via the same expressions used
     * inside {@link #eCircular} (the ISCO-helper code path).
     *
     * <p>Package-private; exposed solely so the test
     * {@code equatorialDerivativesAgreeAcrossPaths} can verify that the
     * ISCO and photon-orbit code paths agree at machine precision. See
     * {@code docs/phase-3a-gates-5-6.md} §3.4 amendment 4. Production
     * code paths ({@link #eCircular}, {@link #residualPhotonOrbit}) call
     * their own internal inline expressions, NOT this helper.
     */
    double[] equatorialGDerivativesViaIscoPath(double r) {
        // Mirrors the inline expressions inside eCircular line-for-line.
        double h = eps3 * M * M * M / (r * r * r);
        double hPrime = -3.0 * eps3 * M * M * M / (r * r * r * r);
        double onePlusH = 1.0 + h;
        double twoMOverR = 2.0 * M / r;
        double twoMaOverR = 2.0 * M * a / r;
        double gtt = -(1.0 - twoMOverR) * onePlusH;
        double gtp = -twoMaOverR * onePlusH;
        double gpp = (r * r + a * a) + (2.0 * M * a * a / r) * onePlusH;
        double dgtt = -(2.0 * M / (r * r)) * onePlusH - (1.0 - twoMOverR) * hPrime;
        double dgtp = (2.0 * M * a / (r * r)) * onePlusH - twoMaOverR * hPrime;
        double dgpp = 2.0 * r + 2.0 * M * a * a * (hPrime / r - onePlusH / (r * r));
        return new double[] { gtt, gtp, gpp, dgtt, dgtp, dgpp };
    }

    /**
     * Returns {@code {g_tt, g_tφ, g_φφ, ∂_r g_tt, ∂_r g_tφ, ∂_r g_φφ}}
     * on the equator (θ = π/2) computed via the same expressions used
     * inside {@link #residualPhotonOrbit} (the photon-orbit code path).
     *
     * <p>Package-private; same role as
     * {@link #equatorialGDerivativesViaIscoPath} but for the photon-
     * orbit code path. Used by
     * {@code equatorialDerivativesAgreeAcrossPaths} for the cross-path
     * comparison.
     */
    double[] equatorialGDerivativesViaPhotonPath(double r) {
        // Mirrors the inline expressions inside residualPhotonOrbit.
        final double r2 = r * r;
        final double r3 = r2 * r;
        final double r4 = r3 * r;
        final double M3 = M * M * M;
        final double Ma = M * a;
        final double Ma2 = Ma * a;
        final double e = eps3;

        final double oneMinusTwoMOverR = 1.0 - 2.0 * M / r;
        final double twoMOverR2 = 2.0 * M / r2;
        final double minus3eM3OverR4 = -3.0 * e * M3 / r4;
        final double eM3OverR3 = e * M3 / r3;
        final double onePlusH_ = 1.0 + eM3OverR3;

        final double gtt_ = -oneMinusTwoMOverR * onePlusH_;
        final double gtp_ = -(2.0 * Ma / r) * onePlusH_;
        final double gpp_ = (r2 + a * a) + (2.0 * Ma2 / r) * onePlusH_;

        final double dgtt_ = -twoMOverR2 * onePlusH_ - oneMinusTwoMOverR * minus3eM3OverR4;
        final double dgtp_ = (2.0 * Ma / r2) * onePlusH_ - (2.0 * Ma / r) * minus3eM3OverR4;
        final double dgpp_ = 2.0 * r + 2.0 * Ma2 * (minus3eM3OverR4 / r - onePlusH_ / r2);

        return new double[] { gtt_, gtp_, gpp_, dgtt_, dgtp_, dgpp_ };
    }

    // ------------------------------------------------------------------
    // Metric tensor and inverse
    // ------------------------------------------------------------------

    @Override
    public double[][] g(double[] x) {
        double r = x[1];
        double th = x[2];
        double s = Math.sin(th);
        double c = Math.cos(th);
        double s2 = s * s;
        double c2 = c * c;

        double r2 = r * r;
        double a2 = a * a;
        double Sigma = r2 + a2 * c2;
        double Delta = r2 - 2.0 * M * r + a2;
        double h = eps3 * M * M * M * r / (Sigma * Sigma);
        double onePlusH = 1.0 + h;

        double[][] gmn = new double[4][4];
        gmn[0][0] = -(1.0 - 2.0 * M * r / Sigma) * onePlusH;
        gmn[1][1] = Sigma * onePlusH / (Delta + a2 * h * s2);
        gmn[2][2] = Sigma;
        gmn[3][3] = s2 * ((r2 + a2) + (2.0 * M * a2 * r * s2 / Sigma) * onePlusH);
        double gtphi = -(2.0 * M * a * r * s2 / Sigma) * onePlusH;
        gmn[0][3] = gtphi;
        gmn[3][0] = gtphi;
        return gmn;
    }

    @Override
    public double[][] gInv(double[] x) {
        double r = x[1];
        double th = x[2];
        double s = Math.sin(th);
        double c = Math.cos(th);
        double s2 = s * s;
        double c2 = c * c;

        double r2 = r * r;
        double a2 = a * a;
        double Sigma = r2 + a2 * c2;
        double Delta = r2 - 2.0 * M * r + a2;
        double h = eps3 * M * M * M * r / (Sigma * Sigma);
        double onePlusH = 1.0 + h;

        double gtt = -(1.0 - 2.0 * M * r / Sigma) * onePlusH;
        double gtp = -(2.0 * M * a * r * s2 / Sigma) * onePlusH;
        double gpp = s2 * ((r2 + a2) + (2.0 * M * a2 * r * s2 / Sigma) * onePlusH);
        double det = gtt * gpp - gtp * gtp;

        double[][] gi = new double[4][4];
        gi[0][0] =  gpp / det;
        gi[1][1] =  (Delta + a2 * h * s2) / (Sigma * onePlusH);
        gi[2][2] =  1.0 / Sigma;
        gi[3][3] =  gtt / det;
        double gitp = -gtp / det;
        gi[0][3] = gitp;
        gi[3][0] = gitp;
        return gi;
    }

    // ------------------------------------------------------------------
    // Christoffel symbols (Sympy-generated CSE block)
    // ------------------------------------------------------------------

    /**
     * Fill a 20-element scratch array with the non-zero Christoffels of the
     * JP metric in the order defined by the {@code G_*} index constants.
     * The body is pasted verbatim from
     * {@code scripts/derive_jp_christoffels.py} (CSE-optimized output); do
     * not hand-edit — rerun the script instead. Kerr reduction at
     * {@code ε₃ = 0} is exact to machine precision (verified in the script
     * and in the test suite).
     *
     * <p>Complexity: O(1), ~88 arithmetic ops per call (68 CSE temporaries
     * + 20 Christoffels). No heap allocation.
     */
    private void fillChristoffelScratch(double r, double theta, double[] out) {
        final double M = this.M;
        final double a = this.a;
        final double eps3 = this.eps3;

        // --- generated by scripts/derive_jp_christoffels.py ---
        final double x0 = Math.pow(r, 2);
        final double x1 = Math.pow(a, 2);
        final double x2 = Math.cos(theta);
        final double x3 = x0 + x1 * Math.pow(x2, 2);
        final double x4 = Math.pow(x3, -2);
        final double x5 = Math.pow(M, 3) * eps3;
        final double x6 = x4 * x5;
        final double x7 = r * x6;
        final double x8 = -x7;
        final double x9 = x7 + 1;
        final double x10 = 1.0 / x3;
        final double x11 = x0 * x10;
        final double x12 = 2 * x11;
        final double x13 = x12 * x9;
        final double x14 = 4 * x11 - 1;
        final double x15 = x13 + x14 * x7 + x8 - 1;
        final double x16 = -x15;
        final double x17 = 2 * r;
        final double x18 = M * x17;
        final double x19 = x10 * x18;
        final double x20 = Math.sin(theta);
        final double x21 = Math.pow(x20, 2);
        final double x22 = x1 * x21;
        final double x23 = x19 * x22;
        final double x24 = x22 * x9;
        final double x25 = x0 + x1;
        final double x26 = x19 * x24 + x25;
        final double x27 = 1.0 / x9;
        final double x28 = x26 * x27;
        final double x29 = Math.pow(M, 2);
        final double x30 = x19 - 1;
        final double x31 = eps3 * x10 * x29 * x30;
        final double x32 = x14 * x31 + 2 * x9 * (x12 - 1);
        final double x33 = (1.0 / 2.0) * x32;
        final double x34 = 4 * x0;
        final double x35 = x29 * x4;
        final double x36 = x34 * x35;
        final double x37 = 1.0 / (x24 * x36 - x26 * x30);
        final double x38 = M * x10 * x37;
        final double x39 = x1 * x10;
        final double x40 = x21 * x39;
        final double x41 = Math.pow(x3, -3);
        final double x42 = x22 * x41;
        final double x43 = x42 * x5;
        final double x44 = x17 * x43 + x40 * x9 + x9;
        final double x45 = x18 * x44;
        final double x46 = x31 + x9;
        final double x47 = x2 * x37;
        final double x48 = x1 * x20;
        final double x49 = x4 * x48;
        final double x50 = 2 * M;
        final double x51 = Math.pow(M, 4) * eps3 * r * x14 * x42 - M * x1 * x10 * x21 * x9 - r + x0 * x24 * x4 * x50;
        final double x52 = -x13 - x14 * x7 + x9;
        final double x53 = a * x38;
        final double x54 = x23 * x44 + x26;
        final double x55 = a * x19 * x47;
        final double x56 = -x18 + x22 * x7 + x25;
        final double x57 = x27 * x56;
        final double x58 = M * x4 * x57;
        final double x59 = 1.0 / x56;
        final double x60 = x59 * x9;
        final double x61 = x10 * x27;
        final double x62 = r * x10;
        final double x63 = x5 * x60 * x62 * (2 * x40 + 1) + x8 + 1;
        final double x64 = x2 * x20;
        final double x65 = x39 * x64;
        final double x66 = x18 * x46;
        final double x67 = 1.0 / x20;

        out[G_T_TR]    = x38 * (x16 * x23 + x28 * x33);
        out[G_T_TTH]   = x18 * x47 * x49 * (-x28 * x46 + x45);
        out[G_T_RPH]   = x21 * x53 * (x17 * x51 + x26 * x27 * x52);
        out[G_T_THPH]  = x20 * x55 * (x26 * x27 * x44 - x54);
        out[G_R_TT]    = x33 * x58;
        out[G_R_TPH]   = -a * x15 * x21 * x58;
        out[G_R_RR]    = (1.0 / 2.0) * x61 * (-x10 * x14 * x5 + x17 * x9 + x3 * x60 * (-x17 - x22 * x6 + x34 * x43 + x50));
        out[G_R_RTH]   = -x27 * x63 * x65;
        out[G_R_THTH]  = -x57 * x62;
        out[G_R_PHPH]  = x21 * x51 * x56 * x61;
        out[G_TH_TT]   = -x1 * x41 * x64 * x66;
        out[G_TH_TPH]  = a * x4 * x45 * x64;
        out[G_TH_RR]   = x59 * x63 * x65;
        out[G_TH_RTH]  = x62;
        out[G_TH_THTH] = -x65;
        out[G_TH_PHPH] = -x10 * x54 * x64;
        out[G_PH_TR]   = x53 * (M * x32 * x62 + x30 * x52);
        out[G_PH_TTH]  = x55 * (x30 * x44 * x67 - x49 * x66);
        out[G_PH_RPH]  = x37 * (x16 * x17 * x22 * x35 + x30 * x51);
        out[G_PH_THPH] = x47 * (-x30 * x54 * x67 + x36 * x44 * x48);
        // --- end generated ---
    }

    @Override
    public double[][][] christoffel(double[] x) {
        double[] G = new double[CHRISTOFFEL_SCRATCH_SIZE];
        fillChristoffelScratch(x[1], x[2], G);

        double[][][] tensor = new double[4][4][4];

        // α = t
        tensor[0][0][1] = G[G_T_TR];   tensor[0][1][0] = G[G_T_TR];
        tensor[0][0][2] = G[G_T_TTH];  tensor[0][2][0] = G[G_T_TTH];
        tensor[0][1][3] = G[G_T_RPH];  tensor[0][3][1] = G[G_T_RPH];
        tensor[0][2][3] = G[G_T_THPH]; tensor[0][3][2] = G[G_T_THPH];

        // α = r
        tensor[1][0][0] = G[G_R_TT];
        tensor[1][0][3] = G[G_R_TPH];  tensor[1][3][0] = G[G_R_TPH];
        tensor[1][1][1] = G[G_R_RR];
        tensor[1][1][2] = G[G_R_RTH];  tensor[1][2][1] = G[G_R_RTH];
        tensor[1][2][2] = G[G_R_THTH];
        tensor[1][3][3] = G[G_R_PHPH];

        // α = θ
        tensor[2][0][0] = G[G_TH_TT];
        tensor[2][0][3] = G[G_TH_TPH]; tensor[2][3][0] = G[G_TH_TPH];
        tensor[2][1][1] = G[G_TH_RR];
        tensor[2][1][2] = G[G_TH_RTH]; tensor[2][2][1] = G[G_TH_RTH];
        tensor[2][2][2] = G[G_TH_THTH];
        tensor[2][3][3] = G[G_TH_PHPH];

        // α = φ
        tensor[3][0][1] = G[G_PH_TR];   tensor[3][1][0] = G[G_PH_TR];
        tensor[3][0][2] = G[G_PH_TTH];  tensor[3][2][0] = G[G_PH_TTH];
        tensor[3][1][3] = G[G_PH_RPH];  tensor[3][3][1] = G[G_PH_RPH];
        tensor[3][2][3] = G[G_PH_THPH]; tensor[3][3][2] = G[G_PH_THPH];

        return tensor;
    }

    // ------------------------------------------------------------------
    // Optimized geodesic acceleration
    // ------------------------------------------------------------------

    /**
     * Non-allocating geodesic acceleration. Uses the same Sympy-generated
     * Christoffel block as {@link #christoffel} via
     * {@link #fillChristoffelScratch}, then fuses the bilinear contraction
     * {@code a^μ = -Γ^μ_{αβ} k^α k^β} into the output array directly. No
     * 4×4×4 tensor is built; the 20-element scratch is local to the call
     * and escape-analysis eligible.
     *
     * <p>Structure (mirrors {@link KerrMetric#geodesicAcceleration}):
     * <ul>
     *   <li>{@code a^t} and {@code a^φ} are each a sum of four bilinear
     *       cross-terms (no diagonal Christoffels for those indices).</li>
     *   <li>{@code a^r} and {@code a^θ} each have four diagonal Christoffels
     *       plus (r,θ) and (t,φ) cross terms.</li>
     * </ul>
     *
     * <p>Complexity: O(1), ~140 arithmetic ops per call (88 to fill scratch,
     * ~52 for the contraction). Thread-safe: no mutable state on the
     * instance beyond the final fields {@code M, a, eps3}.
     */
    @Override
    public void geodesicAcceleration(double[] x, double[] k, double[] out) {
        double[] G = new double[CHRISTOFFEL_SCRATCH_SIZE];
        fillChristoffelScratch(x[1], x[2], G);

        double kt = k[0];
        double kr = k[1];
        double kh = k[2];
        double kf = k[3];

        out[0] = -2.0 * (G[G_T_TR]   * kt * kr
                      + G[G_T_TTH]  * kt * kh
                      + G[G_T_RPH]  * kr * kf
                      + G[G_T_THPH] * kh * kf);

        out[1] = -(G[G_R_TT]   * kt * kt
                  + G[G_R_RR]   * kr * kr
                  + G[G_R_THTH] * kh * kh
                  + G[G_R_PHPH] * kf * kf
                  + 2.0 * G[G_R_RTH] * kr * kh
                  + 2.0 * G[G_R_TPH] * kt * kf);

        out[2] = -(G[G_TH_TT]   * kt * kt
                  + G[G_TH_RR]   * kr * kr
                  + G[G_TH_THTH] * kh * kh
                  + G[G_TH_PHPH] * kf * kf
                  + 2.0 * G[G_TH_RTH] * kr * kh
                  + 2.0 * G[G_TH_TPH] * kt * kf);

        out[3] = -2.0 * (G[G_PH_TR]   * kt * kr
                      + G[G_PH_TTH]  * kt * kh
                      + G[G_PH_RPH]  * kr * kf
                      + G[G_PH_THPH] * kh * kf);
    }
}
