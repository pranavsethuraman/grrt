package com.pranav.grrt.disk;

import com.pranav.grrt.metric.Metric;

/**
 * Novikov-Thorne 1973 thin accretion disk with Page-Thorne 1974 radial
 * flux profile, in the equatorial plane of a Kerr or Kerr-like
 * axisymmetric metric.
 *
 * <h2>Flux formula</h2>
 *
 * Page-Thorne 1974, eq. (15n); cf. Krolik 1999 §7.5:
 *
 * <pre>
 *   F(r, a) = -1 / (4π r √Δ)
 *             · (dΩ/dr) / (E - Ω L)²
 *             · ∫_{r_ISCO}^{r} (E - Ω L) (dL/dr') dr'
 * </pre>
 *
 * with Bardeen-Press-Teukolsky 1972 closed forms for the prograde
 * Kerr equatorial circular-orbit quantities (M-explicit form):
 *
 * <pre>
 *   Ω(r, a) = √M / (r^{3/2} + a √M)
 *   E(r, a) = (r² − 2Mr + a √(Mr)) / (r √(r² − 3Mr + 2a √(Mr)))
 *   L(r, a) = √(Mr) (r² − 2a √(Mr) + a²) / (r √(r² − 3Mr + 2a √(Mr)))
 *   Δ(r)   = r² − 2Mr + a²
 * </pre>
 *
 * Geometrized G = c = 1; the accretion rate Ṁ is normalized to 1 so F
 * is dimensionless.
 *
 * <h2>Reference values</h2>
 *
 * Validated against {@code scripts/page_thorne_reference.py}: 12
 * (a, r) test points at a ∈ {0, 0.9} and r ∈ [r_ISCO + ε, 20 M].
 * Java/Python agreement at 1e-6 or better is Phase 3B gate 1
 * ({@code NovikovThorneDiskTest.surfaceFluxMatchesPageThorneReference}).
 *
 * <h2>Spin parameter</h2>
 *
 * The Bardeen formulas require knowing the Kerr spin a. The current
 * {@link Metric} interface does not expose spin, so it is passed
 * explicitly to the constructor. For Kerr at spin {@code a}, the
 * caller supplies {@code a}. For Johannsen-Psaltis at (a, ε₃), the
 * Page-Thorne formula is treated as the dominant disk-emission
 * contribution, with deformation effects entering only via the host
 * {@link Metric} (used in {@link #keplerianFourVelocity} for the
 * redshift in the shader).
 *
 * <h2>Inner-integral computation</h2>
 *
 * Adaptive Simpson recursion with Richardson extrapolation, tolerance
 * 1e-12, max recursion depth 25. Smooth integrand on
 * [r_ISCO + 1e-12, r]; convergence is typically O(20–50) function
 * evaluations.
 *
 * <h2>Thread-safety</h2>
 *
 * All fields are {@code final}. Methods read state only. Safe for
 * concurrent reads.
 *
 * <h2>Time/space complexity</h2>
 *
 * <ul>
 *   <li>Constructor: O(1) plus one ISCO computation via the host metric.</li>
 *   <li>{@link #surfaceFlux}: O(N) integrand evaluations where N is
 *       the adaptive Simpson depth (~20–50 typical).</li>
 *   <li>{@link #temperature}: O(N) (delegates to {@link #surfaceFlux}).</li>
 *   <li>{@link #keplerianFourVelocity}: O(1) plus one
 *       {@link Metric#g} call.</li>
 *   <li>{@link #crossedEquator}: O(1).</li>
 * </ul>
 *
 * <p>For 3B.1 the surfaceFlux is recomputed per call. If 3B.2/3C
 * profiling shows this dominates frame time, a tabulated cache can
 * be added in the constructor without changing the public API.
 */
public final class NovikovThorneDisk implements Disk {

    /** Lower-bound bump to keep the integration off the ISCO endpoint. */
    private static final double EPS_ISCO = 1.0e-12;
    /** Adaptive Simpson global tolerance. */
    private static final double SIMPSON_TOL = 1.0e-12;
    /** Adaptive Simpson recursion depth ceiling. */
    private static final int SIMPSON_MAX_DEPTH = 25;
    /** σ_SB in geometrized units (normalized to 1). */
    private static final double STEFAN_BOLTZMANN = 1.0;
    private static final double FOUR_PI = 4.0 * Math.PI;

    private final Metric metric;
    private final double a;
    private final double mass;
    private final double rIscoBare;
    private final double rIsco;
    private final double rOuter;

    /**
     * @param metric host spacetime; used for {@link Metric#iscoRadius}
     *               and for {@link Metric#g} at the disk-hit position
     *               in {@link #keplerianFourVelocity}.
     * @param spin   Kerr spin parameter a, in {@code [0, mass)}.
     *               Prograde-only for the current Page-Thorne form.
     * @param rOuter outer disk radius (must exceed inner edge).
     * @param iscoCushion small radial offset added to {@code r_ISCO}
     *               so the disk's inner edge sits just outside the
     *               numerically delicate ISCO. Typical: 1e-3 M.
     * @throws IllegalArgumentException on invalid arguments
     */
    public NovikovThorneDisk(Metric metric, double spin,
                             double rOuter, double iscoCushion) {
        if (metric == null) {
            throw new IllegalArgumentException("metric must not be null");
        }
        double m = metric.mass();
        if (!(spin >= 0.0) || !(spin < m)) {
            throw new IllegalArgumentException(
                    "require 0 <= spin < mass; got spin=" + spin + ", mass=" + m);
        }
        if (!(iscoCushion >= 0.0) || !Double.isFinite(iscoCushion)) {
            throw new IllegalArgumentException(
                    "iscoCushion must be finite and >= 0; got " + iscoCushion);
        }
        if (!(rOuter > 0.0) || !Double.isFinite(rOuter)) {
            throw new IllegalArgumentException(
                    "rOuter must be finite and > 0; got " + rOuter);
        }
        this.metric = metric;
        this.a = spin;
        this.mass = m;
        this.rIscoBare = metric.iscoRadius(true);
        this.rIsco = this.rIscoBare + iscoCushion;
        if (rOuter <= this.rIsco) {
            throw new IllegalArgumentException(
                    "rOuter (" + rOuter + ") must exceed rIsco ("
                            + this.rIsco + ")");
        }
        this.rOuter = rOuter;
    }

    @Override
    public boolean crossedEquator(double[] xPrev, double[] xCurr) {
        if (xPrev == null || xCurr == null
                || xPrev.length != 8 || xCurr.length != 8) {
            throw new IllegalArgumentException("states must have length 8");
        }
        double dPrev = xPrev[2] - 0.5 * Math.PI;
        double dCurr = xCurr[2] - 0.5 * Math.PI;
        if (dPrev * dCurr > 0.0) return false;
        double rMin = Math.min(xPrev[1], xCurr[1]);
        double rMax = Math.max(xPrev[1], xCurr[1]);
        return rMax >= rIsco && rMin <= rOuter;
    }

    @Override public double rIsco() { return rIsco; }
    @Override public double rOuter() { return rOuter; }

    @Override
    public double surfaceFlux(double r) {
        if (!(r >= rIsco) || !(r <= rOuter)) {
            throw new IllegalArgumentException(
                    "r=" + r + " outside disk [rIsco=" + rIsco
                            + ", rOuter=" + rOuter + "]");
        }
        double integral = adaptiveSimpson(rIscoBare + EPS_ISCO, r);

        double delta = r * r - 2.0 * mass * r + a * a;
        if (!(delta > 0.0)) {
            throw new IllegalStateException(
                    "Δ <= 0 at r=" + r + ", a=" + a + ": " + delta);
        }
        double sqrtDelta = Math.sqrt(delta);
        double e = energyE(r, a, mass);
        double l = angularL(r, a, mass);
        double omega = angularVelocity(r, a, mass);
        double dOmega = dAngularVelocity(r, a, mass);
        double diff = e - omega * l;

        return -dOmega / (FOUR_PI * r * sqrtDelta * diff * diff) * integral;
    }

    @Override
    public double temperature(double r) {
        double f = surfaceFlux(r);
        if (f < 0.0) {
            throw new IllegalStateException(
                    "surfaceFlux < 0 at r=" + r + ": " + f);
        }
        return Math.pow(f / STEFAN_BOLTZMANN, 0.25);
    }

    @Override
    public double[] keplerianFourVelocity(double[] xOnDisk) {
        if (xOnDisk == null || xOnDisk.length != 4) {
            throw new IllegalArgumentException("position must have length 4");
        }
        double r = xOnDisk[1];
        if (!(r >= rIsco) || !(r <= rOuter)) {
            throw new IllegalArgumentException(
                    "r=" + r + " outside disk [rIsco=" + rIsco
                            + ", rOuter=" + rOuter + "]");
        }
        double omega = angularVelocity(r, a, mass);
        double[][] gmn = metric.g(xOnDisk);
        double gtt = gmn[0][0];
        double gtp = gmn[0][3];
        double gpp = gmn[3][3];
        double quad = gtt + 2.0 * gtp * omega + gpp * omega * omega;
        if (!(quad < 0.0)) {
            throw new IllegalStateException(
                    "non-timelike circular orbit at r=" + r
                            + " (g_tt + 2 g_tφ Ω + g_φφ Ω² = " + quad + ")");
        }
        double uT = 1.0 / Math.sqrt(-quad);
        double uPhi = omega * uT;
        return new double[] { uT, 0.0, 0.0, uPhi };
    }

    // ---- Page-Thorne integrand and adaptive Simpson ----------------------

    private double integrand(double rp) {
        double e = energyE(rp, a, mass);
        double l = angularL(rp, a, mass);
        double omega = angularVelocity(rp, a, mass);
        double dl = dAngularL(rp, a, mass);
        return (e - omega * l) * dl;
    }

    private double simpson(double aLo, double bHi) {
        double m = 0.5 * (aLo + bHi);
        return ((bHi - aLo) / 6.0)
                * (integrand(aLo) + 4.0 * integrand(m) + integrand(bHi));
    }

    private double adaptiveSimpson(double aLo, double bHi) {
        return adaptiveSimpsonHelper(aLo, bHi, simpson(aLo, bHi),
                                     SIMPSON_TOL, SIMPSON_MAX_DEPTH);
    }

    private double adaptiveSimpsonHelper(double aLo, double bHi, double whole,
                                         double tol, int depthLeft) {
        double m = 0.5 * (aLo + bHi);
        double left  = simpson(aLo, m);
        double right = simpson(m, bHi);
        double diff = left + right - whole;
        if (depthLeft <= 0 || Math.abs(diff) < 15.0 * tol) {
            return left + right + diff / 15.0;
        }
        return adaptiveSimpsonHelper(aLo, m, left,  0.5 * tol, depthLeft - 1)
             + adaptiveSimpsonHelper(m, bHi, right, 0.5 * tol, depthLeft - 1);
    }

    // ---- Bardeen-Press-Teukolsky 1972 closed forms (M-explicit) ----------

    /**
     * Ω(r, a, M) = √M / (r^{3/2} + a √M). Prograde branch.
     *
     * <p>Package-private for cross-test verification against
     * scripts/page_thorne_reference.py at known (r, a, M) points.
     */
    static double angularVelocity(double r, double a, double M) {
        double rootM = Math.sqrt(M);
        double rootR = Math.sqrt(r);
        return rootM / (r * rootR + a * rootM);
    }

    /**
     * dΩ/dr = -(3/2) √M √r / (r^{3/2} + a √M)².
     */
    static double dAngularVelocity(double r, double a, double M) {
        double rootM = Math.sqrt(M);
        double rootR = Math.sqrt(r);
        double denom = r * rootR + a * rootM;
        return -1.5 * rootM * rootR / (denom * denom);
    }

    /** E(r, a, M) per Bardeen 1972 eq. (2.13) prograde. */
    static double energyE(double r, double a, double M) {
        double rootMr = Math.sqrt(M * r);
        double inner = r * r - 3.0 * M * r + 2.0 * a * rootMr;
        if (!(inner > 0.0)) {
            throw new IllegalArgumentException(
                    "circular-orbit denominator nonpositive at r=" + r
                            + ", a=" + a + ", M=" + M + ": " + inner);
        }
        double num = r * r - 2.0 * M * r + a * rootMr;
        return num / (r * Math.sqrt(inner));
    }

    /** L(r, a, M) per Bardeen 1972 eq. (2.14) prograde. */
    static double angularL(double r, double a, double M) {
        double rootMr = Math.sqrt(M * r);
        double inner = r * r - 3.0 * M * r + 2.0 * a * rootMr;
        if (!(inner > 0.0)) {
            throw new IllegalArgumentException(
                    "circular-orbit denominator nonpositive at r=" + r
                            + ", a=" + a + ", M=" + M + ": " + inner);
        }
        double num = rootMr * (r * r - 2.0 * a * rootMr + a * a);
        return num / (r * Math.sqrt(inner));
    }

    /**
     * dL/dr derived via the quotient rule from L = N/D with
     * <pre>
     *   N(r) = √(Mr) (r² − 2a√(Mr) + a²)
     *   D(r) = r √(r² − 3Mr + 2a√(Mr))
     * </pre>
     * Verified at the Schwarzschild ISCO (a=0, r=6M): dL/dr = 0
     * exactly, since L has its minimum at the ISCO for a=0.
     */
    static double dAngularL(double r, double a, double M) {
        double rootM = Math.sqrt(M);
        double rootR = Math.sqrt(r);
        double rootMr = rootM * rootR;
        double inner = r * r - 3.0 * M * r + 2.0 * a * rootMr;
        if (!(inner > 0.0)) {
            throw new IllegalArgumentException(
                    "circular-orbit denominator nonpositive at r=" + r
                            + ", a=" + a + ", M=" + M + ": " + inner);
        }
        double sqrtInner = Math.sqrt(inner);

        double nPolyTerm = r * r - 2.0 * a * rootMr + a * a;
        double N = rootMr * nPolyTerm;
        // d√(Mr)/dr = √M / (2 √r); d(r² − 2a√(Mr) + a²)/dr = 2r − a√M / √r
        double dN = rootM / (2.0 * rootR) * nPolyTerm
                + rootMr * (2.0 * r - a * rootM / rootR);

        double D = r * sqrtInner;
        // d(inner)/dr = 2r − 3M + a√M / √r
        double dInner = 2.0 * r - 3.0 * M + a * rootM / rootR;
        double dD = sqrtInner + r * dInner / (2.0 * sqrtInner);

        return (dN * D - N * dD) / (D * D);
    }
}
