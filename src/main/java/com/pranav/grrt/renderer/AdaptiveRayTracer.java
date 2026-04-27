package com.pranav.grrt.renderer;

import com.pranav.grrt.disk.Disk;
import com.pranav.grrt.integrator.AdaptiveIntegrator;
import com.pranav.grrt.integrator.DormandPrince45;
import com.pranav.grrt.metric.Metric;

/**
 * Adaptive ray tracer wrapping an {@link AdaptiveIntegrator} with
 * per-call tolerances, polar-axis reflection, optional disk-crossing
 * detection, and position-dependent horizon termination.
 *
 * <p>Per step:
 * <ol>
 *   <li>Propose a step of size {@code h}; get back a {@link AdaptiveIntegrator.StepStatus}.</li>
 *   <li>If accepted, and a {@link Disk} is configured, test for an
 *       equator-crossing inside the disk extent. If hit, bisect the
 *       step via {@link DormandPrince45#interpolate(double, double[])}
 *       to find the crossing point and return
 *       {@link RayOutcome#HIT_DISK}.</li>
 *   <li>Otherwise apply the polar-axis reflection (see
 *       {@link #reflectAcrossAxisIfNeeded}) and commit the (possibly
 *       reflected) state to {@code y}.</li>
 *   <li>Test horizon (via {@link Metric#isInsideHorizon}), escape, NaN.</li>
 *   <li>Update {@code h} from {@code status.hNext()} for the next attempt
 *       regardless of accept/reject.</li>
 * </ol>
 *
 * <p>The disk-crossing test runs <i>before</i> axis reflection because:
 * (a) reflection invalidates the integrator's dense-output buffers via
 * {@link AdaptiveIntegrator#resetState}; and (b) crossing the equator
 * (θ = π/2) and crossing the polar axis (θ = 0 or π) are mutually
 * exclusive in a single accepted step at any sane step size.
 *
 * <p>Tracks per-instance accepted and rejected step counts; retrieve via
 * {@link #acceptedSteps()} / {@link #rejectedSteps()}. Since the tracer
 * is thread-local, these are local counts — aggregate across worker
 * threads externally.
 *
 * <p>Not thread-safe: holds scratch buffers, step counters, and wraps a
 * non-thread-safe {@link AdaptiveIntegrator}.
 */
public final class AdaptiveRayTracer implements RayTracer {

    private static final double TWO_PI = 2.0 * Math.PI;
    private static final double HALF_PI = 0.5 * Math.PI;

    /** Bisection tolerance on |θ − π/2| at the disk crossing. */
    private static final double EQUATOR_TOL = 1.0e-9;
    /** Bisection iteration cap (typical convergence: ~30 iters). */
    private static final int    BISECT_MAX_ITERS = 50;

    private final AdaptiveIntegrator integrator;
    private final DormandPrince45 dp45;   // null iff disk == null
    private final Metric metric;
    private final Disk disk;              // null for binary-shader / shadow-only renders

    private final double[] next = new double[8];
    private final double[] hitScratch = new double[8];

    private final double horizonCushion;
    private final double rEscape;
    private final double atol;
    private final double rtol;
    private final double hInitial;
    private final int maxSteps;

    private long acceptedSteps = 0;
    private long rejectedSteps = 0;

    /**
     * Phase-2-compatible constructor (no disk). Equivalent to
     * {@link #AdaptiveRayTracer(AdaptiveIntegrator, Metric, RenderConfig, double, double, double, Disk)}
     * with {@code disk = null}.
     *
     * @param integrator adaptive integrator (e.g. {@link DormandPrince45})
     * @param metric     spacetime metric used by the integrator
     * @param config     termination + shader config; only
     *                   {@code horizonCushion}, {@code escapeRadius},
     *                   and {@code maxSteps} are consulted
     * @param atol       per-call absolute tolerance, must be &gt; 0
     * @param rtol       per-call relative tolerance, must be &gt; 0
     * @param hInitial   first proposed step size (M units); must be ≠ 0
     */
    public AdaptiveRayTracer(AdaptiveIntegrator integrator, Metric metric, RenderConfig config,
                             double atol, double rtol, double hInitial) {
        this(integrator, metric, config, atol, rtol, hInitial, null);
    }

    /**
     * Disk-aware constructor (Phase 3B.2). Disk-crossing detection is
     * enabled when {@code disk != null}; this requires the integrator
     * to be a {@link DormandPrince45} (only DP45 currently exposes the
     * dense-output {@code interpolate} needed for bisection).
     *
     * @param disk optional disk; pass {@code null} for shadow-only renders
     * @throws IllegalArgumentException if {@code disk != null} and
     *                                  {@code integrator} is not
     *                                  {@link DormandPrince45}
     */
    public AdaptiveRayTracer(AdaptiveIntegrator integrator, Metric metric, RenderConfig config,
                             double atol, double rtol, double hInitial, Disk disk) {
        if (integrator == null || metric == null || config == null) {
            throw new IllegalArgumentException("null argument");
        }
        if (!(atol > 0.0) || !(rtol > 0.0)) {
            throw new IllegalArgumentException(
                    "atol, rtol must be positive: atol=" + atol + ", rtol=" + rtol);
        }
        if (hInitial == 0.0 || !Double.isFinite(hInitial)) {
            throw new IllegalArgumentException(
                    "hInitial must be finite and nonzero: " + hInitial);
        }
        if (disk != null && !(integrator instanceof DormandPrince45)) {
            throw new IllegalArgumentException(
                    "disk-aware tracing requires DormandPrince45 integrator;"
                            + " got " + integrator.getClass().getSimpleName());
        }
        this.integrator = integrator;
        this.dp45 = (integrator instanceof DormandPrince45 dp) ? dp : null;
        this.metric = metric;
        this.disk = disk;
        this.horizonCushion = config.horizonCushion();
        this.rEscape = config.escapeRadius();
        this.atol = atol;
        this.rtol = rtol;
        this.hInitial = hInitial;
        this.maxSteps = config.maxSteps();
    }

    @Override
    public RayOutcome trace(double[] y) {
        double h = hInitial;
        integrator.resetState();   // fresh ray — discard any prior FSAL k7
        for (int step = 0; step < maxSteps; step++) {
            AdaptiveIntegrator.StepStatus s =
                    integrator.adaptiveStep(metric, y, h, atol, rtol, next);
            if (s.accepted()) {
                acceptedSteps++;

                // ---- Disk-crossing detection (3B.2) ----
                // Test BEFORE axis reflection. Equator (θ=π/2) and pole
                // (θ ∈ {0, π}) are mutually exclusive in a single step;
                // reflection's resetState would invalidate dense output.
                if (disk != null && disk.crossedEquator(y, next)) {
                    if (findDiskCrossingHit(y[2], next[2])) {
                        System.arraycopy(hitScratch, 0, y, 0, 8);
                        return RayOutcome.HIT_DISK;
                    }
                    // Crossing landed outside disk radial extent; fall
                    // through to the standard accept-step processing.
                }

                boolean reflected = reflectAcrossAxisIfNeeded(next);
                System.arraycopy(next, 0, y, 0, 8);
                if (reflected) integrator.resetState();

                if (Double.isNaN(y[1]) || Double.isNaN(y[4])) return RayOutcome.NAN;
                if (metric.isInsideHorizon(y, horizonCushion)) return RayOutcome.CAPTURED;
                if (y[1] > rEscape) return RayOutcome.ESCAPED;
            } else {
                rejectedSteps++;
            }
            h = s.hNext();
        }
        return RayOutcome.MAX_STEPS;
    }

    /**
     * Bisect θ ∈ [0, 1] over the most recently accepted step to find
     * the equatorial crossing, using
     * {@link DormandPrince45#interpolate(double, double[])}. Writes the
     * bisected 8-state into {@link #hitScratch} (allocation-free).
     *
     * @param thetaPrev value of x[2] at the start of the step (θ-coord)
     * @param thetaCurr value of x[2] at the end of the step (unused
     *                  past the sign-check already done by
     *                  {@link Disk#crossedEquator})
     * @return true iff the bisected crossing radius lies in
     *         {@code [disk.rIsco(), disk.rOuter()]}
     */
    private boolean findDiskCrossingHit(double thetaPrev, double thetaCurr) {
        double thetaLo = 0.0, thetaHi = 1.0;
        double fLo = thetaPrev - HALF_PI;

        for (int iter = 0; iter < BISECT_MAX_ITERS; iter++) {
            double thetaMid = 0.5 * (thetaLo + thetaHi);
            dp45.interpolate(thetaMid, hitScratch);
            double fMid = hitScratch[2] - HALF_PI;
            if (Math.abs(fMid) < EQUATOR_TOL) {
                break;
            }
            if (fLo * fMid < 0.0) {
                thetaHi = thetaMid;
            } else {
                thetaLo = thetaMid;
                fLo = fMid;
            }
        }

        double rEmit = hitScratch[1];
        return rEmit >= disk.rIsco() && rEmit <= disk.rOuter();
    }

    /**
     * Coordinate-singularity fix for Boyer-Lindquist polar axis. Not
     * physics — a ray that crosses the axis at θ = 0 or π sees its
     * (θ, φ) coordinates fold back; we mirror θ and rotate φ by π,
     * flipping k^θ to preserve the geodesic. Invariants (null norm,
     * Killing charges) are unchanged.
     *
     * @return true iff a reflection was applied (caller should reset
     *         any integrator state tied to the pre-reflection point)
     */
    private static boolean reflectAcrossAxisIfNeeded(double[] s) {
        double theta = s[2];
        boolean reflected = false;
        if (theta < 0.0) {
            s[2] = -theta;
            s[3] = s[3] + Math.PI;
            s[6] = -s[6];
            reflected = true;
        } else if (theta > Math.PI) {
            s[2] = 2.0 * Math.PI - theta;
            s[3] = s[3] + Math.PI;
            s[6] = -s[6];
            reflected = true;
        }
        // Normalize φ to [0, 2π) unconditionally — cheap and harmless.
        double phi = s[3];
        s[3] = ((phi % TWO_PI) + TWO_PI) % TWO_PI;
        return reflected;
    }

    public long acceptedSteps() { return acceptedSteps; }
    public long rejectedSteps() { return rejectedSteps; }
}
