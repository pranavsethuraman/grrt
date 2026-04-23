package com.pranav.grrt.renderer;

import com.pranav.grrt.integrator.AdaptiveIntegrator;
import com.pranav.grrt.metric.Metric;

/**
 * Adaptive ray tracer wrapping an {@link AdaptiveIntegrator} with
 * per-call tolerances and a polar-axis reflection fix.
 *
 * <p>Per step:
 * <ol>
 *   <li>Propose a step of size {@code h}; get back a {@link AdaptiveIntegrator.StepStatus}.</li>
 *   <li>If accepted, apply the axis reflection below to the proposed
 *       state (coordinate-singularity fix; see {@link #reflectAcrossAxisIfNeeded}).</li>
 *   <li>Commit the reflected state to {@code y}, check termination.</li>
 *   <li>Update {@code h} from {@code status.hNext()} for the next attempt
 *       regardless of accept/reject.</li>
 * </ol>
 *
 * <p>On acceptance with a reflection, {@link AdaptiveIntegrator#resetState()}
 * is called to invalidate any FSAL buffer the integrator may be holding
 * at the un-reflected state.
 *
 * <p>Tracks per-instance accepted and rejected step counts; retrieve via
 * {@link #acceptedSteps()} / {@link #rejectedSteps()}. Since the tracer
 * is thread-local, these are local counts — aggregate across worker
 * threads externally.
 *
 * <p>Not thread-safe: holds scratch buffer, step counters, and wraps a
 * non-thread-safe {@link AdaptiveIntegrator}.
 */
public final class AdaptiveRayTracer implements RayTracer {

    private static final double TWO_PI = 2.0 * Math.PI;

    private final AdaptiveIntegrator integrator;
    private final Metric metric;
    private final double[] next = new double[8];

    private final double rCapture;
    private final double rEscape;
    private final double atol;
    private final double rtol;
    private final double hInitial;
    private final int maxSteps;

    private long acceptedSteps = 0;
    private long rejectedSteps = 0;

    /**
     * @param integrator adaptive integrator (e.g. {@link com.pranav.grrt.integrator.DormandPrince45})
     * @param metric     spacetime metric used by the integrator
     * @param config     termination + shader config; only
     *                   {@code horizonCushion}, {@code escapeRadius},
     *                   and {@code maxSteps} are consulted — the
     *                   fixed-step zone fields are ignored
     * @param atol       per-call absolute tolerance passed to the
     *                   integrator; must be &gt; 0
     * @param rtol       per-call relative tolerance; must be &gt; 0
     * @param hInitial   first proposed step size (M units); must be &ne; 0
     */
    public AdaptiveRayTracer(AdaptiveIntegrator integrator, Metric metric, RenderConfig config,
                             double atol, double rtol, double hInitial) {
        if (integrator == null || metric == null || config == null) {
            throw new IllegalArgumentException("null argument");
        }
        if (!(atol > 0.0) || !(rtol > 0.0)) {
            throw new IllegalArgumentException(
                    "atol, rtol must be positive: atol=" + atol + ", rtol=" + rtol);
        }
        if (hInitial == 0.0 || !Double.isFinite(hInitial)) {
            throw new IllegalArgumentException("hInitial must be finite and nonzero: " + hInitial);
        }
        this.integrator = integrator;
        this.metric = metric;
        this.rCapture = metric.horizonRadius() + config.horizonCushion();
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
                boolean reflected = reflectAcrossAxisIfNeeded(next);
                System.arraycopy(next, 0, y, 0, 8);
                if (reflected) integrator.resetState();

                double rNow = y[1];
                if (Double.isNaN(rNow) || Double.isNaN(y[4])) return RayOutcome.NAN;
                if (rNow < rCapture) return RayOutcome.CAPTURED;
                if (rNow > rEscape)  return RayOutcome.ESCAPED;
            } else {
                rejectedSteps++;
            }
            h = s.hNext();
        }
        return RayOutcome.MAX_STEPS;
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
