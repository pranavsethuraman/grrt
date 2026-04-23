package com.pranav.grrt.integrator;

import com.pranav.grrt.metric.Metric;

/**
 * Adaptive ODE integrator for null geodesics with per-step error control.
 *
 * <p>Sibling of {@link Integrator}, not an extension: adaptive stepping
 * returns accept/reject status plus a suggested next step size, which do
 * not fit the fixed-step {@code step(metric, y, h, out)} signature. The
 * {@link com.pranav.grrt.renderer.Renderer} selects one interface or the
 * other at construction; no attempt is made to unify them.
 *
 * <p>State vector layout matches {@link Integrator}:
 * <pre>
 *   y[0..3] = x^μ = (t, r, θ, φ)
 *   y[4..7] = k^μ = dx^μ/dλ
 * </pre>
 *
 * <p>Tolerances are supplied per call, not held on the instance; this
 * keeps the integrator instance minimal (scratch buffers plus any
 * controller state) and lets the caller adapt atol/rtol across an
 * integration (e.g. tightening near the photon sphere).
 *
 * <p>Typical loop:
 * <pre>
 *   AdaptiveIntegrator dp = new DormandPrince45();
 *   double h = 0.1;
 *   while (!done) {
 *       StepStatus s = dp.adaptiveStep(metric, y, h, atol, rtol, yNext);
 *       if (s.accepted()) System.arraycopy(yNext, 0, y, 0, 8);
 *       h = s.hNext();
 *   }
 * </pre>
 */
public interface AdaptiveIntegrator {

    /**
     * Outcome of one attempted adaptive step.
     *
     * @param accepted   true iff the step passed the per-component
     *                   error test (err &le; 1). When false, the caller
     *                   must NOT advance {@code y}; it should retry
     *                   with {@code hNext}.
     * @param hNext      suggested next step size (same sign as the input
     *                   {@code hTry}). On acceptance this is the
     *                   controller's growth/shrink proposal; on rejection
     *                   it is a shrink toward the error tolerance.
     * @param errorNorm  RMS error norm actually observed on this attempt;
     *                   useful for diagnostics and test assertions.
     */
    record StepStatus(boolean accepted, double hNext, double errorNorm) {}

    /**
     * Attempt one adaptive step of size {@code hTry}.
     *
     * <p>On acceptance, {@code out} holds the advanced state. On
     * rejection, {@code out} is not required to hold meaningful values
     * and callers must not use it.
     *
     * @param metric spacetime metric providing the geodesic RHS
     * @param y      current state, length 8; not modified
     * @param hTry   attempted affine-parameter step; may be negative
     * @param atol   absolute error tolerance, must be positive
     * @param rtol   relative error tolerance, must be positive
     * @param out    output state, length 8; overwritten only on
     *               acceptance
     * @return the step outcome
     */
    StepStatus adaptiveStep(Metric metric, double[] y, double hTry,
                            double atol, double rtol, double[] out);

    /**
     * Invalidate any controller memory (FSAL buffers, previous-error
     * term, etc.) so the next {@code adaptiveStep} call behaves as if
     * starting from scratch. Call this when the caller introduces a
     * discrete jump into {@code y} that did not come from the previously
     * accepted step — e.g. after an axis reflection in a renderer loop.
     * Default: no-op for stateless integrators.
     */
    default void resetState() {}
}
