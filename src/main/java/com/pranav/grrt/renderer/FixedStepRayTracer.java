package com.pranav.grrt.renderer;

import com.pranav.grrt.integrator.Integrator;
import com.pranav.grrt.metric.Metric;

/**
 * Fixed-step ray tracer with the Phase 1 two-zone step policy.
 *
 * <p>Uses {@link RenderConfig#hNear()} when r &le; {@link RenderConfig#rTransition()}
 * and {@link RenderConfig#hFar()} otherwise. This is the same loop as the
 * Phase 1 {@code Renderer.traceRay} body, extracted for the strategy split.
 *
 * <p>Not thread-safe: holds a scratch buffer and wraps a fixed-step
 * {@link Integrator} which itself is not thread-safe.
 */
public final class FixedStepRayTracer implements RayTracer {

    private final Integrator integrator;
    private final Metric metric;
    private final double[] next = new double[8];

    private final double rCapture;
    private final double rEscape;
    private final double rTransition;
    private final double hNear;
    private final double hFar;
    private final int maxSteps;

    public FixedStepRayTracer(Integrator integrator, Metric metric, RenderConfig config) {
        if (integrator == null || metric == null || config == null) {
            throw new IllegalArgumentException("null argument");
        }
        this.integrator = integrator;
        this.metric = metric;
        this.rCapture = metric.horizonRadius() + config.horizonCushion();
        this.rEscape = config.escapeRadius();
        this.rTransition = config.rTransition();
        this.hNear = config.hNear();
        this.hFar  = config.hFar();
        this.maxSteps = config.maxSteps();
    }

    @Override
    public RayOutcome trace(double[] y) {
        for (int step = 0; step < maxSteps; step++) {
            double r = y[1];
            double h = (r > rTransition) ? hFar : hNear;
            integrator.step(metric, y, h, next);
            System.arraycopy(next, 0, y, 0, 8);

            double rNow = y[1];
            if (Double.isNaN(rNow) || Double.isNaN(y[4])) return RayOutcome.NAN;
            if (rNow < rCapture) return RayOutcome.CAPTURED;
            if (rNow > rEscape)  return RayOutcome.ESCAPED;
        }
        return RayOutcome.MAX_STEPS;
    }
}
