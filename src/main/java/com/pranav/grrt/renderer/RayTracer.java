package com.pranav.grrt.renderer;

/**
 * Strategy for integrating one ray from its initial 8-state to termination.
 *
 * <p>Implementations own the integrator and scratch buffers. The
 * {@link Renderer} holds one {@link RayTracer} per worker thread via
 * {@link ThreadLocal}, so implementations do not need to be thread-safe.
 *
 * <p>{@link FixedStepRayTracer} wraps a fixed-step
 * {@link com.pranav.grrt.integrator.Integrator} with the two-zone step
 * policy from {@link RenderConfig}. {@link AdaptiveRayTracer} wraps a
 * {@link com.pranav.grrt.integrator.AdaptiveIntegrator} with per-call
 * tolerances and an axis-reflection fix for Boyer-Lindquist's polar
 * coordinate singularity.
 */
public interface RayTracer {

    /**
     * Trace a ray from its initial 8-state until termination. The caller
     * passes the initial state in {@code y}; the ray tracer updates
     * {@code y} in place to hold the terminal state. Termination
     * conditions (capture, escape, NaN, step budget) are defined by the
     * {@link RenderConfig} passed at construction.
     *
     * @param y 8-state on entry and exit; overwritten in place
     * @return the terminal outcome of the ray
     */
    RayOutcome trace(double[] y);
}
