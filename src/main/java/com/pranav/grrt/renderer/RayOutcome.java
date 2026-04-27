package com.pranav.grrt.renderer;

/**
 * Terminal state of a traced null geodesic. Produced by {@link Renderer} and
 * consumed by {@link RayShader} to assign a pixel value.
 */
public enum RayOutcome {

    /** Ray reached r &lt; horizonRadius + horizonCushion. Captured by the BH. */
    CAPTURED,

    /** Ray reached r &gt; escapeRadius. Escaped to the celestial sphere. */
    ESCAPED,

    /**
     * Ray exceeded maxSteps without terminating. Practically only happens for
     * rays tangent to an unstable photon orbit; Phase 1 classifies these as
     * {@link #CAPTURED} at the shader level.
     */
    MAX_STEPS,

    /**
     * RK4 produced a NaN. Should not occur with a sensible horizon cushion;
     * treated as {@link #CAPTURED} by the shader. Logged by Renderer.
     */
    NAN,

    /**
     * Ray crossed the equatorial plane inside the disk's radial extent
     * {@code [rIsco, rOuter]}. Produced only by
     * {@code AdaptiveRayTracer} when constructed with a {@code Disk}.
     * The terminal state is the bisected (interpolated) crossing point,
     * suitable for shader redshift / temperature evaluation.
     */
    HIT_DISK
}
