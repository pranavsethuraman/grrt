package com.pranav.grrt.metric;

/**
 * Spacetime metric in spherical-type coordinates x = (t, r, θ, φ).
 *
 * <p>Conventions:
 * <ul>
 *   <li>Signature (-, +, +, +)</li>
 *   <li>Geometrized units: G = c = 1. Mass has dimensions of length.</li>
 *   <li>Coordinate index ordering: 0 = t, 1 = r, 2 = θ, 3 = φ</li>
 *   <li>Geodesic state vector layout: [t, r, θ, φ, k^t, k^r, k^θ, k^φ]</li>
 * </ul>
 */
public interface Metric {

    /** @return gravitational mass M in geometrized units */
    double mass();

    /** @return outer event horizon radius in geometrized units */
    double horizonRadius();

    /**
     * Innermost stable circular orbit (ISCO) radius for equatorial
     * timelike orbits in geometrized units.
     *
     * <p>Default: throws {@link UnsupportedOperationException}. Concrete
     * metrics that admit an equatorial ISCO ({@link KerrMetric},
     * {@link JohannsenPsaltisMetric}) override. Metrics without one
     * ({@link MinkowskiMetric}, and {@link SchwarzschildMetric} which
     * does not currently expose a closed-form ISCO) inherit the
     * throwing default.
     *
     * <p>Used by the disk model in {@code com.pranav.grrt.disk} to
     * compute the inner edge of a Novikov-Thorne disk.
     *
     * @param prograde true for the prograde branch (co-rotating with
     *                 spin); false for retrograde
     * @return the ISCO radius in geometrized units
     * @throws UnsupportedOperationException if the metric does not
     *                                       implement an ISCO formula
     */
    default double iscoRadius(boolean prograde) {
        throw new UnsupportedOperationException(
                "metric does not support iscoRadius");
    }

    /**
     * Test whether the position {@code x} sits at or inside the metric's
     * event horizon, within tolerance {@code tol}.
     *
     * <p>Default implementation: {@code x[1] - horizonRadius() < tol}.
     * Valid for any metric whose horizon locus is purely a function of r
     * (Schwarzschild, Kerr, Minkowski). Position-dependent metrics —
     * notably {@link JohannsenPsaltisMetric} for {@code ε₃ ≠ 0}, where
     * the horizon is the locus {@code Δ + a² h sin²θ = 0} — must
     * override.
     *
     * <p>Used by {@link com.pranav.grrt.renderer.AdaptiveRayTracer} as
     * the per-step ray-termination predicate. The default keeps Phase 2
     * behavior bit-exact: the existing {@code rNow < horizonRadius() +
     * cushion} test reduces to the same comparison.
     *
     * @param x   8-state vector (only positional components 0..3 are
     *            consulted) or 4-position vector
     * @param tol non-negative tolerance / cushion in geometrized units
     * @return true iff x is at or inside the horizon
     */
    default boolean isInsideHorizon(double[] x, double tol) {
        return x[1] - horizonRadius() < tol;
    }

    /**
     * Metric tensor g_{μν} at position x.
     *
     * @param x position 4-vector (t, r, θ, φ)
     * @return new 4x4 array with [μ][ν] = g_{μν}
     */
    double[][] g(double[] x);

    /**
     * Inverse metric tensor g^{μν} at position x.
     *
     * @param x position 4-vector (t, r, θ, φ)
     * @return new 4x4 array with [μ][ν] = g^{μν}
     */
    double[][] gInv(double[] x);

    /**
     * Christoffel symbols of the second kind at position x.
     * Symmetric in the lower indices: Γ^α_{μν} = Γ^α_{νμ}.
     *
     * @param x position 4-vector (t, r, θ, φ)
     * @return new 4x4x4 array with [α][μ][ν] = Γ^α_{μν}
     */
    double[][][] christoffel(double[] x);

    /**
     * Geodesic acceleration: d²x^μ/dλ² = -Γ^μ_{αβ} k^α k^β.
     * Implementations should override for performance (avoid allocating the
     * full Christoffel tensor in the integrator hot path).
     *
     * @param x   position 4-vector
     * @param k   momentum 4-vector (dx/dλ)
     * @param out output acceleration 4-vector; length 4; overwritten
     */
    default void geodesicAcceleration(double[] x, double[] k, double[] out) {
        double[][][] G = christoffel(x);
        for (int mu = 0; mu < 4; mu++) {
            double a = 0.0;
            for (int alpha = 0; alpha < 4; alpha++) {
                double kAlpha = k[alpha];
                double[] row = G[mu][alpha];
                for (int beta = 0; beta < 4; beta++) {
                    a -= row[beta] * kAlpha * k[beta];
                }
            }
            out[mu] = a;
        }
    }

    /**
     * Norm g_{μν} k^μ k^ν. For a null geodesic this is zero; the integrator
     * uses this to detect numerical drift off the light cone.
     *
     * @param x position 4-vector
     * @param k momentum 4-vector
     * @return g_{μν} k^μ k^ν
     */
    default double nullNorm(double[] x, double[] k) {
        double[][] gmn = g(x);
        double sum = 0.0;
        for (int mu = 0; mu < 4; mu++) {
            double[] row = gmn[mu];
            for (int nu = 0; nu < 4; nu++) {
                sum += row[nu] * k[mu] * k[nu];
            }
        }
        return sum;
    }
}
