package com.pranav.grrt.disk;

import com.pranav.grrt.metric.Metric;
import com.pranav.grrt.renderer.RayOutcome;
import com.pranav.grrt.renderer.RayShader;

/**
 * Bolometric Planck emission from a Novikov-Thorne disk, observed at
 * infinity by a static observer aligned with the asymptotic Lorentz
 * frame. Pixel intensity is the redshifted bolometric Planck integral
 *
 * <pre>
 *   I_obs = g⁴ · σ_SB T(r_emit)⁴ / π
 * </pre>
 *
 * where the redshift factor is
 *
 * <pre>
 *   g = -k_μ u_obs^μ / -k_ν u_emit^ν
 * </pre>
 *
 * Both inner products are evaluated at the disk-hit point: the
 * observer term reduces to {@code -k_t} (k_t is conserved along the
 * geodesic via the Killing vector ∂_t, and {@code u_obs^μ = (1,0,0,0)}
 * at infinity), and the emitter term is computed from the Keplerian
 * 4-velocity at the disk.
 *
 * <h2>Phase 3B.1 status</h2>
 *
 * The renderer wiring (returning {@link RayOutcome#HIT_DISK} from
 * {@code AdaptiveRayTracer}) lands in 3B.2. For 3B.1, the
 * {@link #shade} entry point always returns 0 — testing exercises
 * {@link #redshift} and {@link #intensity} directly with synthesized
 * disk-hit states.
 *
 * <h2>Thread-safety</h2>
 *
 * Stateless apart from final fields. Safe for concurrent reads.
 */
public final class DiskEmissionShader implements RayShader {

    private static final double STEFAN_BOLTZMANN = 1.0;

    private final Disk disk;
    private final Metric metric;

    public DiskEmissionShader(Disk disk, Metric metric) {
        if (disk == null) {
            throw new IllegalArgumentException("disk must not be null");
        }
        if (metric == null) {
            throw new IllegalArgumentException("metric must not be null");
        }
        this.disk = disk;
        this.metric = metric;
    }

    /**
     * Dispatch on {@link RayOutcome#HIT_DISK} (Phase 3B.2 wiring):
     * disk-hit rays return the bolometric Planck intensity at the
     * crossing point; all other outcomes return 0.
     */
    @Override
    public float shade(RayOutcome outcome, double[] endState) {
        if (outcome == RayOutcome.HIT_DISK) {
            return intensity(endState);
        }
        return 0.0f;
    }

    /**
     * Redshift factor g = -k_μ u_obs^μ / -k_ν u_emit^ν evaluated at
     * the disk-hit point.
     *
     * @param state 8-state at the disk hit
     *              [t, r, π/2, φ, k^t, k^r, k^θ, k^φ]
     * @return redshift factor g (dimensionless; &lt; 1 for red-shifted
     *         emission, &gt; 1 for blue-shifted)
     * @throws IllegalArgumentException if r is outside the disk extent
     */
    public double redshift(double[] state) {
        if (state == null || state.length != 8) {
            throw new IllegalArgumentException("state must have length 8");
        }
        double[] x = { state[0], state[1], state[2], state[3] };
        double[] k = { state[4], state[5], state[6], state[7] };
        double[] uEmit = disk.keplerianFourVelocity(x);
        double[][] gAtDisk = metric.g(x);

        // k_t at disk = g_{0,μ} k^μ. Conserved along the geodesic, so this
        // equals k_t at the asymptotic observer.
        double kT = 0.0;
        for (int mu = 0; mu < 4; mu++) {
            kT += gAtDisk[0][mu] * k[mu];
        }

        // -k_ν u_emit^ν = -g_{μ,ν} k^μ u_emit^ν (full sum; off-diagonal
        // g_tφ contributes via both (μ,ν) = (t,φ) and (φ,t)).
        double dotEmit = 0.0;
        for (int mu = 0; mu < 4; mu++) {
            for (int nu = 0; nu < 4; nu++) {
                dotEmit += gAtDisk[mu][nu] * k[mu] * uEmit[nu];
            }
        }

        // g = -kT / -dotEmit = kT / dotEmit.
        return kT / dotEmit;
    }

    /**
     * Bolometric Planck intensity at the observer for a disk-hit state:
     * {@code I_obs = g^4 · σ_SB T(r)^4 / π}. Returns a {@code float}
     * to match the rendering pipeline's pixel format.
     *
     * @param state 8-state at the disk hit; r = state[1] must be in
     *              [rIsco, rOuter] of the underlying disk
     * @return observer-frame intensity (geometrized, σ_SB = 1)
     * @throws IllegalArgumentException if r is outside the disk extent
     */
    public float intensity(double[] state) {
        double r = state[1];
        double t = disk.temperature(r);
        double t2 = t * t;
        double t4 = t2 * t2;
        double g = redshift(state);
        double g2 = g * g;
        double g4 = g2 * g2;
        return (float) (g4 * STEFAN_BOLTZMANN * t4 / Math.PI);
    }
}
