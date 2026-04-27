package com.pranav.grrt.disk;

/**
 * Geometrically thin, optically thick accretion disk in the equatorial
 * plane (θ = π/2) of an axisymmetric host metric.
 *
 * <p>Single instance per render frame. The instance caches the inner
 * edge {@code r_ISCO} via the host metric at construction; the renderer
 * reuses the same Disk across all pixels.
 *
 * <p><b>Thread-safety.</b> Implementations must be safe for concurrent
 * reads. The renderer's pixel loop calls {@link #crossedEquator},
 * {@link #temperature}, and {@link #keplerianFourVelocity} from worker
 * threads.
 *
 * <p>See {@code docs/phase-3b-plan.md} §3.1 for the role of this
 * interface in Phase 3B.
 */
public interface Disk {

    /**
     * Coarse pre-bisection equator-crossing test between two
     * consecutive accepted integrator states. True only if
     * (θ − π/2) changes sign between the segment endpoints AND the
     * segment's radial range overlaps the disk extent
     * {@code [rIsco, rOuter]}. The renderer follows up with a
     * {@code DormandPrince45.interpolate}-driven bisection to find the
     * exact crossing radius.
     *
     * @param xPrev 8-state at the start of the segment, length 8
     * @param xCurr 8-state at the end of the segment, length 8
     */
    boolean crossedEquator(double[] xPrev, double[] xCurr);

    /** Disk inner edge in coordinate r (typically {@code r_ISCO + cushion}). */
    double rIsco();

    /** Disk outer edge in coordinate r. */
    double rOuter();

    /**
     * Surface radial flux F(r) per the disk's emission model. For
     * {@link NovikovThorneDisk} this is the Page-Thorne 1974 eq. (15n)
     * form in geometrized {@code M = Ṁ = 1} units; F is dimensionless.
     *
     * @throws IllegalArgumentException if r &lt; rIsco or r &gt; rOuter
     */
    double surfaceFlux(double r);

    /**
     * Local effective temperature: T(r) = (F(r) / σ_SB)^(1/4). σ_SB is
     * normalized to 1 in geometrized units.
     *
     * @throws IllegalArgumentException if r &lt; rIsco or r &gt; rOuter
     */
    double temperature(double r);

    /**
     * Keplerian (prograde equatorial circular) 4-velocity u^μ of the
     * disk material at the equatorial position {@code xOnDisk}.
     * Returns a new {@code double[4]}; not allocation-free (cold path,
     * called only on disk-crossing hits).
     *
     * @param xOnDisk equatorial 4-position [t, r, π/2, φ], length 4
     * @throws IllegalArgumentException if r is outside [rIsco, rOuter]
     */
    double[] keplerianFourVelocity(double[] xOnDisk);
}
