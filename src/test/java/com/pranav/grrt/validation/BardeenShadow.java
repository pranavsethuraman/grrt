package com.pranav.grrt.validation;

/**
 * Test oracle for the Bardeen (1973) equatorial Kerr shadow.
 *
 * <p>Not a library class. Scope is validation-test only; package-private.
 *
 * <p>Computes the analytic horizontal shadow diameter D_H at inclination
 * i by parametrising the photon-orbit boundary by the spherical photon
 * radius r_ph ∈ [r_ph_pro, r_ph_retro]:
 * <pre>
 *   xi(r_ph)  = −(r_ph³ − 3 M r_ph² + a² r_ph + a² M) / (a (r_ph − M))
 *   eta(r_ph) =  r_ph³ (4 a² M − r_ph (r_ph − 3M)²) / (a² (r_ph − M)²)
 * </pre>
 * and taking α(r_ph) = −xi(r_ph) / sin(i), β(r_ph) = √(eta(r_ph) +
 * a² cos²i − xi² cot²i). At i = π/2 this reduces to α = −xi,
 * β = √(eta), so the shadow touches the β = 0 axis at the prograde and
 * retrograde equatorial photon orbits. D_H is the α extent of the
 * visible portion of the curve.
 *
 * <p>Reference: Bardeen, J. M. (1973), in <i>Black Holes</i> (eds.
 * DeWitt & DeWitt), Gordon & Breach, pp. 215–239.
 */
final class BardeenShadow {

    private BardeenShadow() {}

    /**
     * Prograde and retrograde equatorial photon-sphere radii as roots
     * of r (r − 3M)² = 4 M a². Closed-form trig solution (Cardano):
     * <pre>
     *   r_ph(a, σ) = 2 M { 1 + cos[ (2/3) arccos(−σ a / M) ] }
     *   σ = +1 → prograde (smaller r at a > 0)
     *   σ = −1 → retrograde (larger r at a > 0)
     * </pre>
     * At a = 0 both roots collapse to 3 M.
     *
     * @return {r_ph_pro, r_ph_retro}
     */
    static double[] photonOrbitRadii(double M, double a) {
        if (!(M > 0.0)) throw new IllegalArgumentException("M > 0 required");
        if (!(Math.abs(a) < M)) throw new IllegalArgumentException("|a| < M required");
        double argPro  = Math.acos(-a / M);
        double argRetr = Math.acos( a / M);
        double rPro   = 2.0 * M * (1.0 + Math.cos((2.0 / 3.0) * argPro));
        double rRetro = 2.0 * M * (1.0 + Math.cos((2.0 / 3.0) * argRetr));
        return new double[] { rPro, rRetro };
    }

    /** Bardeen's xi = L/E for a spherical photon orbit at radius r. */
    static double xi(double M, double a, double r) {
        return -(r * r * r - 3.0 * M * r * r + a * a * r + a * a * M)
             / (a * (r - M));
    }

    /** Bardeen's eta = Q/E² for a spherical photon orbit at radius r. */
    static double eta(double M, double a, double r) {
        double rm3 = r - 3.0 * M;
        return r * r * r * (4.0 * a * a * M - r * rm3 * rm3)
             / (a * a * (r - M) * (r - M));
    }

    /**
     * Observer-plane α under the Camera.java sign convention.
     *
     * <p>Sign convention matches {@link com.pranav.grrt.camera.Camera}:
     * positive screen-α corresponds to k^φ > 0 at the observer, which
     * (in Kerr with a > 0 equatorial) is a prograde photon. Hence
     * prograde sits on screen-right.
     *
     * <p>This inverts Bardeen 1973's α = −ξ/sin(i) to α = +ξ/sin(i).
     * Use this when comparing against a rendered image produced through
     * {@link com.pranav.grrt.camera.Camera}. D_H (a difference) is
     * convention-insensitive, but which side is flat is not.
     */
    static double alphaCameraConvention(double M, double a, double r,
                                        double inclination) {
        return xi(M, a, r) / Math.sin(inclination);
    }

    /** Same scan as {@link #horizontalExtent} but under the Camera sign
     *  convention. Prograde α is positive. D_H value equals the Bardeen
     *  result (diameter is sign-independent). */
    static ShadowExtent horizontalExtentCameraConvention(
            double M, double a, double inclination, int samples) {
        if (samples < 2) throw new IllegalArgumentException("samples >= 2");
        double[] r = photonOrbitRadii(M, a);
        double rLo = r[0];
        double rHi = r[1];
        double sinI = Math.sin(inclination);

        double aMin = Double.POSITIVE_INFINITY;
        double aMax = Double.NEGATIVE_INFINITY;
        for (int k = 0; k < samples; k++) {
            double t = (double) k / (samples - 1);
            double rPh = rLo + t * (rHi - rLo);
            double alpha = xi(M, a, rPh) / sinI;   // Camera convention: +ξ/sin(i)
            if (alpha < aMin) aMin = alpha;
            if (alpha > aMax) aMax = alpha;
        }
        return new ShadowExtent(aMin, aMax, rLo, rHi);
    }

    /**
     * Scan r_ph uniformly across [r_ph_pro, r_ph_retro] with
     * {@code samples} points and record α_min, α_max at inclination
     * {@code inclination}. At i = π/2 (edge-on) this gives D_H = max α
     * − min α directly; at other inclinations the scan still produces
     * the α extent of the shadow boundary but the geometric D_H may
     * need additional handling of the β condition.
     *
     * @param samples number of uniform r_ph samples; must be ≥ 2
     */
    static ShadowExtent horizontalExtent(double M, double a,
                                         double inclination, int samples) {
        if (samples < 2) throw new IllegalArgumentException("samples >= 2");
        double[] r = photonOrbitRadii(M, a);
        double rLo = r[0];
        double rHi = r[1];
        double sinI = Math.sin(inclination);

        double aMin = Double.POSITIVE_INFINITY;
        double aMax = Double.NEGATIVE_INFINITY;
        for (int k = 0; k < samples; k++) {
            double t = (double) k / (samples - 1);
            double rPh = rLo + t * (rHi - rLo);
            double alpha = -xi(M, a, rPh) / sinI;
            if (alpha < aMin) aMin = alpha;
            if (alpha > aMax) aMax = alpha;
        }
        return new ShadowExtent(aMin, aMax, rLo, rHi);
    }

    /** Result of a horizontal-extent scan. */
    record ShadowExtent(double alphaMin, double alphaMax,
                        double rPhPro, double rPhRetro) {
        double dH() { return alphaMax - alphaMin; }
    }
}
