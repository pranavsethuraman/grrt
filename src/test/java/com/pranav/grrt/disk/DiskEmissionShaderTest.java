package com.pranav.grrt.disk;

import com.pranav.grrt.metric.KerrMetric;
import com.pranav.grrt.metric.Metric;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

/**
 * Phase 3B unit tests for {@link DiskEmissionShader} redshift and
 * intensity. Gate 3 (full face-on render) lands in 3B.2 once the
 * renderer wiring is complete; for 3B.1 we test the shader formulas
 * directly with synthesized disk-hit states.
 */
class DiskEmissionShaderTest {

    private static final double M = 1.0;
    private static final double R_OUTER = 25.0;

    // ------------------------------------------------------------------
    // Schwarzschild ISCO face-on radial photon: redshift = 1/√2
    // ------------------------------------------------------------------

    /**
     * For Schwarzschild ({@code a = 0}) at {@code r_ISCO = 6 M}, a
     * face-on radial photon (k^φ = k^θ = 0) emitted from a Keplerian
     * disk has redshift {@code g = 1/√2} ≈ 0.7071. Derivation:
     * <pre>
     *   g_tt = -(1 - 2/r) = -2/3
     *   g_φφ = r² sin²θ = 36
     *   Ω    = r^{-3/2} = 1/(6√6)
     *   u^t  = 1 / √(-(g_tt + g_φφ Ω²)) = 1/√(1/2) = √2
     *   For radial k (k^φ = 0): g = g_tt / (u^t · g_tt) = 1/u^t = 1/√2
     * </pre>
     */
    @Test
    void redshiftAtSchwarzschildIscoFaceOnRadialPhoton() {
        Metric schwarz = new KerrMetric(M, 0.0);
        Disk disk = new NovikovThorneDisk(schwarz, 0.0, R_OUTER, 0.0);
        DiskEmissionShader shader = new DiskEmissionShader(disk, schwarz);

        double[] state = { 0.0, 6.0, Math.PI / 2.0, 0.0,
                           1.0, -1.0, 0.0, 0.0 };
        double redshift = shader.redshift(state);
        double expected = 1.0 / Math.sqrt(2.0);
        assertEquals(expected, redshift, 1e-12,
                "Schwarzschild ISCO face-on radial: expected 1/√2 ≈ "
                        + expected + ", got " + redshift);
    }

    // ------------------------------------------------------------------
    // Kerr ISCO face-on radial photon: shader matches alternative form
    // ------------------------------------------------------------------

    /**
     * Independent face-on-radial reduction of the redshift formula:
     * <pre>
     *   g = g_tt / (u^t (g_tt + g_tφ Ω))
     * </pre>
     * for k^φ = k^θ = 0. We compute u^t and the (g_tt, g_tφ, g_φφ) at
     * the disk hit using {@code metric.g(x)} (independent of the
     * shader's full {@code -k·u_obs / -k·u_emit} contraction), and
     * assert the two forms agree at machine precision.
     */
    @Test
    void redshiftMatchesFaceOnRadialReductionForKerrA09() {
        Metric kerr = new KerrMetric(M, 0.9);
        Disk disk = new NovikovThorneDisk(kerr, 0.9, R_OUTER, 0.0);
        DiskEmissionShader shader = new DiskEmissionShader(disk, kerr);

        double r = disk.rIsco();   // ≈ 2.32088 M
        double[] state = { 0.0, r, Math.PI / 2.0, 0.0,
                           1.0, -1.0, 0.0, 0.0 };

        double[] x = { 0.0, r, Math.PI / 2.0, 0.0 };
        double[][] gmn = kerr.g(x);
        double gtt = gmn[0][0];
        double gtp = gmn[0][3];
        double gpp = gmn[3][3];
        double rootM = Math.sqrt(M);
        double rootR = Math.sqrt(r);
        double omega = rootM / (r * rootR + 0.9 * rootM);
        double quad = -(gtt + 2.0 * gtp * omega + gpp * omega * omega);
        assertTrue(quad > 0.0, "expected timelike Keplerian normalization");
        double uT = 1.0 / Math.sqrt(quad);
        double expected = gtt / (uT * (gtt + gtp * omega));

        double redshift = shader.redshift(state);
        assertEquals(expected, redshift, 1e-12,
                "Kerr a=0.9 ISCO face-on radial: alternative-form expected "
                        + expected + ", shader returned " + redshift);
        assertTrue(redshift > 0.0 && redshift < 1.0,
                "physical redshift should be in (0,1): " + redshift);
    }

    // ------------------------------------------------------------------
    // Intensity: g^4 sigma T^4 / pi at a known disk-hit configuration
    // ------------------------------------------------------------------

    @Test
    void intensityIsPositiveAndConsistentWithRedshiftAndTemperature() {
        Metric kerr = new KerrMetric(M, 0.9);
        Disk disk = new NovikovThorneDisk(kerr, 0.9, R_OUTER, 1e-3);
        DiskEmissionShader shader = new DiskEmissionShader(disk, kerr);

        double r = 4.0;
        double[] state = { 0.0, r, Math.PI / 2.0, 0.0,
                           1.0, -1.0, 0.0, 0.0 };

        double redshift = shader.redshift(state);
        double t = disk.temperature(r);
        double g2 = redshift * redshift;
        double g4 = g2 * g2;
        double t2 = t * t;
        double t4 = t2 * t2;
        double expected = g4 * t4 / Math.PI;
        float intensity = shader.intensity(state);

        assertTrue(intensity > 0.0f,
                "intensity should be > 0 for a disk hit (got " + intensity + ")");
        assertEquals((float) expected, intensity, (float)(Math.abs(expected) * 1e-6),
                "intensity should equal g^4 σ T^4 / π");
    }

    // ------------------------------------------------------------------
    // Constructor input validation
    // ------------------------------------------------------------------

    @Test
    void constructorRejectsNullArgs() {
        Metric kerr = new KerrMetric(M, 0.9);
        Disk disk = new NovikovThorneDisk(kerr, 0.9, R_OUTER, 0.0);
        assertThrows(IllegalArgumentException.class,
                () -> new DiskEmissionShader(null, kerr));
        assertThrows(IllegalArgumentException.class,
                () -> new DiskEmissionShader(disk, null));
    }

    @Test
    void redshiftRejectsBadStateLength() {
        Metric kerr = new KerrMetric(M, 0.9);
        Disk disk = new NovikovThorneDisk(kerr, 0.9, R_OUTER, 0.0);
        DiskEmissionShader shader = new DiskEmissionShader(disk, kerr);
        assertThrows(IllegalArgumentException.class,
                () -> shader.redshift(new double[7]));
        assertThrows(IllegalArgumentException.class,
                () -> shader.redshift(null));
    }
}
