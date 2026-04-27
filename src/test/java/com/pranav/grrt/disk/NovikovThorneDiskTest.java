package com.pranav.grrt.disk;

import com.pranav.grrt.metric.JohannsenPsaltisMetric;
import com.pranav.grrt.metric.KerrMetric;
import com.pranav.grrt.metric.Metric;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

/**
 * Phase 3B gate 1, gate 2 sanity floor, and JP-smoke tests for
 * {@link NovikovThorneDisk}.
 *
 * <p>Reference values for gate 1 are produced by
 * {@code scripts/page_thorne_reference.py} (Sympy + scipy.integrate.quad)
 * and pasted as a {@code double[][]} literal below. See
 * {@code docs/phase-3b-plan.md} §2.2 for the regeneration policy.
 */
class NovikovThorneDiskTest {

    private static final double M = 1.0;
    private static final double R_OUTER = 25.0;          // > all reference radii
    private static final double ISCO_CUSHION_NONE = 0.0; // match Python script

    /**
     * Page-Thorne 1974 eq. (15n) reference table.
     *
     * <p>Each row: {@code {a, r, F_pageThorne}} from
     * {@code scripts/page_thorne_reference.py} stdout. Values frozen
     * at user approval; regenerate the script and update this table
     * if {@link NovikovThorneDisk#surfaceFlux} ever changes form.
     *
     * <p>Provenance: 12-row table covering 6 radii × 2 spins (a ∈ {0, 0.9}),
     * radii from just outside ISCO to 20 M.
     */
    private static final double[][] PAGE_THORNE_REFERENCE = {
            // {a, r,  F_pageThorne}
            { 0.0,  6.500, 4.369295066852117e-07 },
            { 0.0,  7.000, 1.033625390075104e-06 },
            { 0.0,  8.000, 1.652303013402170e-06 },
            { 0.0, 10.000, 1.517032887034427e-06 },
            { 0.0, 15.000, 5.970509005254869e-07 },
            { 0.0, 20.000, 2.490369586702499e-07 },
            { 0.9,  2.500, 5.136489610260793e-05 },
            { 0.9,  3.000, 1.523412395335334e-04 },
            { 0.9,  4.000, 1.044880487966974e-04 },
            { 0.9,  6.000, 3.144018602870300e-05 },
            { 0.9, 10.000, 5.344236799923510e-06 },
            { 0.9, 20.000, 4.194631402203248e-07 },
    };

    // ------------------------------------------------------------------
    // Gate 1: NT radial flux F(r) vs Page-Thorne 1974 eq. (15n) reference
    // ------------------------------------------------------------------

    @Test
    void surfaceFluxMatchesPageThorneReference() {
        Metric kerrA00 = new KerrMetric(M, 0.0);
        Metric kerrA09 = new KerrMetric(M, 0.9);
        NovikovThorneDisk diskA00 =
                new NovikovThorneDisk(kerrA00, 0.0, R_OUTER, ISCO_CUSHION_NONE);
        NovikovThorneDisk diskA09 =
                new NovikovThorneDisk(kerrA09, 0.9, R_OUTER, ISCO_CUSHION_NONE);

        double maxRelErr = 0.0;
        for (double[] row : PAGE_THORNE_REFERENCE) {
            double a = row[0];
            double r = row[1];
            double fRef = row[2];
            NovikovThorneDisk disk = (a == 0.0) ? diskA00 : diskA09;
            double fJava = disk.surfaceFlux(r);
            double absErr = Math.abs(fJava - fRef);
            double relErr = absErr / Math.abs(fRef);
            maxRelErr = Math.max(maxRelErr, relErr);
            assertEquals(fRef, fJava, Math.abs(fRef) * 1e-6,
                    "Page-Thorne flux at (a=" + a + ", r=" + r
                            + "): Java=" + fJava + ", ref=" + fRef
                            + ", relErr=" + relErr);
        }
        System.out.println("[gate1] max relative error vs Page-Thorne reference = "
                + maxRelErr);
    }

    // ------------------------------------------------------------------
    // Gate 2 sanity floor: r_ISCO via existing Bardeen formula
    // ------------------------------------------------------------------

    @Test
    void iscoMatchesBardeenForKerrA09() {
        Metric kerr = new KerrMetric(M, 0.9);
        NovikovThorneDisk disk =
                new NovikovThorneDisk(kerr, 0.9, R_OUTER, ISCO_CUSHION_NONE);
        assertEquals(2.32088, disk.rIsco(), 1e-4,
                "ISCO at a=0.9 prograde should be 2.32088 M (Bardeen 1972)");
        Metric schwarz = new KerrMetric(M, 0.0);
        NovikovThorneDisk diskSchwarz =
                new NovikovThorneDisk(schwarz, 0.0, R_OUTER, ISCO_CUSHION_NONE);
        assertEquals(6.0, diskSchwarz.rIsco(), 1e-10,
                "ISCO at a=0 should be 6 M");
    }

    // ------------------------------------------------------------------
    // JP smoke: constructor accepts smooth-regime ε₃ values without throw
    // ------------------------------------------------------------------

    @Test
    void constructorWorksForJpInSmoothRegime() {
        // Phase 3-plan §6.5 sweep grid endpoints + interior in the
        // smooth-deformation regime ε₃ ∈ (-2.97, +0.12) at a=0.9.
        // Constructor must compute r_ISCO via JP eCircular without
        // throwing. We don't call surfaceFlux because the Page-Thorne
        // formula uses Kerr Bardeen quantities (the standard
        // approximation; see NovikovThorneDisk Javadoc "Spin parameter").
        double[] eps3Smooth = { -2.5, -1.0, -0.5, +0.05, +0.10 };
        for (double eps3 : eps3Smooth) {
            JohannsenPsaltisMetric jp =
                    new JohannsenPsaltisMetric(M, 0.9, eps3);
            NovikovThorneDisk disk =
                    new NovikovThorneDisk(jp, 0.9, R_OUTER, 1.0e-3);
            assertTrue(Double.isFinite(disk.rIsco()),
                    "rIsco must be finite for ε₃=" + eps3);
            assertTrue(disk.rIsco() > 1.0,
                    "rIsco must exceed Kerr horizon for ε₃=" + eps3
                            + " (got " + disk.rIsco() + ")");
            assertTrue(disk.rIsco() < R_OUTER,
                    "rIsco must be inside disk extent for ε₃=" + eps3);
        }
    }

    // ------------------------------------------------------------------
    // Disk-extent and equator-crossing unit tests
    // ------------------------------------------------------------------

    @Test
    void crossedEquatorReturnsTrueOnSignChangeInsideExtent() {
        Metric kerr = new KerrMetric(M, 0.9);
        NovikovThorneDisk disk =
                new NovikovThorneDisk(kerr, 0.9, R_OUTER, 1.0e-3);
        double[] xPrev = { 0, 10.0, Math.PI / 2 + 0.1, 0,  0, 0, -0.1, 0 };
        double[] xCurr = { 1, 10.0, Math.PI / 2 - 0.1, 0,  0, 0, -0.1, 0 };
        assertTrue(disk.crossedEquator(xPrev, xCurr));
    }

    @Test
    void crossedEquatorReturnsFalseOutsideExtent() {
        Metric kerr = new KerrMetric(M, 0.9);
        NovikovThorneDisk disk =
                new NovikovThorneDisk(kerr, 0.9, R_OUTER, 1.0e-3);
        double[] xPrev = { 0, 100.0, Math.PI / 2 + 0.1, 0,  0, 0, -0.1, 0 };
        double[] xCurr = { 1, 100.0, Math.PI / 2 - 0.1, 0,  0, 0, -0.1, 0 };
        assertFalse(disk.crossedEquator(xPrev, xCurr));
    }

    @Test
    void crossedEquatorReturnsFalseOnSameSide() {
        Metric kerr = new KerrMetric(M, 0.9);
        NovikovThorneDisk disk =
                new NovikovThorneDisk(kerr, 0.9, R_OUTER, 1.0e-3);
        double[] xPrev = { 0, 10.0, Math.PI / 2 + 0.1, 0,  0, 0, -0.1, 0 };
        double[] xCurr = { 1, 10.0, Math.PI / 2 + 0.05, 0, 0, 0, -0.1, 0 };
        assertFalse(disk.crossedEquator(xPrev, xCurr));
    }

    // ------------------------------------------------------------------
    // Keplerian 4-velocity normalization: u^μ u_μ = -1
    // ------------------------------------------------------------------

    @Test
    void keplerianFourVelocityIsTimelikeNormalized() {
        Metric kerr = new KerrMetric(M, 0.9);
        NovikovThorneDisk disk =
                new NovikovThorneDisk(kerr, 0.9, R_OUTER, 0.0);
        for (double r : new double[] { 2.5, 5.0, 10.0, 20.0 }) {
            double[] x = { 0.0, r, Math.PI / 2.0, 0.0 };
            double[] u = disk.keplerianFourVelocity(x);
            double[][] gmn = kerr.g(x);
            double norm = 0.0;
            for (int mu = 0; mu < 4; mu++) {
                for (int nu = 0; nu < 4; nu++) {
                    norm += gmn[mu][nu] * u[mu] * u[nu];
                }
            }
            assertEquals(-1.0, norm, 1e-12,
                    "u·u != -1 at r=" + r + " (got " + norm + ")");
        }
    }

    // ------------------------------------------------------------------
    // Constructor input validation
    // ------------------------------------------------------------------

    @Test
    void constructorRejectsBadInputs() {
        Metric kerr = new KerrMetric(M, 0.9);
        assertThrows(IllegalArgumentException.class,
                () -> new NovikovThorneDisk(null, 0.9, R_OUTER, 0.0));
        assertThrows(IllegalArgumentException.class,
                () -> new NovikovThorneDisk(kerr, -0.1, R_OUTER, 0.0));
        assertThrows(IllegalArgumentException.class,
                () -> new NovikovThorneDisk(kerr, 1.0, R_OUTER, 0.0));
        assertThrows(IllegalArgumentException.class,
                () -> new NovikovThorneDisk(kerr, 0.9, -1.0, 0.0));
        assertThrows(IllegalArgumentException.class,
                () -> new NovikovThorneDisk(kerr, 0.9, 1.0, 0.0));
        assertThrows(IllegalArgumentException.class,
                () -> new NovikovThorneDisk(kerr, 0.9, R_OUTER, -1.0));
    }

    @Test
    void surfaceFluxRejectsOutOfBoundsR() {
        Metric kerr = new KerrMetric(M, 0.9);
        NovikovThorneDisk disk =
                new NovikovThorneDisk(kerr, 0.9, R_OUTER, 0.0);
        assertThrows(IllegalArgumentException.class,
                () -> disk.surfaceFlux(1.5));
        assertThrows(IllegalArgumentException.class,
                () -> disk.surfaceFlux(R_OUTER + 1.0));
    }
}
