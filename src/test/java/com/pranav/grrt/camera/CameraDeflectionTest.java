package com.pranav.grrt.camera;

import com.pranav.grrt.integrator.RK4;
import com.pranav.grrt.metric.SchwarzschildMetric;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

/**
 * Camera + RK4 integration tests. Slower than {@link CameraTest} because they
 * integrate null geodesics through curved spacetime.
 */
class CameraDeflectionTest {

    // ---------------------------------------------------------------
    // A. Weak-field deflection: Δφ_total − (π − 2α) → 4M/b for a ray
    //    at large impact parameter.
    //
    //    b is computed from the conserved Killing charges:
    //      E = −g_tt k^t = f · k^t
    //      L = g_{φφ} k^φ  (at θ = π/2: r² k^φ)
    //      b = L / E
    //    (not from r_obs · sin α directly.)
    //
    //    In flat space the straight-line chord sweeps Δφ = π − 2α. The
    //    Schwarzschild deflection adds to that.
    //
    //    The well-known 4M/b formula is the asymptotic (∞ → ∞) deflection.
    //    For b = 100 M the second-order term (15π/4 − 4)(M/b)² ≈ 7.8×10⁻⁴
    //    is ~2% of 4M/b and cannot be neglected if we want 1%-level accuracy
    //    at the prescribed (b, r_obs). We therefore compare the measured
    //    chord deflection to the analytic series through O((M/b)²); the
    //    residual from O((M/b)³) + finite-r_obs effects is ≲ 1% of 4M/b
    //    at (b = 100 M, r_obs = 1000 M). See Hartle §9.
    // ---------------------------------------------------------------

    @Test
    void weakFieldDeflectionMatchesFourMOverB() {
        double M = 1.0;
        double rObs = 1000.0;
        double b = 100.0;

        SchwarzschildMetric metric = new SchwarzschildMetric(M);
        Camera cam = new Camera(metric, rObs, Math.PI / 2, 64, 64, 0.4);

        // Choose α so that L/E = b exactly.
        double fObs = 1.0 - 2.0 * M / rObs;
        double alpha = Math.asin(b * Math.sqrt(fObs) / rObs);

        double[] y    = new double[8];
        double[] next = new double[8];
        cam.initialStateFromAngles(alpha, 0.0, y);

        // Verify b from conserved quantities (not from r_obs · sin α).
        double E = fObs * y[4];
        double L = rObs * rObs * y[7];   // sin²(π/2) = 1
        assertEquals(b, L / E, 1e-10, "b from L/E");

        RK4 rk = new RK4();
        double h = 0.1;
        int nMax = 200_000;
        boolean turned = false;
        int steps = 0;
        for (int i = 0; i < nMax; i++) {
            rk.step(metric, y, h, next);
            System.arraycopy(next, 0, y, 0, 8);
            steps++;
            if (y[5] > 0.0) turned = true;
            if (turned && y[1] >= rObs) break;
            assertTrue(y[1] > 3.0 * M, "plunged at step " + i + ", r = " + y[1]);
        }
        assertTrue(steps < nMax, "did not return to r_obs within " + nMax + " steps");

        // Subtract flat-space chord Δφ = π − 2α.
        double deltaPhiFlat = Math.PI - 2.0 * alpha;
        double deflection   = y[3] - deltaPhiFlat;

        double leading      = 4.0 * M / b;
        double secondOrder  = (15.0 * Math.PI / 4.0 - 4.0) * (M / b) * (M / b);
        double expected     = leading + secondOrder;

        assertEquals(expected, deflection, 0.01 * leading,
                "deflection " + deflection + " vs analytic " + expected
                        + " (leading " + leading + ")");
    }

    // ---------------------------------------------------------------
    // C. Conservation of E and L along an equatorial ray, integrated
    //    for ~10^4 steps outside the horizon. Both should be preserved
    //    to better than 1 part in 10^8 with h = 0.01.
    // ---------------------------------------------------------------

    @Test
    void equatorialRayConservesEnergyAndAngularMomentum() {
        double M = 1.0;
        double rObs = 100.0;
        double b = 20.0;       // comfortably above b_crit = 3√3 M ≈ 5.196

        SchwarzschildMetric metric = new SchwarzschildMetric(M);
        Camera cam = new Camera(metric, rObs, Math.PI / 2, 64, 64, 0.4);

        double fObs = 1.0 - 2.0 * M / rObs;
        double alpha = Math.asin(b * Math.sqrt(fObs) / rObs);

        double[] y    = new double[8];
        double[] next = new double[8];
        cam.initialStateFromAngles(alpha, 0.0, y);

        double E0 = (1.0 - 2.0 * M / y[1]) * y[4];
        double L0 = y[1] * y[1] * Math.sin(y[2]) * Math.sin(y[2]) * y[7];

        RK4 rk = new RK4();
        double h = 0.01;
        int n = 10_000;
        double maxDE = 0.0;
        double maxDL = 0.0;
        for (int i = 0; i < n; i++) {
            rk.step(metric, y, h, next);
            System.arraycopy(next, 0, y, 0, 8);
            assertTrue(y[1] > 2.1 * M, "plunged at step " + i + ", r = " + y[1]);

            double f  = 1.0 - 2.0 * M / y[1];
            double s  = Math.sin(y[2]);
            double E  = f * y[4];
            double L  = y[1] * y[1] * s * s * y[7];
            maxDE = Math.max(maxDE, Math.abs((E - E0) / E0));
            maxDL = Math.max(maxDL, Math.abs((L - L0) / L0));
        }
        assertTrue(maxDE < 1e-8, "|ΔE/E| too large: " + maxDE);
        assertTrue(maxDL < 1e-8, "|ΔL/L| too large: " + maxDL);
    }
}
