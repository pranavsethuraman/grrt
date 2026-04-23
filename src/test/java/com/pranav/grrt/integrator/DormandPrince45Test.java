package com.pranav.grrt.integrator;

import com.pranav.grrt.metric.Metric;
import com.pranav.grrt.metric.SchwarzschildMetric;
import org.junit.jupiter.api.Test;

import java.util.function.DoubleUnaryOperator;

import static org.junit.jupiter.api.Assertions.*;

class DormandPrince45Test {

    private static final double TABLEAU_TOL = 1e-14;

    // ------------------------------------------------------------------
    // Butcher tableau consistency
    // ------------------------------------------------------------------

    @Test
    void tableauRowSumsEqualC() {
        // Σ_j A[i][j] == c[i] for stages 2..7. Each stage evaluates RHS
        // at t_n + c_i * h, and the A-row sums must match for consistency.
        assertEquals(DormandPrince45.C2,
                DormandPrince45.A21, TABLEAU_TOL, "row 2");
        assertEquals(DormandPrince45.C3,
                DormandPrince45.A31 + DormandPrince45.A32, TABLEAU_TOL, "row 3");
        assertEquals(DormandPrince45.C4,
                DormandPrince45.A41 + DormandPrince45.A42 + DormandPrince45.A43,
                TABLEAU_TOL, "row 4");
        assertEquals(DormandPrince45.C5,
                DormandPrince45.A51 + DormandPrince45.A52
                        + DormandPrince45.A53 + DormandPrince45.A54,
                TABLEAU_TOL, "row 5");
        assertEquals(DormandPrince45.C6,
                DormandPrince45.A61 + DormandPrince45.A62 + DormandPrince45.A63
                        + DormandPrince45.A64 + DormandPrince45.A65,
                TABLEAU_TOL, "row 6");
        // Row 7 (A72 = 0 not stored).
        assertEquals(DormandPrince45.C7,
                DormandPrince45.A71 + DormandPrince45.A73 + DormandPrince45.A74
                        + DormandPrince45.A75 + DormandPrince45.A76,
                TABLEAU_TOL, "row 7");
    }

    @Test
    void tableauBSumsEqualOne() {
        // B5_2 = B5_7 = 0 and B4_2 = 0, not stored; the remaining entries
        // must each sum to 1 (order-1 condition).
        double sumB5 = DormandPrince45.B5_1 + DormandPrince45.B5_3
                     + DormandPrince45.B5_4 + DormandPrince45.B5_5
                     + DormandPrince45.B5_6;
        double sumB4 = DormandPrince45.B4_1 + DormandPrince45.B4_3
                     + DormandPrince45.B4_4 + DormandPrince45.B4_5
                     + DormandPrince45.B4_6 + DormandPrince45.B4_7;
        assertEquals(1.0, sumB5, TABLEAU_TOL, "Σ b5");
        assertEquals(1.0, sumB4, TABLEAU_TOL, "Σ b4");
    }

    @Test
    void fsalPropertyStage7EqualsB5() {
        // FSAL: the stage-7 row A[7][0..5] equals b5[0..5] (with b5[6]=0
        // implicit). This ensures k7 = f(y_{n+1}^(5)) so the accepted
        // step's k7 can be reused as k1 of the next step.
        assertEquals(DormandPrince45.B5_1, DormandPrince45.A71, TABLEAU_TOL);
        assertEquals(DormandPrince45.B5_3, DormandPrince45.A73, TABLEAU_TOL);
        assertEquals(DormandPrince45.B5_4, DormandPrince45.A74, TABLEAU_TOL);
        assertEquals(DormandPrince45.B5_5, DormandPrince45.A75, TABLEAU_TOL);
        assertEquals(DormandPrince45.B5_6, DormandPrince45.A76, TABLEAU_TOL);
    }

    // ------------------------------------------------------------------
    // B2 (a): Schwarzschild shadow edge via DP45 bisection
    // ------------------------------------------------------------------

    @Test
    void schwarzschildShadowEdgeMatchesPhase1() {
        // Bisect on screen-plane angle α (β = 0) to find the capture/escape
        // boundary for DP45 at atol=rtol=1e-8. Compare to the analytic
        // prediction α_pred = asin(3√3 M · √f / r_obs) — Phase 1 RK4 at 256²
        // matches this to sub-pixel precision per CLAUDE.md. Tolerance:
        // 0.5 pixels at 256² with Phase 1 fov = 0.02 rad ⇒ 0.5 · 0.02/256
        // ≈ 3.91e-5 rad.
        Metric m = new SchwarzschildMetric(1.0);
        double rObs = 1000.0;
        double fov = 0.02;
        int W = 256;
        double pixel = fov / W;

        double alphaEdge = bisectShadowEdge(m, rObs, 1e-8, 1e-8);
        double f = 1.0 - 2.0 * m.mass() / rObs;
        double alphaPred = Math.asin(3.0 * Math.sqrt(3.0) * Math.sqrt(f) / rObs);

        double dxPixels = (alphaEdge - alphaPred) / pixel;
        System.out.println("[B2a] α_edge = " + alphaEdge
                + ", α_pred = " + alphaPred
                + ", Δ = " + dxPixels + " px (0.5 px allowed)");
        assertTrue(Math.abs(dxPixels) < 0.5,
                "DP45 shadow edge off by " + dxPixels + " px");
    }

    private double bisectShadowEdge(Metric m, double rObs, double atol, double rtol) {
        double alphaLo = 0.0;       // α = 0 is central → captured
        double alphaHi = 0.01;      // 10 mrad → well outside shadow → escapes
        assertTrue(tracePhotonCaptured(m, rObs, alphaLo, atol, rtol),
                "central ray must be captured");
        assertFalse(tracePhotonCaptured(m, rObs, alphaHi, atol, rtol),
                "α=0.01 ray must escape");
        // Bisect to resolution well under 0.5 px (pixel = fov/W ~ 7.8e-5 rad).
        for (int k = 0; k < 40; k++) {
            double mid = 0.5 * (alphaLo + alphaHi);
            if (tracePhotonCaptured(m, rObs, mid, atol, rtol)) {
                alphaLo = mid;
            } else {
                alphaHi = mid;
            }
            if (alphaHi - alphaLo < 1e-9) break;
        }
        return 0.5 * (alphaLo + alphaHi);
    }

    /**
     * Trace a ray with DP45 from the static-observer tetrad at
     * (r=rObs, θ=π/2), heading inward with screen-plane angle α
     * (β = 0, equatorial). Returns true if captured (r dips below
     * horizonRadius + 0.01 M), false if escaped past 2·rObs.
     */
    private static boolean tracePhotonCaptured(
            Metric m, double rObs, double alpha, double atol, double rtol) {
        double[] y = new double[8];
        double f = 1.0 - 2.0 * m.mass() / rObs;
        double sqrtF = Math.sqrt(f);
        // Static-observer tetrad at equatorial θ = π/2, β = 0. See
        // Camera.java §Tetrad derivation for the same construction.
        y[0] = 0.0;
        y[1] = rObs;
        y[2] = Math.PI / 2.0;
        y[3] = 0.0;
        y[4] = 1.0 / sqrtF;                 // k^t
        y[5] = -sqrtF * Math.cos(alpha);    // k^r (inward)
        y[6] = 0.0;                         // k^θ
        y[7] = Math.sin(alpha) / rObs;      // k^φ (equator: sin θ = 1)

        DormandPrince45 dp = new DormandPrince45();
        double[] yNext = new double[8];
        double h = 1.0;
        double rCapture = m.horizonRadius() + 0.01;
        double rEscape = 2.0 * rObs;
        int maxSteps = 200_000;

        for (int step = 0; step < maxSteps; step++) {
            AdaptiveIntegrator.StepStatus s =
                    dp.adaptiveStep(m, y, h, atol, rtol, yNext);
            if (s.accepted()) {
                System.arraycopy(yNext, 0, y, 0, 8);
                double r = y[1];
                if (Double.isNaN(r)) return true;
                if (r < rCapture) return true;
                if (r > rEscape)  return false;
            }
            h = s.hNext();
        }
        // Ran out of steps near the photon sphere — shade as captured
        // (same convention as Renderer for MAX_STEPS).
        return true;
    }

    // ------------------------------------------------------------------
    // B2 (b): accepted step count scales with tolerance
    // ------------------------------------------------------------------

    @Test
    void stepCountScalesWithToleranceAsOrderFive() {
        // Fixed integration: outgoing Schwarzschild photon from r=10M
        // equatorial with L=0 (radial), integrated to affine λ = 100 M.
        // For an order-5 method, accepted step count N ∝ tol^(-1/5);
        // tightening tol by 1000× grows N by 1000^(1/5) = 3.981.
        Metric m = new SchwarzschildMetric(1.0);

        int n6 = countAcceptedSteps(m, 1e-6, 1e-6);
        int n9 = countAcceptedSteps(m, 1e-9, 1e-9);
        double ratio = (double) n9 / n6;
        double expected = Math.pow(10.0, 3.0 / 5.0);
        System.out.println("[B2b] N(1e-6)=" + n6 + ", N(1e-9)=" + n9
                + ", ratio=" + ratio + ", expected≈" + expected);
        // User-specified band [3.2, 4.8] = roughly ±20% of 3.981.
        assertTrue(ratio > 3.2 && ratio < 4.8,
                "step count ratio = " + ratio + " (expected ≈ " + expected + ")");
    }

    private static int countAcceptedSteps(Metric m, double atol, double rtol) {
        // Inbound Schwarzschild equatorial photon from r=50 M with L=5.5
        // (just above b_crit = 3√3 ≈ 5.196). The photon spirals toward a
        // turning point near r ≈ 3.77 M (root of r³ − L²r + 2L² = 0 with
        // L=5.5), loops around the photon sphere, and escapes. This
        // produces rich dynamics that keep the PI controller in its
        // asymptotic regime, so step count scales cleanly as tol^(-1/5).
        double r = 50.0;
        double L = 5.5;
        double f = 1.0 - 2.0 * m.mass() / r;
        double kT  = 1.0 / f;                   // E = 1
        double kPh = L / (r * r);
        // Null condition: g_tt(k^t)² + g_rr(k^r)² + g_φφ(k^φ)² = 0
        //   ⇒ (k^r)² = f [ f(k^t)² − r²(k^φ)² ]
        //            = f − f r²(k^φ)²
        double kRsq = f - f * r * r * kPh * kPh;
        assertTrue(kRsq > 0.0, "inbound photon needs k^r² > 0; got " + kRsq);
        double kR  = -Math.sqrt(kRsq);          // inbound

        double[] y = { 0.0, r, Math.PI / 2.0, 0.0, kT, kR, 0.0, kPh };

        DormandPrince45 dp = new DormandPrince45();
        double[] yNext = new double[8];
        double h = 0.1;
        double lambdaTotal = 200.0;
        double lambda = 0.0;
        double rCapture = m.horizonRadius() + 0.01;
        double rEscape = 2000.0;
        int accepted = 0;
        int guard = 0;

        while (lambda < lambdaTotal && guard < 1_000_000) {
            double hTry = Math.min(h, lambdaTotal - lambda);
            if (hTry <= 0) break;
            AdaptiveIntegrator.StepStatus s =
                    dp.adaptiveStep(m, y, hTry, atol, rtol, yNext);
            if (s.accepted()) {
                System.arraycopy(yNext, 0, y, 0, 8);
                lambda += hTry;
                accepted++;
                if (y[1] < rCapture || y[1] > rEscape) break;  // done early
            }
            h = s.hNext();
            guard++;
        }
        return accepted;
    }

    // ------------------------------------------------------------------
    // B2 (e): scalar ODE regression (thin test-only driver)
    // ------------------------------------------------------------------

    @Test
    void scalarExponentialDecayToExp100() {
        // dy/dt = -y, y(0) = 1, exact y(t) = exp(-t). Integrate to t=100
        // with atol=rtol=1e-10; final |y - exp(-100)| < 1e-8. The scalar
        // driver uses the same Butcher constants as DormandPrince45, so a
        // correct tableau here transfers to the main integrator.
        ScalarDP45 sdp = new ScalarDP45();
        DoubleUnaryOperator neg = y -> -y;

        double t = 0.0, y = 1.0, h = 0.1;
        double atol = 1e-10, rtol = 1e-10;
        int guard = 0;
        while (t < 100.0 && guard < 1_000_000) {
            double hTry = Math.min(h, 100.0 - t);
            if (hTry <= 0) break;
            double[] r = sdp.attempt(neg, y, hTry, atol, rtol);
            boolean accepted = r[0] > 0.5;
            if (accepted) {
                t += hTry;
                y = r[2];
            }
            h = r[1];
            guard++;
        }

        double exact = Math.exp(-100.0);
        double err = Math.abs(y - exact);
        System.out.println("[B2e] y(100)=" + y + ", exp(-100)=" + exact
                + ", |err|=" + err);
        assertTrue(err < 1e-8, "final error " + err);
    }

    /**
     * Scalar Dormand-Prince 5(4) step, reusing the Butcher constants from
     * {@link DormandPrince45}. Test-only driver: verifies the numerical
     * method (tableau + PI controller) at scalar dimension without
     * running the full 8-component state-vector machinery. Not thread-safe.
     */
    private static final class ScalarDP45 {
        private double k7;
        private double errPrev = 1.0;
        private boolean fsalValid = false;

        /** Returns {accepted (1/0), hNext, yNew} packed in a 3-element array. */
        double[] attempt(DoubleUnaryOperator f, double y, double h,
                         double atol, double rtol) {
            double k1 = fsalValid ? this.k7 : f.applyAsDouble(y);
            double k2 = f.applyAsDouble(y + h * DormandPrince45.A21 * k1);
            double k3 = f.applyAsDouble(y + h * (DormandPrince45.A31 * k1
                                                + DormandPrince45.A32 * k2));
            double k4 = f.applyAsDouble(y + h * (DormandPrince45.A41 * k1
                                                + DormandPrince45.A42 * k2
                                                + DormandPrince45.A43 * k3));
            double k5 = f.applyAsDouble(y + h * (DormandPrince45.A51 * k1
                                                + DormandPrince45.A52 * k2
                                                + DormandPrince45.A53 * k3
                                                + DormandPrince45.A54 * k4));
            double k6 = f.applyAsDouble(y + h * (DormandPrince45.A61 * k1
                                                + DormandPrince45.A62 * k2
                                                + DormandPrince45.A63 * k3
                                                + DormandPrince45.A64 * k4
                                                + DormandPrince45.A65 * k5));
            double y5 = y + h * (DormandPrince45.A71 * k1
                                + DormandPrince45.A73 * k3
                                + DormandPrince45.A74 * k4
                                + DormandPrince45.A75 * k5
                                + DormandPrince45.A76 * k6);
            double k7new = f.applyAsDouble(y5);

            double yErr = h * (DormandPrince45.E1 * k1
                             + DormandPrince45.E3 * k3
                             + DormandPrince45.E4 * k4
                             + DormandPrince45.E5 * k5
                             + DormandPrince45.E6 * k6
                             + DormandPrince45.E7 * k7new);
            double sc = atol + rtol * Math.abs(y);
            double err = Math.max(Math.abs(yErr / sc), 1e-10);

            boolean accepted = err <= 1.0;
            double factor;
            if (accepted) {
                factor = DormandPrince45.SAFETY
                        * Math.pow(err,     -DormandPrince45.ALPHA)
                        * Math.pow(errPrev,  DormandPrince45.BETA);
                factor = Math.max(DormandPrince45.MIN_GROWTH,
                          Math.min(DormandPrince45.MAX_GROWTH, factor));
                this.k7 = k7new;
                this.fsalValid = true;
                this.errPrev = err;
            } else {
                factor = Math.min(1.0,
                        DormandPrince45.SAFETY * Math.pow(err, -1.0 / 5.0));
                factor = Math.max(DormandPrince45.MIN_GROWTH, factor);
                this.fsalValid = false;
            }
            return new double[] { accepted ? 1.0 : 0.0, h * factor,
                                  accepted ? y5 : y };
        }
    }
}
