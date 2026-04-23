package com.pranav.grrt.metric;

import com.pranav.grrt.integrator.RK4;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

class KerrMetricTest {

    private static final double ALG_TOL = 1e-12;  // algebraic identities
    private static final double SCH_TOL = 1e-12;  // a=0 reduction to Schwarzschild

    private static double[] pos(double r, double theta) {
        return new double[] { 0.0, r, theta, 0.0 };
    }

    // ------------------------------------------------------------------
    // Constructor validation
    // ------------------------------------------------------------------

    @Test
    void constructorRejectsNonPositiveMass() {
        assertThrows(IllegalArgumentException.class, () -> new KerrMetric(0.0, 0.5));
        assertThrows(IllegalArgumentException.class, () -> new KerrMetric(-1.0, 0.5));
        assertThrows(IllegalArgumentException.class, () -> new KerrMetric(Double.NaN, 0.5));
        assertThrows(IllegalArgumentException.class,
                () -> new KerrMetric(Double.POSITIVE_INFINITY, 0.5));
    }

    @Test
    void constructorRejectsNonFiniteSpin() {
        assertThrows(IllegalArgumentException.class, () -> new KerrMetric(1.0, Double.NaN));
        assertThrows(IllegalArgumentException.class,
                () -> new KerrMetric(1.0, Double.POSITIVE_INFINITY));
        assertThrows(IllegalArgumentException.class,
                () -> new KerrMetric(1.0, Double.NEGATIVE_INFINITY));
    }

    @Test
    void constructorRejectsExtremalAndSuperextremalSpin() {
        // |a| < M is strict; |a| = M is extremal, rejected.
        assertThrows(IllegalArgumentException.class, () -> new KerrMetric(1.0, 1.0));
        assertThrows(IllegalArgumentException.class, () -> new KerrMetric(1.0, -1.0));
        assertThrows(IllegalArgumentException.class, () -> new KerrMetric(1.0, 1.5));
        assertThrows(IllegalArgumentException.class, () -> new KerrMetric(2.0, 2.5));
    }

    @Test
    void constructorAcceptsSubextremal() {
        assertDoesNotThrow(() -> new KerrMetric(1.0, 0.0));
        assertDoesNotThrow(() -> new KerrMetric(1.0, 0.999));
        assertDoesNotThrow(() -> new KerrMetric(1.0, -0.999));
        assertDoesNotThrow(() -> new KerrMetric(2.5, 1.3));
    }

    // ------------------------------------------------------------------
    // Accessors and landmarks
    // ------------------------------------------------------------------

    @Test
    void massSpinAndHorizonAccessors() {
        KerrMetric k = new KerrMetric(2.0, 1.2);
        assertEquals(2.0, k.mass(), ALG_TOL);
        assertEquals(1.2, k.spin(), ALG_TOL);
        // r_+ = M + sqrt(M² - a²) = 2 + sqrt(4 - 1.44) = 2 + sqrt(2.56) = 2 + 1.6 = 3.6
        assertEquals(3.6, k.horizonRadius(), ALG_TOL);
    }

    @Test
    void horizonAtZeroSpinMatchesSchwarzschild() {
        KerrMetric k = new KerrMetric(1.0, 0.0);
        assertEquals(2.0, k.horizonRadius(), ALG_TOL);
    }

    @Test
    void horizonAtHighSpin() {
        KerrMetric k = new KerrMetric(1.0, 0.9);
        // r_+ = 1 + sqrt(1 - 0.81) = 1 + sqrt(0.19)
        assertEquals(1.0 + Math.sqrt(0.19), k.horizonRadius(), ALG_TOL);
    }

    @Test
    void equatorialErgosphereIsTwoM() {
        assertEquals(4.0, new KerrMetric(2.0, 0.5).ergosphereRadiusEquatorial(), ALG_TOL);
        assertEquals(2.0, new KerrMetric(1.0, 0.9).ergosphereRadiusEquatorial(), ALG_TOL);
    }

    @Test
    void iscoAtZeroSpinIsSixM() {
        // Schwarzschild limit: r_isco = 6M regardless of direction.
        KerrMetric k = new KerrMetric(1.0, 0.0);
        assertEquals(6.0, k.iscoRadius(true),  ALG_TOL);
        assertEquals(6.0, k.iscoRadius(false), ALG_TOL);

        KerrMetric k2 = new KerrMetric(2.5, 0.0);
        assertEquals(15.0, k2.iscoRadius(true),  ALG_TOL);
        assertEquals(15.0, k2.iscoRadius(false), ALG_TOL);
    }

    @Test
    void iscoAtHighSpinBardeenPressTeukolsky() {
        // Bardeen, Press, Teukolsky (1972), published values at M = 1:
        //   a = 0.9 M: prograde r_isco ≈ 2.3209 M, retrograde r_isco ≈ 8.7174 M.
        // Tolerance 1e-4 matches the tabulated precision.
        KerrMetric k = new KerrMetric(1.0, 0.9);
        assertEquals(2.3209, k.iscoRadius(true),  1e-4, "prograde ISCO at a=0.9M");
        assertEquals(8.7174, k.iscoRadius(false), 1e-4, "retrograde ISCO at a=0.9M");
    }

    @Test
    void iscoMonotonicInProgradeSpin() {
        // Prograde ISCO decreases monotonically from 6M (a=0) toward M as
        // spin increases, and never drops below M (extremal limit).
        double r0  = new KerrMetric(1.0, 0.0  ).iscoRadius(true);
        double r5  = new KerrMetric(1.0, 0.5  ).iscoRadius(true);
        double r9  = new KerrMetric(1.0, 0.9  ).iscoRadius(true);
        double r99 = new KerrMetric(1.0, 0.99 ).iscoRadius(true);
        double r999= new KerrMetric(1.0, 0.999).iscoRadius(true);
        assertTrue(r0 > r5,  "r_isco(0) > r_isco(0.5)");
        assertTrue(r5 > r9,  "r_isco(0.5) > r_isco(0.9)");
        assertTrue(r9 > r99, "r_isco(0.9) > r_isco(0.99)");
        assertTrue(r99 > r999, "r_isco(0.99) > r_isco(0.999)");
        assertTrue(r999 > 1.0, "r_isco always > M for subextremal Kerr");
    }

    // ------------------------------------------------------------------
    // Test (a): g · gInv = I at multiple (r, θ)
    // ------------------------------------------------------------------

    @Test
    void metricInverseIdentityAtMultiplePoints() {
        // Spec: {(10M, π/4), (4M, π/3), (2.5M, 1.2)} at a = 0.5 M.
        // All three lie outside r_+ = 1 + sqrt(0.75) ≈ 1.866 M.
        KerrMetric k = new KerrMetric(1.0, 0.5);
        double[][] points = {
                { 10.0, Math.PI / 4.0 },
                { 4.0,  Math.PI / 3.0 },
                { 2.5,  1.2           },
        };
        for (double[] p : points) {
            assertGtimesGinvIsIdentity(k, p[0], p[1]);
        }
    }

    @Test
    void metricInverseIdentityAtHighSpin() {
        // Also exercise a = 0.9 M at r = 3 M, θ = π/2 (equator, near photon sphere).
        KerrMetric k = new KerrMetric(1.0, 0.9);
        assertGtimesGinvIsIdentity(k, 3.0, Math.PI / 2.0);
        assertGtimesGinvIsIdentity(k, 10.0, 0.3);
    }

    private static void assertGtimesGinvIsIdentity(Metric m, double r, double theta) {
        double[] x = pos(r, theta);
        double[][] g  = m.g(x);
        double[][] gi = m.gInv(x);
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                double sum = 0.0;
                for (int k = 0; k < 4; k++) sum += g[i][k] * gi[k][j];
                double expected = (i == j) ? 1.0 : 0.0;
                assertEquals(expected, sum, ALG_TOL,
                        "g·gInv[%d][%d] at (r=%g, θ=%g)".formatted(i, j, r, theta));
            }
        }
    }

    // ------------------------------------------------------------------
    // Test (b): Christoffel symmetric in the lower indices
    // ------------------------------------------------------------------

    @Test
    void christoffelSymmetricInLowerIndices() {
        KerrMetric k = new KerrMetric(1.0, 0.7);
        double[][][] G = k.christoffel(pos(5.0, 0.8));
        for (int a = 0; a < 4; a++) {
            for (int i = 0; i < 4; i++) {
                for (int j = 0; j < 4; j++) {
                    assertEquals(G[a][i][j], G[a][j][i], 1e-14,
                            "Γ^%d_%d%d asymmetric".formatted(a, i, j));
                }
            }
        }
    }

    @Test
    void christoffelSymmetricAtAdditionalPoints() {
        // Sample a few more (a, r, θ) to exercise different sparsity patterns.
        double[][] cases = {
                { 0.0,  8.0, Math.PI / 3.0 },  // Schwarzschild limit
                { 0.5,  3.0, Math.PI / 2.0 }, // equatorial
                { 0.9,  2.0, 1.0           }, // near-horizon, off-equator
                { -0.5, 6.0, 2.4           }, // retrograde, θ > π/2
        };
        for (double[] c : cases) {
            KerrMetric k = new KerrMetric(1.0, c[0]);
            double[][][] G = k.christoffel(pos(c[1], c[2]));
            for (int a = 0; a < 4; a++) {
                for (int i = 0; i < 4; i++) {
                    for (int j = 0; j < 4; j++) {
                        assertEquals(G[a][i][j], G[a][j][i], 1e-14,
                                "Γ asymmetric at (a=%g, r=%g, θ=%g) [%d][%d][%d]"
                                        .formatted(c[0], c[1], c[2], a, i, j));
                    }
                }
            }
        }
    }

    // ------------------------------------------------------------------
    // Test (c): a → 0 limit matches Schwarzschild
    // ------------------------------------------------------------------

    @Test
    void zeroSpinMetricMatchesSchwarzschild() {
        // At a = 0 exactly, the Kerr formulas reduce algebraically to
        // Schwarzschild (A → r⁴, Σ → r², Δ → r² − 2Mr, g_tφ → 0). Check
        // elementwise to algebraic-identity precision.
        KerrMetric k = new KerrMetric(1.0, 0.0);
        SchwarzschildMetric s = new SchwarzschildMetric(1.0);

        double[][] points = {
                { 10.0, Math.PI / 4.0 },
                { 4.0,  Math.PI / 3.0 },
                { 2.5,  1.2           },
                { 100.0, 0.5          },
        };
        for (double[] p : points) {
            double[] x = pos(p[0], p[1]);
            double[][] gK = k.g(x);
            double[][] gS = s.g(x);
            double[][] giK = k.gInv(x);
            double[][] giS = s.gInv(x);
            for (int i = 0; i < 4; i++) {
                for (int j = 0; j < 4; j++) {
                    assertEquals(gS[i][j],  gK[i][j],  SCH_TOL,
                            "g[%d][%d] at (%g, %g)".formatted(i, j, p[0], p[1]));
                    assertEquals(giS[i][j], giK[i][j], SCH_TOL,
                            "gInv[%d][%d] at (%g, %g)".formatted(i, j, p[0], p[1]));
                }
            }
        }
    }

    @Test
    void zeroSpinChristoffelMatchesSchwarzschild() {
        KerrMetric k = new KerrMetric(1.0, 0.0);
        SchwarzschildMetric s = new SchwarzschildMetric(1.0);

        double[][] points = {
                { 10.0, Math.PI / 4.0 },
                { 6.0,  Math.PI / 3.0 },
                { 3.0,  Math.PI / 2.0 },
                { 100.0, 1.0          },
        };
        for (double[] p : points) {
            double[] x = pos(p[0], p[1]);
            double[][][] GK = k.christoffel(x);
            double[][][] GS = s.christoffel(x);
            for (int a = 0; a < 4; a++) {
                for (int i = 0; i < 4; i++) {
                    for (int j = 0; j < 4; j++) {
                        assertEquals(GS[a][i][j], GK[a][i][j], SCH_TOL,
                                "Γ^%d_%d%d at (%g, %g)"
                                        .formatted(a, i, j, p[0], p[1]));
                    }
                }
            }
        }
    }

    // ------------------------------------------------------------------
    // Test (e): optimized geodesicAcceleration matches default contraction
    // ------------------------------------------------------------------

    @Test
    void optimizedAccelerationMatchesDefaultContraction() {
        // Default Metric.geodesicAcceleration contracts Γ^μ_{αβ} k^α k^β
        // via the full 4×4×4 christoffel() tensor. The optimized override
        // fuses that contraction directly. They must agree to algebraic
        // precision on generic Kerr states.
        KerrMetric k = new KerrMetric(1.0, 0.6);

        // Wrapper that delegates everything except geodesicAcceleration to k,
        // forcing the interface-default implementation.
        Metric defaultImpl = new Metric() {
            @Override public double mass() { return k.mass(); }
            @Override public double horizonRadius() { return k.horizonRadius(); }
            @Override public double[][] g(double[] x) { return k.g(x); }
            @Override public double[][] gInv(double[] x) { return k.gInv(x); }
            @Override public double[][][] christoffel(double[] x) { return k.christoffel(x); }
            // intentionally NOT overriding geodesicAcceleration
        };

        java.util.Random rng = new java.util.Random(0x600D5EED);
        double[] aOpt = new double[4];
        double[] aDef = new double[4];
        double rPlus = k.horizonRadius();
        double worst = 0.0;

        for (int trial = 0; trial < 100; trial++) {
            // Random state well outside horizon, away from polar axis.
            double r  = rPlus + 0.5 + 20.0 * rng.nextDouble();    // [r_+ + 0.5, r_+ + 20.5]
            double th = 0.1 + (Math.PI - 0.2) * rng.nextDouble();  // (0.1, π − 0.1)
            double[] x = { 0.0, r, th, 2.0 * Math.PI * rng.nextDouble() };
            // Random momentum; magnitude unconstrained (Γ contraction is
            // bilinear in k, so any k tests the contraction correctly).
            double[] km = {
                    2.0 * rng.nextGaussian(),
                    1.0 * rng.nextGaussian(),
                    0.5 * rng.nextGaussian(),
                    0.5 * rng.nextGaussian(),
            };

            k.geodesicAcceleration(x, km, aOpt);
            defaultImpl.geodesicAcceleration(x, km, aDef);

            for (int mu = 0; mu < 4; mu++) {
                double d = Math.abs(aDef[mu] - aOpt[mu]);
                if (d > worst) worst = d;
                assertEquals(aDef[mu], aOpt[mu], 1e-12,
                        ("mismatch at μ=%d, trial %d, (r=%.3f, θ=%.3f)")
                                .formatted(mu, trial, r, th));
            }
        }
        System.out.println("[A2] max |a_opt - a_default| over 100 random states = " + worst);
    }

    @Test
    void optimizedAccelerationMatchesDefaultAtHighSpin() {
        // Spot-check near the photon sphere at high spin.
        KerrMetric k = new KerrMetric(1.0, 0.9);
        Metric defaultImpl = new Metric() {
            @Override public double mass() { return k.mass(); }
            @Override public double horizonRadius() { return k.horizonRadius(); }
            @Override public double[][] g(double[] x) { return k.g(x); }
            @Override public double[][] gInv(double[] x) { return k.gInv(x); }
            @Override public double[][][] christoffel(double[] x) { return k.christoffel(x); }
        };

        double[] x = { 0.0, 2.5, 1.1, 0.3 };
        double[] km = { 1.2, -0.4, 0.05, 0.07 };
        double[] aOpt = new double[4];
        double[] aDef = new double[4];
        k.geodesicAcceleration(x, km, aOpt);
        defaultImpl.geodesicAcceleration(x, km, aDef);
        for (int mu = 0; mu < 4; mu++) {
            assertEquals(aDef[mu], aOpt[mu], 1e-12, "μ=" + mu);
        }
    }

    // ------------------------------------------------------------------
    // Test (f): equatorial photon conserves E and L under RK4
    // ------------------------------------------------------------------

    @Test
    void equatorialPhotonConservesEnergyAndAngularMomentum() {
        // Kerr has two Killing vectors ∂_t and ∂_φ, which make E = -k_t and
        // L = k_φ constants along every geodesic. An equatorial null
        // geodesic (θ = π/2, k^θ = 0 initial) stays equatorial by symmetry,
        // so the integration exercises Γ^t, Γ^r, Γ^φ Christoffels including
        // all frame-dragging cross terms. RK4 at h = 0.01M for 10⁴ steps
        // (100 M affine range) should drift < 1e-8 on both constants.
        KerrMetric k = new KerrMetric(1.0, 0.5);
        RK4 rk = new RK4();

        // Initial: outgoing equatorial null photon with E = 1, L = 3
        // (units of M) at r = 10 M. Outgoing ensures no horizon plunge
        // over 100 M of affine parameter.
        double r0 = 10.0;
        double E = 1.0;
        double L = 3.0;
        double[] x0 = { 0.0, r0, Math.PI / 2.0, 0.0 };
        double[][] gi = k.gInv(x0);

        double kUpT   = gi[0][0] * (-E) + gi[0][3] * L;       // k^t
        double kUpPhi = gi[0][3] * (-E) + gi[3][3] * L;       // k^φ
        // Null condition (equator, k_θ=0): g^{μν} k_μ k_ν = 0 ⇒
        //   k_r² = -(g^{tt} E² − 2 g^{tφ} E L + g^{φφ} L²) / g^{rr}
        double rhs = gi[0][0] * E * E - 2.0 * gi[0][3] * E * L + gi[3][3] * L * L;
        double kLoRSq = -rhs / gi[1][1];
        assertTrue(kLoRSq > 0.0, "null condition needs positive k_r²; got " + kLoRSq);
        double kUpR = gi[1][1] * Math.sqrt(kLoRSq);           // k^r = g^{rr} k_r, outgoing

        double[] y = { x0[0], x0[1], x0[2], x0[3],
                       kUpT, kUpR, 0.0, kUpPhi };

        // Sanity: initial null-norm should be ~0.
        assertEquals(0.0, k.nullNorm(
                new double[] { y[0], y[1], y[2], y[3] },
                new double[] { y[4], y[5], y[6], y[7] }), 1e-12);

        double E0 = computeE(k, y);
        double L0 = computeL(k, y);
        assertEquals(E, E0, 1e-13, "E round-trip check");
        assertEquals(L, L0, 1e-13, "L round-trip check");

        double[] next = new double[8];
        double h = 0.01;
        int steps = 10_000;
        double maxDriftE = 0.0;
        double maxDriftL = 0.0;

        for (int step = 0; step < steps; step++) {
            rk.step(k, y, h, next);
            System.arraycopy(next, 0, y, 0, 8);
            double dE = Math.abs(computeE(k, y) - E0);
            double dL = Math.abs(computeL(k, y) - L0);
            if (dE > maxDriftE) maxDriftE = dE;
            if (dL > maxDriftL) maxDriftL = dL;
        }

        double finalE = computeE(k, y);
        double finalL = computeL(k, y);
        double relE = Math.abs(finalE - E0) / Math.abs(E0);
        double relL = Math.abs(finalL - L0) / Math.abs(L0);
        System.out.println("[A3] final |ΔE|/|E0| = " + relE);
        System.out.println("[A3] final |ΔL|/|L0| = " + relL);

        assertTrue(maxDriftE < 1e-8, "|ΔE|_max = " + maxDriftE);
        assertTrue(maxDriftL < 1e-8, "|ΔL|_max = " + maxDriftL);
    }

    private static double computeE(Metric m, double[] y) {
        // E = -k_t = -(g_tt k^t + g_tφ k^φ)
        double[][] g = m.g(new double[] { y[0], y[1], y[2], y[3] });
        return -(g[0][0] * y[4] + g[0][3] * y[7]);
    }

    private static double computeL(Metric m, double[] y) {
        // L = k_φ = g_tφ k^t + g_φφ k^φ
        double[][] g = m.g(new double[] { y[0], y[1], y[2], y[3] });
        return g[0][3] * y[4] + g[3][3] * y[7];
    }

    @Test
    void smallSpinOffDiagonalScalesLinearly() {
        // At small a, g_tφ = -2Mar s²/Σ should scale linearly in a.
        // This verifies the Kerr-specific (t,φ) coupling is not mistakenly zero.
        double r = 5.0, theta = Math.PI / 3.0;
        double s2 = Math.sin(theta) * Math.sin(theta);
        double c2 = Math.cos(theta) * Math.cos(theta);
        double M = 1.0;

        double a1 = 1e-4;
        double a2 = 2e-4;
        double gTF1 = new KerrMetric(M, a1).g(pos(r, theta))[0][3];
        double gTF2 = new KerrMetric(M, a2).g(pos(r, theta))[0][3];

        // Expected closed form.
        double expected1 = -2.0 * M * a1 * r * s2 / (r * r + a1 * a1 * c2);
        double expected2 = -2.0 * M * a2 * r * s2 / (r * r + a2 * a2 * c2);
        assertEquals(expected1, gTF1, 1e-14);
        assertEquals(expected2, gTF2, 1e-14);

        // Linearity: g_tφ(2a) / g_tφ(a) ≈ 2 at leading order.
        assertEquals(2.0, gTF2 / gTF1, 1e-6, "g_tφ should scale linearly in a");
    }
}
