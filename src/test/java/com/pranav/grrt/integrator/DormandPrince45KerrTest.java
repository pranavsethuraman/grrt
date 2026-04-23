package com.pranav.grrt.integrator;

import com.pranav.grrt.metric.KerrMetric;
import com.pranav.grrt.metric.Metric;
import org.junit.jupiter.api.Test;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import static org.junit.jupiter.api.Assertions.*;

/**
 * B3 tests: DP45 exercised on Kerr geodesics.
 *
 * <p>Kept separate from {@code DormandPrince45Test} so the Schwarzschild-only
 * B2 gate tests stay independent of {@link KerrMetric}.
 */
class DormandPrince45KerrTest {

    // ------------------------------------------------------------------
    // B3 (c): Kerr a=0.9M sub-critical plunge — adaptive h by radial zone
    // ------------------------------------------------------------------

    @Test
    void kerrPlungingPhotonStepSizeByZone() {
        // a=0.9M, equatorial, inbound with b = L/E = 2.5 M. At a=0.9 the
        // prograde critical impact parameter is b_crit_pro ≈ 2.848 M
        // (Bardeen 1973), so b=2.5 is sub-critical and the photon plunges
        // to r_+ = 1.436 M. Verified algebraically:
        //   Σ²(k^r)² = r(r³ − 5.44 r + 5.12)
        // The cubic has no positive real roots — no turning point above r=0.
        KerrMetric m = new KerrMetric(1.0, 0.9);
        double r0 = 100.0;
        double E  = 1.0;
        double L  = 2.5 * E;

        double[] y0 = equatorialInboundNull(m, r0, E, L);
        assertTrue(y0[5] < 0, "k^r must be inbound (negative): " + y0[5]);
        double n0 = m.nullNorm(slice(y0, 0, 4), slice(y0, 4, 4));
        assertEquals(0.0, n0, 1e-10, "initial state must be null");

        DormandPrince45 dp = new DormandPrince45();
        double[] y = y0.clone();
        double[] yNext = new double[8];
        double h = 1.0;
        double atol = 1e-8, rtol = 1e-8;
        double rStop = 1.01 * m.horizonRadius();
        int maxSteps = 100_000;

        List<Double> hFar  = new ArrayList<>();   // r >= 50 M
        List<Double> hMid  = new ArrayList<>();   // 10 <= r < 50 M
        List<Double> hNear = new ArrayList<>();   // r <  5 M
        int accepted = 0, rejected = 0;

        for (int step = 0; step < maxSteps && y[1] > rStop; step++) {
            AdaptiveIntegrator.StepStatus s =
                    dp.adaptiveStep(m, y, h, atol, rtol, yNext);
            if (s.accepted()) {
                double hMag = Math.abs(h);
                double rMid = 0.5 * (y[1] + yNext[1]);
                System.arraycopy(yNext, 0, y, 0, 8);
                accepted++;
                if (rMid >= 50.0)      hFar.add(hMag);
                else if (rMid >= 10.0) hMid.add(hMag);
                else if (rMid < 5.0)   hNear.add(hMag);
                // 5 <= rMid < 10 intentionally unbinned per the spec
            } else {
                rejected++;
            }
            h = s.hNext();
        }

        double medFar  = median(hFar);
        double medMid  = median(hMid);
        double medNear = median(hNear);

        System.out.println("[B3c] r_+ = " + m.horizonRadius());
        System.out.println("[B3c] accepted=" + accepted
                + ", rejected=" + rejected + ", final r=" + y[1]);
        System.out.println("[B3c] far  (r>=50M): n=" + hFar.size()
                + ", median h=" + medFar);
        System.out.println("[B3c] mid  (10<=r<50M): n=" + hMid.size()
                + ", median h=" + medMid);
        System.out.println("[B3c] near (r<5M): n=" + hNear.size()
                + ", median h=" + medNear);

        assertFalse(hFar.isEmpty(),  "far zone has no accepted steps");
        assertFalse(hMid.isEmpty(),  "mid zone has no accepted steps");
        assertFalse(hNear.isEmpty(), "near zone has no accepted steps");
        assertTrue(medFar  > 1.0, "median far h = " + medFar  + " M (need > 1.0)");
        assertTrue(medNear < 0.1, "median near h = " + medNear + " M (need < 0.1)");
    }

    // ------------------------------------------------------------------
    // B3 (d): long deflecting Kerr integration — conservation of E and L
    // ------------------------------------------------------------------

    @Test
    void kerrDeflectingPhotonConservationUnderLongIntegration() {
        // a=0.5M, equatorial, inbound from r=200 M with L = 10 M·E
        // (well above b_crit ≈ 4.23 M prograde). atol=rtol=1e-10.
        // Integrate until r exits [1.5 M, 500 M] or 5e4 accepted steps,
        // whichever first. Verified algebraically:
        //   Σ²(k^r)² = r(r³ − 99.75 r + 180.5)
        // has outer root r ≈ 8.91, so the photon turns around at
        // r_min ≈ 8.91 M and escapes to r > 500.
        //
        // Bounds rationale:
        //   At rtol=atol=1e-10, per-step position error in r at r~200M
        //   is ~1e-8 absolute. L = g_φφ(r)·k^φ reconstructs L from
        //   position, so dL from position drift per step is
        //   (dg_φφ/dr)·k^φ·dr ~ 1e-9 absolute. Over √N random-walk
        //   accumulation for N=302 accepted steps, expected
        //   |dL| ~ 1.7e-8 abs = 1.7e-9 rel; observed ~3.6e-9 is at that
        //   random-walk floor. The 1e-8 bound gives 3× headroom for
        //   run-to-run variation while still catching any bug that
        //   produces systematic (non-random-walk) drift.
        //
        //   E conservation is tighter because g_tt ≈ -1 makes E ≈ -k^t,
        //   and k^t's floating-point error is dominated by rtol·1 = 1e-10
        //   per step with minimal amplification — hence the 1e-10 bound
        //   on |dE/E|.
        KerrMetric m = new KerrMetric(1.0, 0.5);
        double r0 = 200.0;
        double E  = 1.0;
        double L  = 10.0 * E;

        double[] y0 = equatorialInboundNull(m, r0, E, L);
        assertTrue(y0[5] < 0, "k^r must be inbound");
        double n0 = m.nullNorm(slice(y0, 0, 4), slice(y0, 4, 4));
        assertEquals(0.0, n0, 1e-10, "initial state must be null");

        double E0 = killingE(m, y0);
        double L0 = killingL(m, y0);
        assertEquals(E, E0, 1e-13, "E round-trip");
        assertEquals(L, L0, 1e-13, "L round-trip");

        DormandPrince45 dp = new DormandPrince45();
        double[] y = y0.clone();
        double[] yNext = new double[8];
        double h = 1.0;
        double atol = 1e-10, rtol = 1e-10;
        double rLo = 1.5, rHi = 500.0;
        int maxAccepted = 50_000;
        int accepted = 0, rejected = 0;
        double rMin = y[1];
        double lambda = 0.0;

        while (accepted < maxAccepted && y[1] >= rLo && y[1] <= rHi) {
            AdaptiveIntegrator.StepStatus s =
                    dp.adaptiveStep(m, y, h, atol, rtol, yNext);
            if (s.accepted()) {
                lambda += h;
                System.arraycopy(yNext, 0, y, 0, 8);
                accepted++;
                if (y[1] < rMin) rMin = y[1];
            } else {
                rejected++;
            }
            h = s.hNext();
        }

        double Ef   = killingE(m, y);
        double Lf   = killingL(m, y);
        double relE = Math.abs((Ef - E0) / E0);
        double relL = Math.abs((Lf - L0) / L0);

        System.out.println("[B3d] final |dE/E| = " + relE);
        System.out.println("[B3d] final |dL/L| = " + relL);
        System.out.println("[B3d] r_min reached = " + rMin);
        System.out.println("[B3d] accepted = " + accepted
                + ", rejected = " + rejected);
        System.out.println("[B3d] elapsed λ = " + lambda);
        System.out.println("[B3d] final r = " + y[1]);

        assertTrue(relE < 1e-10, "|dE/E| = " + relE);
        assertTrue(relL < 1e-8,  "|dL/L| = " + relL);
    }

    // ------------------------------------------------------------------
    // Helpers
    // ------------------------------------------------------------------

    /**
     * Build an equatorial inbound null 8-state with conserved energy
     * E = -k_t and angular momentum L = k_φ at radius r. Uses the mixed
     * null condition g_μν k^μ k^ν = 0 (k^θ = 0) to solve for (k^r)²,
     * then takes the negative root for inbound motion.
     */
    private static double[] equatorialInboundNull(
            Metric m, double r, double E, double L) {
        double[] x  = { 0.0, r, Math.PI / 2.0, 0.0 };
        double[][] gmn = m.g(x);
        double[][] gi  = m.gInv(x);
        double kT  = gi[0][0] * (-E) + gi[0][3] * L;
        double kPh = gi[0][3] * (-E) + gi[3][3] * L;
        double num = gmn[0][0] * kT * kT
                   + 2.0 * gmn[0][3] * kT * kPh
                   + gmn[3][3] * kPh * kPh;
        double kRsq = -num / gmn[1][1];
        if (!(kRsq > 0.0)) {
            throw new IllegalStateException(
                    "no real k^r at r=" + r + " (L=" + L + "): (k^r)² = " + kRsq);
        }
        double kR = -Math.sqrt(kRsq);  // inbound
        return new double[] { 0.0, r, Math.PI / 2.0, 0.0, kT, kR, 0.0, kPh };
    }

    private static double killingE(Metric m, double[] y) {
        double[][] g = m.g(new double[] { y[0], y[1], y[2], y[3] });
        return -(g[0][0] * y[4] + g[0][3] * y[7]);
    }

    private static double killingL(Metric m, double[] y) {
        double[][] g = m.g(new double[] { y[0], y[1], y[2], y[3] });
        return g[0][3] * y[4] + g[3][3] * y[7];
    }

    private static double median(List<Double> list) {
        if (list.isEmpty()) return Double.NaN;
        List<Double> s = new ArrayList<>(list);
        Collections.sort(s);
        int n = s.size();
        return (n % 2 == 1) ? s.get(n / 2)
                            : 0.5 * (s.get(n / 2 - 1) + s.get(n / 2));
    }

    private static double[] slice(double[] a, int from, int len) {
        double[] out = new double[len];
        System.arraycopy(a, from, out, 0, len);
        return out;
    }
}
