package com.pranav.grrt.validation;

import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

/**
 * Unit tests on the {@link BardeenShadow} oracle itself, run before the
 * oracle is trusted to judge the rendered shadow in C3/C4.
 */
class BardeenShadowTest {

    // ------------------------------------------------------------------
    // C2 oracle sanity: Schwarzschild limit in the Kerr-coord formula
    // ------------------------------------------------------------------

    @Test
    void dHApproachesSchwarzschildLimitAtZeroSpin() {
        // The Kerr formula for xi has division by a, so a = 0 exactly is
        // 0/0. Use a = 1e-5 M as the numerically stable "a → 0" proxy.
        //
        // Error budget at small a:
        //   analytic:    O(a²)           — expansion xi = ±3√3·M − 2a + O(a²)
        //                                   gives D_H = 6√3·M + O(a²)
        //   cancellation: O(ε_mach / a)  — r³ − 3Mr² subtracts two numbers
        //                                   near 27M³, loss of precision
        //                                   scales inversely with a
        // Balance at a ≈ ε_mach^(1/3) ≈ 5e-6. Chose a = 1e-5 M to hit
        // the 1e-9 tolerance requested in the spec (both terms ≈ 1e-10
        // at this a).
        double M = 1.0;
        double a = 1e-5;

        BardeenShadow.ShadowExtent e = BardeenShadow.horizontalExtent(
                M, a, Math.PI / 2.0, 10001);
        double expected = 6.0 * Math.sqrt(3.0);          // 10.39230484541...

        System.out.println("[C2 a→0]  D_H        = " + String.format("%.15f M", e.dH()));
        System.out.println("[C2 a→0]  expected   = " + String.format("%.15f M", expected));
        System.out.println("[C2 a→0]  |Δ|        = " + String.format("%.3e M",
                Math.abs(e.dH() - expected)));

        assertEquals(expected, e.dH(), 1e-9,
                "D_H in the a→0 limit must equal 2·3√3·M to 1e-9 M");
    }

    // ------------------------------------------------------------------
    // C2 print run: a = 0.9 M, i = 90°
    // ------------------------------------------------------------------

    @Test
    void dHAtHighSpinReportValues() {
        // Spec: scan r_ph uniformly with 10001 samples across
        // [r_ph_pro, r_ph_retro] at a = 0.9 M, i = π/2. Print values to
        // 10 sig figs.
        double M = 1.0;
        double a = 0.9;
        double inclination = Math.PI / 2.0;

        BardeenShadow.ShadowExtent e = BardeenShadow.horizontalExtent(
                M, a, inclination, 10001);

        System.out.printf("[C2 a=0.9]  r_ph_pro   = %.10f M%n", e.rPhPro());
        System.out.printf("[C2 a=0.9]  r_ph_retro = %.10f M%n", e.rPhRetro());
        System.out.printf("[C2 a=0.9]  alpha_min  = %+.10f M%n", e.alphaMin());
        System.out.printf("[C2 a=0.9]  alpha_max  = %+.10f M%n", e.alphaMax());
        System.out.printf("[C2 a=0.9]  D_H        = %.10f M%n", e.dH());

        // Sign convention check (spec): prograde at negative alpha.
        assertTrue(e.alphaMin() < 0,
                "alpha_min should be negative (prograde at -α); got " + e.alphaMin());
        assertTrue(e.alphaMax() > 0,
                "alpha_max should be positive (retrograde at +α); got " + e.alphaMax());
        assertTrue(e.dH() > 0, "D_H must be positive");
    }
}
