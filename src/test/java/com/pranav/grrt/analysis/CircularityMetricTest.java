package com.pranav.grrt.analysis;

import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Phase 3C unit tests for {@link CircularityMetric}: closed-form
 * reductions on synthetic per-bin radius arrays. The pure-math
 * subset of {@code docs/phase-3c-plan.md} §7 ({@code CircularityMetric}
 * isolated table) lives here; the end-to-end image → radii → metric
 * round-trip is in {@code RingExtractorTest}.
 */
class CircularityMetricTest {

    private static final int N_BINS = 180;
    private static final double R0 = 5.5;

    // ------------------------------------------------------------------
    // Input validation
    // ------------------------------------------------------------------

    @Test
    void nullArrayThrowsIAE() {
        assertThrows(IllegalArgumentException.class,
                () -> CircularityMetric.compute(null));
    }

    @Test
    void emptyArrayThrowsIAE() {
        assertThrows(IllegalArgumentException.class,
                () -> CircularityMetric.compute(new double[0]));
    }

    @Test
    void infiniteEntryThrowsIAE() {
        double[] r = { 5.0, Double.POSITIVE_INFINITY, 5.0, 5.0 };
        assertThrows(IllegalArgumentException.class,
                () -> CircularityMetric.compute(r));
    }

    @Test
    void negativeInfiniteEntryThrowsIAE() {
        double[] r = { 5.0, 5.0, Double.NEGATIVE_INFINITY, 5.0 };
        assertThrows(IllegalArgumentException.class,
                () -> CircularityMetric.compute(r));
    }

    // ------------------------------------------------------------------
    // Edge cases: all-NaN, single-bin
    // ------------------------------------------------------------------

    @Test
    void allNaNReturnsZeroValidBinsAndNaNStatistics() {
        double[] r = new double[N_BINS];
        java.util.Arrays.fill(r, Double.NaN);

        var result = CircularityMetric.compute(r);

        assertEquals(0, result.validBins());
        assertTrue(Double.isNaN(result.meanR()));
        assertTrue(Double.isNaN(result.deltaRrms()));
        assertTrue(Double.isNaN(result.deltaRp2p()));
        assertTrue(Double.isNaN(result.fourierM1Amplitude()));
        assertTrue(Double.isNaN(result.fourierM1Phase()));
    }

    /**
     * Single-sample DFT is a degenerate but well-defined edge case.
     * {@code c₁ = R₀} for {@code N = 1}; not physically meaningful as
     * an asymmetry diagnostic but the implementation must not crash.
     */
    @Test
    void singleBinInputDoesNotCrashAndReturnsTrivialDFT() {
        double[] r = { R0 };

        var result = CircularityMetric.compute(r);

        assertEquals(1, result.validBins());
        assertEquals(R0, result.meanR(), 1e-15);
        assertEquals(0.0, result.deltaRrms(), 1e-15);
        assertEquals(0.0, result.deltaRp2p(), 1e-15);
        assertEquals(R0, result.fourierM1Amplitude(), 1e-12);
        assertEquals(0.0, result.fourierM1Phase(), 1e-12);
    }

    // ------------------------------------------------------------------
    // Constant input: all asymmetry metrics zero
    // ------------------------------------------------------------------

    @Test
    void constantInputHasZeroAsymmetry() {
        double[] r = new double[N_BINS];
        java.util.Arrays.fill(r, R0);

        var result = CircularityMetric.compute(r);

        assertEquals(N_BINS, result.validBins());
        assertEquals(R0, result.meanR(), 1e-15);
        assertEquals(0.0, result.deltaRrms(), 1e-15);
        assertEquals(0.0, result.deltaRp2p(), 1e-15);
        assertEquals(0.0, result.fourierM1Amplitude(), 1e-12);
        // Phase is mathematically undefined when amplitude is zero; in
        // floating point, summing 180 cos terms produces real and
        // imaginary parts at the ~1e-16 roundoff floor, so atan2 is
        // dominated by sign noise. We assert only that the amplitude
        // is at the roundoff floor; do not assert the phase.
    }

    // ------------------------------------------------------------------
    // Pure m = 1 sinusoid: closed-form RMS, p2p, |c₁|, arg(c₁)
    // ------------------------------------------------------------------

    /**
     * For {@code rᵢ = R₀ (1 + A cos(2π i / N))} (zero phase offset):
     * <pre>
     *   ⟨r⟩      = R₀
     *   δ_r_rms  = R₀ A / √2          (population SD over a full period)
     *   p2p      = 2 R₀ A             (exact: i=0 hits cos = +1, i=N/2 hits cos = −1)
     *   |c₁|     = R₀ A / 2
     *   arg(c₁)  = 0
     * </pre>
     * The phase test with non-zero {@code φ₀} is split out below to
     * keep the closed-form p2p exact (off-grid {@code φ₀} produces a
     * discrete p2p that differs from {@code 2 R₀ A} by a {@code O(1/N²)}
     * sampling error — real, not a bug, and orthogonal to the metric
     * being verified here).
     */
    @Test
    void pureM1ZeroPhaseHasClosedFormReductions() {
        double a1 = 0.05;
        double[] r = new double[N_BINS];
        for (int i = 0; i < N_BINS; i++) {
            r[i] = R0 * (1.0 + a1 * Math.cos(2.0 * Math.PI * i / N_BINS));
        }

        var result = CircularityMetric.compute(r);

        assertEquals(R0, result.meanR(), 1e-12);
        assertEquals(R0 * a1 / Math.sqrt(2.0), result.deltaRrms(), 1e-12);
        assertEquals(2.0 * R0 * a1, result.deltaRp2p(), 1e-12);
        assertEquals(R0 * a1 / 2.0, result.fourierM1Amplitude(), 1e-12);
        assertEquals(0.0, result.fourierM1Phase(), 1e-12);
    }

    /**
     * Independently verify the phase machinery using a non-zero
     * {@code φ₀} that lands on the {@code N = 180} sample grid
     * ({@code φ₀ = −π/3} aligns with sample {@code i = N/6 = 30}).
     * For grid-aligned phases, both p2p and |c₁| stay exact; we
     * additionally assert {@code arg(c₁) = φ₀}.
     */
    @Test
    void pureM1WithGridAlignedPhaseIsolatesPhase() {
        double a1 = 0.05;
        double phi0 = -Math.PI / 3.0;       // grid-aligned at i = 30 for N = 180
        double[] r = new double[N_BINS];
        for (int i = 0; i < N_BINS; i++) {
            r[i] = R0 * (1.0 + a1 * Math.cos(2.0 * Math.PI * i / N_BINS + phi0));
        }

        var result = CircularityMetric.compute(r);

        assertEquals(R0 * a1 / 2.0, result.fourierM1Amplitude(), 1e-12);
        assertEquals(phi0, result.fourierM1Phase(), 1e-12);
        // Grid-aligned phase: discrete extrema land exactly at cos = ±1
        assertEquals(2.0 * R0 * a1, result.deltaRp2p(), 1e-12);
    }

    /**
     * Off-grid phase {@code φ₀ = 0.7} produces a discrete p2p that
     * misses the continuous {@code 2 R₀ A} by an {@code O(1/N²)}
     * sampling error (~1e-6 at N=180). |c₁| and arg(c₁) remain exact
     * to floating-point precision because the DFT projects onto the
     * pure {@code k = 1} carrier regardless of phase alignment.
     */
    @Test
    void pureM1WithOffGridPhasePreservesAmplitudeAndPhase() {
        double a1 = 0.05;
        double phi0 = 0.7;                  // off-grid for N = 180
        double[] r = new double[N_BINS];
        for (int i = 0; i < N_BINS; i++) {
            r[i] = R0 * (1.0 + a1 * Math.cos(2.0 * Math.PI * i / N_BINS + phi0));
        }

        var result = CircularityMetric.compute(r);

        assertEquals(R0, result.meanR(), 1e-12);
        assertEquals(R0 * a1 / Math.sqrt(2.0), result.deltaRrms(), 1e-12);
        assertEquals(R0 * a1 / 2.0, result.fourierM1Amplitude(), 1e-12);
        assertEquals(phi0, result.fourierM1Phase(), 1e-12);
        // Off-grid p2p: bound the discretization error rather than
        // asserting the continuous value
        double p2pSamplingError = Math.abs(2.0 * R0 * a1 - result.deltaRp2p());
        assertTrue(p2pSamplingError < 1.0e-5,
                "off-grid p2p sampling error " + p2pSamplingError
                        + " exceeds expected O(1/N²) bound at N=" + N_BINS);
    }

    // ------------------------------------------------------------------
    // Pure m = 2 sinusoid: m = 1 component vanishes by orthogonality
    // ------------------------------------------------------------------

    @Test
    void pureM2HasZeroM1ComponentButNonZeroRms() {
        double a2 = 0.05;
        double[] r = new double[N_BINS];
        for (int i = 0; i < N_BINS; i++) {
            r[i] = R0 * (1.0 + a2 * Math.cos(4.0 * Math.PI * i / N_BINS));
        }

        var result = CircularityMetric.compute(r);

        assertEquals(R0, result.meanR(), 1e-12);
        assertEquals(R0 * a2 / Math.sqrt(2.0), result.deltaRrms(), 1e-12);
        assertEquals(2.0 * R0 * a2, result.deltaRp2p(), 1e-12);
        assertEquals(0.0, result.fourierM1Amplitude(), 1e-12);
    }

    // ------------------------------------------------------------------
    // Two-mode mix m = 1 + m = 2: Parseval RMS, m = 1 isolation
    // ------------------------------------------------------------------

    @Test
    void twoModeMixIsolatesM1FromM2() {
        double a1 = 0.05;
        double a2 = 0.03;
        double[] r = new double[N_BINS];
        for (int i = 0; i < N_BINS; i++) {
            double phi = 2.0 * Math.PI * i / N_BINS;
            r[i] = R0 * (1.0 + a1 * Math.cos(phi) + a2 * Math.cos(2.0 * phi));
        }

        var result = CircularityMetric.compute(r);

        assertEquals(R0, result.meanR(), 1e-12);
        // Parseval: σ² = R₀² (A₁² + A₂²) / 2  →  σ = R₀ √(A₁² + A₂²) / √2
        double rmsExpected = R0 * Math.hypot(a1, a2) / Math.sqrt(2.0);
        assertEquals(rmsExpected, result.deltaRrms(), 1e-12);
        // m = 1 amplitude isolates A₁ regardless of A₂
        assertEquals(R0 * a1 / 2.0, result.fourierM1Amplitude(), 1e-12);
        // No phase offset on the m = 1 component
        assertEquals(0.0, result.fourierM1Phase(), 1e-12);
    }

    // ------------------------------------------------------------------
    // NaN-bin handling: excluded from mean / RMS / p2p, replaced by
    // mean in the Fourier sum
    // ------------------------------------------------------------------

    @Test
    void nanBinsExcludedFromMeanRmsAndP2P() {
        double[] r = new double[N_BINS];
        java.util.Arrays.fill(r, R0);
        int[] nanPositions = { 3, 17, 42, 51, 78, 99, 120, 134, 156, 178 };
        for (int p : nanPositions) {
            r[p] = Double.NaN;
        }

        var result = CircularityMetric.compute(r);

        assertEquals(N_BINS - nanPositions.length, result.validBins());
        assertEquals(R0, result.meanR(), 1e-15);
        assertEquals(0.0, result.deltaRrms(), 1e-15);
        assertEquals(0.0, result.deltaRp2p(), 1e-15);
        // NaN bins replaced by mean = R₀ → all 180 effective bins are
        // R₀, so the Fourier m = 1 amplitude is identically zero.
        assertEquals(0.0, result.fourierM1Amplitude(), 1e-12);
    }
}
