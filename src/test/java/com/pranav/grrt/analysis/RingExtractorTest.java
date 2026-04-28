package com.pranav.grrt.analysis;

import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Phase 3C unit tests for {@link RingExtractor}: image → per-bin
 * radius extraction, multi-peak diagnostic, and the
 * "outermost peak wins" rule. End-to-end ring + reduction validation
 * against synthetic geometric primitives covers Gates 1, 2, and 6
 * from {@code docs/phase-3c-plan.md} §7. Sweep-orchestrator gates
 * (3, 4, 7) are slow-tagged and live in {@code EpsilonSweepIT}
 * (Phase 3C.2, not yet committed).
 *
 * <p><b>Tolerance note for Gate 1.</b> The plan's Gate 1 originally
 * specified {@code δ_r/⟨r⟩ &lt; 1e-10} for a synthetic perfect
 * circle. With the current naive integer-pixel peak finder (no
 * sub-pixel parabolic interpolation) on a 512² raster, per-bin
 * radius is pixel-quantised at ±0.5 px, capping the relative
 * dispersion at the {@code 1 / (2 √12 · r_px)} ≈ 0.15 % floor at
 * {@code r = 5.5 M, r_px ≈ 94}. Gate 1 here asserts
 * {@code δ_r / ⟨r⟩ &lt; 0.02} (2 % relative) which is the
 * realistic-and-still-meaningful tolerance for V1; sub-pixel
 * interpolation is a deferred refinement (per
 * {@code docs/phase-3c-plan.md} §13).
 */
class RingExtractorTest {

    private static final int W = 512;
    private static final int H = 512;
    private static final double FIELD_HALF_M = 15.0;
    private static final double PIXEL_TO_M = 2.0 * FIELD_HALF_M / W;        // ≈ 0.0586 M / px
    private static final double CX = 256.0;                                 // half-pixel-grid centre of W=512 image
    private static final double CY = 256.0;
    private static final int N_BINS = 180;

    // ------------------------------------------------------------------
    // Constructor and image-input validation
    // ------------------------------------------------------------------

    @Test
    void constructorRejectsBadNumBins() {
        assertThrows(IllegalArgumentException.class,
                () -> new RingExtractor(0, PIXEL_TO_M, CX, CY));
        assertThrows(IllegalArgumentException.class,
                () -> new RingExtractor(-3, PIXEL_TO_M, CX, CY));
    }

    @Test
    void constructorRejectsBadPixelToM() {
        assertThrows(IllegalArgumentException.class,
                () -> new RingExtractor(N_BINS, 0.0, CX, CY));
        assertThrows(IllegalArgumentException.class,
                () -> new RingExtractor(N_BINS, -1.0, CX, CY));
        assertThrows(IllegalArgumentException.class,
                () -> new RingExtractor(N_BINS, Double.POSITIVE_INFINITY, CX, CY));
        assertThrows(IllegalArgumentException.class,
                () -> new RingExtractor(N_BINS, Double.NaN, CX, CY));
    }

    @Test
    void constructorRejectsNonFiniteCentroid() {
        assertThrows(IllegalArgumentException.class,
                () -> new RingExtractor(N_BINS, PIXEL_TO_M, Double.NaN, CY));
        assertThrows(IllegalArgumentException.class,
                () -> new RingExtractor(N_BINS, PIXEL_TO_M, CX, Double.POSITIVE_INFINITY));
    }

    @Test
    void extractRejectsNullOrMisshapedImage() {
        RingExtractor ext = new RingExtractor(N_BINS, PIXEL_TO_M, CX, CY);
        assertThrows(IllegalArgumentException.class, () -> ext.extract(null, W, H));
        assertThrows(IllegalArgumentException.class, () -> ext.extract(new float[H][W], -1, H));
        assertThrows(IllegalArgumentException.class, () -> ext.extract(new float[H][W], W, 0));
        assertThrows(IllegalArgumentException.class, () -> ext.extract(new float[H - 1][W], W, H));
        float[][] mixed = new float[H][W];
        mixed[5] = new float[W - 1];
        assertThrows(IllegalArgumentException.class, () -> ext.extract(mixed, W, H));
    }

    @Test
    void multipeakBinCountRejectsBadThreshold() {
        RingExtractor ext = new RingExtractor(N_BINS, PIXEL_TO_M, CX, CY);
        float[][] img = new float[H][W];
        assertThrows(IllegalArgumentException.class,
                () -> ext.multipeakBinCount(img, W, H, 0.0));
        assertThrows(IllegalArgumentException.class,
                () -> ext.multipeakBinCount(img, W, H, -1.0));
        assertThrows(IllegalArgumentException.class,
                () -> ext.multipeakBinCount(img, W, H, Double.POSITIVE_INFINITY));
    }

    // ------------------------------------------------------------------
    // Empty-image semantics
    // ------------------------------------------------------------------

    @Test
    void allZeroImageProducesAllNaNRadii() {
        float[][] img = new float[H][W];                                    // all zeros
        RingExtractor ext = new RingExtractor(N_BINS, PIXEL_TO_M, CX, CY);

        double[] radii = ext.extract(img, W, H);

        assertEquals(N_BINS, radii.length);
        for (int i = 0; i < N_BINS; i++) {
            assertTrue(Double.isNaN(radii[i]),
                    "bin " + i + " should be NaN on empty image, got " + radii[i]);
        }
        assertEquals(0, ext.multipeakBinCount(img, W, H, 2.0));
    }

    // ------------------------------------------------------------------
    // Gate 1: synthetic perfect circle
    // ------------------------------------------------------------------

    @Test
    void gate1_perfectCircleHasLowRelativeDispersion() {
        double r0 = 5.5;                                                    // M
        float[][] image = synthesizeGaussianRing(r0, 1.0);
        RingExtractor ext = new RingExtractor(N_BINS, PIXEL_TO_M, CX, CY);

        double[] radii = ext.extract(image, W, H);
        var stats = CircularityMetric.compute(radii);

        // Mean radius must land within one pixel of the true r0
        assertEquals(r0, stats.meanR(), PIXEL_TO_M,
                "mean radius off by more than one pixel: got " + stats.meanR());
        assertEquals(N_BINS, stats.validBins());
        double relDispersion = stats.deltaRrms() / stats.meanR();
        assertTrue(relDispersion < 0.02,
                "expected δ_r / ⟨r⟩ < 0.02 on perfect circle (V1 tolerance); got " + relDispersion);
        assertEquals(0, ext.multipeakBinCount(image, W, H, 2.0),
                "single Gaussian ring should produce zero multi-peak bins");
    }

    // ------------------------------------------------------------------
    // Gate 2: synthetic ellipse, asymmetry vs analytic reference
    // ------------------------------------------------------------------

    @Test
    void gate2_ellipseAsymmetryMatchesAnalyticReference() {
        double aM = 5.5;
        double bM = 5.0;
        float[][] image = synthesizeParametricEllipse(aM, bM, 1.0);
        RingExtractor ext = new RingExtractor(N_BINS, PIXEL_TO_M, CX, CY);

        double[] radii = ext.extract(image, W, H);
        var measured = CircularityMetric.compute(radii);

        double[] refRadii = analyticEllipseRadiiAtBinCenters(aM, bM, N_BINS);
        var reference = CircularityMetric.compute(refRadii);

        double measuredAsym = measured.deltaRrms() / measured.meanR();
        double refAsym = reference.deltaRrms() / reference.meanR();

        // Reference asymmetry for a = 5.5 M, b = 5.0 M, 180 bins is ~3.3 %; pixel quantisation contributes ~0.5 %
        // bin-to-bin noise on top, so we tolerate 10 % relative agreement.
        assertEquals(refAsym, measuredAsym, refAsym * 0.10,
                "ellipse asymmetry: expected ~" + refAsym + " got " + measuredAsym);
    }

    // ------------------------------------------------------------------
    // Gate 6: bin-count convergence on a single frame
    // ------------------------------------------------------------------

    @Test
    void gate6_binCountConvergence() {
        double aM = 5.5;
        double bM = 5.0;
        float[][] image = synthesizeParametricEllipse(aM, bM, 1.0);

        double asym90 = relativeAsymmetry(image, 90);
        double asym180 = relativeAsymmetry(image, 180);
        double asym360 = relativeAsymmetry(image, 360);

        // 5 % relative agreement between 90 / 180 / 360 bins on the same frame
        assertEquals(asym180, asym90, asym180 * 0.05,
                "90-bin vs 180-bin asymmetry diverges: " + asym90 + " vs " + asym180);
        assertEquals(asym180, asym360, asym180 * 0.05,
                "360-bin vs 180-bin asymmetry diverges: " + asym360 + " vs " + asym180);
    }

    // ------------------------------------------------------------------
    // "Outermost peak wins" rule
    // ------------------------------------------------------------------

    @Test
    void extractPicksOutermostLocalMaxEvenWhenInnerPeakIsBrighter() {
        // Inner ring at 4 M, intensity 1.0; outer ring at 7 M, intensity 0.5.
        // Outermost-peak rule: every bin returns ~7 M, not 4 M.
        float[][] image = new float[H][W];
        addGaussianRing(image, 4.0, 1.0, 1.0f);
        addGaussianRing(image, 7.0, 1.0, 0.5f);

        RingExtractor ext = new RingExtractor(N_BINS, PIXEL_TO_M, CX, CY);
        double[] radii = ext.extract(image, W, H);
        var stats = CircularityMetric.compute(radii);

        // Tolerance: pixel quantisation at r = 7 M is ±0.5 px ≈ ±0.03 M
        assertEquals(7.0, stats.meanR(), 0.10,
                "outermost-peak rule should select r ≈ 7 M, got " + stats.meanR());
    }

    @Test
    void extractFallsBackToBrightestPixelWhenNoLocalMax() {
        // Synthetic monotonically-decreasing radial profile (proxy for cusp regime
        // disk-edge bin): brightest at r = 4 M, fades smoothly to zero by r = 8 M.
        // No strict local maximum exists. Fallback should pick the brightest pixel
        // (≈ 4 M).
        float[][] image = new float[H][W];
        for (int y = 0; y < H; y++) {
            for (int x = 0; x < W; x++) {
                double dx = x + 0.5 - CX;
                double dy = y + 0.5 - CY;
                if (dx == 0.0 && dy == 0.0) continue;
                double rPx = Math.hypot(dx, dy);
                double rM = rPx * PIXEL_TO_M;
                if (rM >= 4.0 && rM <= 8.0) {
                    image[y][x] = (float) ((8.0 - rM) / 4.0);                // 1 at r=4, 0 at r=8
                }
            }
        }

        RingExtractor ext = new RingExtractor(N_BINS, PIXEL_TO_M, CX, CY);
        double[] radii = ext.extract(image, W, H);
        var stats = CircularityMetric.compute(radii);

        // Brightest pixel sits at the inner edge r = 4 M ± 1 px
        assertEquals(4.0, stats.meanR(), 0.10,
                "fallback should select brightest pixel near r = 4 M, got " + stats.meanR());
    }

    // ------------------------------------------------------------------
    // Multi-peak diagnostic
    // ------------------------------------------------------------------

    @Test
    void multipeakBinCountFlagsMostBinsForTwoSeparatedRings() {
        // Two concentric rings at r = 4 M and r = 7 M; separation in pixels is
        // 3 M / PIXEL_TO_M ≈ 51 px, far above the 2 px threshold.
        float[][] image = new float[H][W];
        addGaussianRing(image, 4.0, 1.0, 1.0f);
        addGaussianRing(image, 7.0, 1.0, 1.0f);

        RingExtractor ext = new RingExtractor(N_BINS, PIXEL_TO_M, CX, CY);
        int multipeak = ext.multipeakBinCount(image, W, H, 2.0);

        // Two well-separated rings should flag essentially every bin
        assertTrue(multipeak >= N_BINS - 5,
                "expected at least " + (N_BINS - 5) + " multi-peak bins, got " + multipeak);
    }

    @Test
    void multipeakBinCountReturnsZeroForSingleRing() {
        float[][] image = synthesizeGaussianRing(5.5, 1.0);

        RingExtractor ext = new RingExtractor(N_BINS, PIXEL_TO_M, CX, CY);
        int multipeak = ext.multipeakBinCount(image, W, H, 2.0);

        assertEquals(0, multipeak,
                "single ring should produce zero multi-peak bins, got " + multipeak);
    }

    // ------------------------------------------------------------------
    // Helpers
    // ------------------------------------------------------------------

    private double relativeAsymmetry(float[][] image, int numBins) {
        RingExtractor ext = new RingExtractor(numBins, PIXEL_TO_M, CX, CY);
        double[] radii = ext.extract(image, W, H);
        var stats = CircularityMetric.compute(radii);
        return stats.deltaRrms() / stats.meanR();
    }

    /** Single Gaussian ring of unit amplitude, width {@code sigmaPx} px, peak at radius {@code rM} M. */
    private float[][] synthesizeGaussianRing(double rM, double sigmaPx) {
        float[][] img = new float[H][W];
        addGaussianRing(img, rM, sigmaPx, 1.0f);
        return img;
    }

    private void addGaussianRing(float[][] img, double rM, double sigmaPx, float amplitude) {
        double r0Px = rM / PIXEL_TO_M;
        double twoSigmaSq = 2.0 * sigmaPx * sigmaPx;
        for (int y = 0; y < H; y++) {
            for (int x = 0; x < W; x++) {
                double dx = x + 0.5 - CX;
                double dy = y + 0.5 - CY;
                double rPx = Math.hypot(dx, dy);
                double diff = rPx - r0Px;
                img[y][x] += amplitude * (float) Math.exp(-diff * diff / twoSigmaSq);
            }
        }
    }

    /**
     * Synthetic Gaussian ring along the parametric ellipse curve
     * {@code r(θ) = a b / sqrt(b² cos²θ + a² sin²θ)}. For each pixel
     * we compute its radial residual against the ellipse at the same
     * azimuth, then apply a Gaussian profile. Per-bin extraction
     * recovers {@code r(θ_bin)} modulo pixel quantisation.
     */
    private float[][] synthesizeParametricEllipse(double aM, double bM, double sigmaPx) {
        float[][] img = new float[H][W];
        double aPx = aM / PIXEL_TO_M;
        double bPx = bM / PIXEL_TO_M;
        double aPx2 = aPx * aPx;
        double bPx2 = bPx * bPx;
        double twoSigmaSq = 2.0 * sigmaPx * sigmaPx;
        for (int y = 0; y < H; y++) {
            for (int x = 0; x < W; x++) {
                double dx = x + 0.5 - CX;
                double dy = y + 0.5 - CY;
                if (dx == 0.0 && dy == 0.0) continue;
                double rPx = Math.hypot(dx, dy);
                double cos = dx / rPx;
                double sin = dy / rPx;
                double rEllipsePx = aPx * bPx / Math.sqrt(bPx2 * cos * cos + aPx2 * sin * sin);
                double diff = rPx - rEllipsePx;
                img[y][x] = (float) Math.exp(-diff * diff / twoSigmaSq);
            }
        }
        return img;
    }

    /** Sample the parametric ellipse at the centre azimuth of each of {@code numBins} bins. */
    private double[] analyticEllipseRadiiAtBinCenters(double aM, double bM, int numBins) {
        double[] r = new double[numBins];
        for (int i = 0; i < numBins; i++) {
            double theta = 2.0 * Math.PI * (i + 0.5) / numBins;
            double cos = Math.cos(theta);
            double sin = Math.sin(theta);
            r[i] = aM * bM / Math.sqrt(bM * bM * cos * cos + aM * aM * sin * sin);
        }
        return r;
    }
}
