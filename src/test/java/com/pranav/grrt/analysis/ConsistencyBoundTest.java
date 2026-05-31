package com.pranav.grrt.analysis;

import com.pranav.grrt.analysis.ConsistencyBound.Bound;
import com.pranav.grrt.analysis.ConsistencyBound.BoundStatus;
import com.pranav.grrt.analysis.ConsistencyBound.SweepPoint;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

import java.io.IOException;
import java.net.URISyntaxException;
import java.net.URL;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Phase 3D.1 unit tests for {@link ConsistencyBound}, implementing
 * gates 1–3 of {@code docs/phase-3d-plan.md} §5 plus the honest-null
 * regression for the committed sweep.
 *
 * <ul>
 *   <li>Gate 1 — {@code invertFromCsv} on a synthetic linear-monotone
 *       fixture returns the analytic two-sided bound (tol 1e-10).</li>
 *   <li>Gate 2 — {@code toFourierM1} on a pure-cosine ring recovers
 *       {@code A₁ = √2·(δ_r/⟨r⟩)_RMS} (tol 1e-10).</li>
 *   <li>Gate 3 — the {@link BoundStatus} flags fire for the
 *       below-curve, above-curve, and upper-open cases, each with a
 *       non-empty caveat.</li>
 *   <li>Honest null — the real {@code i = 17°} circularity values
 *       (commit {@code 9def498}) invert to
 *       {@link BoundStatus#LOWER_OPEN_BELOW_PATHOLOGY} against the
 *       EHT §7.4 band.</li>
 * </ul>
 *
 * Tests do not depend on {@code output/sweep.csv} (gitignored); the
 * fixture lives at {@code src/test/resources/analysis/synthetic_sweep.csv}.
 */
class ConsistencyBoundTest {

    private static final double TOL = 1.0e-10;
    private static final double SQRT2 = Math.sqrt(2.0);

    // ------------------------------------------------------------------
    // Gate 1: end-to-end CSV inversion -> analytic two-sided bound
    // ------------------------------------------------------------------

    /**
     * The fixture's i=17° circularity is δ_r/⟨r⟩ = 0.072 + 0.01·ε₃ on
     * the domain points {−2.5,−1.5,−1.0,−0.5,−0.2}. The hardcoded EHT
     * band [0.050, 0.060] therefore crosses at ε₃ = −2.2 and −1.2; the
     * exact crossings also pin the band to mid=0.055, σ=0.005.
     */
    @Test
    void invertFromCsvSyntheticLinearReturnsTwoSidedBound() {
        Bound b = ConsistencyBound.invertFromCsv(fixture("synthetic_sweep.csv"));
        assertEquals(BoundStatus.TWO_SIDED, b.status());
        assertEquals(-2.2, b.epsilonLower(), TOL);
        assertEquals(-1.2, b.epsilonUpper(), TOL);
    }

    // ------------------------------------------------------------------
    // Gate 2: toFourierM1 RMS -> m=1 conversion
    // ------------------------------------------------------------------

    @Test
    void toFourierM1PureCosineRingRecoversSqrt2TimesRms() {
        final double r0 = 5.5;
        final double amp = 0.1;
        double[] ring = pureCosineRing(180, r0, amp, 0.0);
        CircularityMetric.Result res = CircularityMetric.compute(ring);

        double fracRms = res.deltaRrms() / res.meanR();          // = A / √2
        double[] a1 = ConsistencyBound.toFourierM1(new double[] {fracRms}, null);
        assertEquals(amp, a1[0], TOL);                            // √2·(A/√2) = A

        // The measured fractional m=1 (2|c₁|/⟨r⟩) agrees with √2·RMS for a pure mode.
        double fracMeasured = 2.0 * res.fourierM1Amplitude() / res.meanR();
        assertEquals(amp, fracMeasured, TOL);
        assertEquals(SQRT2 * fracRms, fracMeasured, TOL);
    }

    @Test
    void toFourierM1PrefersMeasuredColumnAndFallsBackOnNaN() {
        double[] rms = {0.07, 0.08};
        double[] measured = {0.123, Double.NaN};
        double[] out = ConsistencyBound.toFourierM1(rms, measured);
        assertEquals(0.123, out[0], TOL);             // measured value used directly
        assertEquals(SQRT2 * 0.08, out[1], TOL);      // NaN entry falls back to √2·RMS
    }

    @Test
    void toFourierM1NullRmsThrowsIAE() {
        assertThrows(IllegalArgumentException.class,
                () -> ConsistencyBound.toFourierM1(null, null));
    }

    @Test
    void toFourierM1LengthMismatchThrowsIAE() {
        assertThrows(IllegalArgumentException.class,
                () -> ConsistencyBound.toFourierM1(new double[] {0.1}, new double[] {0.1, 0.2}));
    }

    @Test
    void toFourierM1NaNRmsThrowsIAE() {
        assertThrows(IllegalArgumentException.class,
                () -> ConsistencyBound.toFourierM1(new double[] {Double.NaN}, null));
    }

    // ------------------------------------------------------------------
    // Gate 3: status flags for the open / degenerate cases
    // ------------------------------------------------------------------

    @Test
    void invertEhtBandAboveCurveReturnsPlateauConsistent() {
        List<SweepPoint> domain = List.of(
                new SweepPoint(-2.5, 10.0, 0.10, 0.0),   // circularity 0.010
                new SweepPoint(-1.0, 10.0, 0.20, 0.0),   // 0.020
                new SweepPoint(-0.2, 10.0, 0.30, 0.0));  // 0.030
        Bound b = ConsistencyBound.invertCircularity(domain, 0.050, 0.060);
        assertEquals(BoundStatus.PLATEAU_CONSISTENT, b.status());
        assertTrue(Double.isNaN(b.epsilonLower()));
        assertTrue(Double.isNaN(b.epsilonUpper()));
        assertFalse(b.caveat().isEmpty());
    }

    @Test
    void invertUpperEdgeAbovePlateauReturnsUpperOpenPlateau() {
        // cMin=0.045 (ε=−2.5), cMax=0.055 (ε=−0.2); band [0.050,0.060]:
        // lower edge inside (crosses at ε=−1.35), upper edge above.
        List<SweepPoint> domain = List.of(
                new SweepPoint(-2.5, 10.0, 0.45, 0.0),
                new SweepPoint(-0.2, 10.0, 0.55, 0.0));
        Bound b = ConsistencyBound.invertCircularity(domain, 0.050, 0.060);
        assertEquals(BoundStatus.UPPER_OPEN_PLATEAU, b.status());
        assertEquals(-1.35, b.epsilonLower(), TOL);
        assertTrue(Double.isNaN(b.epsilonUpper()));
        assertFalse(b.caveat().isEmpty());
    }

    @Test
    void invertBandBracketsWholeCurveReturnsTwoSidedDomainEdges() {
        // cMin=0.052, cMax=0.058; band [0.050,0.060] brackets the curve.
        List<SweepPoint> domain = List.of(
                new SweepPoint(-2.5, 10.0, 0.52, 0.0),
                new SweepPoint(-0.2, 10.0, 0.58, 0.0));
        Bound b = ConsistencyBound.invertCircularity(domain, 0.050, 0.060);
        assertEquals(BoundStatus.TWO_SIDED, b.status());
        assertEquals(-2.5, b.epsilonLower(), TOL);
        assertEquals(-0.2, b.epsilonUpper(), TOL);
        assertFalse(b.caveat().isEmpty());
    }

    @Test
    void invertInsideInsideReturnsTwoSidedFiniteBound() {
        // circularity 0.072 + 0.01·ε on the real domain points.
        List<SweepPoint> domain = List.of(
                new SweepPoint(-2.5, 10.0, 0.47, 0.0),
                new SweepPoint(-1.5, 10.0, 0.57, 0.0),
                new SweepPoint(-1.0, 10.0, 0.62, 0.0),
                new SweepPoint(-0.5, 10.0, 0.67, 0.0),
                new SweepPoint(-0.2, 10.0, 0.70, 0.0));
        Bound b = ConsistencyBound.invertCircularity(domain, 0.050, 0.060);
        assertEquals(BoundStatus.TWO_SIDED, b.status());
        assertEquals(-2.2, b.epsilonLower(), TOL);
        assertEquals(-1.2, b.epsilonUpper(), TOL);
    }

    // ------------------------------------------------------------------
    // Honest null: real i=17° sweep values invert to LOWER_OPEN
    // ------------------------------------------------------------------

    /**
     * Real {@code i = 17°} circularity at commit {@code 9def498} ranges
     * over {@code [0.121, 0.198]} — wholly above the EHT §7.4 band
     * {@code [0.050, 0.060]} — so the inversion runs off the pathology
     * edge. This is the published "honest null" result of Phase 3D.
     */
    @Test
    void invertRealSweepValuesReturnsLowerOpenBelowPathology() {
        List<SweepPoint> domain = List.of(
                new SweepPoint(-2.5, 6.811198, 0.8248947, 0.2407991),  // δ_r/⟨r⟩ ≈ 0.1211
                new SweepPoint(-1.5, 6.248047, 0.9053509, 0.2395076),  // ≈ 0.1449
                new SweepPoint(-1.0, 5.918945, 0.9528663, 0.2429796),  // ≈ 0.1610
                new SweepPoint(-0.5, 5.542318, 1.0025790, 0.2266303),  // ≈ 0.1809
                new SweepPoint(-0.2, 5.286133, 1.0468330, 0.2260897)); // ≈ 0.1980
        Bound b = ConsistencyBound.invertCircularity(domain, 0.050, 0.060);
        assertEquals(BoundStatus.LOWER_OPEN_BELOW_PATHOLOGY, b.status());
        assertTrue(Double.isNaN(b.epsilonLower()));
        assertTrue(Double.isNaN(b.epsilonUpper()));
        assertFalse(b.caveat().isEmpty());
    }

    // ------------------------------------------------------------------
    // Argument / monotonicity validation
    // ------------------------------------------------------------------

    @Test
    void invertNonMonotoneCircularityThrowsISE() {
        List<SweepPoint> domain = List.of(
                new SweepPoint(-2.5, 10.0, 0.55, 0.0),   // 0.055
                new SweepPoint(-0.2, 10.0, 0.45, 0.0));  // 0.045 (decreasing)
        assertThrows(IllegalStateException.class,
                () -> ConsistencyBound.invertCircularity(domain, 0.050, 0.060));
    }

    @Test
    void invertNonAscendingEpsilonThrowsISE() {
        List<SweepPoint> domain = List.of(
                new SweepPoint(-0.2, 10.0, 0.45, 0.0),
                new SweepPoint(-2.5, 10.0, 0.55, 0.0));  // ε not ascending
        assertThrows(IllegalStateException.class,
                () -> ConsistencyBound.invertCircularity(domain, 0.050, 0.060));
    }

    @Test
    void invertFewerThanTwoPointsThrowsIAE() {
        assertThrows(IllegalArgumentException.class,
                () -> ConsistencyBound.invertCircularity(
                        List.of(new SweepPoint(-1.0, 10.0, 0.55, 0.0)), 0.050, 0.060));
    }

    @Test
    void invertBandLoGreaterThanBandHiThrowsIAE() {
        List<SweepPoint> domain = List.of(
                new SweepPoint(-2.5, 10.0, 0.45, 0.0),
                new SweepPoint(-0.2, 10.0, 0.55, 0.0));
        assertThrows(IllegalArgumentException.class,
                () -> ConsistencyBound.invertCircularity(domain, 0.060, 0.050));
    }

    // ------------------------------------------------------------------
    // CSV parsing
    // ------------------------------------------------------------------

    @Test
    void readInclinationFiltersByInclination() {
        Path csv = fixture("synthetic_sweep.csv");
        List<SweepPoint> at17 = ConsistencyBound.readInclination(csv, 17.0);
        assertEquals(6, at17.size());   // 5 domain rows + 1 out-of-domain (+0.10)
        List<SweepPoint> at60 = ConsistencyBound.readInclination(csv, 60.0);
        assertEquals(1, at60.size());
    }

    @Test
    void readInclinationBadHeaderThrowsIAE(@TempDir Path dir) throws IOException {
        Path bad = dir.resolve("bad.csv");
        Files.writeString(bad, "wrong,header\n1,2\n", StandardCharsets.UTF_8);
        assertThrows(IllegalArgumentException.class,
                () -> ConsistencyBound.readInclination(bad, 17.0));
    }

    @Test
    void readInclinationEmptyFileThrowsIAE(@TempDir Path dir) throws IOException {
        Path empty = dir.resolve("empty.csv");
        Files.writeString(empty, "", StandardCharsets.UTF_8);
        assertThrows(IllegalArgumentException.class,
                () -> ConsistencyBound.readInclination(empty, 17.0));
    }

    @Test
    void invertFromCsvNullPathThrowsIAE() {
        assertThrows(IllegalArgumentException.class,
                () -> ConsistencyBound.invertFromCsv(null));
    }

    @Test
    void invertFromCsvTooFewDomainRowsThrowsISE(@TempDir Path dir) throws IOException {
        Path csv = dir.resolve("one.csv");
        Files.writeString(csv, EpsilonSweep.CSV_HEADER + "\n"
                + "-1.0000,17.0,512,10.0,0.62,3.2,0.47,-2.0,0,1.0,x,2026-01-01T00:00:00Z\n",
                StandardCharsets.UTF_8);
        assertThrows(IllegalStateException.class,
                () -> ConsistencyBound.invertFromCsv(csv));
    }

    // ------------------------------------------------------------------
    // Helpers
    // ------------------------------------------------------------------

    private static double[] pureCosineRing(int n, double r0, double amp, double phase) {
        double[] r = new double[n];
        for (int i = 0; i < n; i++) {
            r[i] = r0 * (1.0 + amp * Math.cos(2.0 * Math.PI * i / n + phase));
        }
        return r;
    }

    private static Path fixture(String name) {
        URL url = ConsistencyBoundTest.class.getResource("/analysis/" + name);
        assertNotNull(url, "fixture not on test classpath: /analysis/" + name);
        try {
            return Path.of(url.toURI());
        } catch (URISyntaxException e) {
            throw new IllegalStateException("bad fixture URI: " + url, e);
        }
    }
}
