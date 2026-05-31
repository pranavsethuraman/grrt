package com.pranav.grrt.analysis;

import java.io.IOException;
import java.io.UncheckedIOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.Locale;

/**
 * Inverts the Phase 3C asymmetry sweep against the EHT M87* circularity
 * constraint to place a bound on the Johannsen-Psaltis deviation
 * parameter {@code ε₃} at {@code a = 0.9}, observer inclination
 * {@code i = 17°}.
 *
 * <p>The pipeline metric is the dimensionless ring circularity
 * {@code δ_r/⟨r⟩ = delta_r_rms / mean_r} (population SD of the per-bin
 * photon-ring radius divided by the mean radius; see
 * {@link CircularityMetric}). On the strictly-monotone sub-interval
 * {@code ε₃ ∈ [−2.5, −0.2]} this is monotone increasing in {@code ε₃},
 * so a target circularity inverts to a unique {@code ε₃} by
 * piecewise-linear interpolation.
 *
 * <p><b>Specification correction (Phase 3D.1).</b> The Phase 3D plan
 * (§2.3, §10.1) and {@code CLAUDE.md} originally cited <i>EHT Paper VI
 * Table 7</i> as the source of the comparison amplitude. Table 7 of
 * Paper VI (arXiv:1906.11243) is the <i>Summary of θ_g and M
 * Measurements</i> — black-hole angular size and mass — and contains
 * <b>no</b> asymmetry or circularity figure. The correct source is
 * Paper VI <b>Section 7.4 "The Circular Shapes of the M87 Images"</b>
 * and <b>Figure 18</b>, which report the fractional spread of the
 * inferred ring diameter along different orientations, peaking at
 * {@code ∼0.05–0.06} for the reconstructed M87 images. The constants
 * below are taken from there.
 *
 * <p><b>Metric mismatch (documented, intentional).</b> Two
 * approximations are baked into this comparison and are reported as
 * caveats in the manuscript rather than corrected here:
 * <ol>
 *   <li>Our {@code δ_r/⟨r⟩} is an all-azimuthal-mode <i>radius</i> RMS;
 *       the EHT §7.4 metric is a <i>diameter</i> spread, in which a
 *       pure {@code m = 1} radius asymmetry cancels
 *       ({@code d(θ) = r(θ) + r(θ+π)} kills odd modes), leaving an
 *       even-{@code m} (ellipticity) measure. The two therefore probe
 *       different multipoles.</li>
 *   <li>Because {@code δ_r/⟨r⟩ = (1/√2)·√(Σ_all aₘ²) ≥
 *       (1/√2)·√(Σ_even aₘ²)} = the diameter-spread, our metric is an
 *       <i>upper</i> envelope: any {@code ε₃} our metric deems
 *       consistent with the EHT band is also consistent under the
 *       stricter diameter-spread test. The bound is conservative in
 *       the correct direction.</li>
 * </ol>
 *
 * <p><b>Expected result for the committed sweep (honest null).</b> At
 * {@code i = 17°} the extractor noise floor (integer-pixel peak
 * finding, no sub-pixel interpolation — see {@link RingExtractor} and
 * {@code docs/phase-3c-plan.md} §13) puts {@code δ_r/⟨r⟩} in the range
 * {@code [0.12, 0.20]} across the whole sweep, well above the EHT band
 * {@code [0.050, 0.060]}. The EHT band thus lies below the
 * monotone-domain minimum {@code δ_r/⟨r⟩(ε₃ = −2.5) ≈ 0.121}, and
 * {@link #invertFromCsv} returns
 * {@link BoundStatus#LOWER_OPEN_BELOW_PATHOLOGY}: this sweep cannot
 * place a two-sided {@code ε₃} bound, and sub-pixel ring extraction is
 * the identified methodological improvement.
 *
 * <p>Reads no external EHT file: the constraint is hardcoded in
 * {@link #EHT_CIRCULARITY_MID}/{@link #EHT_CIRCULARITY_SIGMA}. Updating
 * the EHT input is a one-line source edit + recompile.
 *
 * <p>Time complexity: {@code O(rows)} to parse plus {@code O(D)} to
 * invert over the {@code D} monotone-domain points. No mutation of
 * inputs; all returns are fresh immutable values.
 */
public final class ConsistencyBound {

    /**
     * EHT 2019 Paper VI M87* fractional ring circularity, mid value.
     *
     * @see "EHT Collaboration 2019, ApJL 875:L6 (arXiv:1906.11243),
     *      §7.4 'The Circular Shapes of the M87 Images' and Figure 18:
     *      the distribution of the fractional spread in inferred
     *      diameters peaks at ∼0.05–0.06 (axial ratio ≈ 4:3). Midpoint
     *      of that range. NOT from Table 7 (mass/size); see class
     *      Javadoc."
     */
    private static final double EHT_CIRCULARITY_MID = 0.055;

    /**
     * Half-width of the EHT circularity band. Paper VI quotes no formal
     * 1σ on the §7.4 fractional spread; we take half the published
     * 0.05–0.06 peak range as a conservative band half-width.
     *
     * @see "EHT Collaboration 2019, ApJL 875:L6, §7.4 / Figure 18."
     */
    private static final double EHT_CIRCULARITY_SIGMA = 0.005;

    /**
     * Lower (more negative) edge of the {@code ε₃} sub-interval on
     * which {@code δ_r/⟨r⟩} is strictly monotone at {@code i = 17°}.
     * Hardcoded from sweep inspection at commit {@code 9def498}.
     */
    private static final double EPS3_DOMAIN_LO = -2.5;

    /** Upper (least negative) edge of the monotone {@code ε₃} domain. */
    private static final double EPS3_DOMAIN_HI = -0.2;

    /** Observer inclination (degrees) whose rows drive the inversion. */
    private static final double TARGET_INCLINATION_DEG = 17.0;

    /** Tolerance for matching the inclination column. */
    private static final double INCLINATION_TOL = 1.0e-6;

    /** Tolerance for the domain-membership comparison. */
    private static final double DOMAIN_TOL = 1.0e-9;

    private ConsistencyBound() {
    }

    /**
     * One {@code i = 17°} sweep row, carrying the fields needed for the
     * circularity inversion. {@code meanR}, {@code deltaRrms} and
     * {@code fourierM1Amp} are in gravitational radii {@code M} (the
     * units of the {@link CircularityMetric} output).
     *
     * @param epsilon3     Johannsen-Psaltis deviation parameter ε₃
     * @param meanR        mean ring radius ⟨r⟩ in M
     * @param deltaRrms    RMS radius dispersion δ_r in M
     * @param fourierM1Amp Fourier {@code m = 1} coefficient magnitude
     *                     {@code |c₁|} in M
     */
    public record SweepPoint(double epsilon3, double meanR,
                             double deltaRrms, double fourierM1Amp) {

        /** Dimensionless circularity {@code δ_r/⟨r⟩ = deltaRrms / meanR}. */
        public double circularity() {
            return deltaRrms / meanR;
        }
    }

    /**
     * Outcome class for {@link #invertFromCsv}.
     *
     * <ul>
     *   <li>{@link #TWO_SIDED} — both EHT band edges cross the monotone
     *       curve; a finite {@code (ε₃_low, ε₃_high)} is returned.</li>
     *   <li>{@link #UPPER_OPEN_PLATEAU} — the lower band edge crosses
     *       the monotone curve but the upper edge would cross in the
     *       small-|ε₃| plateau; the upper bound is unconstrained.</li>
     *   <li>{@link #LOWER_OPEN_BELOW_PATHOLOGY} — the EHT band lies
     *       (wholly or partly) below {@code δ_r/⟨r⟩(ε₃ = −2.5)}; the
     *       lower bound runs off the pathology edge of the sweep.</li>
     *   <li>{@link #PLATEAU_CONSISTENT} — the EHT band lies entirely
     *       above the monotone curve / plateau; the model's asymmetry
     *       never reaches the band and every swept {@code ε₃} is
     *       consistent.</li>
     * </ul>
     */
    public enum BoundStatus {
        TWO_SIDED,
        UPPER_OPEN_PLATEAU,
        LOWER_OPEN_BELOW_PATHOLOGY,
        PLATEAU_CONSISTENT
    }

    /**
     * Inversion result. Open edges are {@link Double#NaN}.
     *
     * @param epsilonLower lower (more negative) ε₃ bound, or NaN if open
     * @param epsilonUpper upper (less negative) ε₃ bound, or NaN if open
     * @param status       which case of the inversion fired
     * @param caveat       human-readable explanation; non-empty for
     *                     every non-{@link BoundStatus#TWO_SIDED} status
     */
    public record Bound(double epsilonLower, double epsilonUpper,
                        BoundStatus status, String caveat) {
    }

    /**
     * Read the {@code i = 17°} rows of {@code csv}, restrict to the
     * monotone domain {@code [−2.5, −0.2]}, and invert the hardcoded
     * EHT circularity band {@code [mid − σ, mid + σ]} to an {@code ε₃}
     * bound.
     *
     * @param csv path to a sweep CSV in the {@link EpsilonSweep#CSV_HEADER}
     *            schema
     * @return the bound; never {@code null}
     * @throws IllegalArgumentException if {@code csv} is null, empty, or
     *         has a header that does not match the sweep schema
     * @throws IllegalStateException if fewer than two {@code i = 17°}
     *         rows fall in the monotone domain
     * @throws UncheckedIOException if {@code csv} cannot be read
     */
    public static Bound invertFromCsv(Path csv) {
        List<SweepPoint> domain = filterToMonotoneDomain(
                readInclination(csv, TARGET_INCLINATION_DEG));
        if (domain.size() < 2) {
            throw new IllegalStateException(
                    "need >= 2 i=" + TARGET_INCLINATION_DEG + "° rows in the monotone domain ["
                            + EPS3_DOMAIN_LO + ", " + EPS3_DOMAIN_HI + "], got " + domain.size());
        }
        return invertCircularity(domain,
                EHT_CIRCULARITY_MID - EHT_CIRCULARITY_SIGMA,
                EHT_CIRCULARITY_MID + EHT_CIRCULARITY_SIGMA);
    }

    /**
     * Convert an RMS-dispersion array to the corresponding Fourier
     * {@code m = 1} amplitude {@code A₁ = √2 · rms}, the exact relation
     * for a ring whose radius is a single cosine in azimuth (higher
     * harmonics negligible; cf. {@code docs/phase-3-plan.md} §7.3).
     * When a measured {@code m = 1} value is supplied for an entry
     * (non-null, non-NaN) it is returned directly in preference to the
     * RMS approximation.
     *
     * <p>Inputs must be expressed in a single consistent unit system
     * (both fractional, e.g. {@code rms = deltaRrms/meanR} and
     * {@code fourierColumn = 2|c₁|/meanR}; or both absolute in M). The
     * {@code √2} relation is unit-agnostic.
     *
     * <p>This is an auxiliary used by the figure script's secondary
     * {@code m = 1} axis; {@link #invertFromCsv} does <b>not</b> use it,
     * because the EHT §7.4 value is itself an RMS-type circularity and
     * is compared directly against {@code δ_r/⟨r⟩} with no {@code √2}
     * single-mode conversion.
     *
     * @param rmsArray     per-row RMS dispersion; must be non-null with
     *                     no NaN entries
     * @param fourierColumn measured {@code m = 1} amplitudes, same
     *                     length as {@code rmsArray}, or {@code null} to
     *                     force the RMS approximation for every entry;
     *                     individual NaN entries also fall back
     * @return a new array of {@code m = 1} amplitudes
     * @throws IllegalArgumentException if {@code rmsArray} is null, has
     *         a NaN entry, or {@code fourierColumn} has a different
     *         length
     */
    public static double[] toFourierM1(double[] rmsArray, double[] fourierColumn) {
        if (rmsArray == null) {
            throw new IllegalArgumentException("rmsArray is null");
        }
        if (fourierColumn != null && fourierColumn.length != rmsArray.length) {
            throw new IllegalArgumentException(
                    "fourierColumn length (" + fourierColumn.length
                            + ") != rmsArray length (" + rmsArray.length + ")");
        }
        final double sqrt2 = Math.sqrt(2.0);
        double[] out = new double[rmsArray.length];
        for (int i = 0; i < rmsArray.length; i++) {
            double rms = rmsArray[i];
            if (Double.isNaN(rms)) {
                throw new IllegalArgumentException("rmsArray[" + i + "] is NaN");
            }
            if (fourierColumn != null && !Double.isNaN(fourierColumn[i])) {
                out[i] = fourierColumn[i];
            } else {
                out[i] = sqrt2 * rms;
            }
        }
        return out;
    }

    // ------------------------------------------------------------------
    // Core inversion (package-private for direct unit testing)
    // ------------------------------------------------------------------

    /**
     * Invert a strictly-monotone-increasing circularity curve against
     * the band {@code [bandLo, bandHi]}.
     *
     * @param domain points sorted ascending in {@code ε₃}, with
     *               strictly-increasing {@code circularity()}
     * @param bandLo lower circularity band edge
     * @param bandHi upper circularity band edge ({@code ≥ bandLo})
     * @return the bound
     * @throws IllegalArgumentException on degenerate arguments
     * @throws IllegalStateException if the curve is not strictly
     *         increasing (the domain is assumed monotone)
     */
    static Bound invertCircularity(List<SweepPoint> domain, double bandLo, double bandHi) {
        if (domain == null) {
            throw new IllegalArgumentException("domain is null");
        }
        if (domain.size() < 2) {
            throw new IllegalArgumentException("need >= 2 points, got " + domain.size());
        }
        if (!(bandLo <= bandHi)) {
            throw new IllegalArgumentException("bandLo (" + bandLo + ") must be <= bandHi (" + bandHi + ")");
        }
        final int n = domain.size();
        double[] eps = new double[n];
        double[] c = new double[n];
        for (int i = 0; i < n; i++) {
            eps[i] = domain.get(i).epsilon3();
            c[i] = domain.get(i).circularity();
            if (i > 0 && !(eps[i] > eps[i - 1])) {
                throw new IllegalStateException(
                        "domain not strictly ascending in ε₃ at index " + i);
            }
            if (i > 0 && !(c[i] > c[i - 1])) {
                throw new IllegalStateException(
                        "circularity not strictly increasing at index " + i
                                + " (" + c[i - 1] + " -> " + c[i] + ")");
            }
        }
        final double cMin = c[0];
        final double cMax = c[n - 1];
        Side lo = sideOf(bandLo, cMin, cMax);
        Side hi = sideOf(bandHi, cMin, cMax);
        return combine(eps, c, bandLo, bandHi, lo, hi, cMin, cMax);
    }

    /** Where a target circularity sits relative to the monotone curve. */
    private enum Side { BELOW, INSIDE, ABOVE }

    private static Side sideOf(double target, double cMin, double cMax) {
        if (target < cMin) {
            return Side.BELOW;
        }
        if (target > cMax) {
            return Side.ABOVE;
        }
        return Side.INSIDE;
    }

    private static Bound combine(double[] eps, double[] c,
                                 double bandLo, double bandHi,
                                 Side lo, Side hi, double cMin, double cMax) {
        // bandLo <= bandHi and c is increasing, so the lower band edge
        // can never sit "more above" the curve than the upper edge.
        return switch (lo) {
            case BELOW -> switch (hi) {
                case BELOW -> lowerOpen(Double.NaN, cMin, cMax, bandLo, bandHi,
                        "wholly");
                case INSIDE -> lowerOpen(invert(eps, c, bandHi), cMin, cMax, bandLo, bandHi,
                        "at its lower edge");
                case ABOVE -> bracketsWholeCurve(eps, cMin, cMax, bandLo, bandHi);
            };
            case INSIDE -> switch (hi) {
                case INSIDE -> new Bound(invert(eps, c, bandLo), invert(eps, c, bandHi),
                        BoundStatus.TWO_SIDED, "");
                case ABOVE -> new Bound(invert(eps, c, bandLo), Double.NaN,
                        BoundStatus.UPPER_OPEN_PLATEAU, upperOpenCaveat(cMax, bandHi));
                case BELOW -> throw unreachable(lo, hi);
            };
            case ABOVE -> switch (hi) {
                case ABOVE -> new Bound(Double.NaN, Double.NaN,
                        BoundStatus.PLATEAU_CONSISTENT, plateauCaveat(cMax, bandLo));
                case INSIDE, BELOW -> throw unreachable(lo, hi);
            };
        };
    }

    private static Bound lowerOpen(double epsilonUpper, double cMin, double cMax,
                                   double bandLo, double bandHi, String how) {
        return new Bound(Double.NaN, epsilonUpper,
                BoundStatus.LOWER_OPEN_BELOW_PATHOLOGY,
                String.format(Locale.ROOT,
                        "EHT §7.4 circularity band [%.4f, %.4f] lies %s below the i=%.0f° "
                                + "monotone-domain range δ_r/⟨r⟩ ∈ [%.4f, %.4f] "
                                + "(min at ε₃=%.1f); the extractor noise floor exceeds the EHT "
                                + "circularity, so this sweep places no two-sided ε₃ bound. "
                                + "Sub-pixel ring extraction (phase-3c-plan §13) is the required "
                                + "methodological improvement.",
                        bandLo, bandHi, how, TARGET_INCLINATION_DEG, cMin, cMax, EPS3_DOMAIN_LO));
    }

    private static Bound bracketsWholeCurve(double[] eps, double cMin, double cMax,
                                            double bandLo, double bandHi) {
        return new Bound(eps[0], eps[eps.length - 1], BoundStatus.TWO_SIDED,
                String.format(Locale.ROOT,
                        "EHT band [%.4f, %.4f] brackets the entire swept circularity range "
                                + "[%.4f, %.4f]; every domain ε₃ is consistent, so the bound is "
                                + "the swept-domain edges, not data-driven crossings.",
                        bandLo, bandHi, cMin, cMax));
    }

    private static String upperOpenCaveat(double cMax, double bandHi) {
        return String.format(Locale.ROOT,
                "Upper EHT band edge %.4f exceeds the monotone-domain maximum δ_r/⟨r⟩=%.4f "
                        + "(at ε₃=%.1f); the upper crossing falls in the small-|ε₃| plateau and is "
                        + "unconstrained by this sweep.",
                bandHi, cMax, EPS3_DOMAIN_HI);
    }

    private static String plateauCaveat(double cMax, double bandLo) {
        return String.format(Locale.ROOT,
                "EHT band lower edge %.4f exceeds the monotone-domain maximum δ_r/⟨r⟩=%.4f; "
                        + "the model's asymmetry never reaches the EHT band and every swept ε₃ is "
                        + "consistent (upper bound unconstrained).",
                bandLo, cMax);
    }

    private static IllegalStateException unreachable(Side lo, Side hi) {
        return new IllegalStateException(
                "unreachable band ordering: lo=" + lo + " hi=" + hi
                        + " (bandLo <= bandHi on an increasing curve)");
    }

    /**
     * Piecewise-linear inverse of a strictly-increasing curve:
     * the {@code ε₃} at which {@code circularity = target}, for
     * {@code cMin ≤ target ≤ cMax}.
     */
    private static double invert(double[] eps, double[] c, double target) {
        for (int i = 0; i < c.length - 1; i++) {
            if (target >= c[i] && target <= c[i + 1]) {
                double frac = (target - c[i]) / (c[i + 1] - c[i]);
                return eps[i] + frac * (eps[i + 1] - eps[i]);
            }
        }
        // target == cMax to within rounding: pin to the last node.
        return eps[eps.length - 1];
    }

    // ------------------------------------------------------------------
    // CSV parsing
    // ------------------------------------------------------------------

    /**
     * Parse the sweep CSV and return the rows whose inclination column
     * equals {@code inclDeg} (within {@link #INCLINATION_TOL}).
     *
     * <p>Validates that the header equals {@link EpsilonSweep#CSV_HEADER}
     * so a schema drift fails loudly rather than mis-indexing columns.
     * Malformed (non-numeric) data rows are skipped.
     *
     * @param csv     path to the sweep CSV
     * @param inclDeg inclination to select, in degrees
     * @return the matching rows in file order; never {@code null}
     * @throws IllegalArgumentException if {@code csv} is null, empty, or
     *         header-mismatched
     * @throws UncheckedIOException if the file cannot be read
     */
    static List<SweepPoint> readInclination(Path csv, double inclDeg) {
        if (csv == null) {
            throw new IllegalArgumentException("csv is null");
        }
        final List<String> lines;
        try {
            lines = Files.readAllLines(csv, StandardCharsets.UTF_8);
        } catch (IOException e) {
            throw new UncheckedIOException("failed to read sweep CSV: " + csv, e);
        }
        if (lines.isEmpty()) {
            throw new IllegalArgumentException("empty CSV: " + csv);
        }
        if (!lines.get(0).equals(EpsilonSweep.CSV_HEADER)) {
            throw new IllegalArgumentException(
                    "CSV header mismatch in " + csv + ": expected the EpsilonSweep schema, got '"
                            + lines.get(0) + "'");
        }
        // Frozen-schema column indices (see EpsilonSweep.CSV_HEADER).
        final int colEps3 = 0;
        final int colIncl = 1;
        final int colMeanR = 3;
        final int colDeltaRrms = 4;
        final int colFourierM1 = 6;
        List<SweepPoint> out = new ArrayList<>();
        for (int i = 1; i < lines.size(); i++) {
            String line = lines.get(i).trim();
            if (line.isEmpty()) {
                continue;
            }
            String[] p = line.split(",", -1);
            if (p.length <= colFourierM1) {
                continue;
            }
            try {
                double incl = Double.parseDouble(p[colIncl]);
                if (Math.abs(incl - inclDeg) > INCLINATION_TOL) {
                    continue;
                }
                out.add(new SweepPoint(
                        Double.parseDouble(p[colEps3]),
                        Double.parseDouble(p[colMeanR]),
                        Double.parseDouble(p[colDeltaRrms]),
                        Double.parseDouble(p[colFourierM1])));
            } catch (NumberFormatException ignore) {
                // Malformed row — skip, consistent with EpsilonSweep's reader.
            }
        }
        return out;
    }

    /**
     * Restrict to the monotone {@code ε₃} domain {@code [−2.5, −0.2]}
     * and sort ascending in {@code ε₃}.
     */
    private static List<SweepPoint> filterToMonotoneDomain(List<SweepPoint> points) {
        return points.stream()
                .filter(p -> p.epsilon3() >= EPS3_DOMAIN_LO - DOMAIN_TOL
                        && p.epsilon3() <= EPS3_DOMAIN_HI + DOMAIN_TOL)
                .sorted(Comparator.comparingDouble(SweepPoint::epsilon3))
                .toList();
    }
}
