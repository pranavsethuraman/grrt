package com.pranav.grrt.analysis;

import java.util.ArrayList;
import java.util.List;

/**
 * Extracts a per-azimuthal-bin ring-radius array from a 2-D image
 * produced by the renderer. The returned {@code double[N]} is the
 * input to {@link CircularityMetric#compute(double[])}, completing
 * the image → radii → asymmetry pipeline that {@code EpsilonSweep}
 * (Phase 3C, not yet committed) drives once per rendered frame.
 *
 * <p>See {@code docs/phase-3c-plan.md} §2.1 for the contract,
 * {@code docs/phase-3-plan.md} §6.4 for the "dominant bright ring
 * (n = 0 + n = 1 unresolved), 180 bins, per-bin peak" definition
 * (locked decision D-C3), and {@code docs/phase-3c-plan.md} §9.2
 * for the cusp-regime fallback semantics.
 *
 * <p><b>Algorithm.</b> A single image scan builds a per-bin 1-D
 * radial profile {@code radialProfile[bin][rPx] =} maximum image
 * intensity over all pixels in the {@code (bin, rPx)} cell.
 * {@link #extract} then walks each bin's 1-D profile, finds all
 * strict local maxima (cells where the intensity exceeds both
 * radial neighbours), and returns the radius of the
 * radially-outermost local maximum — the "outermost peak wins"
 * rule ({@code docs/phase-3b-status.md} §6.3 decision 3,
 * {@code docs/phase-3c-plan.md} §2.4). When a bin has no local
 * maximum (typical of the cusp regime past
 * {@code ε₃ = ε₃_crit ≈ 0.1212}: the prograde bins show a monotone
 * radial profile bounded by the disk inner edge instead of a
 * photon-ring peak), the brightest pixel in the bin is used as a
 * fallback (recovery path 3 in {@code docs/phase-3c-plan.md} §9.2).
 *
 * <p><b>Sub-pixel precision.</b> The current implementation snaps
 * each per-bin peak to its integer-pixel radius (no parabolic
 * interpolation). For a 512² raster with field-of-view ±15 M
 * (pixel scale ≈ 0.0586 M / px), this caps the per-bin radius
 * resolution at about 0.03 M and the synthetic perfect-circle
 * dispersion at about 0.3 % relative — see Gate 1 in
 * {@code RingExtractorTest}. Sub-pixel parabolic interpolation is
 * a deferred enhancement (out of scope for 3C.1; see
 * {@code docs/phase-3c-plan.md} §13).
 *
 * <p><b>Multi-peak diagnostic.</b> {@link #multipeakBinCount} counts
 * bins whose radial profile contains two or more local maxima
 * separated by more than {@code pxThreshold} radial pixels. This is
 * a diagnostic for cusp-regime structure: in the smooth-deformation
 * regime we expect zero such bins, while past the cusp the prograde
 * bins may show coexisting disk-edge and faint-ring peaks. The
 * first flagged bin in a frame triggers a single
 * {@code [ring-warn]}-prefixed {@link System#out} message.
 *
 * <p><b>Image convention.</b> Row-major: {@code image[y][x]} is the
 * intensity at pixel column {@code x}, row {@code y}, with
 * {@code 0 ≤ x &lt; width} and {@code 0 ≤ y &lt; height}. Pixel
 * centres are at {@code (x + 0.5, y + 0.5)}; centroid coordinates
 * {@code (cxPx, cyPx)} use the same half-pixel convention so a
 * 512² image has {@code cxPx = cyPx = 255.5} for an origin-centred
 * camera. Time complexity: {@code O(W · H + N · maxRadius)}.
 * Allocation per call: one {@code float[N][maxRadius]} radial
 * profile and one {@code double[N]} result; no instance-level
 * scratch.
 */
public final class RingExtractor {

    private final int numBins;
    private final double pixelToM;
    private final double cxPx;
    private final double cyPx;

    /**
     * @param numBins  number of azimuthal bins; bin {@code i} covers
     *                 {@code θ ∈ [2π i / N, 2π (i + 1) / N)} where
     *                 {@code θ = atan2(y − cyPx, x − cxPx)} mapped
     *                 to {@code [0, 2π)}; must be {@code ≥ 1}
     * @param pixelToM image-plane resolution in {@code M / px}; must
     *                 be strictly positive and finite
     * @param cxPx     image-plane centroid x in pixels (centre is at
     *                 {@code (cxPx, cyPx)} on the half-pixel grid);
     *                 must be finite
     * @param cyPx     image-plane centroid y in pixels; must be
     *                 finite
     * @throws IllegalArgumentException on any of the above
     */
    public RingExtractor(int numBins, double pixelToM, double cxPx, double cyPx) {
        if (numBins < 1) {
            throw new IllegalArgumentException("numBins must be >= 1, got " + numBins);
        }
        if (!(pixelToM > 0.0) || Double.isInfinite(pixelToM)) {
            throw new IllegalArgumentException(
                    "pixelToM must be a positive finite double, got " + pixelToM);
        }
        if (!Double.isFinite(cxPx) || !Double.isFinite(cyPx)) {
            throw new IllegalArgumentException(
                    "centroid (cxPx, cyPx) must be finite, got (" + cxPx + ", " + cyPx + ")");
        }
        this.numBins = numBins;
        this.pixelToM = pixelToM;
        this.cxPx = cxPx;
        this.cyPx = cyPx;
    }

    /** Number of azimuthal bins the extractor was constructed with. */
    public int numBins() {
        return numBins;
    }

    /** Image-plane resolution in {@code M / px}. */
    public double pixelToM() {
        return pixelToM;
    }

    /**
     * Extract per-bin peak-intensity radius. See class Javadoc for
     * algorithm and fallback semantics.
     *
     * @param image  2-D row-major intensity buffer; {@code image[y][x]}
     * @param width  image width in pixels; must equal each
     *               {@code image[y].length}
     * @param height image height in pixels; must equal
     *               {@code image.length}
     * @return new {@code double[numBins]} of per-bin peak radius in
     *         M; {@code Double.NaN} for bins containing no positive
     *         intensity
     * @throws IllegalArgumentException if the image is null, has
     *         non-positive dimensions, or has rows of inconsistent
     *         length
     */
    public double[] extract(float[][] image, int width, int height) {
        validateImage(image, width, height);
        int maxRPx = computeMaxRPx(width, height);
        float[][] profile = buildRadialProfile(image, width, height, maxRPx);
        double[] result = new double[numBins];
        for (int b = 0; b < numBins; b++) {
            int rPx = peakRadiusIndexForBin(profile[b], maxRPx);
            result[b] = (rPx < 0) ? Double.NaN : rPx * pixelToM;
        }
        return result;
    }

    /**
     * Count bins whose radial profile shows two or more local maxima
     * separated by more than {@code pxThreshold} radial pixels. See
     * class Javadoc for diagnostic semantics. The first flagged bin
     * in a call emits one {@code [ring-warn]}-prefixed
     * {@link System#out} line; subsequent flagged bins in the same
     * call are silent (they are still counted).
     *
     * @param image       2-D intensity buffer (same convention as
     *                    {@link #extract})
     * @param width       image width
     * @param height      image height
     * @param pxThreshold radial-pixel threshold above which two
     *                    peaks are treated as distinct; must be
     *                    strictly positive and finite
     * @return count of multi-peak bins, in {@code [0, numBins]}
     * @throws IllegalArgumentException on bad image arguments or
     *         non-positive {@code pxThreshold}
     */
    public int multipeakBinCount(float[][] image, int width, int height, double pxThreshold) {
        validateImage(image, width, height);
        if (!(pxThreshold > 0.0) || Double.isInfinite(pxThreshold)) {
            throw new IllegalArgumentException(
                    "pxThreshold must be a positive finite double, got " + pxThreshold);
        }
        int maxRPx = computeMaxRPx(width, height);
        float[][] profile = buildRadialProfile(image, width, height, maxRPx);
        int count = 0;
        boolean warned = false;
        for (int b = 0; b < numBins; b++) {
            List<Integer> peakRPx = localMaximaInProfile(profile[b], maxRPx);
            if (peakRPx.size() < 2) {
                continue;
            }
            int clusters = countClustersAboveThreshold(peakRPx, pxThreshold);
            if (clusters >= 2) {
                count++;
                if (!warned) {
                    System.out.println("[ring-warn] bin " + b + " of " + numBins
                            + " has " + clusters + " distinct peaks ("
                            + peakRPx.size() + " local maxima); outermost will be used");
                    warned = true;
                }
            }
        }
        return count;
    }

    // ------------------------------------------------------------------
    // Internals
    // ------------------------------------------------------------------

    /** Upper bound on the radial-pixel index any pixel can take. */
    private int computeMaxRPx(int width, int height) {
        return (int) Math.ceil(Math.hypot(width, height)) + 1;
    }

    /**
     * Build {@code profile[bin][rIdx]} = max image intensity over
     * all pixels in the {@code (bin, integer-radial-pixel)} cell.
     * NaN and non-positive intensities are skipped.
     */
    private float[][] buildRadialProfile(float[][] image, int width, int height, int maxRPx) {
        float[][] profile = new float[numBins][maxRPx];
        final double twoPi = 2.0 * Math.PI;
        final double binsOverTwoPi = numBins / twoPi;
        for (int y = 0; y < height; y++) {
            float[] row = image[y];
            double dy = y + 0.5 - cyPx;
            for (int x = 0; x < width; x++) {
                double dx = x + 0.5 - cxPx;
                if (dx == 0.0 && dy == 0.0) {
                    continue;
                }
                float v = row[x];
                if (Float.isNaN(v) || v <= 0.0f) {
                    continue;
                }
                double rPx = Math.hypot(dx, dy);
                int rIdx = (int) Math.floor(rPx);
                if (rIdx < 0 || rIdx >= maxRPx) {
                    continue;
                }
                double theta = Math.atan2(dy, dx);
                if (theta < 0.0) {
                    theta += twoPi;
                }
                int bin = (int) Math.floor(theta * binsOverTwoPi);
                if (bin >= numBins) {
                    bin = numBins - 1;
                }
                if (v > profile[bin][rIdx]) {
                    profile[bin][rIdx] = v;
                }
            }
        }
        return profile;
    }

    /**
     * Return the radial-pixel index for a bin per the
     * "outermost peak wins" rule with brightest-pixel fallback.
     *
     * @return the chosen rPx index, or {@code -1} if the bin has no
     *         positive intensity
     */
    private int peakRadiusIndexForBin(float[] binProfile, int maxRPx) {
        int outermostPeakRPx = -1;
        float profileMaxIntensity = 0.0f;
        int profileMaxRPx = -1;
        for (int r = 1; r < maxRPx - 1; r++) {
            float v = binProfile[r];
            if (v <= 0.0f) {
                continue;
            }
            if (v > binProfile[r - 1] && v > binProfile[r + 1]) {
                outermostPeakRPx = r;
            }
            if (v > profileMaxIntensity) {
                profileMaxIntensity = v;
                profileMaxRPx = r;
            }
        }
        if (outermostPeakRPx >= 0) {
            return outermostPeakRPx;
        }
        return profileMaxRPx;
    }

    /** Strict local maxima ({@code v[r-1] < v[r] > v[r+1]}, {@code v[r] > 0}). */
    private List<Integer> localMaximaInProfile(float[] binProfile, int maxRPx) {
        List<Integer> peaks = new ArrayList<>();
        for (int r = 1; r < maxRPx - 1; r++) {
            float v = binProfile[r];
            if (v <= 0.0f) {
                continue;
            }
            if (v > binProfile[r - 1] && v > binProfile[r + 1]) {
                peaks.add(r);
            }
        }
        return peaks;
    }

    /**
     * Greedy left-to-right cluster count: peaks within
     * {@code pxThreshold} of the previous cluster's last member
     * are absorbed into that cluster; otherwise a new cluster
     * starts. The input {@code peakRPx} is assumed sorted ascending
     * (the caller iterates radii in order).
     */
    private int countClustersAboveThreshold(List<Integer> peakRPx, double pxThreshold) {
        int clusters = 1;
        int prev = peakRPx.get(0);
        for (int j = 1; j < peakRPx.size(); j++) {
            int cur = peakRPx.get(j);
            if (cur - prev > pxThreshold) {
                clusters++;
            }
            prev = cur;
        }
        return clusters;
    }

    private void validateImage(float[][] image, int width, int height) {
        if (image == null) {
            throw new IllegalArgumentException("image is null");
        }
        if (width <= 0 || height <= 0) {
            throw new IllegalArgumentException(
                    "image dimensions must be positive, got width=" + width + " height=" + height);
        }
        if (image.length != height) {
            throw new IllegalArgumentException(
                    "image.length (" + image.length + ") != height (" + height + ")");
        }
        for (int y = 0; y < height; y++) {
            if (image[y] == null) {
                throw new IllegalArgumentException("image[" + y + "] is null");
            }
            if (image[y].length != width) {
                throw new IllegalArgumentException(
                        "image[" + y + "].length (" + image[y].length + ") != width (" + width + ")");
            }
        }
    }
}
