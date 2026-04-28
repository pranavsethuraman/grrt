package com.pranav.grrt.analysis;

/**
 * Stateless reduction from a per-azimuthal-bin ring-radius array
 * (length {@code N}) to four scalar asymmetry diagnostics: mean
 * radius, RMS dispersion, peak-to-peak spread, and the absolute value
 * and argument of the discrete Fourier {@code m = 1} coefficient.
 *
 * <p>Used by {@code EpsilonSweep} (Phase 3C, not yet committed) to
 * produce a row of {@code output/sweep.csv} per rendered frame. See
 * {@code docs/phase-3c-plan.md} §4 for the column definitions and §7
 * for the validation gates.
 *
 * <p><b>Convention.</b> The Fourier coefficient is defined as
 * {@code c_k = (1/N) Σᵢ rᵢ exp(−2π i k i / N)} for {@code i ∈ [0, N)};
 * the {@code k = 1} mode is returned. With this convention, an input
 * {@code rᵢ = R₀ (1 + A cos(2π i / N + φ₀))} yields
 * {@code c₁ = (R₀ A / 2) · exp(i φ₀)}, so {@code |c₁| = R₀ A / 2} and
 * {@code arg(c₁) = φ₀}.
 *
 * <p><b>NaN bins.</b> Bins whose radius is {@code Double.NaN}
 * (interpretation: "no peak found in this bin") are excluded from the
 * mean, RMS dispersion, peak-to-peak spread, and {@code validBins}
 * count, but are replaced by the mean before the Fourier sum so the
 * {@code k = 1} integral is computed over the full {@code N}-bin
 * index. A NaN bin replaced by the mean contributes only to the
 * {@code k = 0} (DC) mode, so the {@code k = 1} amplitude is unbiased
 * to first order. An all-NaN input returns a {@code Result} with
 * {@code validBins = 0} and every scalar set to {@code Double.NaN}.
 *
 * <p>Time complexity: {@code O(N)} in two passes. Allocation: one
 * {@code Result} per call; no scratch arrays.
 */
public final class CircularityMetric {

    /**
     * Reduction summary. Radial quantities ({@code meanR},
     * {@code deltaRrms}, {@code deltaRp2p}, {@code fourierM1Amplitude})
     * carry the same units as the input array — typically
     * gravitational radii {@code M}.
     *
     * @param meanR              mean radius over valid bins; {@code NaN} if no bins are valid
     * @param deltaRrms          RMS dispersion {@code sqrt((1/V) Σ (rᵢ − ⟨r⟩)²)} where the sum runs over the {@code V = validBins} non-NaN entries; population SD (divisor V, not V−1)
     * @param deltaRp2p          peak-to-peak spread {@code max(rᵢ) − min(rᵢ)} over valid bins
     * @param fourierM1Amplitude {@code |c₁|} per the convention in the class Javadoc
     * @param fourierM1Phase     {@code arg(c₁) ∈ (−π, π]}; radians; defaults to {@code 0} when the amplitude is identically zero (Java's {@code Math.atan2(0, 0)} convention)
     * @param validBins          number of non-NaN entries in the input
     */
    public record Result(double meanR,
                         double deltaRrms,
                         double deltaRp2p,
                         double fourierM1Amplitude,
                         double fourierM1Phase,
                         int validBins) {
    }

    private CircularityMetric() {
    }

    /**
     * Reduce the per-bin radius array to a {@link Result}.
     *
     * @param perBinRadii per-bin radius array, length {@code N ≥ 1};
     *                    {@code Double.NaN} entries are interpreted as
     *                    "no peak found in this bin" (see class
     *                    Javadoc for handling)
     * @return the reduction; never {@code null}
     * @throws IllegalArgumentException if {@code perBinRadii} is
     *         {@code null}, empty, or contains a non-NaN infinite
     *         value
     */
    public static Result compute(double[] perBinRadii) {
        if (perBinRadii == null) {
            throw new IllegalArgumentException("perBinRadii is null");
        }
        if (perBinRadii.length == 0) {
            throw new IllegalArgumentException("perBinRadii is empty");
        }
        final int n = perBinRadii.length;

        // Pass 1: count valid bins, mean, min/max
        int validBins = 0;
        double sum = 0.0;
        double minR = Double.POSITIVE_INFINITY;
        double maxR = Double.NEGATIVE_INFINITY;
        for (int i = 0; i < n; i++) {
            double r = perBinRadii[i];
            if (Double.isNaN(r)) {
                continue;
            }
            if (Double.isInfinite(r)) {
                throw new IllegalArgumentException(
                        "perBinRadii[" + i + "] is infinite (" + r + ")");
            }
            validBins++;
            sum += r;
            if (r < minR) minR = r;
            if (r > maxR) maxR = r;
        }
        if (validBins == 0) {
            return new Result(Double.NaN, Double.NaN, Double.NaN,
                    Double.NaN, Double.NaN, 0);
        }
        final double mean = sum / validBins;

        // Pass 2: RMS dispersion (over valid bins only) and Fourier m=1
        // coefficient (over all N bins, NaN bins replaced by mean so
        // they contribute only to the DC mode and not to k=1).
        double sse = 0.0;
        double cReal = 0.0;
        double cImag = 0.0;
        final double twoPiOverN = 2.0 * Math.PI / n;
        for (int i = 0; i < n; i++) {
            double r = perBinRadii[i];
            double rEff;
            if (Double.isNaN(r)) {
                rEff = mean;
            } else {
                rEff = r;
                double dev = r - mean;
                sse += dev * dev;
            }
            double angle = -twoPiOverN * i;
            cReal += rEff * Math.cos(angle);
            cImag += rEff * Math.sin(angle);
        }
        cReal /= n;
        cImag /= n;

        double fourierM1Amplitude = Math.hypot(cReal, cImag);
        double fourierM1Phase = Math.atan2(cImag, cReal);
        double deltaRrms = Math.sqrt(sse / validBins);
        double deltaRp2p = maxR - minR;
        return new Result(mean, deltaRrms, deltaRp2p,
                fourierM1Amplitude, fourierM1Phase, validBins);
    }
}
