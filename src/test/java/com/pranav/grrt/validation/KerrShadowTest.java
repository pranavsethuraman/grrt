package com.pranav.grrt.validation;

import com.pranav.grrt.camera.Camera;
import com.pranav.grrt.integrator.DormandPrince45;
import com.pranav.grrt.metric.KerrMetric;
import com.pranav.grrt.renderer.AdaptiveRayTracer;
import com.pranav.grrt.renderer.RayOutcome;
import com.pranav.grrt.renderer.RayShader;
import com.pranav.grrt.renderer.RayTracer;
import com.pranav.grrt.renderer.RenderConfig;
import com.pranav.grrt.renderer.Renderer;
import org.junit.jupiter.api.Test;

import javax.imageio.ImageIO;
import java.awt.image.BufferedImage;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.function.Supplier;

import static org.junit.jupiter.api.Assertions.*;

/**
 * C3/C4: render the Kerr shadow at a = 0.9 M, i = 90° with DP45 and
 * validate against the Bardeen (1973) analytic D_H under the
 * <b>Camera.java sign convention</b> (prograde on screen-right; see
 * {@link BardeenShadow#horizontalExtentCameraConvention}).
 *
 * <p>Two tests:
 * <ul>
 *   <li>{@link #kerrShadowSmoke256} — 256×256 fast smoke test with
 *       tolerance 1.0 · pixel_scale_M = 0.078125 M.</li>
 *   <li>{@link #kerrShadowPrimary1024} — 1024×1024 primary validation
 *       with tolerance 1.0 · pixel_scale_M = 0.01953125 M.</li>
 * </ul>
 * The 1.0 · pixel_scale floor reflects the worst-case quantization
 * error of integer-pixel band-width measurement when both edges fall
 * in the unfavourable half of their containing pixels — a 0.5 · pixel
 * bound cannot be satisfied without sub-pixel edge methods.
 *
 * <p>Convention derivation from Camera.java:
 * <ul>
 *   <li>{@code e_(φ̂)^μ = (0, 0, 0, 1/(r sinθ))} — points in +φ direction.</li>
 *   <li>{@code n̂^φ̂ = sin α · cos β} — positive α ⇒ positive n̂^φ̂
 *       ⇒ {@code k^φ > 0} at the observer.</li>
 *   <li>In Kerr equatorial with a &gt; 0, {@code k^φ > 0} carries the
 *       photon with the BH rotation direction — prograde.</li>
 * </ul>
 * Camera convention: positive screen-alpha ⇒ prograde on screen-right.
 * Bardeen 1973 uses α = −ξ/sin(i), placing prograde on screen-left
 * (opposite). {@link BardeenShadow#horizontalExtentCameraConvention}
 * inverts Bardeen's sign so the oracle matches the rendered image.
 *
 * <h2>Secondary asymmetry check — centroid offset, not bbox-edge</h2>
 *
 * <p>The Kerr shadow at inclination i = π/2 has algebraic <b>cusps</b>
 * at its horizontal tips: near each tip, β ∼ (δα)^(1/4) rather than the
 * parabolic (δα)^(1/2) expected at a generic turning point. Integer-
 * pixel measurements of vertical extent at the bbox-extreme columns are
 * therefore dominated by the sub-pixel position of the cusp inside the
 * bounding pixel — at W = 256 the two tips happen to land at similar
 * fractional offsets from their columns and the asymmetry (L−R)/H reads
 * −0.016; at W = 1024 the offsets differ by 13× and the measurement
 * <i>flips sign</i> to +0.027. The min/max edge-extent ratio is
 * unstable for the same reason (0.85 at 256, 0.36 at 1024). Neither
 * measurement can be asserted without sub-pixel edge extraction.
 *
 * <p>The asserted secondary is therefore the <b>centroid offset</b> of
 * the captured region, (c_x − W/2) / W. This quantity integrates over
 * every captured pixel and is dominated by the intrinsic asymmetry of
 * [α_min, α_max] = [−6.83 M, +2.84 M] (at a = 0.9 M, i = π/2) about
 * zero. It reads −0.0965 to four decimals at both W = 256 and W = 1024
 * — resolution-invariant, physically meaningful, and it flips sign
 * (+0.09×) under a g_tφ sign error or a Camera handedness inversion.
 * The required threshold is centroid_offset_frac &lt; −0.05 (observed
 * margin ≈ 2×).
 *
 * <p>The bbox-edge extents, the edge ratio, and the full bbox are
 * still <i>printed</i> as diagnostics to help eyeball the shadow
 * shape, but they are not asserted.
 */
class KerrShadowTest {

    private static final float SHADE_CAPTURED   = 0.0f;
    private static final float SHADE_ESCAPED    = 1.0f;
    private static final float SHADE_UNRESOLVED = 0.5f;

    // Binding setup — same for both tests.
    private static final double M    = 1.0;
    private static final double A    = 0.9;
    private static final double R_OBS = 1000.0;
    private static final double FOV   = 0.02;
    private static final double ATOL  = 1e-8;
    private static final double RTOL  = 1e-8;
    private static final double H_INIT = 1.0;

    /** Outcome-encoding shader so C3 can distinguish captured / escaped / unresolved. */
    private static final class OutcomeShader implements RayShader {
        @Override
        public float shade(RayOutcome outcome, double[] y) {
            return switch (outcome) {
                case CAPTURED  -> SHADE_CAPTURED;
                case ESCAPED   -> SHADE_ESCAPED;
                case MAX_STEPS -> SHADE_UNRESOLVED;
                case NAN       -> Float.NaN;
                case HIT_DISK  -> Float.NaN;   // not expected: this test has no Disk
            };
        }
    }

    // ------------------------------------------------------------------
    // Tests
    // ------------------------------------------------------------------

    @Test
    void kerrShadowSmoke256() throws IOException {
        MeasurementResult r = renderAndMeasure(256, 256, "C3-256",
                "kerr_shadow_a0p9_i90_256.png");
        assertAuxiliaryChecks(r, "256");
        double tolPx = 1.0;
        assertTrue(Math.abs(r.diffM) < tolPx * r.pixelScaleM,
                "256: |D_H_measured − D_H_analytic| = " + r.diffM
                        + " M exceeds ±" + (tolPx * r.pixelScaleM)
                        + " M (" + tolPx + " px)");
    }

    @Test
    void kerrShadowPrimary1024() throws IOException {
        MeasurementResult r = renderAndMeasure(1024, 1024, "C3-1024",
                "kerr_shadow_a0p9_i90_1024.png");
        assertAuxiliaryChecks(r, "1024");
        double tolPx = 1.0;
        assertTrue(Math.abs(r.diffM) < tolPx * r.pixelScaleM,
                "1024: |D_H_measured − D_H_analytic| = " + r.diffM
                        + " M exceeds ±" + (tolPx * r.pixelScaleM)
                        + " M (" + tolPx + " px)");
    }

    // ------------------------------------------------------------------
    // Helper: render + measure, parameterised by resolution
    // ------------------------------------------------------------------

    private static MeasurementResult renderAndMeasure(
            int W, int H, String label, String pngFilename) throws IOException {

        KerrMetric metric = new KerrMetric(M, A);
        Camera cam = new Camera(metric, R_OBS, Math.PI / 2.0, W, H, FOV);

        double rPlus = metric.horizonRadius();
        double horizonCushion = 0.01 * rPlus;         // rCapture = 1.01·r_+
        double escapeRadius   = 1.1 * R_OBS;          // 1100 M
        int    maxSteps       = 200_000;              // attempts; ≫ 10⁵ accepted

        RenderConfig config = new RenderConfig(
                horizonCushion, escapeRadius,
                0.01, 0.5, 10.0,                      // unused by AdaptiveRayTracer
                maxSteps,
                new OutcomeShader());

        var tracers = new ConcurrentLinkedQueue<AdaptiveRayTracer>();
        Supplier<RayTracer> factory = () -> {
            var dp = new DormandPrince45();
            var tracer = new AdaptiveRayTracer(dp, metric, config, ATOL, RTOL, H_INIT);
            tracers.add(tracer);
            return tracer;
        };

        Renderer renderer = Renderer.withRayTracer(cam, metric, factory, config);

        long t0 = System.nanoTime();
        float[][] image = renderer.render();
        double elapsedSec = (System.nanoTime() - t0) / 1e9;

        // --- Classify pixels ---
        int captured = 0, escaped = 0, unresolved = 0, nanCount = 0;
        for (int j = 0; j < H; j++) {
            for (int i = 0; i < W; i++) {
                float v = image[j][i];
                if (Float.isNaN(v))                nanCount++;
                else if (v == SHADE_CAPTURED)      captured++;
                else if (v == SHADE_ESCAPED)       escaped++;
                else if (v == SHADE_UNRESOLVED)    unresolved++;
            }
        }
        double unresolvedFrac = (double) unresolved / (W * H);

        long totalAccepted = 0, totalRejected = 0;
        for (AdaptiveRayTracer t : tracers) {
            totalAccepted += t.acceptedSteps();
            totalRejected += t.rejectedSteps();
        }

        // --- D_H via row-scan on captured pixels + shadow bounding box ---
        int dHpx = 0;
        int bboxXMin = W, bboxXMax = -1, bboxYMin = H, bboxYMax = -1;
        for (int j = 0; j < H; j++) {
            int xLeft = -1, xRight = -1;
            for (int i = 0; i < W; i++) {
                if (image[j][i] == SHADE_CAPTURED) {
                    if (xLeft < 0) xLeft = i;
                    xRight = i;
                    if (i < bboxXMin) bboxXMin = i;
                    if (i > bboxXMax) bboxXMax = i;
                    if (j < bboxYMin) bboxYMin = j;
                    if (j > bboxYMax) bboxYMax = j;
                }
            }
            if (xLeft >= 0) {
                int width = xRight - xLeft + 1;
                if (width > dHpx) dHpx = width;
            }
        }
        double pixelScaleM  = FOV * R_OBS / W;
        double dHMeasuredM  = dHpx * pixelScaleM;
        double dHAnalyticM  = BardeenShadow.horizontalExtentCameraConvention(
                M, A, Math.PI / 2.0, 10001).dH();
        double diffM        = dHMeasuredM - dHAnalyticM;
        double diffPx       = diffM / pixelScaleM;

        // --- Vertical extent at bbox edges (asymmetry + ratio) ---
        int vertLeft = 0, vertRight = 0, vertCenter = 0, maxColExtent = 0;
        if (bboxXMax >= 0) {
            for (int j = 0; j < H; j++) {
                if (image[j][bboxXMin] == SHADE_CAPTURED) vertLeft++;
                if (image[j][bboxXMax] == SHADE_CAPTURED) vertRight++;
            }
            int xMid = (bboxXMin + bboxXMax) / 2;
            for (int j = 0; j < H; j++) {
                if (image[j][xMid] == SHADE_CAPTURED) vertCenter++;
            }
            for (int i = bboxXMin; i <= bboxXMax; i++) {
                int n = 0;
                for (int j = 0; j < H; j++) if (image[j][i] == SHADE_CAPTURED) n++;
                if (n > maxColExtent) maxColExtent = n;
            }
        }
        double asymmetry = (double) (vertLeft - vertRight) / H;
        int minEdge = Math.min(vertLeft, vertRight);
        int maxEdge = Math.max(vertLeft, vertRight);
        double ratio = (maxEdge > 0) ? (double) minEdge / maxEdge : Double.NaN;

        // --- Centroid offset (robust resolution-invariant asymmetry) ---
        // Sum x-coordinates of captured pixels; compute centroid offset
        // from image center as a fraction of W. Under Camera convention
        // with the D-flat on the prograde (right) side, the shadow
        // extends further on the retrograde side, so the centroid is
        // left of center ⇒ centroidOffsetFrac < 0.
        double cx = 0.0;
        long capCount = 0;
        for (int j = 0; j < H; j++) {
            for (int i = 0; i < W; i++) {
                if (image[j][i] == SHADE_CAPTURED) {
                    cx += (i + 0.5);
                    capCount++;
                }
            }
        }
        double centroidX = (capCount > 0) ? (cx / capCount) : Double.NaN;
        double centroidOffsetFrac = (centroidX - 0.5 * W) / W;

        // --- Save PNG ---
        Path outDir = Path.of("output");
        Files.createDirectories(outDir);
        Path pngPath = outDir.resolve(pngFilename);
        savePng(image, pngPath);

        MeasurementResult result = new MeasurementResult(
                W, H, elapsedSec,
                captured, escaped, unresolved, nanCount, unresolvedFrac,
                totalAccepted, totalRejected,
                pixelScaleM, dHpx, dHMeasuredM, dHAnalyticM, diffM, diffPx,
                bboxXMin, bboxXMax, bboxYMin, bboxYMax,
                vertLeft, vertRight, vertCenter, maxColExtent,
                asymmetry, ratio, centroidOffsetFrac, pngPath);

        printReport(label, result);
        return result;
    }

    // ------------------------------------------------------------------
    // Auxiliary checks, same for 256 and 1024
    // ------------------------------------------------------------------

    private static void assertAuxiliaryChecks(MeasurementResult r, String label) {
        assertEquals(0, r.nanCount, label + ": unexpected NaN pixels");
        assertTrue(r.elapsedSec < 30.0,
                label + ": render wall time " + r.elapsedSec + " s exceeds 30 s hard-stop");
        assertTrue(r.unresolvedFrac < 0.001,
                label + ": unresolved fraction " + r.unresolvedFrac + " exceeds 0.1%");
        // Sign check for g_tφ correctness and Camera handedness. See class
        // Javadoc for why this supersedes the bbox-edge asymmetry / ratio
        // checks as the asserted secondary.
        assertTrue(r.centroidOffsetFrac < -0.05,
                label + ": centroid_offset_frac = " + r.centroidOffsetFrac
                        + " — under Camera convention must be < −0.05 (shadow"
                        + " shifted retrograde-ward, i.e. left of image center);"
                        + " ≥ 0 would indicate a g_tφ sign error or a Camera"
                        + " handedness inversion");
    }

    // ------------------------------------------------------------------
    // Report
    // ------------------------------------------------------------------

    private static void printReport(String label, MeasurementResult r) {
        System.out.println("[" + label + "] --- " + r.W + "×" + r.H + " ---");
        System.out.println("[" + label + "] render time           = "
                + String.format("%.2f s", r.elapsedSec));
        System.out.println("[" + label + "] captured_count        = " + r.captured);
        System.out.println("[" + label + "] escaped_count         = " + r.escaped);
        System.out.println("[" + label + "] unresolved_count      = " + r.unresolved
                + "   (" + String.format("%.4f%%", r.unresolvedFrac * 100) + ")");
        if (r.nanCount > 0)
            System.out.println("[" + label + "] NaN pixels            = " + r.nanCount);
        System.out.println("[" + label + "] total accepted DP45   = " + r.totalAccepted);
        System.out.println("[" + label + "] total rejected DP45   = " + r.totalRejected);
        System.out.println("[" + label + "] pixel_scale_M         = " + r.pixelScaleM);
        System.out.println("[" + label + "] D_H_px                = " + r.dHpx);
        System.out.println("[" + label + "] D_H_analytic_M        = "
                + String.format("%.10f M", r.dHAnalyticM));
        System.out.println("[" + label + "] D_H_measured_M        = "
                + String.format("%.10f M", r.dHMeasuredM));
        System.out.println("[" + label + "] difference            = "
                + String.format("%.6f M  (%.4f px)", r.diffM, r.diffPx));
        System.out.println("[" + label + "] tolerance (1.0 px)    = "
                + String.format("%.6f M", r.pixelScaleM));
        System.out.println("[" + label + "] bbox                  = x∈[" + r.bboxXMin
                + "," + r.bboxXMax + "] y∈[" + r.bboxYMin + "," + r.bboxYMax + "]");
        System.out.println("[" + label + "] vertical_extent_left  = " + r.vertLeft
                + " px  (bbox x_min)");
        System.out.println("[" + label + "] vertical_extent_right = " + r.vertRight
                + " px  (bbox x_max)");
        System.out.println("[" + label + "] vertical_center       = " + r.vertCenter
                + " px  (bbox midpoint)");
        System.out.println("[" + label + "] max_col_extent        = " + r.maxColExtent + " px");
        System.out.println("[" + label + "] asymmetry (L−R)/H     = "
                + String.format("%+.4f", r.asymmetry));
        System.out.println("[" + label + "] ratio min/max edge    = "
                + String.format("%.4f", r.ratio));
        System.out.println("[" + label + "] centroid_offset_frac  = "
                + String.format("%+.4f  (captured pixels, (cx-W/2)/W)", r.centroidOffsetFrac));
        System.out.println("[" + label + "] PNG saved to          = "
                + r.pngPath.toAbsolutePath());
    }

    private static void savePng(float[][] img, Path out) throws IOException {
        int H = img.length;
        int W = img[0].length;
        BufferedImage bi = new BufferedImage(W, H, BufferedImage.TYPE_INT_RGB);
        for (int j = 0; j < H; j++) {
            for (int i = 0; i < W; i++) {
                float v = img[j][i];
                int rgb;
                if (Float.isNaN(v)) {
                    rgb = 0xFF0000;                             // red for NaN
                } else if (v == SHADE_UNRESOLVED) {
                    rgb = 0x00C8FF;                             // cyan for unresolved
                } else {
                    int g = (int) Math.round(Math.max(0.0f, Math.min(1.0f, v)) * 255.0f);
                    rgb = (g << 16) | (g << 8) | g;
                }
                bi.setRGB(i, j, rgb);
            }
        }
        ImageIO.write(bi, "png", out.toFile());
    }

    // ------------------------------------------------------------------
    // Measurement result bundle
    // ------------------------------------------------------------------

    /** All quantities produced by one render + measurement pass. */
    private record MeasurementResult(
            int W, int H,
            double elapsedSec,
            int captured, int escaped, int unresolved, int nanCount,
            double unresolvedFrac,
            long totalAccepted, long totalRejected,
            double pixelScaleM,
            int dHpx, double dHMeasuredM, double dHAnalyticM,
            double diffM, double diffPx,
            int bboxXMin, int bboxXMax, int bboxYMin, int bboxYMax,
            int vertLeft, int vertRight, int vertCenter, int maxColExtent,
            double asymmetry, double ratio,
            double centroidOffsetFrac,
            Path pngPath) {}
}
