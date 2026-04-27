package com.pranav.grrt.renderer;

import com.pranav.grrt.camera.Camera;
import com.pranav.grrt.disk.Disk;
import com.pranav.grrt.disk.DiskEmissionShader;
import com.pranav.grrt.disk.NovikovThorneDisk;
import com.pranav.grrt.integrator.DormandPrince45;
import com.pranav.grrt.integrator.RK4;
import com.pranav.grrt.metric.KerrMetric;
import com.pranav.grrt.metric.SchwarzschildMetric;
import org.junit.jupiter.api.Test;

import java.util.function.Supplier;

import static org.junit.jupiter.api.Assertions.*;

class RendererTest {

    // ---------------------------------------------------------------
    // a. The centermost pixel should be captured (it points directly
    //    at the black hole). Run at small resolution for speed.
    // ---------------------------------------------------------------

    @Test
    void centerPixelIsCaptured() {
        double M = 1.0;
        double rObs = 100.0;
        int W = 32, H = 32;
        double fov = 0.3;  // wide enough that center pixel is well inside shadow

        SchwarzschildMetric metric = new SchwarzschildMetric(M);
        Camera cam = new Camera(metric, rObs, Math.PI / 2, W, H, fov);
        Renderer r = new Renderer(cam, metric, RK4::new, RenderConfig.defaults(rObs));

        float[][] img = r.render();

        // "Central" pixels — the true optical center lies between pixels for
        // even resolution, so check the 2x2 block straddling the center.
        for (int dj = 0; dj < 2; dj++) {
            for (int di = 0; di < 2; di++) {
                float v = img[H / 2 - 1 + dj][W / 2 - 1 + di];
                assertEquals(0.0f, v, 0.0f,
                        "center pixel (" + (W/2 - 1 + di) + "," + (H/2 - 1 + dj)
                                + ") should be captured: " + v);
            }
        }
    }

    // ---------------------------------------------------------------
    // b. Corner pixels (largest |α|, |β|) should escape cleanly at any
    //    fov wide enough that they're well outside b_crit.
    // ---------------------------------------------------------------

    @Test
    void cornerPixelsEscape() {
        double M = 1.0;
        double rObs = 1000.0;
        int W = 32, H = 32;
        double fov = 0.05;  // shadow radius ≈ 5.2e-3 rad; corner pixels are at ~3.5e-2 rad

        SchwarzschildMetric metric = new SchwarzschildMetric(M);
        Camera cam = new Camera(metric, rObs, Math.PI / 2, W, H, fov);
        Renderer r = new Renderer(cam, metric, RK4::new, RenderConfig.defaults(rObs));

        float[][] img = r.render();

        assertEquals(1.0f, img[0][0],             0.0f, "top-left");
        assertEquals(1.0f, img[0][W - 1],         0.0f, "top-right");
        assertEquals(1.0f, img[H - 1][0],         0.0f, "bottom-left");
        assertEquals(1.0f, img[H - 1][W - 1],     0.0f, "bottom-right");
    }

    // ---------------------------------------------------------------
    // c. Phase 1 validation milestone: at 256×256, the shadow boundary
    //    along the middle row should match b_crit = 3√3 M within ±2
    //    pixels (edge-bracketing, per R3).
    //
    //    Rationale for edge-bracketing: the region at |α| ≈ b_crit/r_obs
    //    is numerically chaotic — rays with impact parameter b ≈ b_crit
    //    spend many orbits on the unstable photon sphere, and tiny
    //    round-off differences flip the CAPTURED/ESCAPED classification.
    //    Demanding a single monotonic 0→1 transition is brittle. Instead
    //    we locate the outermost CAPTURED pixel and the innermost
    //    ESCAPED pixel past it, and assert the midpoint is within ±2
    //    pixels of the analytic prediction.
    // ---------------------------------------------------------------

    @Test
    void shadowRadiusAt256MatchesCriticalImpactParameter() {
        double M = 1.0;
        double rObs = 1000.0;
        int W = 256, H = 256;
        double fov = 0.02;  // shadow diameter ≈ 0.01 rad → 2× margin

        SchwarzschildMetric metric = new SchwarzschildMetric(M);
        Camera cam = new Camera(metric, rObs, Math.PI / 2, W, H, fov);
        Renderer r = new Renderer(cam, metric, RK4::new, RenderConfig.defaults(rObs));

        float[][] img = r.render();

        // Scan the middle row to the right of center.
        int jMid = H / 2;
        int lastCaptured = -1;
        for (int i = W / 2; i < W; i++) {
            if (img[jMid][i] < 0.5f) lastCaptured = i;
        }
        assertTrue(lastCaptured >= 0, "no captured pixel found on midrow");

        int firstEscaped = -1;
        for (int i = lastCaptured + 1; i < W; i++) {
            if (img[jMid][i] >= 0.5f) {
                firstEscaped = i;
                break;
            }
        }
        assertTrue(firstEscaped > 0, "no escaped pixel past captured region on midrow");

        double midpointPx = 0.5 * (lastCaptured + firstEscaped);

        // Analytic: b_crit = 3√3 M corresponds to
        //   α_edge = asin(b_crit · √f / r_obs)
        // Camera maps α → pixel via i = W/2 − 0.5 + α · W / fov.
        double bCrit   = 3.0 * Math.sqrt(3.0) * M;
        double f       = 1.0 - 2.0 * M / rObs;
        double alpha   = Math.asin(bCrit * Math.sqrt(f) / rObs);
        double predPx  = 0.5 * W - 0.5 + alpha * W / fov;

        assertEquals(predPx, midpointPx, 2.0,
                "shadow edge midpoint " + midpointPx + " vs predicted " + predPx
                        + " (last captured " + lastCaptured
                        + ", first escaped " + firstEscaped + ")");
    }

    // ---------------------------------------------------------------
    // e. Determinism: two renders with identical configuration produce
    //    byte-identical arrays despite running on the ForkJoin pool.
    // ---------------------------------------------------------------

    @Test
    void determinismTwoRendersAreIdentical() {
        SchwarzschildMetric metric = new SchwarzschildMetric(1.0);
        Camera cam = new Camera(metric, 500.0, Math.PI / 2, 64, 64, 0.03);
        Renderer r1 = new Renderer(cam, metric, RK4::new, RenderConfig.defaults(500.0));
        Renderer r2 = new Renderer(cam, metric, RK4::new, RenderConfig.defaults(500.0));

        float[][] a = r1.render();
        float[][] b = r2.render();

        assertEquals(a.length, b.length);
        for (int j = 0; j < a.length; j++) {
            assertArrayEquals(a[j], b[j], 0.0f, "row " + j);
        }
    }

    // ---------------------------------------------------------------
    // f. Smoke test: all pixels are 0.0 or 1.0 (binary shader), no NaNs.
    // ---------------------------------------------------------------

    @Test
    void smokeTestBinaryValuesAndFiniteness() {
        SchwarzschildMetric metric = new SchwarzschildMetric(1.0);
        Camera cam = new Camera(metric, 200.0, Math.PI / 2, 16, 16, 0.08);
        Renderer r = new Renderer(cam, metric, RK4::new, RenderConfig.defaults(200.0));

        float[][] img = r.render();
        int zeros = 0, ones = 0;
        for (float[] row : img) {
            for (float v : row) {
                assertTrue(Float.isFinite(v), "non-finite pixel: " + v);
                if (v == 0.0f) zeros++;
                else if (v == 1.0f) ones++;
                else fail("non-binary pixel: " + v);
            }
        }
        // Both classes should be represented.
        assertTrue(zeros > 0, "no captured pixels");
        assertTrue(ones > 0,  "no escaped pixels");
    }

    // ---------------------------------------------------------------
    // Phase 3B gate 3: face-on Schwarzschild disk at 256² should
    // produce a centroid-symmetric annulus.
    //
    // Schwarzschild via KerrMetric(M, 0). Inclination is set to 1e-4
    // rad (effectively face-on); θ = 0 exactly is a Boyer-Lindquist
    // coordinate singularity in Camera (1/sin(θ) tetrad factor) and
    // cannot be used directly. The 1e-4 rad offset shifts the centroid
    // by ~r_avg · sin(i) ≈ 13 · 1e-4 = 1.3e-3 M, which is ~0.007 px at
    // our pixel scale — well below the 0.5 px tolerance.
    //
    // Disk r ∈ [r_ISCO + 1e-3, 20 M]; r_ISCO = 6 M for Schwarzschild.
    // ---------------------------------------------------------------

    @Test
    void faceOnSchwarzschildDiskCentroidIsNearOrigin() {
        double M = 1.0;
        double rObs = 100.0;
        double inclination = 1.0e-4;
        int W = 256, H = 256;
        double fov = 0.5;

        KerrMetric metric = new KerrMetric(M, 0.0);   // Schwarzschild via Kerr a=0
        Disk disk = new NovikovThorneDisk(metric, 0.0, 20.0, 1.0e-3);
        DiskEmissionShader shader = new DiskEmissionShader(disk, metric);

        // RenderConfig fields hNear/hFar/rTransition are ignored by
        // AdaptiveRayTracer; they only affect FixedStepRayTracer.
        RenderConfig config = new RenderConfig(
                0.01, 2.0 * rObs,
                0.01, 0.5, 10.0,
                1_000_000,
                shader);

        Camera cam = new Camera(metric, rObs, inclination, W, H, fov);

        double atol = 1.0e-8, rtol = 1.0e-8, hInit = 1.0;
        Supplier<RayTracer> factory = () -> new AdaptiveRayTracer(
                new DormandPrince45(), metric, config, atol, rtol, hInit, disk);

        Renderer r = Renderer.withRayTracer(cam, metric, factory, config);
        float[][] img = r.render();

        // Centroid of pixels with positive intensity (= disk hits).
        double sumX = 0.0, sumY = 0.0;
        int hitCount = 0;
        for (int j = 0; j < H; j++) {
            for (int i = 0; i < W; i++) {
                if (img[j][i] > 0.0f) {
                    sumX += i;
                    sumY += j;
                    hitCount++;
                }
            }
        }
        assertTrue(hitCount > 1000,
                "expected substantial disk-hit count for the face-on annulus; got " + hitCount);

        double cx = sumX / hitCount;
        double cy = sumY / hitCount;
        double centerX = (W - 1) / 2.0;   // 127.5 for W=256
        double centerY = (H - 1) / 2.0;
        double dx = cx - centerX;
        double dy = cy - centerY;

        System.out.println("[gate3] face-on Schwarzschild disk render");
        System.out.println("[gate3] hit pixel count = " + hitCount);
        System.out.println("[gate3] centroid = (" + cx + ", " + cy + ")");
        System.out.println("[gate3] |Δ| from image center = (" + dx + ", " + dy + ")");

        assertTrue(Math.abs(dx) < 0.5, "centroid X off-center: " + dx);
        assertTrue(Math.abs(dy) < 0.5, "centroid Y off-center: " + dy);
    }

    // ---------------------------------------------------------------
    // Config validation sanity checks.
    // ---------------------------------------------------------------

    @Test
    void renderConfigValidatesArguments() {
        RayShader shader = new BinaryShader();
        assertThrows(IllegalArgumentException.class,
                () -> new RenderConfig(0, 100, 0.01, 0.5, 10.0, 1000, shader));
        assertThrows(IllegalArgumentException.class,
                () -> new RenderConfig(0.01, 100, 0.01, 0.5, 10.0, 1000, null));
        assertThrows(IllegalArgumentException.class,
                () -> new RenderConfig(0.01, 100, -0.01, 0.5, 10.0, 1000, shader));
        assertThrows(IllegalArgumentException.class,
                () -> new RenderConfig(0.01, 100, 0.01, 0.5, 10.0, 0, shader));
    }
}
