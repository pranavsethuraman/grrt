package com.pranav.grrt.camera;

import com.pranav.grrt.metric.Metric;
import com.pranav.grrt.metric.SchwarzschildMetric;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

/**
 * Algebraic Camera tests. No integration of geodesics — those live in
 * {@link CameraDeflectionTest}.
 */
class CameraTest {

    // ---------------------------------------------------------------
    // 1. Center ray: α = β = 0 gives purely radial inward photon.
    // ---------------------------------------------------------------

    @Test
    void centerRayIsPurelyRadialInward() {
        double M = 1.0;
        double rObs = 100.0;
        double f = 1.0 - 2.0 * M / rObs;
        Camera cam = new Camera(new SchwarzschildMetric(M),
                                rObs, Math.PI / 2, 64, 64, 0.02);

        double[] y = new double[8];
        cam.initialStateFromAngles(0.0, 0.0, y);

        assertEquals(0.0,              y[0], 0.0, "t");
        assertEquals(rObs,             y[1], 0.0, "r");
        assertEquals(Math.PI / 2,      y[2], 0.0, "θ");
        assertEquals(0.0,              y[3], 0.0, "φ");
        assertEquals(1.0 / Math.sqrt(f), y[4], 1e-15, "k^t");
        assertEquals(-Math.sqrt(f),      y[5], 1e-15, "k^r inward");
        assertEquals(0.0,              y[6], 0.0, "k^θ");
        assertEquals(0.0,              y[7], 0.0, "k^φ");
    }

    // ---------------------------------------------------------------
    // 2. Null condition holds for every pixel on a grid sweep, at
    //    several inclinations.
    // ---------------------------------------------------------------

    @Test
    void nullConditionHoldsOnPixelGrid() {
        Metric m = new SchwarzschildMetric(1.0);
        double[] inclinations = {Math.PI / 4, Math.PI / 2, 2.0 * Math.PI / 3};

        double[] y = new double[8];
        double[] x = new double[4];
        double[] k = new double[4];

        for (double inc : inclinations) {
            Camera cam = new Camera(m, 500.0, inc, 32, 32, 0.05);
            double maxNorm = 0.0;
            for (int i = 0; i < cam.width(); i++) {
                for (int j = 0; j < cam.height(); j++) {
                    cam.initialState(i, j, y);
                    System.arraycopy(y, 0, x, 0, 4);
                    System.arraycopy(y, 4, k, 0, 4);
                    double n = Math.abs(m.nullNorm(x, k));
                    if (n > maxNorm) maxNorm = n;
                }
            }
            assertTrue(maxNorm < 1e-12,
                    "null-norm violated at inclination " + inc + ": " + maxNorm);
        }
    }

    // ---------------------------------------------------------------
    // 3. Tetrad orthonormality: extract e_(a)^μ by probing (α, β)
    //    and verify g_μν e_(a)^μ e_(b)^ν = η_(a)(b) to 1e-12.
    // ---------------------------------------------------------------

    @Test
    void tetradIsOrthonormalAtObserver() {
        Metric m = new SchwarzschildMetric(1.0);
        double rObs = 100.0;
        double theta = 1.1;   // deliberately off-equator
        Camera cam = new Camera(m, rObs, theta, 64, 64, 0.02);

        double[] y = new double[8];
        double f = 1.0 - 2.0 / rObs;

        // e_(t̂) is the t-component of any ray (spatial zero).
        double[] et = { 1.0 / Math.sqrt(f), 0, 0, 0 };

        // At (α=0, β=0): k_spatial = −e_(r̂).
        cam.initialStateFromAngles(0.0, 0.0, y);
        double[] er = { 0, -y[5], -y[6], -y[7] };

        // At (α=0, β=π/2): k_spatial = +e_(θ̂).
        cam.initialStateFromAngles(0.0, Math.PI / 2, y);
        double[] eth = { 0, y[5], y[6], y[7] };

        // At (α=π/2, β=0): k_spatial = +e_(φ̂).
        cam.initialStateFromAngles(Math.PI / 2, 0.0, y);
        double[] eph = { 0, y[5], y[6], y[7] };

        double[][] legs = { et, er, eth, eph };
        double[] etaDiag = { -1, 1, 1, 1 };
        double[][] gmn = m.g(new double[]{0.0, rObs, theta, 0.0});

        for (int a = 0; a < 4; a++) {
            for (int b = 0; b < 4; b++) {
                double inner = 0.0;
                for (int mu = 0; mu < 4; mu++) {
                    for (int nu = 0; nu < 4; nu++) {
                        inner += gmn[mu][nu] * legs[a][mu] * legs[b][nu];
                    }
                }
                double expected = (a == b) ? etaDiag[a] : 0.0;
                assertEquals(expected, inner, 1e-12,
                        "g(e_" + a + ", e_" + b + ") = " + inner);
            }
        }
    }

    // ---------------------------------------------------------------
    // 4. Reflection symmetry: pixel (i, j) and (W-1-i, j) produce
    //    4-momenta related by k^φ → −k^φ; other components unchanged.
    // ---------------------------------------------------------------

    @Test
    void reflectionSymmetryAcrossVerticalCenter() {
        Metric m = new SchwarzschildMetric(1.0);
        Camera cam = new Camera(m, 500.0, Math.PI / 2, 64, 32, 0.02);

        double[] left  = new double[8];
        double[] right = new double[8];

        for (int j = 0; j < cam.height(); j++) {
            for (int i = 0; i < cam.width() / 2; i++) {
                int iMirror = cam.width() - 1 - i;
                cam.initialState(i,       j, left);
                cam.initialState(iMirror, j, right);

                assertEquals(left[0], right[0], 0.0, "t");
                assertEquals(left[1], right[1], 0.0, "r");
                assertEquals(left[2], right[2], 0.0, "θ");
                assertEquals(left[3], right[3], 0.0, "φ");
                assertEquals(left[4], right[4], 1e-15, "k^t");
                assertEquals(left[5], right[5], 1e-15, "k^r");
                assertEquals(left[6], right[6], 1e-15, "k^θ");
                assertEquals(left[7], -right[7], 1e-15, "k^φ should flip sign");
            }
        }
    }

    // ---------------------------------------------------------------
    // 5. Input validation.
    // ---------------------------------------------------------------

    @Test
    void rejectsInvalidConstructorArguments() {
        Metric m = new SchwarzschildMetric(1.0);
        assertThrows(IllegalArgumentException.class,
                () -> new Camera(null, 100, Math.PI / 2, 64, 64, 0.02),
                "null metric");
        assertThrows(IllegalArgumentException.class,
                () -> new Camera(m, 2.0, Math.PI / 2, 64, 64, 0.02),
                "rObs at horizon");
        assertThrows(IllegalArgumentException.class,
                () -> new Camera(m, 1.0, Math.PI / 2, 64, 64, 0.02),
                "rObs inside horizon");
        assertThrows(IllegalArgumentException.class,
                () -> new Camera(m, Double.POSITIVE_INFINITY, Math.PI / 2, 64, 64, 0.02),
                "non-finite rObs");
        assertThrows(IllegalArgumentException.class,
                () -> new Camera(m, 100, 0.0, 64, 64, 0.02),
                "inclination at pole");
        assertThrows(IllegalArgumentException.class,
                () -> new Camera(m, 100, Math.PI, 64, 64, 0.02),
                "inclination at south pole");
        assertThrows(IllegalArgumentException.class,
                () -> new Camera(m, 100, Math.PI / 2, 0, 64, 0.02),
                "zero width");
        assertThrows(IllegalArgumentException.class,
                () -> new Camera(m, 100, Math.PI / 2, 64, -1, 0.02),
                "negative height");
        assertThrows(IllegalArgumentException.class,
                () -> new Camera(m, 100, Math.PI / 2, 64, 64, 0.0),
                "zero fov");
        assertThrows(IllegalArgumentException.class,
                () -> new Camera(m, 100, Math.PI / 2, 64, 64, Math.PI),
                "fov at π");
    }

    @Test
    void rejectsInvalidPixelOrBufferArguments() {
        Metric m = new SchwarzschildMetric(1.0);
        Camera cam = new Camera(m, 100, Math.PI / 2, 64, 64, 0.02);

        assertThrows(IllegalArgumentException.class,
                () -> cam.initialState(-1, 0, new double[8]));
        assertThrows(IllegalArgumentException.class,
                () -> cam.initialState(64, 0, new double[8]));
        assertThrows(IllegalArgumentException.class,
                () -> cam.initialState(0, 64, new double[8]));
        assertThrows(IllegalArgumentException.class,
                () -> cam.initialState(0, 0, new double[7]));
    }
}
