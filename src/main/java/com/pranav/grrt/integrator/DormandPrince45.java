package com.pranav.grrt.integrator;

import com.pranav.grrt.metric.Metric;

/**
 * Dormand-Prince 5(4) adaptive integrator with First Same As Last (FSAL)
 * and PI step-size control, for null-geodesic state vectors.
 *
 * <p>Tableau source: Hairer, Nørsett, Wanner,
 * <i>Solving Ordinary Differential Equations I: Nonstiff Problems</i>
 * (2nd ed., 1993), Table 5.2, p. 178. The embedded 4th-order estimator
 * {@code b4} gives a local error of order O(h⁵); the step is propagated
 * with the 5th-order {@code b5}. Local truncation error of the
 * propagated solution is O(h⁶); global error O(h⁵).
 *
 * <h2>Butcher tableau</h2>
 *
 * <pre>
 *   c       | A
 *   ------- | ---------------------------------------------------------------------
 *   0       |
 *   1/5     | 1/5
 *   3/10    | 3/40         9/40
 *   4/5     | 44/45       -56/15        32/9
 *   8/9     | 19372/6561  -25360/2187   64448/6561    -212/729
 *   1       | 9017/3168   -355/33       46732/5247    49/176      -5103/18656
 *   1       | 35/384       0            500/1113      125/192     -2187/6784     11/84
 *   --------+-----------------------------------------------------------------------------
 *   b5      | 35/384       0            500/1113      125/192     -2187/6784     11/84     0
 *   b4      | 5179/57600   0            7571/16695    393/640     -92097/339200  187/2100  1/40
 *   e=b5-b4 | 71/57600     0           -71/16695      71/1920     -17253/339200  22/525   -1/40
 * </pre>
 *
 * <p>Tableau identities (verified in
 * {@link com.pranav.grrt.integrator.DormandPrince45Test}):
 * <ul>
 *   <li>&Sigma;_j A[i][j] == c[i] for stages i = 2..7 (row-sum
 *       consistency: each stage evaluates RHS at a consistent abscissa).</li>
 *   <li>&Sigma;_i b5[i] = 1 and &Sigma;_i b4[i] = 1 (order-1 condition).</li>
 *   <li>b5[7] = 0: the 5th-order solution does not use k7, so that
 *       k7 = f(y_{n+1}^{(5)}) becomes k1 of the next accepted step.
 *       This is the FSAL property.</li>
 *   <li>e[2] = 0: stage 2 contributes equally to b5 and b4.</li>
 * </ul>
 *
 * <h2>Error estimate and PI controller</h2>
 *
 * <p>Per-component error and RMS norm (Gustafsson 1991):
 * <pre>
 *   y_err[i] = h * &Sigma;_j e[j] k_j[i]
 *   sc[i]    = atol + rtol * |y[i]|
 *   err      = sqrt( (1/N) * &Sigma;_i (y_err[i] / sc[i])^2 )
 * </pre>
 * Accept if {@code err <= 1}.
 *
 * <p>PI step-size update:
 * <pre>
 *   h_new = h * safety * err^(-&alpha;) * err_prev^(&beta;)
 * </pre>
 * with {@code safety = 0.9}, {@code &alpha; = 0.7/5}, {@code &beta; = 0.4/5}
 * (standard DP5 PI gains). The step ratio h_new/h is clamped to
 * [0.2, 5]; on rejection the ratio is further clamped at 1 (no growth).
 * Absolute |h_new| is clamped to [hMin, hMax] (constructor arguments).
 *
 * <p>On rejection the controller falls back to a pure I law
 * {@code h_new = h * safety * err^(-1/5)} to shrink more aggressively,
 * and FSAL is invalidated (next k1 is recomputed from {@code y}).
 *
 * <h2>Dense output (continuous extension)</h2>
 *
 * <p>{@link #interpolate(double, double[])} returns the integrator's
 * 5th-order Hermite-style continuous extension over the most recently
 * <i>accepted</i> step, parameterized by θ ∈ [0, 1] where θ = 0
 * corresponds to the start of the step and θ = 1 to its end. The
 * continuous extension is the Dormand-Prince formula given by
 * Hairer-Nørsett-Wanner (1993) and used in Hairer's reference
 * Fortran/C code {@code dopri5.c} (function {@code contd5}). It
 * achieves global O(h⁵) accuracy on the interpolant — matching the
 * discrete-step order — using only the seven existing stage values
 * {@code k_1..k_7}. No additional RHS evaluation is required.
 *
 * <p>The constants {@code D1, D3, D4, D5, D6, D7} below are exact
 * rational coefficients lifted verbatim from Hairer's reference code.
 *
 * <p>Phase 3B.1 wires this method through the
 * {@code DormandPrince45Test.kerrInterpolateMatchesHalfStepDirect}
 * gate at 1e-8 mid-step accuracy. Phase 3B.2 uses it in
 * {@code AdaptiveRayTracer} for disk-crossing bisection.
 *
 * <h2>Thread-safety</h2>
 *
 * <p>This integrator holds per-step scratch buffers plus persistent FSAL
 * {@code k7}, controller {@code err_prev}, and dense-output state
 * ({@code yPrev}, {@code yNext}, {@code hAccepted}, {@code denseValid}),
 * and is <b>not</b> thread-safe. Construct one instance per worker
 * thread.
 *
 * <h2>FSAL contract</h2>
 *
 * <p>FSAL is valid only when the caller loops in the standard pattern
 * (pass the previously-accepted output as the next {@code y}). If the
 * caller restarts from a different state (e.g. between unrelated ray
 * traces), call {@link #resetState()} first to avoid reusing a stale k7.
 */
public final class DormandPrince45 implements AdaptiveIntegrator {

    // --- Butcher tableau constants.
    // Package-private so the scalar test driver can reuse exactly the
    // same coefficients (test (e) does not replicate them).

    static final double A21 =  1.0 / 5.0;

    static final double A31 =  3.0 / 40.0;
    static final double A32 =  9.0 / 40.0;

    static final double A41 =  44.0 / 45.0;
    static final double A42 = -56.0 / 15.0;
    static final double A43 =  32.0 / 9.0;

    static final double A51 =  19372.0 / 6561.0;
    static final double A52 = -25360.0 / 2187.0;
    static final double A53 =  64448.0 / 6561.0;
    static final double A54 =   -212.0 / 729.0;

    static final double A61 =  9017.0 / 3168.0;
    static final double A62 =  -355.0 / 33.0;
    static final double A63 = 46732.0 / 5247.0;
    static final double A64 =    49.0 / 176.0;
    static final double A65 = -5103.0 / 18656.0;

    // Stage-7 row equals b5 (without trailing 0): FSAL.
    static final double A71 =    35.0 / 384.0;
    static final double A73 =   500.0 / 1113.0;
    static final double A74 =   125.0 / 192.0;
    static final double A75 = -2187.0 / 6784.0;
    static final double A76 =    11.0 / 84.0;

    static final double B5_1 =    35.0 / 384.0;
    static final double B5_3 =   500.0 / 1113.0;
    static final double B5_4 =   125.0 / 192.0;
    static final double B5_5 = -2187.0 / 6784.0;
    static final double B5_6 =    11.0 / 84.0;
    // B5_2 = B5_7 = 0 (not stored)

    static final double B4_1 =   5179.0 / 57600.0;
    static final double B4_3 =   7571.0 / 16695.0;
    static final double B4_4 =    393.0 / 640.0;
    static final double B4_5 = -92097.0 / 339200.0;
    static final double B4_6 =    187.0 / 2100.0;
    static final double B4_7 =      1.0 / 40.0;

    // e = b5 - b4 (error vector).
    static final double E1 = B5_1 - B4_1;     //  71/57600
    static final double E3 = B5_3 - B4_3;     // -71/16695
    static final double E4 = B5_4 - B4_4;     //  71/1920
    static final double E5 = B5_5 - B4_5;     // -17253/339200
    static final double E6 = B5_6 - B4_6;     //  22/525
    static final double E7 = 0.0   - B4_7;    // -1/40

    static final double C2 = 1.0 / 5.0;
    static final double C3 = 3.0 / 10.0;
    static final double C4 = 4.0 / 5.0;
    static final double C5 = 8.0 / 9.0;
    static final double C6 = 1.0;
    static final double C7 = 1.0;

    // --- Dense-output ("contd5") coefficients from Hairer-Nørsett-Wanner
    //     (1993), implemented as in Hairer's reference Fortran/C code
    //     (dopri5.f / dopri5.c, function contd5). The formula is
    //
    //       y_dense(θ) = rcont1 + θ ( rcont2 + (1−θ)(rcont3 + θ(rcont4 + (1−θ) rcont5)) )
    //
    //     with
    //
    //       rcont1 = y_n
    //       rcont2 = y_{n+1} − y_n
    //       rcont3 = h k_1 − rcont2
    //       rcont4 = rcont2 − h k_7 − rcont3
    //       rcont5 = h ( D1 k_1 + D3 k_3 + D4 k_4 + D5 k_5 + D6 k_6 + D7 k_7 )

    static final double D1 = -12715105075.0 / 11282082432.0;
    static final double D3 =   87487479700.0 / 32700410799.0;
    static final double D4 =  -10690763975.0 /  1880347072.0;
    static final double D5 =  701980252875.0 / 199316789632.0;
    static final double D6 =   -1453857185.0 /    822651844.0;
    static final double D7 =      69997945.0 /     29380423.0;

    // --- PI controller constants.
    static final double SAFETY     = 0.9;
    static final double ALPHA      = 0.7 / 5.0;    // 0.14
    static final double BETA       = 0.4 / 5.0;    // 0.08
    static final double MIN_GROWTH = 0.2;
    static final double MAX_GROWTH = 5.0;
    private static final double MIN_ERR = 1e-10;   // floor to keep PI stable

    // --- Per-instance scratch buffers and persistent state.

    private final double[] k1 = new double[8];
    private final double[] k2 = new double[8];
    private final double[] k3 = new double[8];
    private final double[] k4 = new double[8];
    private final double[] k5 = new double[8];
    private final double[] k6 = new double[8];
    private final double[] k7 = new double[8];
    private final double[] yTmp = new double[8];
    private final double[] xBuf = new double[4];
    private final double[] kBuf = new double[4];
    private final double[] accel = new double[4];

    // --- Dense-output persistent state.
    //
    // Populated only on accepted steps. yPrev = y_n, yNext = y_{n+1},
    // hAccepted = the propagating h that was accepted. denseValid is
    // false until the first accepted step, after a rejection, or after
    // resetState().

    private final double[] yPrev = new double[8];
    private final double[] yNext = new double[8];
    private double hAccepted = 0.0;
    private boolean denseValid = false;

    private final double hMin;
    private final double hMax;

    private double errPrev = 1.0;
    private boolean fsalValid = false;

    /** Default step caps: hMin = 1e-12, hMax = 1e6 (loose; ample for ray tracing). */
    public DormandPrince45() {
        this(1e-12, 1e6);
    }

    /**
     * @param hMin minimum absolute step size (floor); must be &gt; 0
     * @param hMax maximum absolute step size (ceiling); must be &gt; hMin
     */
    public DormandPrince45(double hMin, double hMax) {
        if (!(hMin > 0.0) || !(hMax > hMin) || !Double.isFinite(hMax)) {
            throw new IllegalArgumentException(
                    "require 0 < hMin < hMax, both finite: hMin="
                            + hMin + ", hMax=" + hMax);
        }
        this.hMin = hMin;
        this.hMax = hMax;
    }

    /**
     * Invalidate the stored FSAL k7, the dense-output buffers, and the
     * PI controller's previous-error memory. Call this when restarting
     * integration from a state unrelated to the previously-accepted
     * output (e.g. between independent ray traces, or after an
     * axis-reflection jump in the renderer loop).
     */
    @Override
    public void resetState() {
        fsalValid = false;
        denseValid = false;
        errPrev = 1.0;
    }

    @Override
    public StepStatus adaptiveStep(Metric metric, double[] y, double hTry,
                                   double atol, double rtol, double[] out) {
        if (y.length != 8 || out.length != 8) {
            throw new IllegalArgumentException("state vectors must have length 8");
        }
        if (!(atol > 0.0) || !(rtol > 0.0)) {
            throw new IllegalArgumentException(
                    "atol, rtol must be positive: atol=" + atol + ", rtol=" + rtol);
        }
        if (hTry == 0.0 || !Double.isFinite(hTry)) {
            throw new IllegalArgumentException("hTry must be finite and nonzero: " + hTry);
        }

        // Snapshot y_n for dense output BEFORE the step is propagated.
        // Cheap (System.arraycopy is allocation-free); harmless if the
        // step is later rejected — denseValid gates whether yPrev is
        // meaningful.
        System.arraycopy(y, 0, yPrev, 0, 8);

        // Clamp input hTry magnitude to [hMin, hMax], preserving sign.
        double h = clampStep(hTry);

        // Stage 1: k1 = f(y_n), reusing FSAL k7 if the previous step
        // was accepted and the caller has not reset.
        if (fsalValid) {
            System.arraycopy(k7, 0, k1, 0, 8);
        } else {
            rhs(metric, y, k1);
        }

        // Stage 2
        for (int i = 0; i < 8; i++) yTmp[i] = y[i] + h * A21 * k1[i];
        rhs(metric, yTmp, k2);

        // Stage 3
        for (int i = 0; i < 8; i++)
            yTmp[i] = y[i] + h * (A31 * k1[i] + A32 * k2[i]);
        rhs(metric, yTmp, k3);

        // Stage 4
        for (int i = 0; i < 8; i++)
            yTmp[i] = y[i] + h * (A41 * k1[i] + A42 * k2[i] + A43 * k3[i]);
        rhs(metric, yTmp, k4);

        // Stage 5
        for (int i = 0; i < 8; i++)
            yTmp[i] = y[i] + h * (A51 * k1[i] + A52 * k2[i] + A53 * k3[i] + A54 * k4[i]);
        rhs(metric, yTmp, k5);

        // Stage 6
        for (int i = 0; i < 8; i++)
            yTmp[i] = y[i] + h * (A61 * k1[i] + A62 * k2[i] + A63 * k3[i]
                                  + A64 * k4[i] + A65 * k5[i]);
        rhs(metric, yTmp, k6);

        // Stage 7: input coincides with y5 (FSAL). yTmp now holds y5.
        for (int i = 0; i < 8; i++)
            yTmp[i] = y[i] + h * (A71 * k1[i] + A73 * k3[i] + A74 * k4[i]
                                  + A75 * k5[i] + A76 * k6[i]);
        rhs(metric, yTmp, k7);

        // Per-component error: y_err = h * Σ e_i k_i. RMS norm on scaled
        // components. N = 8 (full 8-state); no per-component disabling.
        double sqSum = 0.0;
        for (int i = 0; i < 8; i++) {
            double yErr = h * (E1 * k1[i] + E3 * k3[i] + E4 * k4[i]
                              + E5 * k5[i] + E6 * k6[i] + E7 * k7[i]);
            double sc = atol + rtol * Math.abs(y[i]);
            double ratio = yErr / sc;
            sqSum += ratio * ratio;
        }
        double err = Math.sqrt(sqSum / 8.0);
        err = Math.max(err, MIN_ERR);

        if (err <= 1.0) {
            // Accepted: copy y5 to out, keep k7 for FSAL, update controller.
            System.arraycopy(yTmp, 0, out, 0, 8);
            fsalValid = true;

            // Persist dense-output state for this accepted step.
            System.arraycopy(yTmp, 0, yNext, 0, 8);
            hAccepted = h;
            denseValid = true;

            double factor = SAFETY
                    * Math.pow(err,     -ALPHA)
                    * Math.pow(errPrev,  BETA);
            factor = Math.max(MIN_GROWTH, Math.min(MAX_GROWTH, factor));
            double hNew = clampStep(h * factor);

            errPrev = err;
            return new StepStatus(true, hNew, err);
        } else {
            // Rejected: fall back to pure I (no β term) for a firmer shrink,
            // cap growth at 1 (never grow on rejection), invalidate FSAL,
            // and invalidate dense output (the buffered yNext/hAccepted
            // belong to whatever previous accepted step the user holds).
            fsalValid = false;
            denseValid = false;

            double factor = SAFETY * Math.pow(err, -1.0 / 5.0);
            factor = Math.max(MIN_GROWTH, Math.min(1.0, factor));
            double hNew = clampStep(h * factor);

            return new StepStatus(false, hNew, err);
        }
    }

    /**
     * Continuous extension of the most recently accepted step,
     * evaluated at parameter θ ∈ [0, 1]. θ = 0 returns y_n
     * (the start of the step); θ = 1 returns y_{n+1} (the end).
     * Globally O(h⁵) accurate per Hairer-Nørsett-Wanner (1993)
     * Eq. II.6.20 / Hairer's {@code contd5} reference implementation.
     *
     * <p>No additional RHS evaluation is performed; the formula uses
     * only the seven stages {@code k_1..k_7} stored from the last
     * accepted step.
     *
     * @param theta interpolation parameter in [0, 1]
     * @param out   8-vector buffer; overwritten with y(θ)
     * @throws IllegalStateException    if no step has been accepted, or
     *                                  if the most recent step was
     *                                  rejected, or after
     *                                  {@link #resetState()}.
     * @throws IllegalArgumentException if {@code theta} is outside
     *                                  [0, 1] or {@code out.length != 8}
     */
    public void interpolate(double theta, double[] out) {
        if (out == null || out.length != 8) {
            throw new IllegalArgumentException("out must have length 8");
        }
        if (!denseValid) {
            throw new IllegalStateException(
                    "no accepted step available for dense output");
        }
        if (!(theta >= 0.0) || !(theta <= 1.0) || !Double.isFinite(theta)) {
            throw new IllegalArgumentException(
                    "theta must be in [0, 1]: " + theta);
        }
        double h = hAccepted;
        double oneMinusTheta = 1.0 - theta;
        for (int i = 0; i < 8; i++) {
            double rcont1 = yPrev[i];
            double rcont2 = yNext[i] - yPrev[i];
            double rcont3 = h * k1[i] - rcont2;
            double rcont4 = rcont2 - h * k7[i] - rcont3;
            double rcont5 = h * (D1 * k1[i] + D3 * k3[i] + D4 * k4[i]
                               + D5 * k5[i] + D6 * k6[i] + D7 * k7[i]);
            out[i] = rcont1
                    + theta * (rcont2
                    + oneMinusTheta * (rcont3
                    + theta * (rcont4
                    + oneMinusTheta * rcont5)));
        }
    }

    private double clampStep(double h) {
        double mag = Math.abs(h);
        if (mag < hMin) return Math.signum(h) * hMin;
        if (mag > hMax) return Math.signum(h) * hMax;
        return h;
    }

    /**
     * Geodesic RHS: f(y) into {@code dy}. Identical structure to
     * {@link RK4#rhs}: position derivatives are the momenta, momentum
     * derivatives come from {@link Metric#geodesicAcceleration}.
     */
    private void rhs(Metric metric, double[] y, double[] dy) {
        dy[0] = y[4]; dy[1] = y[5]; dy[2] = y[6]; dy[3] = y[7];
        xBuf[0] = y[0]; xBuf[1] = y[1]; xBuf[2] = y[2]; xBuf[3] = y[3];
        kBuf[0] = y[4]; kBuf[1] = y[5]; kBuf[2] = y[6]; kBuf[3] = y[7];
        metric.geodesicAcceleration(xBuf, kBuf, accel);
        dy[4] = accel[0]; dy[5] = accel[1]; dy[6] = accel[2]; dy[7] = accel[3];
    }
}
