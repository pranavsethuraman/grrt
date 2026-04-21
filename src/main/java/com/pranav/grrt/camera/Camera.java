package com.pranav.grrt.camera;

import com.pranav.grrt.metric.Metric;

/**
 * Converts image-plane pixel coordinates to the initial 8-state of a null
 * geodesic traced outward from a static observer in Schwarzschild spacetime.
 *
 * <h2>Physical setup</h2>
 *
 * <p>The observer sits at Schwarzschild coordinates (t=0, r=rObs, θ=inclination,
 * φ=0), at rest with respect to the coordinate system. Its normalized
 * 4-velocity is u^μ = (1/√f, 0, 0, 0), where f ≡ 1 − 2M/r_obs. This static
 * approximation is appropriate for M87*-style imaging: Earth's peculiar
 * velocity relative to M87's rest frame is β ~ 10⁻³, which distorts apparent
 * angles by O(β) — well below EHT's angular resolution and below our
 * RK4 + finite-resolution error budget for Phase 1.
 *
 * <h2>Tetrad derivation</h2>
 *
 * <p>An orthonormal tetrad e_{(a)}^μ attached to the static observer is
 * diagonal in Schwarzschild coordinates:
 *
 * <pre>
 *   e_(t̂)^μ  = ( 1/√f ,    0  ,    0   ,    0           )
 *   e_(r̂)^μ  = (   0  ,  √f   ,    0   ,    0           )
 *   e_(θ̂)^μ  = (   0  ,    0  ,  1/r   ,    0           )
 *   e_(φ̂)^μ  = (   0  ,    0  ,    0   ,  1/(r sin θ)   )
 * </pre>
 *
 * <p>Orthonormality g_μν e_(a)^μ e_(b)^ν = η_(a)(b) (with g diagonal, only
 * (a,a) terms survive):
 * <ul>
 *   <li>g_tt · (1/√f)²        = −f · (1/f)          = −1 ✓</li>
 *   <li>g_rr · (√f)²          = (1/f) · f           = +1 ✓</li>
 *   <li>g_θθ · (1/r)²         = r² · (1/r²)         = +1 ✓</li>
 *   <li>g_φφ · 1/(r² sin²θ)   = r² sin²θ · (...)    = +1 ✓</li>
 * </ul>
 *
 * <p>The tetrad leg e_(r̂) points radially outward; −e_(r̂) points at the
 * black hole.
 *
 * <h2>Pixel → 4-momentum</h2>
 *
 * <p>Pixel (i, j) with resolution (W, H) maps to angular screen coordinates
 * (in radians):
 *
 * <pre>
 *   α = fov · ((i + 0.5)/W − 0.5)
 *   β = fov · (H/W) · (0.5 − (j + 0.5)/H)
 * </pre>
 *
 * <p>where {@code fov} is the full horizontal angular extent. Row j=0 is the
 * top of the image (FITS convention); α increases screen-right (along
 * e_(φ̂)); β increases screen-up (along e_(θ̂)); pixels are square.
 *
 * <p>The outgoing ray's spatial direction in the tetrad frame is:
 *
 * <pre>
 *   n̂^r̂  = − cos α · cos β      (mostly inward toward the BH)
 *   n̂^θ̂  =           sin β      (screen-up)
 *   n̂^φ̂  =   sin α · cos β      (screen-right)
 * </pre>
 *
 * <p>|n̂|² = cos²α cos²β + sin²β + sin²α cos²β = cos²β + sin²β = 1.
 *
 * <p>The locally measured photon energy is normalized to E_tetrad = 1. The
 * conserved Killing energy is then E_conserved = √f · E_tetrad = √f. The
 * tetrad-frame 4-momentum is k^(â) = (1, n̂^r̂, n̂^θ̂, n̂^φ̂) (future-directed,
 * null: η_(â)(b̂) k^(â) k^(b̂) = −1 + |n̂|² = 0).
 *
 * <p>Pushing to the coordinate basis, k^μ = e_(a)^μ k^(â):
 *
 * <pre>
 *   k^t = 1 / √f
 *   k^r = √f · (− cos α · cos β)
 *   k^θ = sin β / r_obs
 *   k^φ = sin α · cos β / (r_obs · sin θ_obs)
 * </pre>
 *
 * <p>This is null in the coordinate basis to round-off precision, because
 * the tetrad push-forward preserves the null condition exactly.
 *
 * <h2>Integration convention</h2>
 *
 * <p>k^μ as defined above is future-directed and points outward from the
 * observer into the scene. The caller integrates with positive affine-
 * parameter step (h &gt; 0). The traced geodesic is the time-reverse of the
 * physical photon's path; intersections with the horizon, disk, or
 * celestial sphere identify the source. This is the standard convention
 * in GRRT codes (Chan, Psaltis, Özel 2013).
 *
 * <h2>Thread-safety</h2>
 *
 * <p>All fields are final; {@link #initialState} writes only to a caller-owned
 * buffer. Concurrent calls from multiple worker threads are safe.
 */
public final class Camera {

    private final double rObs;
    private final double inclination;
    private final int width;
    private final int height;
    private final double fovRadians;

    // Pre-computed from observer position (does not depend on pixel).
    private final double sqrtF;             // √(1 − 2M/r_obs)
    private final double invR;              // 1 / r_obs
    private final double invRSinTheta;      // 1 / (r_obs · sin θ_obs)
    private final double pixelScale;        // fov / W  (square pixels: same for x and y)

    /**
     * @param metric       spacetime metric; must be the same instance used by the
     *                     integrator. Used only in the constructor for {@code mass()}
     *                     and {@code horizonRadius()} reads.
     * @param rObs         observer r-coordinate in geometrized units; must exceed
     *                     {@code metric.horizonRadius()}. For M87*-style imaging,
     *                     r_obs &gt;&gt; M (recommended r_obs ≥ 1000 M) so the static-
     *                     observer and asymptotic-flatness approximations hold.
     * @param inclination  θ_obs in radians, strictly in (0, π). Exposed even though
     *                     Schwarzschild is spherically symmetric so the interface
     *                     is stable for Kerr (Phase 2).
     * @param width        image width in pixels; must be positive.
     * @param height       image height in pixels; must be positive.
     * @param fovRadians   full horizontal angular extent of the image in radians;
     *                     must be in (0, π). At r_obs = 1000 M, fov ≈ 20 M / r_obs
     *                     = 0.02 rad gives ~2× margin around the ~10 M shadow
     *                     diameter (shadow diameter ≈ 6√3 M ≈ 10.4 M for
     *                     Schwarzschild).
     * @throws IllegalArgumentException if any argument is outside its valid domain
     */
    public Camera(Metric metric, double rObs, double inclination,
                  int width, int height, double fovRadians) {
        if (metric == null) {
            throw new IllegalArgumentException("metric must not be null");
        }
        if (!Double.isFinite(rObs) || rObs <= metric.horizonRadius()) {
            throw new IllegalArgumentException(
                    "rObs must be finite and > horizonRadius ("
                            + metric.horizonRadius() + "): " + rObs);
        }
        if (!(inclination > 0.0 && inclination < Math.PI)) {
            throw new IllegalArgumentException(
                    "inclination must be in (0, π): " + inclination);
        }
        if (width <= 0 || height <= 0) {
            throw new IllegalArgumentException(
                    "width and height must be positive: " + width + ", " + height);
        }
        if (!(fovRadians > 0.0 && fovRadians < Math.PI)) {
            throw new IllegalArgumentException(
                    "fovRadians must be in (0, π): " + fovRadians);
        }

        this.rObs = rObs;
        this.inclination = inclination;
        this.width = width;
        this.height = height;
        this.fovRadians = fovRadians;

        double f = 1.0 - 2.0 * metric.mass() / rObs;
        this.sqrtF = Math.sqrt(f);
        this.invR = 1.0 / rObs;
        this.invRSinTheta = 1.0 / (rObs * Math.sin(inclination));
        this.pixelScale = fovRadians / width;    // square pixels
    }

    /**
     * Write the initial 8-state [t, r, θ, φ, k^t, k^r, k^θ, k^φ] for the ray
     * passing through pixel (i, j) into the caller's buffer. No allocation.
     *
     * @param i   pixel column, 0 ≤ i &lt; width (0 = leftmost)
     * @param j   pixel row,    0 ≤ j &lt; height (0 = topmost, FITS convention)
     * @param out output buffer, length 8; overwritten
     * @throws IllegalArgumentException if indices out of range or out.length ≠ 8
     */
    public void initialState(int i, int j, double[] out) {
        if (i < 0 || i >= width || j < 0 || j >= height) {
            throw new IllegalArgumentException(
                    "pixel (" + i + ", " + j + ") out of range "
                            + "[0," + width + ") x [0," + height + ")");
        }
        double alpha = pixelScale * (i + 0.5 - 0.5 * width);
        double beta  = pixelScale * (0.5 * height - j - 0.5);
        writeState(alpha, beta, out);
    }

    /**
     * Test hook: write the initial state for a ray whose screen-plane
     * direction is specified by (α, β) in radians. Bypasses pixel centering
     * so tests can assert exact central-ray values at α = β = 0. Not part
     * of the public API.
     */
    void initialStateFromAngles(double alpha, double beta, double[] out) {
        writeState(alpha, beta, out);
    }

    private void writeState(double alpha, double beta, double[] out) {
        if (out.length != 8) {
            throw new IllegalArgumentException(
                    "out.length must be 8, got " + out.length);
        }
        double cosA = Math.cos(alpha);
        double sinA = Math.sin(alpha);
        double cosB = Math.cos(beta);
        double sinB = Math.sin(beta);

        out[0] = 0.0;                               // t
        out[1] = rObs;                              // r
        out[2] = inclination;                       // θ
        out[3] = 0.0;                               // φ

        out[4] = 1.0 / sqrtF;                       // k^t
        out[5] = sqrtF * (-cosA * cosB);            // k^r  (inward near center)
        out[6] = sinB * invR;                       // k^θ  (screen-up)
        out[7] = sinA * cosB * invRSinTheta;        // k^φ  (screen-right)
    }

    public int width()      { return width; }
    public int height()     { return height; }
    public double fovRadians()   { return fovRadians; }
    public double rObs()         { return rObs; }
    public double inclination()  { return inclination; }
}
