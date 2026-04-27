package com.pranav.grrt.metric;

/**
 * Kerr metric for a rotating, uncharged black hole of mass M and spin
 * parameter a in Boyer-Lindquist coordinates (t, r, θ, φ).
 *
 * <p>Line element (with Σ = r² + a² cos²θ, Δ = r² − 2Mr + a²):
 * <pre>
 *   ds² = -(1 - 2Mr/Σ) dt²
 *         - (4Mar sin²θ / Σ) dt dφ
 *         + (Σ/Δ) dr²
 *         + Σ dθ²
 *         + ((r² + a²) + 2Ma²r sin²θ/Σ) sin²θ dφ²
 * </pre>
 *
 * <p>Signature (-, +, +, +), geometrized units G = c = 1. Spin is
 * restricted to subextremal |a| &lt; M; extremal and superextremal cases
 * are rejected so the outer horizon r₊ = M + √(M² - a²) stays real and
 * separated from the inner horizon.
 *
 * <p>Landmarks:
 * <ul>
 *   <li>Outer event horizon: r₊ = M + √(M² − a²)</li>
 *   <li>Equatorial ergosphere: r_ergo(θ=π/2) = 2M</li>
 *   <li>Prograde/retrograde ISCO: Bardeen, Press, Teukolsky (1972)</li>
 * </ul>
 *
 * <p>Valid for r &gt; r₊ and θ ∈ (0, π). Boyer-Lindquist has coordinate
 * singularities at r = r₊ and on the polar axis; integration must stop
 * before reaching either.
 *
 * <h2>Inverse metric</h2>
 *
 * The (r, θ) block is diagonal (g^{rr} = Δ/Σ, g^{θθ} = 1/Σ). The
 * (t, φ) 2×2 block is inverted in closed form using the identity
 * det((t,φ) block) = -Δ sin²θ (derivation: g_tt g_φφ − g_tφ²
 * = [A(Σ − 2Mr) + 4M²a²r² sin²θ] · (−sin²θ/Σ²), and the bracket simplifies
 * to Σ²Δ via A − 2Mr(r²+a²) = ΣΔ):
 * <pre>
 *   g^{tt} = -((r² + a²)² - a²Δ sin²θ) / (Σ Δ)     [= -A/(ΣΔ)]
 *   g^{φφ} = (Δ - a² sin²θ) / (Σ Δ sin²θ)
 *   g^{tφ} = -2Mar / (Σ Δ)
 * </pre>
 * No 4×4 generic inversion is performed.
 *
 * <h2>Christoffel symbols</h2>
 *
 * {@link #christoffel} uses the textbook formula
 *   Γ^α_{μν} = ½ g^{ασ} (∂_μ g_{σν} + ∂_ν g_{σμ} - ∂_σ g_{μν})
 * with ∂_t = ∂_φ = 0 (Kerr is stationary and axisymmetric). The required
 * partial derivatives are all computed in closed form (derived by hand
 * inside the method). This avoids tabulating every Γ component explicitly
 * and keeps the sign conventions in one place; it costs one 4×4×4
 * allocation per call and is therefore not the hot path.
 *
 * <p>The optimized {@code geodesicAcceleration} override (Phase 2 A2)
 * fuses the non-zero contractions directly into the four output
 * components and is regression-tested against the default contraction
 * that goes through this method.
 *
 * <p>References:
 * <ul>
 *   <li>Chandrasekhar (1983), <i>The Mathematical Theory of Black Holes</i>, Ch. 6 §58.</li>
 *   <li>Chan, Psaltis, Özel (2013), <i>GRay</i>, ApJ 777 13, Appendix A.</li>
 *   <li>Bardeen, Press, Teukolsky (1972), ApJ 178 347 (ISCO).</li>
 * </ul>
 */
public final class KerrMetric implements Metric {

    private final double M;
    private final double a;

    /**
     * @param mass gravitational mass M in geometrized units (G = c = 1);
     *             must be positive and finite
     * @param spin Kerr spin parameter a in geometrized units; must be
     *             finite with |a| &lt; M (subextremal)
     */
    public KerrMetric(double mass, double spin) {
        if (mass <= 0.0 || !Double.isFinite(mass)) {
            throw new IllegalArgumentException("Mass must be positive and finite: " + mass);
        }
        if (!Double.isFinite(spin)) {
            throw new IllegalArgumentException("Spin must be finite: " + spin);
        }
        if (Math.abs(spin) >= mass) {
            throw new IllegalArgumentException(
                    "Spin must be subextremal |a| < M: a = " + spin + ", M = " + mass);
        }
        this.M = mass;
        this.a = spin;
    }

    /** Unit-mass Kerr hole (M = 1) with the given spin. */
    public KerrMetric(double spin) {
        this(1.0, spin);
    }

    @Override public double mass() { return M; }

    /** @return Kerr spin parameter a in geometrized units */
    public double spin() { return a; }

    /** @return outer event horizon r₊ = M + √(M² − a²) in geometrized units */
    @Override
    public double horizonRadius() {
        return M + Math.sqrt(M * M - a * a);
    }

    /** @return equatorial ergosphere radius 2M (θ = π/2); off-equator the
     *  static limit is r = M + √(M² − a² cos²θ), which is not reported
     *  here to keep the interface minimal. */
    public double ergosphereRadiusEquatorial() {
        return 2.0 * M;
    }

    /**
     * Innermost stable circular orbit radius for timelike equatorial
     * geodesics, Bardeen, Press, Teukolsky (1972), ApJ 178 347, eq. 2.21.
     *
     * <p>With χ = a/M:
     * <pre>
     *   Z1 = 1 + (1 - χ²)^(1/3) [(1 + χ)^(1/3) + (1 - χ)^(1/3)]
     *   Z2 = √(3 χ² + Z1²)
     *   r_isco / M = 3 + Z2 ∓ √((3 - Z1)(3 + Z1 + 2 Z2))
     * </pre>
     * with the minus sign for prograde (co-rotating) and plus for
     * retrograde (counter-rotating) orbits.
     *
     * <p>Limits: at a = 0, r_isco = 6M either way. At a → M (prograde),
     * r_isco → M (extremal co-rotating marginal orbit).
     *
     * @param prograde true for co-rotating ISCO, false for counter-rotating
     * @return ISCO radius in geometrized units
     */
    @Override
    public double iscoRadius(boolean prograde) {
        double chi = a / M;
        double chi2 = chi * chi;
        double Z1 = 1.0
                + Math.cbrt(1.0 - chi2)
                * (Math.cbrt(1.0 + chi) + Math.cbrt(1.0 - chi));
        double Z2 = Math.sqrt(3.0 * chi2 + Z1 * Z1);
        double root = Math.sqrt((3.0 - Z1) * (3.0 + Z1 + 2.0 * Z2));
        double rOverM = prograde ? (3.0 + Z2 - root) : (3.0 + Z2 + root);
        return rOverM * M;
    }

    // ------------------------------------------------------------------
    // Metric tensor and inverse
    // ------------------------------------------------------------------

    @Override
    public double[][] g(double[] x) {
        double r = x[1];
        double th = x[2];
        double s = Math.sin(th);
        double c = Math.cos(th);
        double s2 = s * s;
        double c2 = c * c;

        double r2 = r * r;
        double a2 = a * a;
        double Sigma = r2 + a2 * c2;
        double Delta = r2 - 2.0 * M * r + a2;
        double R2 = r2 + a2;
        double A = R2 * R2 - a2 * Delta * s2;

        double[][] gmn = new double[4][4];
        gmn[0][0] = -(1.0 - 2.0 * M * r / Sigma);
        gmn[1][1] = Sigma / Delta;
        gmn[2][2] = Sigma;
        gmn[3][3] = A * s2 / Sigma;
        double gtphi = -2.0 * M * a * r * s2 / Sigma;
        gmn[0][3] = gtphi;
        gmn[3][0] = gtphi;
        return gmn;
    }

    @Override
    public double[][] gInv(double[] x) {
        double r = x[1];
        double th = x[2];
        double s = Math.sin(th);
        double c = Math.cos(th);
        double s2 = s * s;
        double c2 = c * c;

        double r2 = r * r;
        double a2 = a * a;
        double Sigma = r2 + a2 * c2;
        double Delta = r2 - 2.0 * M * r + a2;
        double R2 = r2 + a2;
        double A = R2 * R2 - a2 * Delta * s2;
        double SD = Sigma * Delta;

        double[][] gi = new double[4][4];
        gi[0][0] = -A / SD;
        gi[1][1] = Delta / Sigma;
        gi[2][2] = 1.0 / Sigma;
        gi[3][3] = (Delta - a2 * s2) / (SD * s2);
        double gitp = -2.0 * M * a * r / SD;
        gi[0][3] = gitp;
        gi[3][0] = gitp;
        return gi;
    }

    // ------------------------------------------------------------------
    // Christoffel symbols
    // ------------------------------------------------------------------

    @Override
    public double[][][] christoffel(double[] x) {
        double r = x[1];
        double th = x[2];
        double s = Math.sin(th);
        double c = Math.cos(th);
        double s2 = s * s;
        double c2 = c * c;
        double sc = s * c;

        double r2 = r * r;
        double a2 = a * a;
        double Sigma = r2 + a2 * c2;
        double Delta = r2 - 2.0 * M * r + a2;
        double R2 = r2 + a2;
        double A = R2 * R2 - a2 * Delta * s2;

        // Inverse metric (closed-form, as in gInv).
        double SD = Sigma * Delta;
        double gTT = -A / SD;
        double gRR = Delta / Sigma;
        double gHH = 1.0 / Sigma;
        double gFF = (Delta - a2 * s2) / (SD * s2);
        double gTF = -2.0 * M * a * r / SD;

        // Partial derivatives of g_{μν}. Only ∂_r and ∂_θ are non-zero
        // (Kerr is stationary and axisymmetric). Derivations:
        //   g_tt = -1 + 2Mr/Σ
        //     ∂_r g_tt = 2M(Σ - r Σ,r)/Σ² = -2M(r² - a²c²)/Σ²
        //     ∂_θ g_tt = -2Mr Σ,θ/Σ² = 4Mra²sc/Σ²
        //   g_rr = Σ/Δ
        //     ∂_r g_rr = 2[rΔ - (r-M)Σ]/Δ² = 2[a²rs² - M(r²-a²c²)]/Δ²
        //     ∂_θ g_rr = Σ,θ/Δ = -2a²sc/Δ
        //   g_θθ = Σ
        //     ∂_r g_θθ = 2r
        //     ∂_θ g_θθ = -2a²sc
        //   g_tφ = -2Mars²/Σ
        //     ∂_r g_tφ = -2Mas²(Σ - r Σ,r)/Σ² = 2Mas²(r² - a²c²)/Σ²
        //     ∂_θ g_tφ = -2Mar(2scΣ - s² Σ,θ)/Σ² = -4Marsc(r²+a²)/Σ²
        //   g_φφ = A s²/Σ   (A,r = 4r(r²+a²) - 2a²(r-M)s²; A,θ = -2a²Δsc)
        //     ∂_r g_φφ = s²(A,r Σ - 2rA)/Σ²
        //     ∂_θ g_φφ = 2sc[A(r²+a²) - a²Δs²Σ]/Σ²
        double Sig2 = Sigma * Sigma;
        double Del2 = Delta * Delta;
        double rrMinusA2c2 = r2 - a2 * c2;
        double Ar = 4.0 * r * R2 - 2.0 * a2 * (r - M) * s2;

        double drgTT = -2.0 * M * rrMinusA2c2 / Sig2;
        double dhgTT = 4.0 * M * r * a2 * sc / Sig2;
        double drgRR = 2.0 * (a2 * r * s2 - M * rrMinusA2c2) / Del2;
        double dhgRR = -2.0 * a2 * sc / Delta;
        double drgHH = 2.0 * r;
        double dhgHH = -2.0 * a2 * sc;
        double drgTF = 2.0 * M * a * s2 * rrMinusA2c2 / Sig2;
        double dhgTF = -4.0 * M * a * r * sc * R2 / Sig2;
        double drgFF = s2 * (Ar * Sigma - 2.0 * r * A) / Sig2;
        double dhgFF = 2.0 * sc * (A * R2 - a2 * Delta * s2 * Sigma) / Sig2;

        // Assemble Γ^α_{μν} from the formula. Only non-zero gInv
        // components contribute: g^{tt}, g^{tφ}, g^{rr}, g^{θθ}, g^{φφ}.
        double[][][] G = new double[4][4][4];

        // ---- α = r : Γ^r_{μν} = ½ g^{rr} (∂_μ g_{rν} + ∂_ν g_{rμ} - ∂_r g_{μν}) ----
        // g_{rν} ≠ 0 only for ν = r.
        G[1][1][1] =  0.5 * gRR * drgRR;           // Γ^r_rr
        double Grth = 0.5 * gRR * dhgRR;           // Γ^r_rθ = Γ^r_θr
        G[1][1][2] = Grth; G[1][2][1] = Grth;
        G[1][0][0] = -0.5 * gRR * drgTT;           // Γ^r_tt
        G[1][2][2] = -0.5 * gRR * drgHH;           // Γ^r_θθ  (= -rΔ/Σ)
        G[1][3][3] = -0.5 * gRR * drgFF;           // Γ^r_φφ
        double Grtf = -0.5 * gRR * drgTF;          // Γ^r_tφ = Γ^r_φt
        G[1][0][3] = Grtf; G[1][3][0] = Grtf;

        // ---- α = θ : Γ^θ_{μν} = ½ g^{θθ} (∂_μ g_{θν} + ∂_ν g_{θμ} - ∂_θ g_{μν}) ----
        // g_{θν} ≠ 0 only for ν = θ.
        double Ghrh = 0.5 * gHH * drgHH;           // Γ^θ_rθ = Γ^θ_θr = r/Σ
        G[2][1][2] = Ghrh; G[2][2][1] = Ghrh;
        G[2][2][2] =  0.5 * gHH * dhgHH;           // Γ^θ_θθ = -a²sc/Σ
        G[2][1][1] = -0.5 * gHH * dhgRR;           // Γ^θ_rr = a²sc/(ΣΔ)
        G[2][0][0] = -0.5 * gHH * dhgTT;           // Γ^θ_tt
        G[2][3][3] = -0.5 * gHH * dhgFF;           // Γ^θ_φφ
        double Ghtf = -0.5 * gHH * dhgTF;          // Γ^θ_tφ = Γ^θ_φt
        G[2][0][3] = Ghtf; G[2][3][0] = Ghtf;

        // ---- α = t : Γ^t_{μν} = ½[g^{tt}(∂_μ g_{tν}+∂_ν g_{tμ})
        //                        + g^{tφ}(∂_μ g_{φν}+∂_ν g_{φμ})]  ----
        //   (∂_t = ∂_φ = 0, so only (μ,ν) pairs with at least one of {r, θ}
        //    contribute; g_{tν} ≠ 0 only for ν ∈ {t, φ}; likewise g_{φν}.)
        double Gttr = 0.5 * (gTT * drgTT + gTF * drgTF);     // Γ^t_tr
        G[0][0][1] = Gttr; G[0][1][0] = Gttr;
        double Gtth = 0.5 * (gTT * dhgTT + gTF * dhgTF);     // Γ^t_tθ
        G[0][0][2] = Gtth; G[0][2][0] = Gtth;
        double Gtrf = 0.5 * (gTT * drgTF + gTF * drgFF);     // Γ^t_rφ
        G[0][1][3] = Gtrf; G[0][3][1] = Gtrf;
        double Gthf = 0.5 * (gTT * dhgTF + gTF * dhgFF);     // Γ^t_θφ
        G[0][2][3] = Gthf; G[0][3][2] = Gthf;

        // ---- α = φ : Γ^φ_{μν} = ½[g^{tφ}(∂_μ g_{tν}+∂_ν g_{tμ})
        //                        + g^{φφ}(∂_μ g_{φν}+∂_ν g_{φμ})]  ----
        double Gftr = 0.5 * (gTF * drgTT + gFF * drgTF);     // Γ^φ_tr
        G[3][0][1] = Gftr; G[3][1][0] = Gftr;
        double Gfth = 0.5 * (gTF * dhgTT + gFF * dhgTF);     // Γ^φ_tθ
        G[3][0][2] = Gfth; G[3][2][0] = Gfth;
        double Gfrf = 0.5 * (gTF * drgTF + gFF * drgFF);     // Γ^φ_rφ
        G[3][1][3] = Gfrf; G[3][3][1] = Gfrf;
        double Gfhf = 0.5 * (gTF * dhgTF + gFF * dhgFF);     // Γ^φ_θφ
        G[3][2][3] = Gfhf; G[3][3][2] = Gfhf;

        return G;
    }

    // ------------------------------------------------------------------
    // Optimized geodesic acceleration
    // ------------------------------------------------------------------

    /**
     * Non-allocating geodesic acceleration exploiting Kerr sparsity.
     * Follows the same closed-form partial-derivative derivation as
     * {@link #christoffel}, but skips the 4×4×4 tensor build and fuses
     * the non-zero contractions directly into {@code out[0..3]}.
     *
     * <p>Structure:
     * <ul>
     *   <li>{@code a^t} and {@code a^φ} are each a sum of four bilinear
     *       cross-terms (no diagonal Christoffels for those indices).</li>
     *   <li>{@code a^r} and {@code a^θ} each have four diagonal
     *       Christoffels plus (r,θ) and (t,φ) cross terms.</li>
     * </ul>
     *
     * <p>Complexity: O(1), ~60 FLOPs, no heap allocation, independent of
     * the state components t and φ (Kerr is stationary and axisymmetric).
     */
    @Override
    public void geodesicAcceleration(double[] x, double[] k, double[] out) {
        double r = x[1];
        double th = x[2];
        double s = Math.sin(th);
        double c = Math.cos(th);
        double s2 = s * s;
        double c2 = c * c;
        double sc = s * c;

        double r2 = r * r;
        double a2 = a * a;
        double Sigma = r2 + a2 * c2;
        double Delta = r2 - 2.0 * M * r + a2;
        double R2 = r2 + a2;
        double A = R2 * R2 - a2 * Delta * s2;

        // Inverse metric components (closed-form; see gInv).
        double SD = Sigma * Delta;
        double gTT = -A / SD;
        double gRR = Delta / Sigma;
        double gHH = 1.0 / Sigma;
        double gFF = (Delta - a2 * s2) / (SD * s2);
        double gTF = -2.0 * M * a * r / SD;

        // ∂g partials (derivation in christoffel()).
        double Sig2 = Sigma * Sigma;
        double Del2 = Delta * Delta;
        double rrMinusA2c2 = r2 - a2 * c2;
        double Ar = 4.0 * r * R2 - 2.0 * a2 * (r - M) * s2;

        double drgTT = -2.0 * M * rrMinusA2c2 / Sig2;
        double dhgTT = 4.0 * M * r * a2 * sc / Sig2;
        double drgRR = 2.0 * (a2 * r * s2 - M * rrMinusA2c2) / Del2;
        double dhgRR = -2.0 * a2 * sc / Delta;
        double drgHH = 2.0 * r;
        double dhgHH = -2.0 * a2 * sc;
        double drgTF = 2.0 * M * a * s2 * rrMinusA2c2 / Sig2;
        double dhgTF = -4.0 * M * a * r * sc * R2 / Sig2;
        double drgFF = s2 * (Ar * Sigma - 2.0 * r * A) / Sig2;
        double dhgFF = 2.0 * sc * (A * R2 - a2 * Delta * s2 * Sigma) / Sig2;

        // 20 non-zero Christoffels (each unique; symmetric partners reused
        // implicitly through the 2× factors on cross terms below).
        double Grrr  =  0.5 * gRR * drgRR;
        double Grrth =  0.5 * gRR * dhgRR;
        double Grtt  = -0.5 * gRR * drgTT;
        double Grhh  = -0.5 * gRR * drgHH;
        double Grff  = -0.5 * gRR * drgFF;
        double Grtf  = -0.5 * gRR * drgTF;

        double Ghrh  =  0.5 * gHH * drgHH;
        double Ghhh  =  0.5 * gHH * dhgHH;
        double Ghrr  = -0.5 * gHH * dhgRR;
        double Ghtt  = -0.5 * gHH * dhgTT;
        double Ghff  = -0.5 * gHH * dhgFF;
        double Ghtf  = -0.5 * gHH * dhgTF;

        double Gttr  =  0.5 * (gTT * drgTT + gTF * drgTF);
        double Gtth  =  0.5 * (gTT * dhgTT + gTF * dhgTF);
        double Gtrf  =  0.5 * (gTT * drgTF + gTF * drgFF);
        double Gthf  =  0.5 * (gTT * dhgTF + gTF * dhgFF);

        double Gftr  =  0.5 * (gTF * drgTT + gFF * drgTF);
        double Gfth  =  0.5 * (gTF * dhgTT + gFF * dhgTF);
        double Gfrf  =  0.5 * (gTF * drgTF + gFF * drgFF);
        double Gfhf  =  0.5 * (gTF * dhgTF + gFF * dhgFF);

        // Contract: a^μ = -Γ^μ_{αβ} k^α k^β. Cross terms carry factor 2
        // from the (α,β)+(β,α) sum.
        double kt = k[0];
        double kr = k[1];
        double kh = k[2];
        double kf = k[3];

        out[0] = -2.0 * (Gttr * kt * kr + Gtth * kt * kh
                       + Gtrf * kr * kf + Gthf * kh * kf);

        out[1] = -(Grtt * kt * kt + Grrr * kr * kr
                  + Grhh * kh * kh + Grff * kf * kf
                  + 2.0 * Grrth * kr * kh
                  + 2.0 * Grtf  * kt * kf);

        out[2] = -(Ghtt * kt * kt + Ghrr * kr * kr
                  + Ghhh * kh * kh + Ghff * kf * kf
                  + 2.0 * Ghrh * kr * kh
                  + 2.0 * Ghtf * kt * kf);

        out[3] = -2.0 * (Gftr * kt * kr + Gfth * kt * kh
                       + Gfrf * kr * kf + Gfhf * kh * kf);
    }
}
