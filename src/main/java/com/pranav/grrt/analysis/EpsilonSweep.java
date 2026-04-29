package com.pranav.grrt.analysis;

import com.pranav.grrt.camera.Camera;
import com.pranav.grrt.disk.Disk;
import com.pranav.grrt.disk.DiskEmissionShader;
import com.pranav.grrt.disk.NovikovThorneDisk;
import com.pranav.grrt.integrator.DormandPrince45;
import com.pranav.grrt.io.FitsWriter;
import com.pranav.grrt.metric.JohannsenPsaltisMetric;
import com.pranav.grrt.metric.Metric;
import com.pranav.grrt.renderer.AdaptiveRayTracer;
import com.pranav.grrt.renderer.RayTracer;
import com.pranav.grrt.renderer.RenderConfig;
import com.pranav.grrt.renderer.Renderer;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardOpenOption;
import java.time.Instant;
import java.util.HashSet;
import java.util.List;
import java.util.Locale;
import java.util.Objects;
import java.util.Set;
import java.util.concurrent.TimeUnit;
import java.util.function.Supplier;

/**
 * Phase 3C orchestrator that turns the {@code (a, ε₃, i)} parameter
 * grid of {@code docs/phase-3c-plan.md} §3 into a populated
 * {@code output/sweep.csv} (schema in §4) plus one rendered FITS
 * frame per row. Resumable per §5.2: an interrupted run can be
 * restarted and will skip rows whose
 * {@code (epsilon_3, inclination_deg, resolution)} key already
 * appears in the CSV.
 *
 * <p><b>Per-frame work.</b> Build a fresh
 * {@link JohannsenPsaltisMetric} for the current {@code ε₃}, attach
 * a {@link NovikovThorneDisk} on the equator and a {@link Camera} at
 * the configured inclination, render via
 * {@link AdaptiveRayTracer}/{@link DormandPrince45} +
 * {@link DiskEmissionShader}, write the float image to FITS, then
 * run {@link RingExtractor} + {@link CircularityMetric} and append
 * the resulting CSV row with an {@code fsync} after each row so a
 * mid-sweep crash leaves the file consistent.
 *
 * <p><b>git_sha column.</b> Captured at the start of {@link #runSweep}
 * via {@code git rev-parse --short=7 HEAD}; if the working tree is
 * dirty the suffix {@code -dirty} is appended (e.g.
 * {@code 0bd68bc-dirty}). If git is unavailable the literal
 * {@code unknown} is recorded.
 *
 * <p><b>Concurrency.</b> Single-process: do not run two
 * {@code EpsilonSweep} instances against the same CSV concurrently.
 * Pixel-level parallelism within a frame is delegated to the
 * existing {@link Renderer} {@code parallelStream}.
 *
 * <p>See {@code docs/phase-3c-plan.md} §6.2 for the commit's place
 * in the 3C sub-phase decomposition and §7 / §9 for the validation
 * gates and failure-mode catalogue.
 */
public final class EpsilonSweep {

    /**
     * Frozen 12-column CSV header. Any deviation between the on-disk
     * header and this string causes {@link #runSweep} to fail fast on
     * resume (no schema migration in 3C). Mirrors
     * {@code docs/phase-3c-plan.md} §4.
     */
    public static final String CSV_HEADER =
            "epsilon_3,inclination_deg,resolution,mean_r,delta_r_rms,delta_r_p2p,"
            + "fourier_m1_amp,fourier_m1_phase,multipeak_bins,render_wallclock_s,"
            + "git_sha,timestamp_iso";

    private static final String UNKNOWN_GIT_SHA = "unknown";

    /**
     * Sweep configuration.
     *
     * @param spin               black-hole spin {@code a / M}, e.g. {@code 0.9}
     * @param epsilonGrid        ε₃ values to render (one frame per inclination per ε₃)
     * @param inclinationDegrees observer inclinations in degrees, in {@code (0, 180)}
     * @param resolution         square image side in pixels (production: {@code 512})
     * @param rOuter             disk outer edge in M (production: {@code 20.0})
     * @param csvPath            destination of the resumable CSV (e.g. {@code output/sweep.csv})
     * @param fitsDir            directory to receive per-frame FITS files
     * @param rObs               observer radius in M (production: {@code 1000.0})
     * @param fovRadians         full horizontal FOV in radians; linear extent at the BH is {@code rObs · fovRadians}
     * @param atol               DP45 absolute tolerance (production: {@code 1e-10})
     * @param rtol               DP45 relative tolerance (production: {@code 1e-10})
     * @param hInitial           initial DP45 step size in M
     * @param horizonCushion     {@code ε_h} in M (production: {@code 0.01})
     * @param maxSteps           per-ray integrator step cap
     * @param numBins            azimuthal bin count for {@link RingExtractor} (production: {@code 180})
     * @param pxThreshold        radial-pixel threshold for the multi-peak diagnostic (production: {@code 2.0})
     */
    public record Config(
            double spin,
            double[] epsilonGrid,
            double[] inclinationDegrees,
            int resolution,
            double rOuter,
            Path csvPath,
            Path fitsDir,
            double rObs,
            double fovRadians,
            double atol,
            double rtol,
            double hInitial,
            double horizonCushion,
            int maxSteps,
            int numBins,
            double pxThreshold
    ) {
        public Config {
            Objects.requireNonNull(epsilonGrid, "epsilonGrid");
            Objects.requireNonNull(inclinationDegrees, "inclinationDegrees");
            Objects.requireNonNull(csvPath, "csvPath");
            Objects.requireNonNull(fitsDir, "fitsDir");
            if (epsilonGrid.length == 0) {
                throw new IllegalArgumentException("epsilonGrid must be non-empty");
            }
            if (inclinationDegrees.length == 0) {
                throw new IllegalArgumentException("inclinationDegrees must be non-empty");
            }
            if (resolution < 2) {
                throw new IllegalArgumentException("resolution must be >= 2, got " + resolution);
            }
            if (!(rOuter > 0.0)) {
                throw new IllegalArgumentException("rOuter must be > 0, got " + rOuter);
            }
            if (!(rObs > rOuter)) {
                throw new IllegalArgumentException(
                        "rObs (" + rObs + ") must exceed rOuter (" + rOuter + ")");
            }
            if (!(fovRadians > 0.0 && fovRadians < Math.PI)) {
                throw new IllegalArgumentException("fovRadians must be in (0, π), got " + fovRadians);
            }
            if (!(atol > 0.0) || !(rtol > 0.0)) {
                throw new IllegalArgumentException("atol, rtol must be > 0");
            }
            if (!(hInitial > 0.0)) {
                throw new IllegalArgumentException("hInitial must be > 0");
            }
            if (!(horizonCushion > 0.0)) {
                throw new IllegalArgumentException("horizonCushion must be > 0");
            }
            if (maxSteps < 1) {
                throw new IllegalArgumentException("maxSteps must be >= 1");
            }
            if (numBins < 1) {
                throw new IllegalArgumentException("numBins must be >= 1");
            }
            if (!(pxThreshold > 0.0)) {
                throw new IllegalArgumentException("pxThreshold must be > 0");
            }
            for (double inc : inclinationDegrees) {
                if (!(inc > 0.0 && inc < 180.0)) {
                    throw new IllegalArgumentException(
                            "inclinationDegrees entry must be in (0, 180), got " + inc);
                }
            }
        }

        /** Image-plane resolution {@code rObs · fov / W} in M / px (linear). */
        public double pixelToM() {
            return rObs * fovRadians / resolution;
        }
    }

    /**
     * Summary returned by {@link #runSweep}.
     *
     * @param rendered count of frames freshly written to {@code csvPath} this run
     * @param skipped  count of frames whose key was already in {@code csvPath} on resume
     * @param failed   count of frames that raised an
     *                 {@link IllegalStateException} from the disk / metric
     *                 setup (recorded via the {@code [sweep-skip]} stdout
     *                 line; no CSV row written for these). Typical cause:
     *                 deep post-cusp JP at large positive ε₃ where the
     *                 prograde ISCO doesn't exist and
     *                 {@link NovikovThorneDisk}'s constructor cannot
     *                 bracket {@code dE/dr = 0}.
     */
    public record SweepResult(int rendered, int skipped, int failed) {
    }

    /** Per-frame analysis bundle (package-private; used by tests). */
    record FrameResult(CircularityMetric.Result stats, int multipeakBins, double renderWallclockS) {
    }

    public EpsilonSweep() {
    }

    /**
     * Drive the full {@code (epsilonGrid × inclinationDegrees)}
     * sweep. Skips frames whose key already appears in
     * {@code cfg.csvPath()}; emits one CSV row plus one FITS file
     * per rendered frame.
     *
     * @return summary of frames rendered vs skipped
     * @throws IOException on any file-system error or CSV header
     *                     mismatch on resume
     */
    public SweepResult runSweep(Config cfg) throws IOException {
        Objects.requireNonNull(cfg, "cfg");
        Path csvParent = cfg.csvPath().toAbsolutePath().getParent();
        if (csvParent != null) {
            Files.createDirectories(csvParent);
        }
        Files.createDirectories(cfg.fitsDir());

        Set<String> completedKeys = readCompletedKeys(cfg.csvPath());
        if (!Files.exists(cfg.csvPath())) {
            Files.writeString(cfg.csvPath(), CSV_HEADER + "\n", StandardCharsets.UTF_8,
                    StandardOpenOption.CREATE_NEW, StandardOpenOption.WRITE);
        }

        String gitSha = currentGitSha();
        int rendered = 0;
        int skipped = 0;
        int failed = 0;
        int totalFrames = cfg.epsilonGrid().length * cfg.inclinationDegrees().length;
        int frameIdx = 0;

        for (double eps : cfg.epsilonGrid()) {
            for (double inclDeg : cfg.inclinationDegrees()) {
                frameIdx++;
                String key = formatKey(eps, inclDeg, cfg.resolution());
                if (completedKeys.contains(key)) {
                    skipped++;
                    System.out.printf(Locale.ROOT,
                            "[sweep] frame %d/%d skip: eps3=%+.4f incl=%.1f res=%d (already in CSV)%n",
                            frameIdx, totalFrames, eps, inclDeg, cfg.resolution());
                    continue;
                }
                System.out.printf(Locale.ROOT,
                        "[sweep] frame %d/%d start: eps3=%+.4f incl=%.1f res=%d%n",
                        frameIdx, totalFrames, eps, inclDeg, cfg.resolution());
                Metric metric = new JohannsenPsaltisMetric(cfg.spin(), eps);
                Path fitsPath = fitsPathFor(cfg, eps, inclDeg);
                try {
                    FrameResult fr = renderAndAnalyze(metric, inclDeg, cfg, fitsPath);
                    String row = buildRow(eps, inclDeg, cfg.resolution(), fr, gitSha);
                    appendRow(cfg.csvPath(), row);
                    completedKeys.add(key);
                    rendered++;
                    System.out.printf(Locale.ROOT,
                            "[sweep] frame %d/%d done: wallclock=%.3fs mean_r=%.4f delta_r_rms/mean=%.4f multipeak=%d%n",
                            frameIdx, totalFrames, fr.renderWallclockS(),
                            fr.stats().meanR(),
                            fr.stats().deltaRrms() / fr.stats().meanR(),
                            fr.multipeakBins());
                } catch (IllegalStateException e) {
                    // Typical cause: post-cusp JP at large +ε₃ where the prograde
                    // ISCO doesn't exist and NovikovThorneDisk can't bracket
                    // dE/dr = 0. Skip the frame; do not write a CSV row.
                    failed++;
                    System.out.printf(Locale.ROOT,
                            "[sweep-skip] frame %d/%d: eps3=%+.4f incl=%.1f res=%d reason=%s%n",
                            frameIdx, totalFrames, eps, inclDeg, cfg.resolution(),
                            e.getMessage());
                }
            }
        }
        return new SweepResult(rendered, skipped, failed);
    }

    /**
     * Render one frame for the given metric / inclination, write the
     * FITS file, and run extraction + reduction. Package-private so
     * {@code EpsilonSweepIT} can drive a Kerr render through the
     * same pipeline for the consistency gate without going through
     * {@link #runSweep}'s CSV machinery.
     */
    FrameResult renderAndAnalyze(Metric metric, double inclinationDeg, Config cfg, Path fitsPath)
            throws IOException {
        double inclRad = Math.toRadians(inclinationDeg);
        // 1e-2 (vs the 3B.1 default of 1e-3) keeps the Page-Thorne radial-flux
        // integral comfortably away from the integrand's mild ISCO singularity,
        // where adaptive Simpson can otherwise return tiny noise-floor negatives
        // (~1e-4) for JP at non-zero ε₃. 0.01 M is 0.4–1.7 % of r_ISCO across
        // the sweep grid — negligible vs the ~10 % asymmetry signal.
        double iscoCushion = 1.0e-2;
        Disk disk = new NovikovThorneDisk(metric, cfg.spin(), cfg.rOuter(), iscoCushion);
        DiskEmissionShader baseShader = new DiskEmissionShader(disk, metric);

        // Defensive wrapper: NovikovThorneDisk uses Kerr-only Bardeen-1972
        // closed forms for E, L, Ω, dL/dr (a 3B implementation choice, see
        // NovikovThorneDisk.java §"Bardeen-Press-Teukolsky 1972 closed forms").
        // For JP at non-zero ε₃ the integrand is evaluated against Kerr
        // formulas integrated over [JP-r_ISCO + EPS, r], which can produce
        // noise-floor-negative Page-Thorne flux at a vanishing fraction of
        // pixels, throwing IllegalStateException out of disk.temperature(r).
        // Treat such pixels as zero emission (the failed integral contributes
        // ~1e-4 of the disk's nominal flux range; affects << 1 % of disk-hit
        // pixels in practice). Per docs/phase-3c-plan.md §13, a JP-aware NT
        // disk is deferred.
        com.pranav.grrt.renderer.RayShader shader = (outcome, state) -> {
            try {
                return baseShader.shade(outcome, state);
            } catch (IllegalStateException e) {
                return 0.0f;
            }
        };

        RenderConfig renderConfig = new RenderConfig(
                cfg.horizonCushion(),
                2.0 * cfg.rObs(),
                0.01, 0.5, 10.0,                                            // hNear/hFar/rTransition (FixedStep-only; ignored by AdaptiveRayTracer)
                cfg.maxSteps(),
                shader);

        Camera camera = new Camera(metric, cfg.rObs(), inclRad,
                cfg.resolution(), cfg.resolution(), cfg.fovRadians());

        Supplier<RayTracer> tracerFactory = () -> new AdaptiveRayTracer(
                new DormandPrince45(), metric, renderConfig,
                cfg.atol(), cfg.rtol(), cfg.hInitial(), disk);

        Renderer renderer = Renderer.withRayTracer(camera, metric, tracerFactory, renderConfig);

        long renderStart = System.nanoTime();
        float[][] image = renderer.render();
        long renderEnd = System.nanoTime();
        double renderWallclockS = (renderEnd - renderStart) / 1.0e9;

        FitsWriter.Metadata meta = new FitsWriter.Metadata(
                "sweep-" + metric.getClass().getSimpleName(),
                metric.mass(),
                cfg.rObs(),
                inclRad,
                cfg.fovRadians(),
                cfg.horizonCushion(),
                "DP45 atol=" + cfg.atol() + " rtol=" + cfg.rtol(),
                cfg.maxSteps(),
                FitsWriter.gitHash());
        FitsWriter.write(fitsPath, image, meta);

        double cxPx = 0.5 * cfg.resolution();
        double cyPx = 0.5 * cfg.resolution();
        RingExtractor extractor = new RingExtractor(cfg.numBins(), cfg.pixelToM(), cxPx, cyPx);
        double[] radii = extractor.extract(image, cfg.resolution(), cfg.resolution());
        CircularityMetric.Result stats = CircularityMetric.compute(radii);
        int multipeakBins = extractor.multipeakBinCount(image,
                cfg.resolution(), cfg.resolution(), cfg.pxThreshold());

        return new FrameResult(stats, multipeakBins, renderWallclockS);
    }

    /**
     * Production sweep entry point. Writes
     * {@code output/sweep.csv} + 32 FITS frames at 512² for the
     * {@code a = 0.9}, 16-point ε₃ × {17°, 60°} grid of
     * {@code docs/phase-3c-plan.md} §3.
     */
    public static void main(String[] args) throws IOException {
        double spin = 0.9;
        // Grid reordered for early-feedback execution: ε₃ = 0 first (closest
        // to Kerr, fastest expected), then small-|ε₃| values, then the cusp
        // neighborhood, then mid magnitudes, then the extremes (-2.5, +1.00)
        // last so they don't block all preceding frames if integration in
        // the deep-deformation regime is slow. Resumability preserves
        // progress if a late frame stalls. Sixteen values total — same set
        // as docs/jp-parameter-space-notes.md §5.2.
        double[] epsGrid = {
                 0.00, -0.05, +0.05, +0.08, -0.10, +0.10, +0.11, +0.12,
                +0.13, -0.20, +0.20, -0.50, -1.00, -1.50, -2.50, +1.00
        };
        double[] inclDeg = { 17.0, 60.0 };
        int resolution = 512;
        double rOuter = 20.0;
        Path csvPath = Path.of("output", "sweep.csv");
        Path fitsDir = Path.of("output");

        Config cfg = new Config(
                spin, epsGrid, inclDeg, resolution, rOuter, csvPath, fitsDir,
                1000.0,                                                     // rObs
                30.0 / 1000.0,                                              // ±15 M field of view
                1.0e-10, 1.0e-10,
                1.0,
                0.01,
                // 200k chosen empirically: at the production 1M cap, a few
                // pathological pixels at large |ε₃| (e.g. -0.20 onward)
                // consumed ~12 min CPU each before classifying MAX_STEPS,
                // dominating wall-clock. 200k bounds the per-pixel CPU at
                // ~25 s and lets the full sweep finish in a few hours rather
                // than overnight. Pixels that would have taken 200k–1M steps
                // render as MAX_STEPS (black) rather than capturing ground
                // truth — acceptable for RNAAS-scope shadow-edge measurement;
                // revisit if the prograde bin radii at extreme ε₃ disagree
                // with visual inspection.
                200_000,
                180,
                2.0
        );

        long t0 = System.nanoTime();
        SweepResult result = new EpsilonSweep().runSweep(cfg);
        double wallS = (System.nanoTime() - t0) / 1.0e9;
        System.out.printf(Locale.ROOT,
                "Sweep complete: %d rendered, %d skipped, %d failed (%.1f s wall-clock)%n",
                result.rendered(), result.skipped(), result.failed(), wallS);
    }

    // ------------------------------------------------------------------
    // CSV r/w + key handling
    // ------------------------------------------------------------------

    static String formatKey(double eps, double inclDeg, int resolution) {
        return String.format(Locale.ROOT, "%+.4f|%.1f|%d", eps, inclDeg, resolution);
    }

    private Path fitsPathFor(Config cfg, double eps, double inclDeg) {
        String name = String.format(Locale.ROOT,
                "sweep_a%.1f_eps%+.4f_i%.1f_res%d.fits",
                cfg.spin(), eps, inclDeg, cfg.resolution());
        return cfg.fitsDir().resolve(name);
    }

    private static String buildRow(double eps, double inclDeg, int resolution,
                                   FrameResult fr, String gitSha) {
        CircularityMetric.Result s = fr.stats();
        return String.format(Locale.ROOT,
                "%+.4f,%.1f,%d,%.6e,%.6e,%.6e,%.6e,%.6e,%d,%.3f,%s,%s",
                eps,
                inclDeg,
                resolution,
                s.meanR(),
                s.deltaRrms(),
                s.deltaRp2p(),
                s.fourierM1Amplitude(),
                s.fourierM1Phase(),
                fr.multipeakBins(),
                fr.renderWallclockS(),
                gitSha,
                Instant.now().toString());
    }

    /**
     * Read the CSV at {@code csvPath} (if it exists) and return the
     * set of completed-frame keys. Validates that the on-disk header
     * exactly matches {@link #CSV_HEADER}; throws {@link IOException}
     * on mismatch.
     */
    static Set<String> readCompletedKeys(Path csvPath) throws IOException {
        Set<String> keys = new HashSet<>();
        if (!Files.exists(csvPath)) {
            return keys;
        }
        List<String> lines = Files.readAllLines(csvPath, StandardCharsets.UTF_8);
        if (lines.isEmpty()) {
            return keys;
        }
        if (!lines.get(0).equals(CSV_HEADER)) {
            throw new IOException(
                    "CSV header mismatch in " + csvPath
                            + ": expected '" + CSV_HEADER + "', got '" + lines.get(0) + "'");
        }
        for (int i = 1; i < lines.size(); i++) {
            String line = lines.get(i).trim();
            if (line.isEmpty()) {
                continue;
            }
            String[] parts = line.split(",", -1);
            if (parts.length < 3) {
                continue;
            }
            try {
                double eps = Double.parseDouble(parts[0]);
                double incl = Double.parseDouble(parts[1]);
                int res = Integer.parseInt(parts[2]);
                keys.add(formatKey(eps, incl, res));
            } catch (NumberFormatException ignore) {
                // Malformed row — skip; the next sweep run will overwrite.
            }
        }
        return keys;
    }

    /**
     * Append one row to {@code csvPath} with a sync after the write
     * so a crash between rows does not leave the file in a torn
     * state. The CSV row must NOT include the trailing newline; this
     * method adds one.
     */
    private static void appendRow(Path csvPath, String row) throws IOException {
        Files.writeString(csvPath, row + "\n", StandardCharsets.UTF_8,
                StandardOpenOption.APPEND, StandardOpenOption.WRITE,
                StandardOpenOption.SYNC);
    }

    // ------------------------------------------------------------------
    // git SHA capture (7-char + -dirty)
    // ------------------------------------------------------------------

    /** {@code git rev-parse --short=7 HEAD} + {@code -dirty} if working tree differs from HEAD. */
    static String currentGitSha() {
        String shortSha = runGitOutput("rev-parse", "--short=7", "HEAD");
        if (shortSha.isEmpty()) {
            return UNKNOWN_GIT_SHA;
        }
        boolean dirty = !runGitExitZero("diff-index", "--quiet", "HEAD", "--");
        return dirty ? shortSha + "-dirty" : shortSha;
    }

    private static String runGitOutput(String... args) {
        String[] cmd = new String[args.length + 1];
        cmd[0] = "git";
        System.arraycopy(args, 0, cmd, 1, args.length);
        try {
            Process p = new ProcessBuilder(cmd).redirectErrorStream(true).start();
            String line;
            try (BufferedReader r = new BufferedReader(
                    new InputStreamReader(p.getInputStream(), StandardCharsets.UTF_8))) {
                line = r.readLine();
            }
            if (!p.waitFor(2, TimeUnit.SECONDS)) {
                p.destroyForcibly();
                return "";
            }
            if (p.exitValue() != 0 || line == null) {
                return "";
            }
            return line.trim();
        } catch (IOException | InterruptedException e) {
            if (e instanceof InterruptedException) {
                Thread.currentThread().interrupt();
            }
            return "";
        }
    }

    private static boolean runGitExitZero(String... args) {
        String[] cmd = new String[args.length + 1];
        cmd[0] = "git";
        System.arraycopy(args, 0, cmd, 1, args.length);
        try {
            Process p = new ProcessBuilder(cmd).redirectErrorStream(true).start();
            if (!p.waitFor(2, TimeUnit.SECONDS)) {
                p.destroyForcibly();
                return false;
            }
            return p.exitValue() == 0;
        } catch (IOException | InterruptedException e) {
            if (e instanceof InterruptedException) {
                Thread.currentThread().interrupt();
            }
            return false;
        }
    }
}
