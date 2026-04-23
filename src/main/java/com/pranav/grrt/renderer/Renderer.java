package com.pranav.grrt.renderer;

import com.pranav.grrt.camera.Camera;
import com.pranav.grrt.integrator.Integrator;
import com.pranav.grrt.io.FitsWriter;
import com.pranav.grrt.metric.Metric;

import java.io.IOException;
import java.nio.file.Path;
import java.util.function.Supplier;
import java.util.stream.IntStream;

/**
 * Pixel loop + ray tracer. For each pixel the {@link Camera} provides an
 * initial 8-state, a {@link RayTracer} advances it to termination, and
 * the {@link RayShader} assigns a float value based on the terminal
 * outcome.
 *
 * <h2>Strategy split (Phase 2 C1)</h2>
 *
 * Ray tracing is delegated to a {@link RayTracer}:
 * <ul>
 *   <li>{@link FixedStepRayTracer} — Phase 1 two-zone fixed-step RK4,
 *       used by the backward-compatible constructor that takes a
 *       {@code Supplier<Integrator>}.</li>
 *   <li>{@link AdaptiveRayTracer} — Phase 2 DP45 + polar-axis reflection,
 *       selected via {@link #withRayTracer}.</li>
 * </ul>
 * The strategy is fixed at construction. A {@link Supplier} of
 * {@link RayTracer} is used to hand each worker thread its own instance
 * (both tracer flavours are stateful and not thread-safe).
 *
 * <h2>Termination</h2>
 *
 * Delegated to the {@link RayTracer} implementation, which reads
 * {@link RenderConfig#horizonCushion()}, {@link RenderConfig#escapeRadius()},
 * and {@link RenderConfig#maxSteps()}.
 *
 * <h2>Parallelism</h2>
 *
 * Pixels are processed in parallel via {@code IntStream.range(0, W*H).parallel()}.
 * Each worker uses its own {@link ThreadLocal} {@link RayTracer} plus a
 * thread-local state buffer. {@link Camera} is shared (thread-safe).
 * Because each output pixel is written by exactly one worker and the
 * shader is pure, the result is deterministic regardless of thread
 * scheduling.
 */
public final class Renderer {

    private final Camera camera;
    private final Metric metric;
    private final Supplier<RayTracer> tracerFactory;
    private final RenderConfig config;

    /**
     * Phase 1 backward-compatible constructor. The supplied fixed-step
     * {@link Integrator} is wrapped in a {@link FixedStepRayTracer} per
     * worker thread.
     *
     * @param camera            provides initial states
     * @param metric            spacetime metric
     * @param integratorFactory one fresh fixed-step integrator per worker
     *                          (e.g. {@code RK4::new})
     * @param config            termination + shader settings
     */
    public Renderer(Camera camera, Metric metric,
                    Supplier<Integrator> integratorFactory,
                    RenderConfig config) {
        this(camera, metric, config,
             integratorFactoryToTracerFactory(integratorFactory, metric, config));
    }

    /**
     * Phase 2+ constructor taking a {@link RayTracer} factory directly.
     * Use this for adaptive integration or any other custom strategy.
     *
     * @param camera        provides initial states
     * @param metric        spacetime metric
     * @param tracerFactory one fresh {@link RayTracer} per worker thread
     * @param config        termination + shader settings
     */
    public static Renderer withRayTracer(Camera camera, Metric metric,
                                         Supplier<RayTracer> tracerFactory,
                                         RenderConfig config) {
        return new Renderer(camera, metric, config, tracerFactory);
    }

    /** Shared private constructor. Parameter order differs from the public
     *  Phase 1 constructor to avoid a {@code Supplier<X>} erasure collision. */
    private Renderer(Camera camera, Metric metric, RenderConfig config,
                     Supplier<RayTracer> tracerFactory) {
        if (camera == null || metric == null || tracerFactory == null || config == null) {
            throw new IllegalArgumentException("null argument");
        }
        this.camera = camera;
        this.metric = metric;
        this.tracerFactory = tracerFactory;
        this.config = config;
    }

    private static Supplier<RayTracer> integratorFactoryToTracerFactory(
            Supplier<Integrator> integratorFactory, Metric metric, RenderConfig config) {
        if (integratorFactory == null) {
            throw new IllegalArgumentException("integratorFactory must not be null");
        }
        return () -> new FixedStepRayTracer(integratorFactory.get(), metric, config);
    }

    /**
     * Render the full image. Returns a newly-allocated {@code float[H][W]} with
     * {@code image[0]} as the top row (j = 0 is top, FITS convention).
     * Pure — no internal state is mutated; two calls with identical configuration
     * produce byte-identical arrays.
     */
    public float[][] render() {
        final int W = camera.width();
        final int H = camera.height();
        final float[][] image = new float[H][W];

        final ThreadLocal<RayTracer> tlTracer = ThreadLocal.withInitial(tracerFactory);
        final ThreadLocal<double[]>  tlY      = ThreadLocal.withInitial(() -> new double[8]);

        IntStream.range(0, W * H).parallel().forEach(idx -> {
            int i = idx % W;
            int j = idx / W;
            double[] y = tlY.get();
            RayTracer tracer = tlTracer.get();

            camera.initialState(i, j, y);
            RayOutcome outcome = tracer.trace(y);
            image[j][i] = config.shader().shade(outcome, y);
        });

        return image;
    }

    /** Render and write a FITS file with a populated header. */
    public void renderToFits(Path outPath, String scene) throws IOException {
        float[][] image = render();
        FitsWriter.write(outPath, image, fitsMetadata(scene));
    }

    FitsWriter.Metadata fitsMetadata(String scene) {
        return new FitsWriter.Metadata(
                scene,
                metric.mass(),
                camera.rObs(),
                camera.inclination(),
                camera.fovRadians(),
                config.horizonCushion(),
                config.stepPolicyString(),
                config.maxSteps(),
                FitsWriter.gitHash());
    }
}
