@testitem "Ascertainment scales expected observations by a latent model" begin
    using ComposableTuringIDModels, Distributions, Random
    Random.seed!(51)
    Y = fill(100.0, 14)
    asc = Ascertainment(model = NegativeBinomialError(), latent_model = FixedIntercept(0.1))
    sim = as_turing_model(asc, missing, Y)().y_t
    @test length(sim) == length(Y)
    @test all(>=(0), sim)

    adw = ascertainment_dayofweek(PoissonError())
    @test length(as_turing_model(adw, missing, Y)().y_t) == length(Y)
end

@testitem "Ascertainment accepts a prior model (constant vs time-varying)" begin
    using ComposableTuringIDModels, Distributions, Random
    using ComposableTuringIDModels: AbstractPriorModel
    Random.seed!(53)
    # A bare Distribution is wrapped in an `Intercept`, giving a single constant
    # factor (one shared draw) — not `n` iid values. `as_prior`/`BroadcastPrior`
    # are gone; the constant semantics are preserved by the `Intercept` wrapping.
    asc_const = Ascertainment(PoissonError(), Normal(0.0, 0.1))
    @test asc_const.latent_model isa AbstractPriorModel
    sim = as_turing_model(asc_const, missing, fill(50.0, 8))().y_t
    @test length(sim) == 8
    @test all(>=(0), sim)
    # Pin the "single shared value" semantics: with a multiplicative transform the
    # ascertained expected series is a constant scaling of Y_t (same factor at
    # every t), i.e. one shared draw broadcast — not a per-t iid effect.
    # Use the exp-scale (default-style) transform so a negative shared draw still
    # yields a positive expected series (a plain `Y_t .* x` can go negative for a
    # `Normal` draw and make the inner `PoissonError` throw a DomainError — the
    # constant-scaling invariant we are pinning holds either way).
    exp_series = as_turing_model(
        Ascertainment(PoissonError(), Normal(0.0, 0.1);
            transform = (Y_t, x) -> Y_t .* exp.(x)),
        missing, fill(20.0, 6))().expected
    ratio = exp_series ./ 20.0
    @test all(≈(first(ratio)), ratio)
    # The constant factor is one shared draw: the underlying expected series is a
    # single global exp(θ) scaling. Recover the scaling with FixedIntercept.
    Y = fill(20.0, 6)
    asc_fixed = Ascertainment(PoissonError(), FixedIntercept(0.0);
        transform = (Y_t, x) -> Y_t .* x)
    @test asc_fixed.latent_model isa AbstractPriorModel
    # A latent model still gives a time-varying effect (existing behaviour).
    asc_tv = Ascertainment(PoissonError(), RandomWalk())
    @test asc_tv.latent_model isa AbstractPriorModel
    @test length(as_turing_model(asc_tv, missing, Y)().y_t) == length(Y)
end

@testitem "Aggregate sums expected observations over windows" begin
    using ComposableTuringIDModels, Random
    Random.seed!(52)
    agg = Aggregate(PoissonError(), [0, 0, 0, 0, 7, 0, 0])
    res = as_turing_model(agg, missing, fill(1.0, 28))()
    out = res.y_t
    @test length(out) == 28
    # Only the present (weekly) positions are non-zero.
    @test count(!=(0), out) == 4
    # `Aggregate` also threads through the uniform contract: its `expected` is the
    # scattered per-window expected means (so it can feed a `Split` stream).
    @test length(res.expected) == 28
    @test count(!=(0), res.expected) == 4
end

@testitem "ReportingCDF produces a length-n completeness curve" begin
    using ComposableTuringIDModels, Distributions

    # Role + interface conformance: the correction is a latent-role component.
    c = ReportingCDF([0.2, 0.6, 1.0])
    @test c isa AbstractLatentModel
    @test implements_latent_interface(c)

    # A precomputed curve is padded with ones to the requested length (older
    # reference days are fully reported).
    @test as_turing_model(c, 5)() == [0.2, 0.6, 1.0, 1.0, 1.0]
    # When the curve is at least as long as n, only its head is used.
    @test as_turing_model(ReportingCDF(collect(range(0.05, 1.0; length = 20))),
        5)() == collect(range(0.05, 1.0; length = 20))[1:5]

    # A curve out of [0, 1] is rejected at construction.
    @test_throws Exception ReportingCDF([0.2, 0.6, 1.5])
    @test_throws Exception ReportingCDF([-0.1, 0.6, 1.0])
    # A NON-monotonic curve is allowed: the correction is a free completeness
    # curve, not necessarily a CDF, so over-/under-reporting that recovers works.
    nonmono = ReportingCDF([0.6, 0.2, 0.9])
    @test as_turing_model(nonmono, 3)() == [0.6, 0.2, 0.9]

    # The distribution constructor builds the CDF from the released-CD
    # double-interval-censored PMF (the LatentDelay path): cumulative,
    # non-decreasing, ending at 1.
    cd = ReportingCDF(truncated(Normal(5.0, 2.0), 0.0, Inf))
    F = as_turing_model(cd, 10)()
    @test issorted(F)
    @test isapprox(F[end], 1.0)
    @test all(0 .<= F .<= 1 + 1e-8)
end

@testitem "RightTruncate construction and interface" begin
    using ComposableTuringIDModels, Distributions

    # Role + interface conformance (the nowcasting modifier is an observation model).
    o = RightTruncate(PoissonError(), [0.2, 0.6, 1.0])
    @test o isa AbstractObservationModel
    @test implements_observation_interface(o)
    # The convenience constructors wrap the correction in a ReportingCDF submodel.
    @test o.cdf_model isa ReportingCDF
    @test RightTruncate(PoissonError(),
        truncated(Normal(5.0, 2.0), 0.0, Inf)).cdf_model isa ReportingCDF
    # A correction submodel can be supplied directly.
    @test RightTruncate(PoissonError(),
        ReportingCDF([0.2, 0.6, 1.0])).cdf_model isa ReportingCDF
end

@testitem "RightTruncate scales expected eventual totals by the delay CDF" begin
    using ComposableTuringIDModels, Distributions, Random
    Random.seed!(60)

    # Huge expected means make the Poisson draw ≈ its mean, so the realised
    # observed-so-far reveals the applied scaling. With a 3-bin CDF and a length-6
    # series, reference days older than the delay support are fully reported.
    n = 6
    Y = fill(1.0e6, n)
    F = [0.2, 0.5, 1.0]                       # completeness by age 0, 1, 2
    obs = RightTruncate(PoissonError(), F)
    sim = as_turing_model(obs, missing, Y)().y_t
    # Completeness by reference day t = reverse of [F; ones(n - length(F))]:
    # the most recent day (t = n, age 0) is scaled by F[1], the oldest by 1.
    expected_scale = reverse(vcat(F, ones(n - length(F))))
    @test all(abs.(sim ./ Y .- expected_scale) .< 0.01)

    # The oldest day (t = 1) is fully reported; the most recent (t = n) least so.
    @test isapprox(sim[1] / Y[1], 1.0; atol = 0.01)
    @test isapprox(sim[n] / Y[n], F[1]; atol = 0.01)
end

@testitem "RightTruncate with a complete delay reduces to the inner model" begin
    using ComposableTuringIDModels, Random
    Random.seed!(61)
    n = 8
    Y = fill(50.0, n)

    # A fully-reported delay (CDF all ones, or a single bin F = [1.0] padded to
    # ones) means every reference day is complete, so the modifier is a no-op:
    # it reduces exactly to the wrapped error model.
    for F in (ones(n), [1.0])
        Random.seed!(123)
        scaled = as_turing_model(RightTruncate(PoissonError(), F), missing, Y)().y_t
        Random.seed!(123)
        inner = as_turing_model(PoissonError(), missing, Y)().y_t
        @test scaled == inner
    end
end

@testitem "RightTruncate accepts a custom (non-default) correction submodel" begin
    using ComposableTuringIDModels, Random
    Random.seed!(64)
    # The correction is a submodel, so a user can supply any latent component
    # producing the length-n completeness curve — here a flat 0.5 correction
    # (not a CDF shape at all), exercising the composable-submodel slot.
    custom = TransformLatentModel(FixedIntercept(0.0), x -> fill(0.5, length(x)))
    obs = RightTruncate(PoissonError(), custom)
    @test implements_observation_interface(obs)

    Y = fill(1.0e6, 5)
    sim = as_turing_model(obs, missing, Y)().y_t
    # Every reference day is down-weighted by the flat 0.5 correction.
    @test all(abs.(sim ./ Y .- 0.5) .< 0.01)
end

@testitem "RightTruncate simulate-then-condition and length handling" begin
    using ComposableTuringIDModels, Random
    Random.seed!(62)
    n = 10
    Y = fill(100.0, n)
    F = collect(range(0.1, 1.0; length = n))
    obs = RightTruncate(PoissonError(), F)

    sim = as_turing_model(obs, missing, Y)().y_t
    @test length(sim) == n
    @test all(>=(0), sim)
    # Conditioning on the simulated data returns it.
    @test as_turing_model(obs, sim, Y)().y_t == sim

    # A CDF shorter than the series is fine: older days are taken complete.
    short = RightTruncate(PoissonError(), [0.3, 0.7, 1.0])
    @test length(as_turing_model(short, missing, fill(10.0, 20))().y_t) == 20
    # A CDF longer than the series is also fine (only its head is used).
    long = RightTruncate(PoissonError(), collect(range(0.05, 1.0; length = 30)))
    @test length(as_turing_model(long, missing, fill(10.0, 5))().y_t) == 5
end

@testitem "RightTruncate composes with a renewal model end-to-end" begin
    using ComposableTuringIDModels, Distributions, Random
    Random.seed!(63)
    gen_int = [0.2, 0.3, 0.5]
    n = 25
    model = IDModel(
        Renewal(gen_int; rt = RandomWalk(), initialisation = Normal()),
        RightTruncate(NegativeBinomialError(),
            truncated(Normal(4.0, 1.5), 0.0, Inf)))

    sim = as_turing_model(model, missing, n)()
    y = sim.generated_y_t
    @test length(y) == n
    # The composed posterior builds and evaluates when conditioned on the data,
    # reproducing the conditioned observations as a generated quantity.
    @test as_turing_model(model, y, n)().generated_y_t == y

    # It also composes as a single-stream Split.
    stk = Split((reported = RightTruncate(PoissonError(),
        truncated(Normal(4.0, 1.5), 0.0, Inf)),))
    out = as_turing_model(stk, (reported = missing,), fill(100.0, n))().y_t
    @test length(out) == 1
    @test length(out.reported) == n
end

@testitem "PrefixObservationModel prefixes observation parameters" begin
    using ComposableTuringIDModels, Random
    Random.seed!(53)
    pom = PrefixObservationModel(model = NegativeBinomialError(), prefix = "Test")
    names = string.(collect(keys(rand(as_turing_model(pom, missing, fill(10.0, 5))))))
    @test any(startswith("Test."), names)
end

@testitem "RecordExpectedObs and TransformObservationModel wrap an error model" begin
    using ComposableTuringIDModels, Random
    Random.seed!(54)
    Y = fill(10.0, 30)
    reo = RecordExpectedObs(NegativeBinomialError())
    @test length(as_turing_model(reo, missing, Y)().y_t) == length(Y)

    tom = TransformObservationModel(NegativeBinomialError())
    @test length(as_turing_model(tom, missing, Y)().y_t) == length(Y)
end

@testitem "ReportTriangle constructs the reporting triangle and masks t+d≤now" begin
    using ComposableTuringIDModels, Distributions, Random
    Random.seed!(70)
    obs = ReportTriangle(PoissonError(), [0.5, 0.3, 0.2])     # Dmax = 2

    # Role + interface conformance: it is an observation model whose `y_t` is a
    # reporting triangle (built from a matrix via `define_y_t`).
    @test obs isa AbstractObservationModel
    N = [10 5 2; 12 6 3; 14 7 4; 16 8 5]
    tri = define_y_t(obs, N, fill(20.0, 4); now = 4)
    @test tri isa ReportingTriangle
    @test implements_observation_interface(obs; y_t = tri, Y_t = fill(20.0, 4))

    # The mask is exactly `t + d ≤ now` (delays d = 0,1,2 on the columns).
    expected = BitMatrix([t + d <= 4 for t in 1:4, d in 0:2])
    @test tri.observed == expected
    @test tri.observed == BitMatrix([true true true; true true true;
                     true true false; true false false])
    @test tri.Dmax == 2

    # A triangle whose Dmax does not match the model's delay PMF is rejected.
    @test_throws Exception define_y_t(obs,
        ReportingTriangle(zeros(Int, 4, 5), trues(4, 5), 4), fill(20.0, 4))
    # A count matrix with the wrong number of delay columns is rejected.
    @test_throws Exception define_y_t(obs, zeros(Int, 4, 5), fill(20.0, 4))
end

@testitem "ReportTriangle takes the delay as a composable PMF submodel" begin
    using ComposableTuringIDModels, Distributions, Random
    Random.seed!(76)

    # ReportingPMF is the delay submodel (a latent-role component), mirroring how
    # ReportingCDF is RightTruncate's correction submodel.
    pm = ReportingPMF([0.5, 0.3, 0.2])
    @test pm isa AbstractLatentModel
    @test implements_latent_interface(pm)
    @test as_turing_model(pm, 10)() == [0.5, 0.3, 0.2]   # n is ignored (PMF by delay)
    # The distribution constructor builds the PMF via the released-CD path.
    pmd = ReportingPMF(truncated(Normal(2.0, 1.0), 0.0, Inf))
    @test isapprox(sum(as_turing_model(pmd, 5)()), 1.0)
    # A PMF that does not sum to 1 / is negative is rejected.
    @test_throws Exception ReportingPMF([0.5, 0.3])
    @test_throws Exception ReportingPMF([-0.1, 0.6, 0.5])

    # The three ReportTriangle constructors all wrap a delay submodel.
    @test ReportTriangle(PoissonError(), [0.5, 0.3, 0.2]).delay_model isa ReportingPMF
    @test ReportTriangle(PoissonError(),
        truncated(Normal(2.0, 1.0), 0.0, Inf)).delay_model isa ReportingPMF
    obs = ReportTriangle(PoissonError(), ReportingPMF([0.6, 0.25, 0.15]))
    @test obs.delay_model isa ReportingPMF

    # Dmax is read statically from the delay submodel (before the PMF is sampled),
    # so the triangle sizes correctly and the model simulates.
    tri = define_y_t(obs, fill(0, 5, 3), fill(20.0, 5); now = 5)
    @test tri.Dmax == 2
    @test as_turing_model(obs, missing, fill(20.0, 6))().y_t isa ReportingTriangle
end

@testitem "define_y_t builds a triangle from a matrix and a long-form table" begin
    using ComposableTuringIDModels, Random
    Random.seed!(71)
    obs = ReportTriangle(PoissonError(), [0.6, 0.25, 0.15])   # Dmax = 2
    N = [10 5 2; 12 6 3; 14 7 4]                            # 3 reference days

    from_matrix = define_y_t(obs, N, fill(20.0, 3); now = 3)

    # The same data as long-form `(reference, delay, count)` rows.
    refs = Int[];
    dls = Int[];
    cts = Int[]
    for t in 1:3, d in 0:2

        push!(refs, t)
        push!(dls, d)
        push!(cts, N[t, d + 1])
    end
    tbl = (reference = refs, delay = dls, count = cts)
    from_table = define_y_t(obs, tbl, fill(20.0, 3); now = 3)

    @test from_table.counts == N
    @test from_table.observed == from_matrix.observed
    @test from_table.Dmax == from_matrix.Dmax == 2

    # A `missing` series gives a fully observed triangle (now = n + Dmax) of
    # `missing` cells, sized to `Y_t`.
    mt = define_y_t(obs, missing, fill(20.0, 5))
    @test size(mt.counts) == (5, 3)
    @test all(ismissing, mt.counts)
    @test all(mt.observed)
end

@testitem "ReportTriangle simulates, conditions, and recovers per-cell means" begin
    using ComposableTuringIDModels, Distributions, Random
    Random.seed!(72)
    pmf = [0.5, 0.3, 0.2]
    obs = ReportTriangle(PoissonError(), pmf)
    n = 7
    Y = fill(40.0, n)

    # Simulate: a fully observed triangle of non-negative integer counts.
    sim = as_turing_model(obs, missing, Y)().y_t
    @test sim isa ReportingTriangle
    @test size(sim.counts) == (n, 3)
    @test all(c -> c >= 0, sim.counts[sim.observed])

    # Condition on the simulated triangle: the observed cells round-trip.
    cond = as_turing_model(obs, sim, Y)().y_t
    @test cond.counts[sim.observed] == sim.counts[sim.observed]

    # Each cell mean is `Y_t · p[d+1]`: average a stack of simulated triangles
    # (a plain sum / count, to avoid a Statistics dep in the clean test env).
    Random.seed!(73)
    draws = [Float64.(as_turing_model(obs, missing, fill(100.0, 5))().y_t.counts)
             for _ in 1:4000]
    cell_means = sum(draws) ./ length(draws)
    for d in 0:2
        @test isapprox(cell_means[1, d + 1], 100.0 * pmf[d + 1]; rtol = 0.05)
    end
end

@testitem "ReportTriangle observed row-sums reconcile with RightTruncate marginal" begin
    using ComposableTuringIDModels, Distributions, Random
    Random.seed!(74)
    # #50's consistency check: the observed row-sums of the triangle are the
    # marginal that the right-truncation (CDF-scaling) nowcast conditions on. With
    # `now = n`, reference day `t` has age `a = n - t` and observed delays
    # `d = 0 … a`, so its expected observed row-sum is
    #   Σ_{d=0}^{a} μ_t · p[d+1] = μ_t · F[a+1],
    # exactly the CDF-scaled expected observed-so-far that `RightTruncate` applies.
    # We assert the noise-free triangle row-sums equal the `RightTruncate` scaling
    # (recovered here from the `ReportingCDF` completeness curve) to tolerance.
    pmf = [0.5, 0.3, 0.2]
    Dmax = length(pmf) - 1
    n = 8
    μ = collect(range(20.0, 80.0; length = n))      # expected eventual totals

    # Expected (noise-free) observed row-sums from the masked triangle.
    triangle_rowsums = zeros(n)
    for t in 1:n, d in 0:Dmax

        (t + d <= n) || continue                    # the `t + d ≤ now` mask, now = n
        triangle_rowsums[t] += μ[t] * pmf[d + 1]
    end

    # The `RightTruncate` scaling: build the completeness curve from the SAME
    # delay (via `ReportingCDF`, the released-CD cumsum path) and reverse it onto
    # the reference-day axis, exactly as `RightTruncate` does internally.
    completeness = as_turing_model(ReportingCDF(cumsum(pmf)), n)()
    right_truncate_scaled = μ .* reverse(completeness)

    @test triangle_rowsums≈right_truncate_scaled rtol=1e-10
    # The oldest day is fully reported; the most recent only its delay-0 fraction.
    @test triangle_rowsums[1]≈μ[1] rtol=1e-10
    @test triangle_rowsums[n]≈μ[n] * pmf[1] rtol=1e-10
end

@testitem "ReportTriangle composes with the renewal infection process" begin
    using ComposableTuringIDModels, Distributions, Random
    Random.seed!(75)
    gen_int = [0.2, 0.3, 0.5]
    renewal = Renewal(gen_int; rt = RandomWalk(), initialisation = Normal())

    # Released-CD discretised-delay constructor (the LatentDelay path).
    obs = ReportTriangle(NegativeBinomialError(),
        truncated(Normal(2.0, 1.0), 0.0, Inf))
    model = IDModel(renewal, obs)

    sim = as_turing_model(model, missing, 15)()
    @test sim.generated_y_t isa ReportingTriangle
    @test length(sim.I_t) == 15

    # Conditioning on the simulated triangle builds and evaluates.
    cond = as_turing_model(model, sim.generated_y_t, 15)()
    @test cond.generated_y_t isa ReportingTriangle

    # A triangle stream composes alongside a plain-vector stream under Split.
    stk = Split((
        tri = ReportTriangle(PoissonError(), [0.6, 0.25, 0.15]),
        cases = PoissonError()))
    out = as_turing_model(stk, (tri = missing, cases = missing), fill(20.0, 8))().y_t
    @test keys(out) == (:tri, :cases)
    @test out.tri isa ReportingTriangle
    @test length(out.cases) == 8
end

@testitem "UncertainDelay builds an inferred-delay LatentDelay" begin
    using ComposableTuringIDModels, Distributions
    u = UncertainDelay(LogNormal,
        [Normal(1.5, 0.4), truncated(Normal(0.4, 0.2), 0, Inf)]; D = 20.0)
    @test u isa ComposableTuringIDModels.AbstractPriorModel
    obs = LatentDelay(NegativeBinomialError(), u)
    @test obs isa AbstractObservationModel
    # A fixed horizon `D` is required so the PMF length is constant across draws.
    @test_throws AssertionError UncertainDelay(LogNormal,
        [Normal(1.5, 0.4), truncated(Normal(0.4, 0.2), 0, Inf)]; D = nothing)
    # The fixed-PMF constructors are unaffected (no regression).
    @test LatentDelay(PoissonError(), [0.3, 0.4, 0.3]) isa AbstractObservationModel
    @test LatentDelay(PoissonError(), truncated(Normal(5, 2), 0, Inf)) isa
          AbstractObservationModel
end

@testitem "UncertainDelay samples a valid delay PMF per draw" begin
    using ComposableTuringIDModels, Distributions, Random
    Random.seed!(101)
    u = UncertainDelay(LogNormal,
        [Normal(1.5, 0.4), truncated(Normal(0.4, 0.2), 0, Inf)]; D = 15.0)
    # Each prior draw builds a normalised, non-negative PMF of constant length.
    lens = Int[]
    for _ in 1:20
        pmf = as_turing_model(u)()
        @test isapprox(sum(pmf), 1.0)
        @test all(>=(0), pmf)
        push!(lens, length(pmf))
    end
    @test all(==(first(lens)), lens)

    obs = LatentDelay(PoissonError(), u)
    Y = fill(100.0, 40)
    sim = as_turing_model(obs, missing, Y)().y_t
    @test length(sim) == length(Y)
    @test all(>=(0), filter(!ismissing, sim))
end

@testitem "LatentDelay recovers an uncertain delay's parameters" tags=[:sample] begin
    using ComposableTuringIDModels, Distributions, Turing, Random, Statistics
    Random.seed!(20240716)
    # Simulate reports from a KNOWN LogNormal reporting delay convolving a smooth
    # epidemic bump, then infer the delay parameters through the priors seam.
    n = 60
    t = 1:n
    Y_t = 800.0 .* exp.(-((t .- 32) ./ 10.0) .^ 2) .+ 5.0
    μ0, σ0 = 1.5, 0.40
    truth = LatentDelay(PoissonError(), LogNormal(μ0, σ0); D = 20.0)
    y = as_turing_model(truth, missing, Y_t)().y_t

    fit = LatentDelay(PoissonError(),
        UncertainDelay(LogNormal,
            [Normal(1.1, 0.4), truncated(Normal(0.6, 0.3), 0, Inf)]; D = 20.0))
    chain = sample(as_turing_model(fit, y, Y_t),
        NUTS(0.8; adtype = Turing.AutoForwardDiff()), 1000; progress = false)

    # The delay parameters are namespaced under the `delay` prior slot.
    dvec = vec(chain[@varname(delay.θ)])
    meanlog = mean(getindex.(dvec, 1))
    sdlog = mean(getindex.(dvec, 2))
    @test isapprox(meanlog, μ0; atol = 0.15)
    @test isapprox(sdlog, σ0; atol = 0.15)
end
