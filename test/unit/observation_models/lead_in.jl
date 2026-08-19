@testitem "observation_lead_in sums nested LatentDelay lead-ins" begin
    using ComposableTuringIDModels, Distributions

    # A bare error family drops nothing.
    @test observation_lead_in(PoissonError()) == 0

    # One delay costs `length(pmf) - 1`.
    d1 = fill(1 / 15, 15)
    @test observation_lead_in(LatentDelay(PoissonError(), d1)) == 14

    # Nested delays sum.
    d2 = fill(1 / 30, 30)
    chain = LatentDelay(LatentDelay(PoissonError(), d1), d2)
    @test observation_lead_in(chain) == 14 + 29
end

@testitem "observation_lead_in walks every observation nesting shape" begin
    using ComposableTuringIDModels, Distributions

    pmf = fill(1 / 5, 5)                       # lead-in 4
    delayed = LatentDelay(PoissonError(), pmf)

    # Modifiers that keep the series length pass the lead-in through.
    @test observation_lead_in(Ascertainment(delayed, FixedIntercept(0.0))) == 4
    @test observation_lead_in(Aggregate(delayed, [0, 0, 0, 0, 0, 0, 7])) == 4
    @test observation_lead_in(RightTruncate(delayed, fill(1.0, 10))) == 4
    @test observation_lead_in(TransformObservationModel(delayed)) == 4
    @test observation_lead_in(PrefixObservationModel(delayed, "a")) == 4
    @test observation_lead_in(RecordExpectedObs(delayed)) == 4

    # A ReportTriangle's delay is a nowcasting kernel, not a series-shortening
    # convolution, so it contributes nothing of its own.
    @test observation_lead_in(ReportTriangle(PoissonError(), fill(1 / 3, 3))) == 0

    # Streams of a Split run in parallel: the shared lead-in, not their sum.
    split = Split((cases = delayed, deaths = LatentDelay(PoissonError(), pmf)))
    @test observation_lead_in(split) == 4
    # A Split placed under a delay adds that delay's lead-in.
    @test observation_lead_in(LatentDelay(split, fill(1 / 3, 3))) == 4 + 2
    # A strata-template Split reports its template's lead-in.
    @test observation_lead_in(Split(delayed, [1.0 1.0])) == 4

    # Streams that drop different amounts report one lead-in each, and a
    # surrounding delay adds to every stream.
    mismatched = Split(
        (
            cases = delayed, deaths = LatentDelay(PoissonError(), fill(1 / 9, 9)),
        )
    )
    @test observation_lead_in(mismatched) == (cases = 4, deaths = 8)
    @test observation_lead_in(LatentDelay(mismatched, fill(1 / 3, 3))) ==
        (cases = 6, deaths = 10)
    @test observation_lead_in(Ascertainment(mismatched, FixedIntercept(0.0))) ==
        (cases = 4, deaths = 8)

    # A stacked chain sums each delay it passes through.
    stacked = LatentDelay(
        Ascertainment(
            LatentDelay(RecordExpectedObs(NegativeBinomialError()), fill(1 / 4, 4)),
            FixedIntercept(0.0)
        ),
        pmf
    )
    @test observation_lead_in(stacked) == 4 + 3

    # An IDModel and an IDProblem delegate to their observation model.
    infection = DirectInfections(; Z = RandomWalk(), initialisation = Normal())
    @test observation_lead_in(IDModel(infection, stacked)) == 7
    problem = IDProblem(
        infection = infection, observation_model = stacked, tspan = (1, 30)
    )
    @test observation_lead_in(problem) == 7
end

@testitem "observation_lead_in covers uncertain and time-varying delays" begin
    using ComposableTuringIDModels, Distributions, Random
    Random.seed!(1064)

    # An `UncertainDelay` has a fixed horizon, so its pmf length is known
    # before the parameters are drawn.
    uncertain = UncertainDelay(
        LogNormal, [Normal(1.5, 0.4), truncated(Normal(0.4, 0.2), 0, Inf)]; D = 8.0
    )
    obs_uncertain = LatentDelay(NegativeBinomialError(), uncertain)
    @test observation_lead_in(obs_uncertain) == 7
    sim = as_turing_model(obs_uncertain, missing, fill(100.0, 40))()
    @test length(sim.expected) == 40 - observation_lead_in(obs_uncertain)

    # A deterministic per-time sequence of pmfs: all pmfs share a length.
    n = 20
    pmfs = [fill(1 / 4, 4) for _ in 1:n]
    obs_tv = LatentDelay(PoissonError(), pmfs)
    @test observation_lead_in(obs_tv) == 3
    sim_tv = as_turing_model(obs_tv, missing, fill(100.0, n))()
    @test length(sim_tv.expected) == n - observation_lead_in(obs_tv)

    # A delay whose length cannot be read off the component is an error, not a
    # silent zero.
    opaque = LatentDelay(PoissonError(), RandomWalk())
    @test_throws ArgumentError observation_lead_in(opaque)
end

@testitem "observation_lead_in predicts the expected-series length" begin
    using ComposableTuringIDModels, Distributions, Random
    Random.seed!(1065)

    n = 40
    pmf = fill(1 / 6, 6)
    chains = [
        LatentDelay(PoissonError(), pmf),
        LatentDelay(Ascertainment(PoissonError(), FixedIntercept(0.0)), pmf),
        LatentDelay(LatentDelay(PoissonError(), pmf), fill(1 / 4, 4)),
        Ascertainment(LatentDelay(PoissonError(), pmf), FixedIntercept(0.0)),
        TransformObservationModel(LatentDelay(PoissonError(), pmf); transform = identity),
        RightTruncate(LatentDelay(PoissonError(), pmf), fill(1.0, n)),
    ]
    for obs in chains
        sim = as_turing_model(obs, missing, fill(100.0, n))()
        @test length(sim.expected) == n - observation_lead_in(obs)
    end
end

@testitem "observation_coverage reports the observations a chain scores" begin
    using ComposableTuringIDModels, Distributions

    obs = LatentDelay(LatentDelay(PoissonError(), fill(1 / 15, 15)), fill(1 / 30, 30))
    y = fill(10, 60)

    # The naive call: `n = length(y)` silently drops the lead-in.
    naive = observation_coverage(obs, y, length(y))
    @test naive.lead_in == 43
    @test naive.n_observations == 60
    @test naive.n_scored == 17
    @test naive.n_unscored == 43

    # The documented idiom scores every observation.
    n = length(y) + observation_lead_in(obs)
    full = observation_coverage(obs, y, n)
    @test full.n_scored == 60
    @test full.n_unscored == 0

    # Simulating from the prior scores whatever the chain produces.
    @test observation_coverage(obs, missing, 60).n_unscored == 0

    # An `IDProblem` reads its series length from `tspan`.
    problem = IDProblem(
        infection = DirectInfections(; Z = RandomWalk(), initialisation = Normal()),
        observation_model = obs, tspan = (1, 60)
    )
    @test observation_coverage(problem, (; y_t = y)).n_unscored == 43
end

@testitem "observation_coverage reports a Split stream by stream" begin
    using ComposableTuringIDModels, Distributions, Random
    Random.seed!(1067)

    # Two streams with different delays, as a cases/deaths model has.
    streams = Split(
        (
            cases = LatentDelay(PoissonError(), fill(1 / 5, 5)),
            deaths = LatentDelay(PoissonError(), fill(1 / 20, 20)),
        )
    )
    lead_in = observation_lead_in(streams)
    @test lead_in == (cases = 4, deaths = 19)

    # Each stream scores the last `n - lead_in[stream]` of its own series, so
    # the two are fully scored at different data lengths.
    n = 40
    y = (cases = fill(10, n - 4), deaths = fill(1, n - 19))
    cover = observation_coverage(streams, y, n)
    @test cover.cases.n_unscored == 0
    @test cover.deaths.n_unscored == 0
    # Equal-length streams leave the longer-delayed one short.
    equal_length = observation_coverage(
        streams, (cases = fill(10, n), deaths = fill(1, n)), n
    )
    @test equal_length.cases.n_unscored == 4
    @test equal_length.deaths.n_unscored == 19

    # The reported lead-ins match the expected series each stream is scored
    # against.
    sim = as_turing_model(streams, (cases = missing, deaths = missing), fill(100.0, n))()
    @test length(sim.expected.cases) == n - 4
    @test length(sim.expected.deaths) == n - 19

    # A shared lead-in still reports each stream's data separately.
    shared = Split(
        (
            cases = LatentDelay(PoissonError(), fill(1 / 5, 5)),
            deaths = LatentDelay(PoissonError(), fill(1 / 5, 5)),
        )
    )
    @test observation_lead_in(shared) == 4
    shared_cover = observation_coverage(
        shared, (cases = fill(10, n), deaths = fill(1, n - 4)), n
    )
    @test shared_cover.cases.n_unscored == 4
    @test shared_cover.deaths.n_unscored == 0
end

@testitem "a multi-delay chain scores only observations above the lead-in" begin
    using ComposableTuringIDModels, Distributions, Random
    using DynamicPPL: VarInfo, logjoint
    Random.seed!(1066)

    d1 = fill(1 / 15, 15)
    d2 = fill(1 / 30, 30)
    obs = LatentDelay(LatentDelay(PoissonError(), d1), d2)
    # A tight random walk keeps the prior draw's expected series near 100, so
    # the log-joint stays well scaled and a perturbation of a scored entry is
    # not lost to floating-point rounding.
    infection = Renewal(;
        generation_time = [0.2, 0.3, 0.3, 0.2],
        rt = RandomWalk(; init = Normal(0.0, 0.05), ϵ_t = Normal(0.0, 0.01)),
        initialisation = Normal(log(100.0), 0.1)
    )
    model = IDModel(infection, obs)
    lead_in = observation_lead_in(model)
    @test lead_in == 43

    # The exact arithmetic the issue pins: 60 observations, 17 expected.
    n = 60
    sim = as_turing_model(model, missing, n)()
    @test length(sim.generated_y_t) == n
    @test length(sim.expected_y_t) == n - lead_in == 17

    # With `n = length(y)` the leading `lead_in` observations are unscored:
    # perturbing them leaves the log-joint bit-identical.
    y = fill(100, n)
    vi = VarInfo(as_turing_model(model, y, n))
    lj(yv) = logjoint(as_turing_model(model, yv, n), vi)
    base = lj(y)
    for i in (1, lead_in ÷ 2, lead_in)
        perturbed = copy(y)
        perturbed[i] = 10^6
        @test lj(perturbed) == base
    end
    # The first scored index is `lead_in + 1`.
    scored = copy(y)
    scored[lead_in + 1] = 10^6
    @test lj(scored) != base
    @test observation_coverage(model, y, n).n_unscored == lead_in

    # With the lead-in added back, every observation is scored.
    m = length(y) + lead_in
    vi_full = VarInfo(as_turing_model(model, y, m))
    lj_full(yv) = logjoint(as_turing_model(model, yv, m), vi_full)
    base_full = lj_full(y)
    for i in (1, lead_in, lead_in + 1, n)
        perturbed = copy(y)
        perturbed[i] = 10^6
        @test lj_full(perturbed) != base_full
    end
    @test observation_coverage(model, y, m).n_unscored == 0
end

@testitem "observation_coverage counts each data shape's observations" begin
    using ComposableTuringIDModels, Distributions

    pmf = fill(1 / 5, 5)                       # lead-in 4
    n = 30

    # A plain vector: the series itself carries the observations.
    plain = LatentDelay(PoissonError(), pmf)
    @test observation_coverage(plain, fill(10, n), n) ==
        (n_observations = 30, lead_in = 4, n_scored = 26, n_unscored = 4)

    # A `MissingObservations` carrier reports the series it stands in for.
    carrier = ComposableTuringIDModels.MissingObservations(
        fill(10.0, n), fill(true, n)
    )
    @test observation_coverage(plain, carrier, n).n_observations == 30

    # A stratified model's data is `strata x time`: the time axis counts.
    @test observation_coverage(plain, fill(10, 2, n), (2, n)).n_observations == 30

    # `BinomialError` takes `(; y, N)`. The trials are known covariates, not
    # observations, so only the `y` field is counted.
    binom = LatentDelay(BinomialError(), pmf)
    scalar_N = observation_coverage(binom, (y = fill(3, n), N = 20), n)
    @test scalar_N ==
        (n_observations = 30, lead_in = 4, n_scored = 26, n_unscored = 4)
    vector_N = observation_coverage(binom, (y = fill(3, n), N = fill(20, n)), n)
    @test vector_N == scalar_N
    # Adding the lead-in back scores every trial-weighted observation.
    @test observation_coverage(
        binom, (y = fill(3, n), N = 20), n + observation_lead_in(binom)
    ).n_unscored == 0
    # Simulating scores whatever the chain produces.
    @test observation_coverage(binom, (y = missing, N = 20), n).n_unscored == 0
end

@testitem "observation_coverage counts a ReportTriangle's reference days" begin
    using ComposableTuringIDModels, Distributions

    n = 30
    tri = ReportTriangle(PoissonError(), fill(1 / 3, 3))    # Dmax = 2

    # The reference days run down the rows; the delay columns are not
    # observations of their own.
    counts = fill(2, n, 3)
    @test observation_coverage(tri, counts, n) ==
        (n_observations = 30, lead_in = 0, n_scored = 30, n_unscored = 0)

    # An already-built `ReportingTriangle` reports the same.
    built = define_y_t(tri, counts, fill(20.0, n))
    @test observation_coverage(tri, built, n).n_observations == 30

    # Simulating sizes the triangle from the expected series.
    @test observation_coverage(tri, missing, n) ==
        (n_observations = 30, lead_in = 0, n_scored = 30, n_unscored = 0)

    # Under a delay the triangle needs `n - lead_in` reference days, and a
    # full-length one is over-supplied.
    delayed = LatentDelay(tri, fill(1 / 5, 5))
    @test observation_lead_in(delayed) == 4
    @test observation_coverage(delayed, counts, n).n_unscored == 4
    @test observation_coverage(delayed, fill(2, n - 4, 3), n).n_unscored == 0
end

@testitem "observation_coverage keys a Split by stream, not by data field" begin
    using ComposableTuringIDModels, Distributions

    n = 40
    streams = Split(
        (
            cases = LatentDelay(PoissonError(), fill(1 / 5, 5)),
            deaths = LatentDelay(PoissonError(), fill(1 / 20, 20)),
        )
    )

    # One series shared by both streams still reports per stream.
    shared = observation_coverage(streams, fill(10, n), n)
    @test shared.cases.n_unscored == 4
    @test shared.deaths.n_unscored == 19

    # A stream taking `(; y, N)` data is counted by its own contract, not
    # unpacked as if `N` were another stream.
    mixed = Split((cases = PoissonError(), positivity = BinomialError()))
    cover = observation_coverage(
        mixed,
        (cases = fill(10, n), positivity = (y = fill(3, n), N = 50)),
        n
    )
    @test cover.cases == (
        n_observations = 40, lead_in = 0, n_scored = 40, n_unscored = 0,
    )
    @test cover.positivity == cover.cases
end

@testitem "forecast follows the fitted series length" tags = [:sample] begin
    using ComposableTuringIDModels, Distributions, Turing, Random
    Random.seed!(1068)

    # A delay chain fitted with the documented idiom runs its infection series
    # over `length(y) + lead_in` steps, so the forecast must be rebuilt from
    # that length rather than from `length(y)`.
    obs = LatentDelay(PoissonError(), fill(1 / 6, 6))
    model = IDModel(
        DirectInfections(; Z = RandomWalk(), initialisation = Normal(1.0, 0.5)),
        obs
    )
    y = fill(5, 20)
    n = length(y) + observation_lead_in(model)
    @test observation_coverage(model, y, n).n_unscored == 0

    chain = sample(
        as_turing_model(model, y, n), Prior(), 20; progress = false
    )
    h = 3

    # Rebuilding at `length(y) + h` asks for a shorter latent stream than the
    # chain holds; say so rather than failing on a dimension mismatch.
    err = try
        forecast(model, y, chain, h)
        nothing
    catch e
        e
    end
    @test err isa ErrorException
    @test occursin("observation_lead_in", err.msg)

    fc = forecast(model, y, chain, h; n = n)
    @test size(fc, 1) == 20

    # An `IDProblem` records the fitted length in its `tspan`, so it needs no
    # keyword.
    problem = IDProblem(
        infection = model.infection_model, observation_model = obs,
        tspan = (1, n)
    )
    @test observation_coverage(problem, (; y_t = y)).n_unscored == 0
    pchain = sample(
        as_turing_model(problem, (; y_t = y)), Prior(), 20; progress = false
    )
    pfc = forecast(problem, y, pchain, h)
    @test size(pfc, 1) == 20
end

@testitem "observation_coverage counts an Aggregate's windows" begin
    using ComposableTuringIDModels, Distributions

    n = 28
    weekly = Aggregate(PoissonError(), [0, 0, 0, 0, 0, 0, 7])

    # Only the day closing each week reaches the error model, so a full daily
    # series carries four observations rather than 28.
    @test observation_coverage(weekly, fill(10, n), n) ==
        (n_observations = 4, lead_in = 0, n_scored = 4, n_unscored = 0)

    # Simulating is sized from the expected series, in windows.
    @test observation_coverage(weekly, missing, n).n_observations == 4

    # A delay before the aggregation drops time points, and the windows are
    # counted in what is left of the series.
    delayed = LatentDelay(weekly, fill(1 / 8, 8))
    @test observation_lead_in(delayed) == 7
    @test observation_coverage(delayed, fill(10, n - 7), n) ==
        (n_observations = 3, lead_in = 7, n_scored = 3, n_unscored = 0)

    # An aggregation does not right-align either: a full-length series over
    # the same chain is a call that fails, not one that drops its head.
    @test observation_coverage(delayed, fill(10, n), n).n_unscored == 1
end
