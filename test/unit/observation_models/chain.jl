@testitem "Chain composes ordered streams with name prefixes" begin
    using ComposableTuringIDModels, Distributions, Random
    Random.seed!(70)
    chain = Chain((cases = PoissonError(), deaths = NegativeBinomialError()))
    @test chain isa AbstractObservationModel
    @test implements_observation_interface(chain)

    yt = (cases = missing, deaths = missing)
    cm = as_turing_model(chain, yt, fill(10.0, 10))
    # Each stream's variables are prefixed by its name.
    names = string.(collect(keys(rand(cm))))
    @test any(startswith("cases."), names)
    @test any(startswith("deaths."), names)

    # The uniform return contract: per-stream `y_t` / `expected` NamedTuples in
    # chain order.
    out = cm()
    @test keys(out) == (:y_t, :expected)
    @test keys(out.y_t) == (:cases, :deaths)
    @test length(out.y_t.cases) == 10
    @test length(out.y_t.deaths) == 10
end

@testitem "Chain threads each stream's expected output into the next" begin
    using ComposableTuringIDModels, Distributions, Random
    Random.seed!(71)
    # cases: bare error on the incoming infections (expected == input).
    # deaths: observed DOWNSTREAM of the expected cases, delayed (10 → 8) and
    # ascertained ×0.1. A parallel Split would feed deaths the raw infections.
    chain = Chain((
        cases = PoissonError(),
        deaths = LatentDelay(
            Ascertainment(PoissonError(), FixedIntercept(log(0.1))),
            [0.2, 0.3, 0.5])))
    out = as_turing_model(
        chain, (cases = missing, deaths = missing), fill(100.0, 10))()

    # Stream 1 sees the raw incoming series.
    @test length(out.expected.cases) == 10
    @test all(≈(100.0), out.expected.cases)
    # Stream 2 sees stream 1's expected output, delayed (10 → 8) and scaled ×0.1.
    @test length(out.expected.deaths) == 8
    @test all(≈(10.0), out.expected.deaths)
end

@testitem "Chain propagates an upstream stream's own ascertainment downstream" begin
    using ComposableTuringIDModels, Distributions, Random
    Random.seed!(72)
    # The distinguishing case vs a Split branch: cases carry their OWN
    # ascertainment (×0.6), and deaths are observed off the ASCERTAINED expected
    # cases (a further ×0.5), so the death expectation is 100 × 0.6 × 0.5 = 30.
    chain = Chain((
        cases = Ascertainment(PoissonError(), FixedIntercept(log(0.6))),
        deaths = Ascertainment(PoissonError(), FixedIntercept(log(0.5)))))
    out = as_turing_model(
        chain, (cases = missing, deaths = missing), fill(100.0, 8))()

    @test all(≈(60.0), out.expected.cases)
    @test all(≈(30.0), out.expected.deaths)

    # A parallel Split of the same streams feeds BOTH the raw infections, so its
    # deaths expectation is 100 × 0.5 = 50, not 30 — the behaviours differ.
    split = Split((
        cases = Ascertainment(PoissonError(), FixedIntercept(log(0.6))),
        deaths = Ascertainment(PoissonError(), FixedIntercept(log(0.5)))))
    osplit = as_turing_model(
        split, (cases = missing, deaths = missing), fill(100.0, 8))()
    @test all(≈(50.0), osplit.expected.deaths)
end

@testitem "Chain shortens the series step by step down successive delays" begin
    using ComposableTuringIDModels, Distributions, Random
    Random.seed!(73)
    # Three streams, each adding a delay: 12 → 10 → 8 → 6 as the chain descends.
    chain = Chain((
        a = LatentDelay(PoissonError(), [0.5, 0.3, 0.2]),
        b = LatentDelay(PoissonError(), [0.5, 0.3, 0.2]),
        c = LatentDelay(PoissonError(), [0.5, 0.3, 0.2])))
    out = as_turing_model(
        chain, (a = missing, b = missing, c = missing), fill(50.0, 12))()
    @test length(out.expected.a) == 10
    @test length(out.expected.b) == 8
    @test length(out.expected.c) == 6
    @test all(≈(50.0), out.expected.c)
end

@testitem "Chain nests inside a Split via the shared contract" begin
    using ComposableTuringIDModels, Distributions, Random
    Random.seed!(74)
    # A Split branch whose stream is itself a Chain: the contracts compose, so the
    # nested per-stream NamedTuple appears under the branch name.
    split = Split((
        region = Chain((cases = PoissonError(),
        deaths = LatentDelay(PoissonError(), [0.5, 0.3, 0.2]))),))
    out = as_turing_model(split, (region = missing,), fill(20.0, 9))()
    @test keys(out.expected.region) == (:cases, :deaths)
    @test length(out.expected.region.cases) == 9
    @test length(out.expected.region.deaths) == 7
end

@testitem "Chain requires at least one stream" begin
    using ComposableTuringIDModels
    @test_throws AssertionError Chain(NamedTuple())
end

@testitem "Chain simulate-then-condition round-trips" begin
    using ComposableTuringIDModels, Random
    Random.seed!(75)
    chain = Chain((cases = PoissonError(), deaths = PoissonError()))
    sim = as_turing_model(
        chain, (cases = missing, deaths = missing), fill(15.0, 7))().y_t
    cond = as_turing_model(chain,
        (cases = sim.cases, deaths = sim.deaths), fill(15.0, 7))().y_t
    @test cond.cases == sim.cases
    @test cond.deaths == sim.deaths
end

@testitem "Chain composes with a renewal model and samples under NUTS" tags=[:sample] begin
    using ComposableTuringIDModels, Distributions, Random, Turing
    Random.seed!(76)
    data = IDData([0.2, 0.3, 0.5], exp)
    n = 25

    # Deaths observed downstream of the ascertained, delayed expected cases off a
    # shared renewal infection trajectory.
    obs = Chain((
        cases = LatentDelay(
            Ascertainment(NegativeBinomialError(), FixedIntercept(log(0.6))),
            truncated(Normal(3.0, 1.0), 0.0, Inf)),
        deaths = LatentDelay(
            Ascertainment(NegativeBinomialError(), FixedIntercept(log(0.05))),
            truncated(Normal(7.0, 2.0), 0.0, Inf))))
    model = IDModel(
        Renewal(data; rt = RandomWalk(), initialisation_prior = Normal()), obs)

    sim = as_turing_model(model, (cases = missing, deaths = missing), n)()
    @test keys(sim.generated_y_t) == (:cases, :deaths)
    ydata = (cases = sim.generated_y_t.cases, deaths = sim.generated_y_t.deaths)

    # The composed posterior conditions on the two-stream data and reproduces it.
    # `isequal` because the delayed streams carry a leading `missing` head.
    cond = as_turing_model(model, ydata, n)()
    @test isequal(cond.generated_y_t.cases, ydata.cases)

    # It samples under NUTS (default ForwardDiff): a few steps exercises the
    # gradient path on the composed renewal → chain posterior.
    chain = sample(as_turing_model(model, ydata, n), NUTS(), 5; progress = false)
    @test size(chain, 1) == 5
end
