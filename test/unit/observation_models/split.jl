@testitem "Split composes several streams in parallel with name prefixes" begin
    using ComposableTuringIDModels, Distributions, Random
    Random.seed!(80)
    split = Split((cases = PoissonError(), deaths = NegativeBinomialError()))
    @test split isa AbstractObservationModel
    @test implements_observation_interface(split)

    yt = (cases = missing, deaths = missing)
    sm = as_turing_model(split, yt, fill(10.0, 10))
    # Each stream's variables are prefixed by its name.
    names = string.(collect(keys(rand(sm))))
    @test any(startswith("cases."), names)
    @test any(startswith("deaths."), names)

    # The uniform return contract: per-stream `y_t` / `expected` NamedTuples.
    out = sm()
    @test keys(out) == (:y_t, :expected)
    @test keys(out.y_t) == (:cases, :deaths)
    @test length(out.y_t.cases) == 10
    @test length(out.y_t.deaths) == 10
    # In parallel every stream sees the SAME expected input.
    @test out.expected.cases == fill(10.0, 10)
    @test out.expected.deaths == fill(10.0, 10)
end

@testitem "Split accepts a per-stream NamedTuple of expected series" begin
    using ComposableTuringIDModels, Random
    Random.seed!(81)
    split = Split((cases = PoissonError(), deaths = PoissonError()))
    Y = (cases = fill(5.0, 6), deaths = fill(50.0, 6))
    out = as_turing_model(split, (cases = missing, deaths = missing), Y)()
    @test out.expected.cases == fill(5.0, 6)
    @test out.expected.deaths == fill(50.0, 6)
end

@testitem "Split cascades a stream downstream by nesting after a shared delay" begin
    using ComposableTuringIDModels, Distributions, Random
    Random.seed!(82)
    # Share the case delay, then split: cases apply their error to the delayed
    # expectation and deaths sit DOWNSTREAM, delayed again and scaled ×0.1.
    cascade = LatentDelay(
        Split((
            cases = PoissonError(),
            deaths = LatentDelay(
                Ascertainment(PoissonError(), FixedIntercept(log(0.1))),
                [0.2, 0.3, 0.5]))),
        [0.5, 0.3, 0.2])

    out = as_turing_model(
        cascade, (cases = missing, deaths = missing), fill(100.0, 12))()

    # cases see the shared delayed expectation (12 → 10, mean 100).
    @test length(out.expected.cases) == 10
    @test all(≈(100.0), out.expected.cases)
    # deaths are delayed again (10 → 8) and scaled ×0.1.
    @test length(out.expected.deaths) == 8
    @test all(≈(10.0), out.expected.deaths)
end

@testitem "Placing Split high vs low chooses parallel vs cascade" begin
    using ComposableTuringIDModels, Distributions, Random
    Random.seed!(83)
    # High (on infections): deaths fork off the raw input alongside a delayed cases
    # stream — parallel. Low (under the case delay): deaths see the shortened
    # expectation — cascade. Only the placement differs.
    high = Split((cases = LatentDelay(PoissonError(), [0.5, 0.3, 0.2]),
        deaths = PoissonError()))
    low = LatentDelay(Split((cases = PoissonError(), deaths = PoissonError())),
        [0.5, 0.3, 0.2])

    Yin = fill(20.0, 12)
    ohigh = as_turing_model(high, (cases = missing, deaths = missing), Yin)()
    olow = as_turing_model(low, (cases = missing, deaths = missing), Yin)()

    # High: deaths fed the raw infections (length 12).
    @test length(ohigh.expected.deaths) == 12
    # Low: deaths fed the shared shortened expected (length 10).
    @test length(olow.expected.deaths) == 10
end

@testitem "Nesting Splits builds an arbitrary branch DAG" begin
    using ComposableTuringIDModels, Random
    Random.seed!(84)
    # b and c both branch off a's delayed expectation by nesting a Split under a's
    # delay: a general DAG built purely by placement, no source wiring.
    split = Split((
        a = LatentDelay(
        Split((leaf = PoissonError(), b = PoissonError(), c = PoissonError())),
        [0.5, 0.3, 0.2]),))
    out = as_turing_model(split, (a = missing,), fill(30.0, 8))()
    # a's streams share a's delayed expectation (8 → 6, mean 30).
    @test length(out.expected.a.b) == 6
    @test out.expected.a.b == out.expected.a.c
    @test all(≈(30.0), out.expected.a.leaf)
end

@testitem "Split builds data-driven strata from a single template" begin
    using ComposableTuringIDModels, Random
    Random.seed!(85)
    # One template, replicated per data stream: the number and names of streams
    # come from the `y_t` data, not the struct.
    template = Split(PoissonError())
    yt = (young = missing, old = missing)
    # A 2-strata × 6-time expected matrix: one stream per row (one-to-one).
    M = [10.0 10 10 10 10 10; 30.0 30 30 30 30 30]
    sm = as_turing_model(template, yt, M)
    out = sm()
    @test keys(out.y_t) == (:young, :old)
    @test out.expected.young == fill(10.0, 6)
    @test out.expected.old == fill(30.0, 6)
    names = string.(collect(keys(rand(sm))))
    @test any(startswith("young."), names)
    @test any(startswith("old."), names)

    # A three-stream template also works — the count is purely data-driven.
    out3 = as_turing_model(template,
        (a = missing, b = missing, c = missing),
        [1.0 1; 2.0 2; 3.0 3])()
    @test keys(out3.y_t) == (:a, :b, :c)
end

@testitem "StrataMap maps infection strata to observation streams (1:1, m:1, m:m)" begin
    using ComposableTuringIDModels, Random
    Random.seed!(86)
    template = Split(PoissonError())
    # 2 infection strata over 6 days.
    M = [10.0 10 10 10 10 10; 30.0 30 30 30 30 30]

    # one-to-one (identity map).
    id = StrataMap(M, [1.0 0.0; 0.0 1.0])
    o1 = as_turing_model(template, (a = missing, b = missing), id)()
    @test o1.expected.a == fill(10.0, 6)
    @test o1.expected.b == fill(30.0, 6)

    # many-to-one and many-to-many in one 3 × 2 map: stratum 1, stratum 2, and the
    # aggregate of both.
    W = [1.0 0.0; 0.0 1.0; 1.0 1.0]
    mm = StrataMap(M, W)
    o2 = as_turing_model(template,
        (a = missing, b = missing, total = missing), mm)()
    @test o2.expected.total == fill(40.0, 6)   # many-to-one aggregate

    # A map whose columns disagree with the strata rows is rejected.
    @test_throws Exception StrataMap(M, [1.0 0.0 0.0])
end

@testitem "Split(template, map) projects infection strata inside a composed model" begin
    using ComposableTuringIDModels, Distributions, Random
    Random.seed!(89)
    # The weight map lives in the observation model, so the strata come from the
    # composed infection process at run time (not a hand-built matrix).
    W = reshape([0.7, 0.3, 1.0], 3, 1)          # young, old, and their total
    weighted = Split(PoissonError(), W)
    model = IDModel(
        DirectInfections(; Z = RandomWalk(), initialisation = Normal()),
        weighted)
    out = as_turing_model(
        model, (young = missing, old = missing, total = missing), 12)()
    @test keys(out.generated_y_t) == (:young, :old, :total)
    e = out.expected_y_t
    # young/old are weighted fractions of the one infection stratum; total sums.
    @test e.young ≈ 0.7 .* out.I_t
    @test e.old ≈ 0.3 .* out.I_t
    @test e.total ≈ e.young .+ e.old
end

@testitem "Split simulate-then-condition round-trips" begin
    using ComposableTuringIDModels, Random
    Random.seed!(87)
    split = Split((cases = PoissonError(), deaths = PoissonError()))
    sim = as_turing_model(
        split, (cases = missing, deaths = missing), fill(15.0, 7))().y_t
    cond = as_turing_model(split,
        (cases = sim.cases, deaths = sim.deaths), fill(15.0, 7))().y_t
    @test cond.cases == sim.cases
    @test cond.deaths == sim.deaths
end

@testitem "Split composes with a renewal model and samples under NUTS" tags=[:sample] begin
    using ComposableTuringIDModels, Distributions, Random, Turing
    Random.seed!(88)
    gen_int = [0.2, 0.3, 0.5]
    n = 25

    # Parallel cases + deaths off one shared renewal infection trajectory.
    obs = Split((
        cases = LatentDelay(NegativeBinomialError(),
            truncated(Normal(3.0, 1.0), 0.0, Inf)),
        deaths = LatentDelay(
            Ascertainment(NegativeBinomialError(), FixedIntercept(log(0.05))),
            truncated(Normal(7.0, 2.0), 0.0, Inf))))
    model = IDModel(
        Renewal(gen_int; rt = RandomWalk(), initialisation = Normal()), obs)

    sim = as_turing_model(model, missing, n)()
    @test keys(sim.generated_y_t) == (:cases, :deaths)
    ydata = (cases = sim.generated_y_t.cases, deaths = sim.generated_y_t.deaths)

    # The composed posterior conditions on the two-stream data and reproduces it.
    # `isequal` because the delayed streams carry a leading `missing` head.
    cond = as_turing_model(model, ydata, n)()
    @test isequal(cond.generated_y_t.cases, ydata.cases)

    # It samples under NUTS (default ForwardDiff): a few steps is enough to
    # exercise the gradient path on the composed renewal → split posterior.
    chain = sample(as_turing_model(model, ydata, n), NUTS(), 5; progress = false)
    @test size(chain, 1) == 5
end

@testitem "Split cascade composes with a renewal model and samples under NUTS" tags=[:sample] begin
    using ComposableTuringIDModels, Distributions, Random, Turing
    Random.seed!(90)
    gen_int = [0.2, 0.3, 0.5]
    n = 25

    # Cascade: deaths downstream of the delayed expected cases (placement, no flag).
    cascade = LatentDelay(
        Split((
            cases = NegativeBinomialError(),
            deaths = LatentDelay(
                Ascertainment(NegativeBinomialError(), FixedIntercept(log(0.02))),
                truncated(Normal(5.0, 1.5), 0.0, Inf)))),
        truncated(Normal(3.0, 1.0), 0.0, Inf))
    model = IDModel(
        Renewal(gen_int; rt = RandomWalk(), initialisation = Normal()), cascade)

    sim = as_turing_model(model, (cases = missing, deaths = missing), n)()
    @test keys(sim.generated_y_t) == (:cases, :deaths)
    ydata = (cases = sim.generated_y_t.cases, deaths = sim.generated_y_t.deaths)

    chain = sample(as_turing_model(model, ydata, n), NUTS(), 5; progress = false)
    @test size(chain, 1) == 5
end

@testitem "Split strata composes with a renewal model and samples under NUTS" tags=[:sample] begin
    using ComposableTuringIDModels, Distributions, Random, Turing
    Random.seed!(91)
    gen_int = [0.2, 0.3, 0.5]
    n = 25

    # Strata: one full stream per band off the shared renewal infections.
    strata = Split((
        young = LatentDelay(
            Ascertainment(NegativeBinomialError(), FixedIntercept(log(0.7))),
            truncated(Normal(2.0, 1.0), 0.0, Inf)),
        old = LatentDelay(
            Ascertainment(NegativeBinomialError(), FixedIntercept(log(0.4))),
            truncated(Normal(3.0, 1.0), 0.0, Inf))))
    model = IDModel(
        Renewal(gen_int; rt = RandomWalk(), initialisation = Normal()), strata)

    sim = as_turing_model(model, (young = missing, old = missing), n)()
    @test keys(sim.generated_y_t) == (:young, :old)
    ydata = (young = sim.generated_y_t.young, old = sim.generated_y_t.old)

    chain = sample(as_turing_model(model, ydata, n), NUTS(), 5; progress = false)
    @test size(chain, 1) == 5
end
