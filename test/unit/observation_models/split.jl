@testitem "Split composes several streams in parallel with name prefixes" begin
    using EpiAwarePrototype, Distributions, Random
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
    using EpiAwarePrototype, Random
    Random.seed!(81)
    split = Split((cases = PoissonError(), deaths = PoissonError()))
    Y = (cases = fill(5.0, 6), deaths = fill(50.0, 6))
    out = as_turing_model(split, (cases = missing, deaths = missing), Y)()
    @test out.expected.cases == fill(5.0, 6)
    @test out.expected.deaths == fill(50.0, 6)
end

@testitem "Split threads a stream downstream of another (sequential)" begin
    using EpiAwarePrototype, Distributions, Random
    Random.seed!(82)
    # cases: a delayed fraction of infections; deaths: a 0.1x delayed fraction of
    # the EXPECTED cases (downstream), not of infections directly.
    split = Split(
        (
            cases = LatentDelay(PoissonError(), [0.5, 0.3, 0.2]),
            deaths = LatentDelay(
                Ascertainment(PoissonError(), FixedIntercept(log(0.1))),
                [0.2, 0.3, 0.5]));
        sequential = true)

    out = as_turing_model(
        split, (cases = missing, deaths = missing), fill(100.0, 12))()

    # cases' expected is the post-delay mean of the length-12 input, shortened by
    # the 3-bin delay to length 10 (mean 100).
    @test length(out.expected.cases) == 10
    @test all(≈(100.0), out.expected.cases)
    # deaths is fed cases' expected (len 10), delayed again (→ 8) and scaled ×0.1.
    @test length(out.expected.deaths) == 8
    @test all(≈(10.0), out.expected.deaths)
end

@testitem "Split sequential differs from parallel on the first hop" begin
    using EpiAwarePrototype, Distributions, Random
    Random.seed!(83)
    # A length-shortening delay on the first stream: under `sequential` the second
    # stream sees the SHORTENED expected (a genuine cascade), under parallel it
    # sees the raw input (a fork off the same infections).
    streams = (cases = LatentDelay(PoissonError(), [0.5, 0.3, 0.2]),
        deaths = PoissonError())
    seq = Split(streams; sequential = true)
    par = Split(streams)

    Yin = fill(20.0, 12)
    oseq = as_turing_model(seq, (cases = missing, deaths = missing), Yin)()
    opar = as_turing_model(par, (cases = missing, deaths = missing), Yin)()

    # Sequential: deaths fed cases' shortened expected (length 10).
    @test length(oseq.expected.deaths) == 10
    # Parallel: deaths fed the raw infections (length 12).
    @test length(opar.expected.deaths) == 12
end

@testitem "Split wires an explicit source DAG and rejects bad sources" begin
    using EpiAwarePrototype, Random
    Random.seed!(84)
    # c is sourced from a (not the immediately-preceding b): a branch, not a chain.
    split = Split((a = PoissonError(), b = PoissonError(), c = PoissonError());
        sources = (c = :a,))
    out = as_turing_model(split,
        (a = missing, b = missing, c = missing), fill(30.0, 5))()
    @test out.expected.a == fill(30.0, 5)   # a: root
    @test out.expected.b == fill(30.0, 5)   # b: root (default)
    @test out.expected.c == out.expected.a  # c: sourced from a's expected

    # A source pointing at a LATER stream is rejected at construction.
    @test_throws Exception Split(
        (a = PoissonError(), b = PoissonError()); sources = (a = :b,))
    # A source naming an unknown stream is rejected.
    @test_throws Exception Split(
        (a = PoissonError(),); sources = (a = :ghost,))
    # `sequential` and `sources` cannot both be given.
    @test_throws Exception Split(
        (a = PoissonError(), b = PoissonError()); sequential = true,
        sources = (b = :a,))
end

@testitem "Split builds data-driven strata from a single template" begin
    using EpiAwarePrototype, Random
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
    using EpiAwarePrototype, Random
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

@testitem "Split simulate-then-condition round-trips" begin
    using EpiAwarePrototype, Random
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
    using EpiAwarePrototype, Distributions, Random, Turing
    Random.seed!(88)
    data = EpiData([0.2, 0.3, 0.5], exp)
    n = 25

    # Parallel cases + deaths off one shared renewal infection trajectory.
    obs = Split((
        cases = LatentDelay(NegativeBinomialError(),
            truncated(Normal(3.0, 1.0), 0.0, Inf)),
        deaths = LatentDelay(
            Ascertainment(NegativeBinomialError(), FixedIntercept(log(0.05))),
            truncated(Normal(7.0, 2.0), 0.0, Inf))))
    model = EpiAwareModel(
        Renewal(data; rt = RandomWalk(), initialisation_prior = Normal()), obs)

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
