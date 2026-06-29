@testitem "SequentialObservationModels: ordered construction and name-prefixing" begin
    using EpiAwarePrototype, Distributions, Random
    Random.seed!(61)

    # Order is significant: the NamedTuple key order is preserved as the cascade
    # order, and each stream is prefix-wrapped by its name (like the parallel
    # stack).
    seq = SequentialObservationModels((
        cases = LatentDelay(PoissonError(), [0.4, 0.3, 0.2, 0.1]),
        deaths = (inner -> Ascertainment(inner, FixedIntercept(log(0.1)))) =>
            PoissonError()))
    @test seq isa AbstractObservationModel
    @test seq.model_names == ["cases", "deaths"]
    # Each stream is wrapped in a PrefixObservationModel keyed by its name.
    @test all(m -> m isa PrefixObservationModel, seq.models)
    @test [m.prefix for m in seq.models] == ["cases", "deaths"]

    yt = (cases = missing, deaths = missing)
    mdl = as_turing_model(seq, yt, fill(100.0, 12))
    names = string.(collect(keys(rand(mdl))))
    @test any(startswith("cases."), names)
    @test any(startswith("deaths."), names)

    # The reversed order is a different (still valid) construct — order matters.
    rev = SequentialObservationModels((
        deaths = PoissonError(),
        cases = LatentDelay(PoissonError(), [0.4, 0.3, 0.2, 0.1])))
    @test rev.model_names == ["deaths", "cases"]
end

@testitem "SequentialObservationModels: interface conformance" begin
    using EpiAwarePrototype, Distributions
    seq = SequentialObservationModels((
        cases = PoissonError(),
        deaths = (inner -> Ascertainment(inner, FixedIntercept(log(0.1)))) =>
            PoissonError()))
    @test implements_observation_interface(seq;
        y_t = (cases = missing, deaths = missing), Y_t = fill(50.0, 10))
end

@testitem "SequentialObservationModels: recorder threads the expected output" begin
    using EpiAwarePrototype, Distributions, Random
    Random.seed!(62)

    # The internal recorder returns (; obs, expected) with `expected` the
    # error-leaf input it was handed — the threaded (pre-error) series, NOT the
    # sampled output.
    rec = EpiAwarePrototype._SeqExpectedRecorder(PoissonError())
    out = as_turing_model(rec, missing, fill(42.0, 5))()
    @test out.expected == fill(42.0, 5)
    @test length(out.obs) == 5
    @test out.obs != out.expected   # sampled counts differ from the expected mean
end

@testitem "SequentialObservationModels: cascade threads genuinely from stream 1" begin
    using EpiAwarePrototype, Distributions, Random
    Random.seed!(63)

    # The FIRST hop must thread: stream 2 is fed stream 1's expected (post-its-own
    # transform) output, not the raw `I_t`. With a length-shortening delay on
    # stream 1, that shows up as a length difference — if the raw `I_t` were
    # threaded (as the parallel stack does), stream 2 would see the full length.
    pmf = [0.4, 0.3, 0.2, 0.1]
    ramp = collect(1.0:12.0)
    seq = SequentialObservationModels((
        cases = LatentDelay(PoissonError(), pmf),     # shortens by length(pmf)-1
        sink = PoissonError()))                       # pure leaf on the threaded series
    out = as_turing_model(seq, (cases = missing, sink = missing), ramp)()
    @test length(out) == 2
    # Stream 2 sees stream 1's delay-SHORTENED expected (12 - 3 = 9), proving the
    # first hop threads rather than forking off the raw 12-length `I_t`.
    @test length(out[2]) == length(ramp) - (length(pmf) - 1)
    @test length(out[2]) == 9
end

@testitem "SequentialObservationModels: downstream consumes upstream .expected" begin
    using EpiAwarePrototype, Distributions, Random
    Random.seed!(631)

    # cases  : delayed Poisson on I_t = 100
    # deaths : ascertained (×0.1) Poisson leaf -> its error-leaf input ≈ 10
    # tert   : a pure Poisson leaf -> the stack feeds it the DEATHS stream's
    #          expected (≈ 10), not the original 100. So its sampled counts must
    #          cluster near 10, demonstrating the cascade threads `.expected`.
    seq = SequentialObservationModels((
        cases = LatentDelay(PoissonError(), [0.4, 0.3, 0.2, 0.1]),
        deaths = (inner -> Ascertainment(inner, FixedIntercept(log(0.1)))) =>
            PoissonError(),
        tert = PoissonError()))
    yt = (cases = missing, deaths = missing, tert = missing)
    out = as_turing_model(seq, yt, fill(100.0, 12))()
    @test length(out) == 3

    tert = collect(skipmissing(out[3]))
    @test !isempty(tert)
    # Fed ≈10 (the deaths post-ascertainment mean), the tertiary counts are far
    # below 100 — if the raw input had been threaded they would sit near 100.
    @test sum(tert) / length(tert) < 40
end

@testitem "SequentialObservationModels: a bare nested stream is peeled at its leaf" begin
    using EpiAwarePrototype, Distributions, Random
    Random.seed!(632)

    # A stream written as a plain nested observation model (no explicit
    # `transform_chain => error_leaf`) is peeled so the recorder lands at the error
    # leaf: the threaded expected series is then the post-transform leaf input.
    # Here stream 1 = LatentDelay(Ascertainment(×0.1, Poisson), pmf): its expected
    # output is BOTH ascertained (≈10) AND delay-shortened, so the bare-leaf second
    # stream sees ≈10 over the shortened length (not 100 over 12).
    pmf = [0.4, 0.3, 0.2, 0.1]
    nested = LatentDelay(
        Ascertainment(PoissonError(), FixedIntercept(log(0.1))), pmf)
    seq = SequentialObservationModels((cases = nested, sink = PoissonError()))
    out = as_turing_model(seq, (cases = missing, sink = missing), fill(100.0, 12))()
    @test length(out[2]) == 12 - (length(pmf) - 1)        # shortened by the delay
    sink = collect(skipmissing(out[2]))
    @test !isempty(sink)
    @test sum(sink) / length(sink) < 40                   # ascertained ≈10, not 100

    # An un-peelable wrapper (no registered transform-chain rule) errors clearly,
    # steering the user to the explicit `transform_chain => error_leaf` form.
    @test_throws Exception SequentialObservationModels((
        a = RecordExpectedObs(PoissonError()), b = PoissonError()))
end

@testitem "SequentialObservationModels: simulate then condition" begin
    using EpiAwarePrototype, Distributions, Random
    Random.seed!(64)

    seq = SequentialObservationModels((
        cases = PoissonError(),
        deaths = (inner -> Ascertainment(inner, FixedIntercept(log(0.2)))) =>
            PoissonError()))
    yt = (cases = missing, deaths = missing)
    Y = fill(80.0, 14)

    sim = as_turing_model(seq, yt, Y)()
    @test length(sim) == 2
    @test length(sim[1]) == 14
    @test length(sim[2]) == 14

    # Conditioning on the simulated data rebuilds and evaluates to the same draw.
    data = (cases = sim[1], deaths = sim[2])
    cond = as_turing_model(seq, data, Y)
    @test cond() == sim
end

@testitem "SequentialObservationModels: NamedTuple expected seeds the first stream" begin
    using EpiAwarePrototype, Distributions, Random
    Random.seed!(65)

    # A per-stream NamedTuple `Y_t` is accepted; only its first entry seeds the
    # cascade (later streams are fed by the threading), mirroring the data
    # contract of the parallel stack.
    seq = SequentialObservationModels((
        cases = PoissonError(),
        deaths = PoissonError()))
    yt = (cases = missing, deaths = missing)
    Yt = (cases = fill(30.0, 8), deaths = fill(999.0, 8))   # deaths entry ignored
    out = as_turing_model(seq, yt, Yt)()
    @test length(out) == 2
    @test length(out[1]) == 8

    # The whole key set is validated (names AND order), so a misnamed or
    # misordered `Y_t` errors loudly rather than silently using the wrong column.
    @test_throws Exception as_turing_model(
        seq, yt, (cases = fill(30.0, 8), WRONG = fill(1.0, 8)))()
    @test_throws Exception as_turing_model(
        seq, yt, (deaths = fill(1.0, 8), cases = fill(30.0, 8)))()
end

@testitem "SequentialObservationModels: end-to-end renewal compose, simulate, NUTS" tags=[:sample] begin
    using EpiAwarePrototype, Distributions, Turing, Random
    using ADTypes: AutoForwardDiff
    Random.seed!(66)

    data = EpiData([0.2, 0.3, 0.5], exp)
    n = 12
    seq = SequentialObservationModels((
        cases = LatentDelay(PoissonError(), [0.4, 0.3, 0.2, 0.1]),
        deaths = (inner -> Ascertainment(inner, FixedIntercept(log(0.1)))) =>
            PoissonError()))
    model = EpiAwareModel(
        Renewal(data; rt = RandomWalk(), initialisation_prior = Normal()), seq)

    # Simulate the joint cases/deaths cascade off the renewal infections.
    sim = as_turing_model(model, (cases = missing, deaths = missing), n)()
    y = sim.generated_y_t
    @test length(y) == 2

    # Condition and run a short NUTS sample (ForwardDiff). The Mooncake gradient
    # for this composed model is exercised in the AD fixture registry
    # (`test/ADFixtures`), which the per-backend AD CI runs under Mooncake.
    cond = as_turing_model(model, (cases = y[1], deaths = y[2]), n)
    chn = sample(cond, NUTS(; adtype = AutoForwardDiff()), 20; progress = false)
    @test chn !== nothing
end
