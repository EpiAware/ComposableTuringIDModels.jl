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

@testitem "SequentialObservationModels: downstream consumes upstream .expected" begin
    using EpiAwarePrototype, Distributions, Random
    Random.seed!(63)

    # cases  : delayed Poisson on I_t = 100  (a full first stream on I_t)
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
