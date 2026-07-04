@testitem "EpiProblem assembles and simulates a composed model" begin
    using EpiAwarePrototype, Distributions, Random
    Random.seed!(71)
    problem = EpiProblem(
        epi_model = DirectInfections(; Z = RandomWalk(), initialisation_prior = Normal()),
        observation_model = PoissonError(),
        tspan = (1, 20))
    m = as_turing_model(problem, (; y_t = missing))
    sim = m()
    @test length(sim.generated_y_t) == 20
    @test length(sim.Z_t) == 20
end

@testitem "apply_method runs a NUTSampler over an EpiProblem" tags=[:sample] begin
    using EpiAwarePrototype, Distributions, Random
    Random.seed!(72)
    problem = EpiProblem(
        epi_model = DirectInfections(; Z = RandomWalk(), initialisation_prior = Normal()),
        observation_model = PoissonError(),
        tspan = (1, 20))
    ydata = as_turing_model(problem, (; y_t = missing))().generated_y_t
    res = apply_method(problem, NUTSampler(; ndraws = 40, nchains = 1), (; y_t = ydata))
    @test res isa EpiAwareObservables
    @test res.samples !== nothing
    # Generated quantities are now recovered from the chain (previously always
    # `missing`); `returned` yields the model's `(; generated_y_t, I_t, Z_t)`
    # per draw.
    @test res.generated !== missing
end

@testitem "EpiMethod threads a Pathfinder pre-step into NUTS" tags=[:sample] begin
    using EpiAwarePrototype, Distributions, Random
    Random.seed!(73)
    problem = EpiProblem(
        epi_model = DirectInfections(; Z = RandomWalk(), initialisation_prior = Normal()),
        observation_model = PoissonError(),
        tspan = (1, 20))
    ydata = as_turing_model(problem, (; y_t = missing))().generated_y_t
    method = EpiMethod(
        pre_sampler_steps = [ManyPathfinder(; ndraws = 10, nruns = 2)],
        sampler = NUTSampler(; ndraws = 40, nchains = 1))
    res = apply_method(problem, method, (; y_t = ydata))
    @test res isa EpiAwareObservables
    @test res.samples !== nothing
end

@testitem "spread_draws produces tidy draw/chain/iteration columns" tags=[:sample] begin
    using EpiAwarePrototype, Distributions, Turing, MCMCChains, Random
    Random.seed!(74)
    @model f() = (x ~ Normal())
    chn = MCMCChains.Chains(sample(f(), NUTS(), 30; progress = false))
    df = spread_draws(chn)
    @test all(c -> c in names(df), ["draw", "chain", "iteration"])
    @test size(df, 1) == 30
end

@testitem "generated_observables wraps model, data, and solution" begin
    using EpiAwarePrototype, Distributions, Random
    Random.seed!(75)
    m = as_turing_model(
        EpiAwareModel(
            DirectInfections(; Z = RandomWalk(), initialisation_prior = Normal()),
            PoissonError()), missing, 10)
    obs = generated_observables(m, (; y_t = missing), rand(m))
    @test obs isa EpiAwareObservables
    @test obs.model === m
    # `returned` recovers the model's generated quantities for any draw it can
    # consume (a chain or, as here, a single draw); only solutions it cannot
    # consume (e.g. an optimiser result) leave `generated` as `missing`.
    @test obs.generated !== missing
end

@testitem "generated_observables leaves non-chain solutions missing" begin
    using EpiAwarePrototype, Distributions, Random
    Random.seed!(76)
    m = as_turing_model(
        EpiAwareModel(
            DirectInfections(; Z = RandomWalk(), initialisation_prior = Normal()),
            PoissonError()), missing, 10)
    # A solution `returned` cannot consume (here a bare marker) has no generated
    # quantities, so the field stays `missing`.
    obs = generated_observables(m, (; y_t = missing), :no_solution)
    @test obs.generated === missing
end
