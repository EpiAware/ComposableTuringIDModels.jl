@testitem "IDProblem assembles and simulates a composed model" begin
    using ComposableTuringIDModels, Distributions, Random
    Random.seed!(71)
    problem = IDProblem(
        infection = DirectInfections(; Z = RandomWalk(), initialisation_prior = Normal()),
        observation_model = PoissonError(),
        tspan = (1, 20))
    m = as_turing_model(problem, (; y_t = missing))
    sim = m()
    @test length(sim.generated_y_t) == 20
    @test length(sim.Z_t) == 20
end

@testitem "apply_method runs a NUTSampler over an IDProblem" tags=[:sample] begin
    using ComposableTuringIDModels, Distributions, Random
    Random.seed!(72)
    problem = IDProblem(
        infection = DirectInfections(; Z = RandomWalk(), initialisation_prior = Normal()),
        observation_model = PoissonError(),
        tspan = (1, 20))
    ydata = as_turing_model(problem, (; y_t = missing))().generated_y_t
    res = apply_method(problem, NUTSampler(; ndraws = 40, nchains = 1), (; y_t = ydata))
    @test res isa IDObservables
    @test res.samples !== nothing
end

@testitem "IDMethod threads a Pathfinder pre-step into NUTS" tags=[:sample] begin
    using ComposableTuringIDModels, Distributions, Random
    Random.seed!(73)
    problem = IDProblem(
        infection = DirectInfections(; Z = RandomWalk(), initialisation_prior = Normal()),
        observation_model = PoissonError(),
        tspan = (1, 20))
    ydata = as_turing_model(problem, (; y_t = missing))().generated_y_t
    method = IDMethod(
        pre_sampler_steps = [ManyPathfinder(; ndraws = 10, nruns = 2)],
        sampler = NUTSampler(; ndraws = 40, nchains = 1))
    res = apply_method(problem, method, (; y_t = ydata))
    @test res isa IDObservables
    @test res.samples !== nothing
end

@testitem "spread_draws produces tidy draw/chain/iteration columns" tags=[:sample] begin
    using ComposableTuringIDModels, Distributions, Turing, MCMCChains, Random
    Random.seed!(74)
    @model f() = (x ~ Normal())
    chn = MCMCChains.Chains(sample(f(), NUTS(), 30; progress = false))
    df = spread_draws(chn)
    @test all(c -> c in names(df), ["draw", "chain", "iteration"])
    @test size(df, 1) == 30
end

@testitem "generated_observables wraps model, data, and solution" begin
    using ComposableTuringIDModels, Distributions, Random
    Random.seed!(75)
    m = as_turing_model(
        IDModel(
            DirectInfections(; Z = RandomWalk(), initialisation_prior = Normal()),
            PoissonError()), missing, 10)
    obs = generated_observables(m, (; y_t = missing), rand(m))
    @test obs isa IDObservables
    @test obs.model === m
end
