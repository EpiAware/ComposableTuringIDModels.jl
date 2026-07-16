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
    # Generated quantities are now recovered from the chain (previously always
    # `missing`); `returned` yields the model's `(; generated_y_t, I_t, Z_t)`
    # per draw.
    @test res.generated !== missing
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
    # `returned` recovers the model's generated quantities for any draw it can
    # consume (a chain or, as here, a single draw); only solutions it cannot
    # consume (e.g. an optimiser result) leave `generated` as `missing`.
    @test obs.generated !== missing
end

@testitem "generated_observables leaves non-chain solutions missing" begin
    using ComposableTuringIDModels, Distributions, Random
    Random.seed!(76)
    m = as_turing_model(
        IDModel(
            DirectInfections(; Z = RandomWalk(), initialisation_prior = Normal()),
            PoissonError()), missing, 10)
    # A solution `returned` cannot consume (here a bare marker) has no generated
    # quantities, so the field stays `missing`.
    obs = generated_observables(m, (; y_t = missing), :no_solution)
    @test obs.generated === missing
    # And the untyped fallback: a non-model `model` also yields `missing`.
    @test generated_observables(:not_a_model, (; y_t = missing), :no_solution).generated ===
          missing
end
