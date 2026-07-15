@testitem "IDProblem assembles and simulates a composed model" begin
    using ComposableTuringIDModels, Distributions, Random
    Random.seed!(71)
    problem = IDProblem(
        infection = DirectInfections(; Z = RandomWalk(), initialisation = Normal()),
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
        infection = DirectInfections(; Z = RandomWalk(), initialisation = Normal()),
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

@testitem "IDMethod threads a Pathfinder pre-step into NUTS" tags=[:sample] begin
    using ComposableTuringIDModels, Distributions, Random
    Random.seed!(73)
    problem = IDProblem(
        infection = DirectInfections(; Z = RandomWalk(), initialisation = Normal()),
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
            DirectInfections(; Z = RandomWalk(), initialisation = Normal()),
            PoissonError()), missing, 10)
    obs = generated_observables(m, (; y_t = missing), rand(m))
    @test obs isa IDObservables
    @test obs.model === m
    # `returned` recovers the model's generated quantities for any draw it can
    # consume (a chain or, as here, a single draw); only solutions it cannot
    # consume (e.g. an optimiser result) leave `generated` as `missing`.
    @test obs.generated !== missing
end

@testitem "forecast rejects a non-positive horizon" begin
    using ComposableTuringIDModels, Distributions
    model = IDModel(
        DirectInfections(; Z = RandomWalk(), initialisation = Normal()),
        PoissonError())
    @test_throws ArgumentError forecast(model, fill(5, 10), :chain, 0)
end

@testitem "forecast extends a RandomWalk model over the horizon" tags=[:sample] begin
    using ComposableTuringIDModels, Distributions, Turing, Random
    Random.seed!(101)
    model = IDModel(
        DirectInfections(; Z = RandomWalk(), initialisation = Normal(1.0, 0.5)),
        PoissonError())
    T, h = 15, 6
    y = as_turing_model(model, fill(missing, T), T)().generated_y_t
    chain = sample(as_turing_model(model, y, T), Prior(), 40; progress = false)
    fc = forecast(model, y, chain, h)
    # The horizon points are predicted and integer-valued counts.
    @test size(fc, 1) == 40
    fut = vec(fc[@varname(y_t[T + 1])])
    @test length(fut) == 40
    @test all(x -> x isa Integer && x ≥ 0, fut)
    @test vec(fc[@varname(y_t[T + h])]) |> length == 40
    # The extended latent path continues the fitted trajectory rather than
    # overwriting it: in-sample Zₜ is unchanged and the path reaches T + h.
    gens_fit = vec(returned(as_turing_model(model, y, T), chain))
    gens_fc = vec(returned(as_turing_model(model, vcat(y, fill(missing, h)),
            T + h), fc))
    @test length(gens_fc[1].Z_t) == T + h
    @test all(d -> isapprox(gens_fit[d].Z_t, gens_fc[d].Z_t[1:T]; atol = 1e-8),
        1:10)
end

@testitem "forecast refuses a correlated (non-iid) latent stream" tags=[:sample] begin
    using ComposableTuringIDModels, Distributions, Turing, Random, LinearAlgebra
    Random.seed!(104)
    # A deliberately correlated latent: its stored stream is a smooth MvNormal,
    # so it does NOT factorise across time. Independent-tail extension would be
    # statistically wrong, and forecast must refuse rather than mis-forecast.
    struct CorrLatent <: AbstractPriorModel end
    @model function ComposableTuringIDModels.as_turing_model(::CorrLatent, n)
        Σ = [exp(-abs(i - j) / 5) for i in 1:n, j in 1:n] + 1.0e-6 * I
        z ~ MvNormal(zeros(n), Σ)
        return z
    end
    model = IDModel(
        DirectInfections(; Z = CorrLatent(), initialisation = Normal()),
        PoissonError())
    T, h = 15, 5
    y = as_turing_model(model, fill(missing, T), T)().generated_y_t
    chain = sample(as_turing_model(model, y, T), Prior(), 30; progress = false)
    @test_throws ErrorException forecast(model, y, chain, h)
end

@testitem "forecast works through an IDProblem and an AR latent" tags=[:sample] begin
    using ComposableTuringIDModels, Distributions, Turing, Random
    Random.seed!(102)
    problem = IDProblem(
        infection = DirectInfections(; Z = AR(), initialisation = Normal()),
        observation_model = PoissonError(),
        tspan = (1, 18))
    T, h = 18, 5
    y = as_turing_model(problem, (; y_t = missing))().generated_y_t
    chain = sample(as_turing_model(problem, (; y_t = y)), Prior(), 30;
        progress = false)
    fc = forecast(problem, y, chain, h)
    @test size(fc, 1) == 30
    @test length(vec(fc[@varname(y_t[T + h])])) == 30
end

@testitem "generated_observables leaves non-chain solutions missing" begin
    using ComposableTuringIDModels, Distributions, Random
    Random.seed!(76)
    m = as_turing_model(
        IDModel(
            DirectInfections(; Z = RandomWalk(), initialisation = Normal()),
            PoissonError()), missing, 10)
    # A solution `returned` cannot consume (here a bare marker) has no generated
    # quantities, so the field stays `missing`.
    obs = generated_observables(m, (; y_t = missing), :no_solution)
    @test obs.generated === missing
    # And the untyped fallback: a non-model `model` also yields `missing`.
    @test generated_observables(:not_a_model, (; y_t = missing), :no_solution).generated ===
          missing
end
