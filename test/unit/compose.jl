@testitem "composed model: prior simulation and generated quantities" begin
    using ComposableTuringIDModels, Distributions, Random
    Random.seed!(21)
    model = IDModel(
        DirectInfections(; Z = RandomWalk(), initialisation = Normal()),
        PoissonError())
    n = 20
    mdl = as_turing_model(model, missing, n)

    draw = rand(mdl)
    @test !isempty(keys(draw))

    gen = mdl()
    @test length(gen.generated_y_t) == n
    @test length(gen.I_t) == n
    @test length(gen.Z_t) == n
    @test all(>=(0), gen.I_t)
end

@testitem "composed model: fix and condition" begin
    using ComposableTuringIDModels, Distributions, Random
    using DynamicPPL: fix, condition
    Random.seed!(22)
    model = IDModel(
        DirectInfections(; Z = RandomWalk(), initialisation = Normal()),
        PoissonError())
    n = 20
    mdl = as_turing_model(model, missing, n)

    # fix removes a parameter from the sampled set.
    fixed = fix(mdl, (init_incidence = 0.0,))
    fixed_names = string.(collect(keys(rand(fixed))))
    @test !("init_incidence" in fixed_names)

    # condition and the | syntax produce equivalent conditioned models.
    c1 = condition(mdl, (std = 0.1,))
    c2 = mdl | (std = 0.1,)
    @test typeof(c1) == typeof(c2)
end

@testitem "composed model: short NUTS sample runs" tags=[:sample] begin
    using ComposableTuringIDModels, Distributions, Turing, Random
    Random.seed!(23)
    model = IDModel(
        DirectInfections(; Z = RandomWalk(), initialisation = Normal()),
        NegativeBinomialError())
    n = 20
    y = as_turing_model(model, missing, n)()
    cond_model = as_turing_model(model, y.generated_y_t, n)
    chn = sample(cond_model, NUTS(), 50; progress = false)
    @test chn !== nothing
end

# --- grouped/panel IDModel: the shared constructor from issue #180 ---------

@testitem "IDModel(inf, group_effect, obs) builds a grouped panel" begin
    using ComposableTuringIDModels, Distributions
    inf = DirectInfections(; Z = RandomWalk(), initialisation = Normal())
    h = Hierarchy(; mean = Normal(0.0, 0.5), across = IID(Normal(0.0, 0.5)))
    model = IDModel(inf, h, PoissonError())
    @test model isa IDModel
    @test model isa AbstractComposableModel
    @test model.infection_model isa GroupedInfections
    @test model.infection_model.infection_model isa DirectInfections
    @test model.infection_model.group_effect isa Hierarchy
    @test model.observation_model isa Split
    # Lifting a plain IDModel to a panel keeps its infection/observation parts.
    idmodel = IDModel(inf, PoissonError())
    lifted = IDModel(idmodel, h)
    @test lifted.infection_model.infection_model === idmodel.infection_model
    @test lifted.observation_model.streams === idmodel.observation_model
    # A bare Distribution is a valid group effect (independent per-group levels).
    @test IDModel(inf, Normal(0.0, 0.5), PoissonError()) isa IDModel
end

@testitem "grouped IDModel simulates a panel with the grouping dim from the data" begin
    using ComposableTuringIDModels, Distributions, Random
    Random.seed!(451)
    model = IDModel(
        DirectInfections(; Z = RandomWalk(), initialisation = Normal(log(50.0), 0.2)),
        Hierarchy(; mean = Normal(0.0, 0.5), across = IID(Normal(0.0, 0.5))),
        PoissonError())
    for (n_time, n_groups) in ((24, 4), (12, 6))
        Ymiss = Matrix{Union{Missing, Float64}}(missing, n_groups, n_time)
        sim = as_turing_model(model, Ymiss)()
        @test size(sim.I_t) == (n_groups, n_time)
        @test length(sim.Z_t.Z_t) == n_time
        @test length(sim.Z_t.group_levels) == n_groups
        @test length(sim.generated_y_t) == n_groups
        @test all(length(sim.generated_y_t[g]) == n_time for g in 1:n_groups)
        @test all(>=(0), reduce(vcat, collect(sim.generated_y_t)))
    end
end

@testitem "grouped IDModel: variable names are namespaced per group" begin
    using ComposableTuringIDModels, Distributions
    model = IDModel(
        DirectInfections(; Z = RandomWalk(), initialisation = Normal()),
        Hierarchy(; across = IID(Normal(0.0, 0.5))),
        PoissonError())
    Ymiss = Matrix{Union{Missing, Float64}}(missing, 3, 8)
    mdl = as_turing_model(model, Ymiss)
    names = string.(collect(keys(rand(mdl))))
    @test any(startswith(n, "group1.") for n in names)
    @test any(startswith(n, "group2.") for n in names)
    @test any(startswith(n, "group_levels.") for n in names)
end

@testitem "grouped IDModel recovers per-group levels under NUTS" tags=[:sample] begin
    using ComposableTuringIDModels, Distributions, Turing, Random, Statistics
    using Turing: returned
    Random.seed!(388)
    model = IDModel(
        DirectInfections(; Z = RandomWalk(), initialisation = Normal(log(50.0), 0.2)),
        Hierarchy(; mean = Normal(0.0, 0.5), across = IID(Normal(0.0, 0.5))),
        PoissonError())

    n_time, n_groups = 20, 8
    Ymiss = Matrix{Union{Missing, Float64}}(missing, n_groups, n_time)
    sim = as_turing_model(model, Ymiss)()
    Ydata = Float64.(reduce(vcat, [permutedims(sim.generated_y_t[g])
                                   for g in 1:n_groups]))

    posterior = as_turing_model(model, Ydata)
    chain = sample(posterior, NUTS(0.8; adtype = Turing.AutoForwardDiff()), 80;
        progress = false)
    @test size(chain, 1) == 80
    draws = reduce(hcat, [g.Z_t.group_levels for g in vec(returned(posterior, chain))])
    @test size(draws) == (n_groups, 80)
    post_mean = vec(mean(draws; dims = 2))
    @test cor(sim.Z_t.group_levels, post_mean) > 0.5
end
