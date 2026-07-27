# Tests for the `GroupedIDModel` composed panel model: one shared infection
# process observed by several groups, each reporting it at its own partially
# pooled level. The grouping dimension is read from the data (columns of `Y`),
# never stored on the struct, and the whole panel composes from components with no
# hand-written `@model`.

@testitem "GroupedIDModel constructs and is an AbstractComposableModel" begin
    using ComposableTuringIDModels, Distributions
    inf = DirectInfections(; Z = RandomWalk(), initialisation = Normal())
    h = Hierarchy(; mean = Normal(0.0, 0.5), across = IID(Normal(0.0, 0.5)))
    model = GroupedIDModel(inf, h, PoissonError())
    @test model isa GroupedIDModel
    @test model isa AbstractComposableModel
    @test model.infection_model isa DirectInfections
    @test model.group_effect isa Hierarchy
    @test model.observation_model isa PoissonError
    @test model.combiner isa Function
    # The grouping dimension is NOT a field: it is read from the data.
    @test !(:n_groups in fieldnames(typeof(model)))
    # Lifting an IDModel to a panel keeps its infection and observation parts.
    idmodel = IDModel(inf, PoissonError())
    lifted = GroupedIDModel(idmodel, h)
    @test lifted.infection_model === idmodel.infection_model
    @test lifted.observation_model === idmodel.observation_model
    # A bare Distribution is a valid group effect (independent per-group levels).
    @test GroupedIDModel(inf, Normal(0.0, 0.5), PoissonError()) isa GroupedIDModel
end

@testitem "GroupedIDModel simulates a panel with the grouping dim from the data" begin
    using ComposableTuringIDModels, Distributions, Random
    Random.seed!(451)
    model = GroupedIDModel(
        DirectInfections(; Z = RandomWalk(), initialisation = Normal(log(50.0), 0.2)),
        Hierarchy(; mean = Normal(0.0, 0.5), across = IID(Normal(0.0, 0.5))),
        PoissonError())
    for (n_time, n_groups) in ((24, 4), (12, 6))
        Ymiss = Matrix{Union{Missing, Float64}}(missing, n_time, n_groups)
        sim = as_turing_model(model, Ymiss)()
        # The shared infection curve is drawn once; the group axis threads from Y.
        @test length(sim.I_t) == n_time
        @test length(sim.Z_t) == n_time
        @test length(sim.group_levels) == n_groups
        @test size(sim.generated_y_t) == (n_time, n_groups)
        @test size(sim.expected_y_t) == (n_time, n_groups)
        @test length(sim.y) == n_groups
        @test length(sim.y[1].y_t) == n_time
        @test all(>=(0), sim.generated_y_t)
    end
end

@testitem "GroupedIDModel conforms to as_turing_model and prints as a tree" begin
    using ComposableTuringIDModels, Distributions
    model = GroupedIDModel(
        DirectInfections(; Z = RandomWalk(), initialisation = Normal()),
        Hierarchy(; across = IID(Normal(0.0, 0.5))),
        PoissonError())
    Ymiss = Matrix{Union{Missing, Float64}}(missing, 8, 3)
    @test as_turing_model(model, Ymiss) isa ComposableTuringIDModels.DynamicPPL.Model
    # The tree display recurses through the component slots.
    s = sprint(show, MIME("text/plain"), model)
    @test occursin("GroupedIDModel", s)
    @test occursin("infection", s)
    @test occursin("observation", s)
end

@testitem "GroupedIDModel: a swappable combiner changes the group mapping" begin
    using ComposableTuringIDModels, Distributions, Random
    Random.seed!(452)
    inf = DirectInfections(; Z = RandomWalk(), initialisation = Normal(log(50.0), 0.2))
    h = Hierarchy(; mean = Normal(0.0, 0.5), across = IID(Normal(0.0, 0.5)))
    # An additive combiner on the natural scale instead of the multiplicative
    # default: proves the mapping to the observation scale is a swappable field.
    additive = GroupedIDModel(inf, h, PoissonError();
        combiner = (I_t, level) -> max.(I_t .+ level, 0.0))
    Ymiss = Matrix{Union{Missing, Float64}}(missing, 10, 3)
    sim = as_turing_model(additive, Ymiss)()
    @test size(sim.generated_y_t) == (10, 3)
end

@testitem "GroupedIDModel recovers per-group levels under NUTS" tags=[:sample] begin
    using ComposableTuringIDModels, Distributions, Turing, Random, Statistics
    using Turing: returned
    Random.seed!(388)
    # A full composed panel: shared infection + observation, per-group reporting
    # level partially pooled by a Hierarchy. 8 groups for a stable correlation.
    model = GroupedIDModel(
        IDModel(
            DirectInfections(;
                Z = RandomWalk(), initialisation = Normal(log(50.0), 0.2)),
            PoissonError()),
        Hierarchy(; mean = Normal(0.0, 0.5), across = IID(Normal(0.0, 0.5))))

    n_time, n_groups = 20, 8
    Ymiss = Matrix{Union{Missing, Float64}}(missing, n_time, n_groups)
    sim = as_turing_model(model, Ymiss)()
    Ydata = reduce(hcat, [Int.(sim.y[g].y_t) for g in 1:n_groups])
    @test size(Ydata) == (n_time, n_groups)

    posterior = as_turing_model(model, Float64.(Ydata))
    chain = sample(posterior, NUTS(0.8; adtype = Turing.AutoForwardDiff()), 80;
        progress = false)
    @test size(chain, 1) == 80
    draws = reduce(hcat, [g.group_levels for g in vec(returned(posterior, chain))])
    @test size(draws) == (n_groups, 80)
    post_mean = vec(mean(draws; dims = 2))
    @test cor(sim.group_levels, post_mean) > 0.5
end
