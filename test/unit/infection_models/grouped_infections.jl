# Tests for `GroupedInfections`: one shared infection process, replicated
# across groups by a per-group effect, into one `n_groups x n_time` `I_t`
# matrix — the "same many" many-to-many infection-side mapping, and the
# direct successor of `GroupedIDModel`'s core mechanism (issue #180).

@testitem "GroupedInfections constructs and is an AbstractInfectionModel" begin
    using ComposableTuringIDModels, Distributions
    inf = DirectInfections(; Z = RandomWalk(), initialisation = Normal())
    h = Hierarchy(; mean = Normal(0.0, 0.5), across = IID(Normal(0.0, 0.5)))
    model = GroupedInfections(inf, h)
    @test model isa GroupedInfections
    @test model isa AbstractInfectionModel
    @test model.infection_model isa DirectInfections
    @test model.group_effect isa Hierarchy
    @test model.combiner isa Function
    @test !(:n_groups in fieldnames(typeof(model)))
    @test GroupedInfections(inf, Normal(0.0, 0.5)) isa GroupedInfections
end

@testitem "GroupedInfections stacks a shared curve across groups" begin
    using ComposableTuringIDModels, Distributions, Random
    Random.seed!(902)
    inf = DirectInfections(; Z = RandomWalk(), initialisation = Normal(log(50.0), 0.2))
    h = Hierarchy(; mean = Normal(0.0, 0.5), across = IID(Normal(0.0, 0.5)))
    model = GroupedInfections(inf, h)
    for (n_time, n_groups) in ((24, 4), (12, 6))
        sim = as_turing_model(model, (n_time = n_time, n_groups = n_groups))()
        @test size(sim.I_t) == (n_groups, n_time)
        @test all(>=(0), sim.I_t)
        @test length(sim.Z_t.Z_t) == n_time
        @test length(sim.Z_t.group_levels) == n_groups
    end
end

@testitem "GroupedInfections: a swappable combiner changes the group mapping" begin
    using ComposableTuringIDModels, Distributions, Random
    Random.seed!(903)
    inf = DirectInfections(; Z = RandomWalk(), initialisation = Normal(log(50.0), 0.2))
    h = Hierarchy(; mean = Normal(0.0, 0.5), across = IID(Normal(0.0, 0.5)))
    additive = GroupedInfections(inf, h; combiner = (I_t, level) -> max.(I_t .+ level, 0.0))
    sim = as_turing_model(additive, (n_time = 10, n_groups = 3))()
    @test size(sim.I_t) == (3, 10)
end

@testitem "GroupedInfections namespaces the group prior away from the infection latent" begin
    using ComposableTuringIDModels, Distributions
    inf = DirectInfections(; Z = RandomWalk(), initialisation = Normal())
    h = Hierarchy(; across = IID(Normal(0.0, 0.5)))
    model = GroupedInfections(inf, h)
    mdl = as_turing_model(model, (n_time = 8, n_groups = 3))
    draw = rand(mdl)
    names = string.(collect(keys(draw)))
    # The infection process's own innovations are flat (prefix off, unchanged).
    @test any(n -> !occursin('.', n) && occursin("ϵ_t", n), names)
    # The group prior is namespaced under its prior-slot name.
    @test any(startswith(n, "group_levels.") for n in names)
end

@testitem "GroupedInfections conforms to as_turing_model and prints as a tree" begin
    using ComposableTuringIDModels, Distributions
    inf = DirectInfections(; Z = RandomWalk(), initialisation = Normal())
    h = Hierarchy(; across = IID(Normal(0.0, 0.5)))
    model = GroupedInfections(inf, h)
    @test implements_infection_interface(model; n = (n_time = 8, n_groups = 3))
    s = sprint(show, MIME("text/plain"), model)
    @test occursin("GroupedInfections", s)
    @test occursin("DirectInfections", s)
    @test occursin("Hierarchy", s)
end
