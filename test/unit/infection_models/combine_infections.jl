# Tests for `CombineInfections`: several (possibly different) infection
# processes stacked into one `n_strata x n_time` `I_t` matrix, the "different
# many" many-to-many infection-side mapping (issue #180).

@testitem "CombineInfections constructs and is an AbstractInfectionModel" begin
    using ComposableTuringIDModels, Distributions
    a = DirectInfections(; Z = RandomWalk(), initialisation = Normal())
    b = DirectInfections(; Z = RandomWalk(), initialisation = Normal())
    model = CombineInfections([a, b])
    @test model isa CombineInfections
    @test model isa AbstractInfectionModel
    @test model.names == ["inf1", "inf2"]

    named = CombineInfections([a, b], ["north", "south"])
    @test named.names == ["north", "south"]
end

@testitem "CombineInfections requires at least one model and matching names" begin
    using ComposableTuringIDModels, Distributions
    a = DirectInfections(; Z = RandomWalk(), initialisation = Normal())
    @test_throws AssertionError CombineInfections(typeof(a)[])
    @test_throws AssertionError CombineInfections([a], ["x", "y"])
end

@testitem "CombineInfections stacks independent infection curves" begin
    using ComposableTuringIDModels, Distributions, Random
    Random.seed!(901)
    a = DirectInfections(; Z = RandomWalk(), initialisation = Normal(log(50.0), 0.2))
    b = DirectInfections(; Z = RandomWalk(), initialisation = Normal(log(20.0), 0.2))
    model = CombineInfections([a, b], ["north", "south"])
    n = 12
    sim = as_turing_model(model, n)()
    @test size(sim.I_t) == (2, n)
    @test all(>=(0), sim.I_t)
    @test keys(sim.Z_t) == (:north, :south)
    @test length(sim.Z_t.north) == n
    @test length(sim.Z_t.south) == n
    # The two curves are drawn independently, so they are not deterministically
    # equal (different `initialisation` priors, checked via distinct means).
    @test sim.I_t[1, :] != sim.I_t[2, :]
end

@testitem "CombineInfections prefixes each model so variables never collide" begin
    using ComposableTuringIDModels, Distributions
    using DynamicPPL: VarName
    a = DirectInfections(; Z = RandomWalk(), initialisation = Normal())
    b = DirectInfections(; Z = RandomWalk(), initialisation = Normal())
    model = CombineInfections([a, b], ["north", "south"])
    mdl = as_turing_model(model, 8)
    draw = rand(mdl)
    names = string.(collect(keys(draw)))
    @test any(startswith(n, "north.") for n in names)
    @test any(startswith(n, "south.") for n in names)
end

@testitem "CombineInfections conforms to as_turing_model and prints as a tree" begin
    using ComposableTuringIDModels, Distributions
    a = DirectInfections(; Z = RandomWalk(), initialisation = Normal())
    b = DirectInfections(; Z = RandomWalk(), initialisation = Normal())
    model = CombineInfections([a, b], ["north", "south"])
    @test implements_infection_interface(model)
    s = sprint(show, MIME("text/plain"), model)
    @test occursin("CombineInfections", s)
    @test occursin("DirectInfections", s)
end
