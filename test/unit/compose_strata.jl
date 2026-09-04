# End-to-end coverage of `Split`'s weight-matrix mapping through a composed
# `IDModel`, for every infection <-> observation cardinality, and the
# two-argument `as_turing_model(model, Y)` shape resolution that reads the
# infection strata count from the observation model and the data alone.

@testitem "infection_strata reads a weighted Split's map column count" begin
    using ComposableTuringIDModels, Distributions
    using ComposableTuringIDModels: infection_strata
    @test infection_strata(PoissonError(), 4) == 4
    @test infection_strata(Split(PoissonError()), 4) == 4
    @test infection_strata(Split(PoissonError(), [1.0 1.0 1.0]), 1) == 3
end

@testitem "Split(streams::NamedTuple, map) checks the row count" begin
    using ComposableTuringIDModels, Distributions
    ok = Split((a = PoissonError(), b = PoissonError()), [1.0 0.0; 0.0 1.0])
    @test ok isa Split
    @test ok.map == [1.0 0.0; 0.0 1.0]
    @test_throws AssertionError Split(
        (a = PoissonError(), b = PoissonError()), [1.0 0.0 0.0]
    )
end

@testitem "IDModel one-to-one: strata resolved from the data alone" begin
    using ComposableTuringIDModels, Distributions, Random
    Random.seed!(801)
    model = IDModel(
        DirectInfections(;
            Z = Stratify(
                RandomWalk(), Hierarchy(; across = IID(Normal(0, 0.5)))
            ),
            initialisation = Normal(log(50), 0.2)
        ),
        Split(PoissonError())
    )
    Ymiss = Matrix{Union{Missing, Float64}}(missing, 3, 12)
    sim = as_turing_model(model, Ymiss)()
    @test size(sim.I_t) == (3, 12)
    @test keys(sim.generated_y_t) == (:group1, :group2, :group3)
    e = sim.expected_y_t
    @test e.group1 ≈ sim.I_t[1, :]
    @test e.group3 ≈ sim.I_t[3, :]
end

@testitem "IDModel many-to-one: three strata aggregate into one stream" begin
    using ComposableTuringIDModels, Distributions, Random
    Random.seed!(802)
    model = IDModel(
        DirectInfections(;
            Z = Stratify(
                RandomWalk(), Hierarchy(; across = IID(Normal(0, 0.5)))
            ),
            initialisation = Normal(log(30), 0.2)
        ),
        Split(NegativeBinomialError(), [1.0 1.0 1.0])
    )
    Ymiss = Matrix{Union{Missing, Float64}}(missing, 1, 14)   # 1 obs stream
    sim = as_turing_model(model, Ymiss)()
    @test size(sim.I_t) == (3, 14)   # the map's column count builds 3 strata
    @test keys(sim.generated_y_t) == (:group1,)
    @test sim.expected_y_t.group1 ≈ vec(sum(sim.I_t; dims = 1))
end

@testitem "IDModel many-to-many: a weight matrix maps strata onto streams" begin
    using ComposableTuringIDModels, Distributions, Random
    Random.seed!(803)
    W = [1.0 0.0; 0.0 1.0; 1.0 1.0]   # stratum 1, stratum 2, and their sum
    model = IDModel(
        DirectInfections(;
            Z = Stratify(
                RandomWalk(), Hierarchy(; across = IID(Normal(0, 0.5)))
            ),
            initialisation = Normal(log(40), 0.2)
        ),
        Split(
            (
                a = PoissonError(), b = PoissonError(),
                total = PoissonError(),
            ), W
        )
    )
    Ymiss = Matrix{Union{Missing, Float64}}(missing, 3, 10)
    sim = as_turing_model(model, Ymiss)()
    @test size(sim.I_t) == (2, 10)   # size(W, 2) infection strata
    @test keys(sim.generated_y_t) == (:a, :b, :total)
    e = sim.expected_y_t
    @test e.a ≈ sim.I_t[1, :]
    @test e.b ≈ sim.I_t[2, :]
    @test e.total ≈ e.a .+ e.b
end

@testitem "IDModel finer than infections: one stratum into two streams" begin
    using ComposableTuringIDModels, Distributions, Random
    Random.seed!(804)
    W = reshape([0.7, 0.3], 2, 1)   # young, old fractions of one stratum
    model = IDModel(
        DirectInfections(;
            Z = Stratify(RandomWalk(), FixedIntercept(0.0)),
            initialisation = Normal(log(50), 0.2)
        ),
        Split((young = PoissonError(), old = PoissonError()), W)
    )
    Ymiss = Matrix{Union{Missing, Float64}}(missing, 2, 12)
    sim = as_turing_model(model, Ymiss)()
    @test size(sim.I_t) == (1, 12)   # size(W, 2) == 1 infection stratum
    e = sim.expected_y_t
    @test e.young ≈ 0.7 .* vec(sim.I_t)
    @test e.old ≈ 0.3 .* vec(sim.I_t)
end

@testitem "two-argument as_turing_model reads the shape from any data form" begin
    using ComposableTuringIDModels, Distributions, DynamicPPL, Random
    Random.seed!(806)
    infection = DirectInfections(; Z = RandomWalk(), initialisation = Normal())

    # A vector of observations gives its own length, so nothing is restated.
    model = IDModel(infection, PoissonError())
    y = fill(5, 14)
    @test keys(VarInfo(as_turing_model(model, y))) ==
        keys(VarInfo(as_turing_model(model, y, length(y))))
    @test length(as_turing_model(model, Vector{Missing}(missing, 14))().Z_t) == 14

    # A `NamedTuple` of streams gives the stream count and the shared time
    # length.
    split = IDModel(
        DirectInfections(;
            Z = Stratify(RandomWalk(), FixedIntercept(0.0)),
            initialisation = IID(Normal(log(20), 0.2))
        ),
        Split((a = PoissonError(), b = PoissonError()))
    )
    sim = as_turing_model(split, (a = fill(5, 12), b = fill(3, 12)))()
    @test size(sim.I_t) == (2, 12)

    # A scalar `missing` carries no length, and says so.
    @test_throws ArgumentError as_turing_model(model, missing)
end
