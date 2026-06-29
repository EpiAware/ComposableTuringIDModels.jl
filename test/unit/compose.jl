@testitem "composed model: prior simulation and generated quantities" begin
    using EpiAwarePrototype, Distributions, Random
    Random.seed!(21)
    model = EpiAwareModel(
        DirectInfections(; Z = RandomWalk(), initialisation_prior = Normal()),
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
    using EpiAwarePrototype, Distributions, Random
    using DynamicPPL: fix, condition
    Random.seed!(22)
    model = EpiAwareModel(
        DirectInfections(; Z = RandomWalk(), initialisation_prior = Normal()),
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
    using EpiAwarePrototype, Distributions, Turing, Random
    Random.seed!(23)
    model = EpiAwareModel(
        DirectInfections(; Z = RandomWalk(), initialisation_prior = Normal()),
        NegativeBinomialError())
    n = 20
    y = as_turing_model(model, missing, n)()
    cond_model = as_turing_model(model, y.generated_y_t, n)
    chn = sample(cond_model, NUTS(), 50; progress = false)
    @test chn !== nothing
end
