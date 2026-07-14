@testitem "Hierarchy is an observation model that fans out per data stream" begin
    using ComposableTuringIDModels, Distributions, Accessors, Random
    Random.seed!(38)
    base = Ascertainment(PoissonError(), FixedIntercept(0.0); latent_prefix = "")
    h = Hierarchy(base, (@optic _.latent_model.intercept), HierarchicalNormal())
    @test h isa AbstractObservationModel
    @test implements_observation_interface(h)

    yt = (cases = missing, deaths = missing)
    sm = as_turing_model(h, yt, fill(100.0, 12))
    draw = rand(sm)
    names = string.(collect(keys(draw)))
    # The cross-group relationship (HierarchicalNormal) is sampled once, flat.
    @test any(==("std"), names)
    # Each stream's error variables are prefixed by its name.
    @test any(startswith("cases."), names)
    @test any(startswith("deaths."), names)
end

@testitem "Hierarchy returns the Split contract with per-stream series" begin
    using ComposableTuringIDModels, Distributions, Accessors, Random
    Random.seed!(39)
    base = Ascertainment(PoissonError(), FixedIntercept(0.0); latent_prefix = "")
    h = Hierarchy(base, (@optic _.latent_model.intercept), HierarchicalNormal())
    out = as_turing_model(h, (cases = missing, deaths = missing), fill(50.0, 8))()
    @test keys(out) == (:y_t, :expected)
    @test keys(out.y_t) == (:cases, :deaths)
    @test length(out.y_t.cases) == 8
    @test length(out.y_t.deaths) == 8
end

@testitem "Hierarchy pools a field buried inside a LatentDelay wrapper" begin
    using ComposableTuringIDModels, Distributions, Accessors, Random
    Random.seed!(40)
    # The pooled ascertainment effect sits two levels deep, under the delay.
    base = LatentDelay(
        Ascertainment(PoissonError(), FixedIntercept(0.0); latent_prefix = ""),
        [0.4, 0.3, 0.2, 0.1])
    lens = @optic _.model.latent_model.intercept
    @test lens(base) == 0.0
    h = Hierarchy(base, lens, HierarchicalNormal())
    out = as_turing_model(h, (cases = missing, deaths = missing), fill(50.0, 12))()
    # The shared delay shortens each stream's expected series (12 → 9), proving
    # the lens reached the pooled field inside the LatentDelay wrapper.
    @test keys(out.expected) == (:cases, :deaths)
    @test length(out.expected.cases) == 9
    @test length(out.expected.deaths) == 9
end

@testitem "Hierarchy needs a NamedTuple y_t to derive its groups" begin
    using ComposableTuringIDModels, Distributions, Accessors
    base = Ascertainment(PoissonError(), FixedIntercept(0.0); latent_prefix = "")
    h = Hierarchy(base, (@optic _.latent_model.intercept), HierarchicalNormal())
    @test_throws Exception as_turing_model(h, missing, fill(10.0, 5))()
end

@testitem "Hierarchy is gradient-safe under NUTS (value-threaded pooling)" begin
    using ComposableTuringIDModels, Distributions, Accessors, Turing, Random
    Random.seed!(41)
    base = Ascertainment(PoissonError(), FixedIntercept(0.0); latent_prefix = "")
    h = Hierarchy(base, (@optic _.latent_model.intercept), HierarchicalNormal())
    ydata = (cases = fill(7, 10), deaths = fill(3, 10))
    mdl = as_turing_model(h, ydata, fill(10.0, 10))
    chn = sample(mdl, NUTS(), 20; progress = false)
    @test size(chn, 1) == 20
end
