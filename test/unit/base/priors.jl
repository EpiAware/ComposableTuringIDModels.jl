@testitem "AbstractPriorModel role: latent-as-prior, not bare distributions" begin
    using ComposableTuringIDModels, Distributions
    @test AbstractLatentModel <: AbstractPriorModel
    # Latent models satisfy the prior contract, so they are prior models.
    @test RandomWalk() isa AbstractPriorModel
    @test AR() isa AbstractPriorModel
    # A bare Distribution is not a prior model (it composes through the seam).
    @test !(Normal() isa AbstractPriorModel)
    @test implements_prior_interface(RandomWalk())
    @test !implements_prior_interface(Normal())
end

@testitem "as_turing_submodel composes a component or a raw prior" begin
    using ComposableTuringIDModels, Distributions, Turing, Random
    Random.seed!(100)
    # The seam threads a component or a raw prior through a submodel.
    @model function use_prior(prior, n)
        θ ~ as_turing_submodel(prior, n)
        return θ
    end
    # A bare Distribution flows through identically to a model.
    @test length(use_prior(Normal(), 4)()) == 4
    @test length(use_prior(RandomWalk(), 6)()) == 6
    @test length(use_prior([Normal(), Normal()], 2)()) == 2
end

@testitem "as_turing_model over a Distribution draws n i.i.d. values" begin
    using ComposableTuringIDModels, Distributions, Random
    Random.seed!(101)
    v = as_turing_model(Normal(), 5)()
    @test length(v) == 5
    # Independent draws, not a single shared value repeated.
    @test !all(==(first(v)), v)
    # Length 1 is a single draw (the scalar case).
    @test length(as_turing_model(Normal(), 1)()) == 1
end

@testitem "as_turing_model over a vector gives one draw per element" begin
    using ComposableTuringIDModels, Distributions, Random
    Random.seed!(102)
    # Homogeneous vector: filldist path.
    vh = as_turing_model([Normal(), Normal()], 2)()
    @test length(vh) == 2
    # Heterogeneous vector: arraydist, element i drawn from distribution i.
    vt = as_turing_model([Normal(0, 1), Normal(5, 0.1)], 2)()
    @test length(vt) == 2
    @test vt[2] > vt[1]                  # second is tight around 5, first around 0
    # Length must match the vector length.
    @test_throws Exception as_turing_model([Normal(), Normal()], 3)()
end

@testitem "bare latent-model-as-prior threads under a linked log-density (#80)" begin
    using ComposableTuringIDModels, Distributions, Turing, Random
    using DynamicPPL: LogDensityFunction, VarInfo, link, getlogjoint
    import LogDensityProblems as LDP
    Random.seed!(180)
    # The exact #80 repro: a bare `AR(damp = RandomWalk())` used to sample via
    # `rand` but ERROR as a linked log-density, because the damping RandomWalk's
    # inner `std`/`ϵ_t`/`rw_init` collided with the AR innovation's. The prior slot
    # now namespaces the whole submodel (prefix-on `as_turing_submodel`), so they
    # cannot collide.
    m = as_turing_model(AR(; damp = RandomWalk()), 8)
    @test rand(m) !== nothing                        # sampled fine before too
    vi = link(VarInfo(m), m)
    ldf = LogDensityFunction(m, getlogjoint, vi)
    val = LDP.logdensity(ldf, zeros(LDP.dimension(ldf)))   # previously threw
    @test isfinite(val)
    # And it samples under NUTS end-to-end.
    chn = sample(m, NUTS(0.8; adtype = Turing.AutoForwardDiff()), 40;
        progress = false)
    @test size(chn, 1) == 40
end

@testitem "order inference: single distribution ⇒ 1, vector ⇒ length" begin
    using ComposableTuringIDModels, Distributions
    using ComposableTuringIDModels: _prior_order
    # A single distribution (or richer prior) implies order 1, a vector fixes it
    # to the vector length.
    @test _prior_order(Normal()) == 1
    @test _prior_order([Normal(), Normal()]) == 2
    @test _prior_order(RandomWalk()) == 1
    # A single (non-vector) damping distribution coerces to an order-1 AR.
    @test AR(; damp = Normal()).p == 1
    @test AR(; damp = [truncated(Normal(0, 0.05), 0, 1)]).p == 1
    @test MA(; θ = [Normal(), Normal()]).q == 2
    # A vector-of-distributions init sets the differencing order.
    @test DiffLatentModel(; model = RandomWalk(), init = [Normal(), Normal()]).d == 2
end

@testitem "raw priors are stored unwrapped in component fields" begin
    using ComposableTuringIDModels, Distributions
    # Scalar prior slots hold the bare distribution (drawn with a native tilde).
    @test NormalError(truncated(Normal(0, 1), 0, Inf)).std isa Distribution
    @test NegativeBinomialError(HalfNormal(0.1)).cluster_factor isa Distribution
    @test Intercept(; intercept = Normal()).intercept isa Distribution
    @test HierarchicalNormal(truncated(Normal(0, 1), 0, Inf)).std isa Distribution
    # Vector prior slots hold the raw vector.
    @test AR(; damp = [Normal(), Normal()],
        init = [Normal(), Normal()]).damp isa AbstractVector
    # A process prior slot holds the latent model unchanged.
    @test AR(; damp = RandomWalk()).damp isa RandomWalk
    # Positional bare-`Distribution` constructor forms build valid latents.
    for m in (AR(Normal(), Normal()),
        HierarchicalNormal(truncated(Normal(0, 1), 0, Inf)),
        HierarchicalNormal(0.5, truncated(Normal(0, 1), 0, Inf)),
        Intercept(; intercept = Normal()),
        MA(truncated(Normal(0.0, 0.05), -1, 1)),
        DiffLatentModel(RandomWalk(), Normal(); d = 2))
        @test m isa AbstractLatentModel
    end
end

@testitem "priors compose as submodels (distribution and latent)" begin
    using ComposableTuringIDModels, Distributions, Turing, Random
    Random.seed!(103)
    # A RandomWalk used directly as a (time-varying) prior returns a length-n path.
    @test length(as_turing_model(RandomWalk(), 8)()) == 8
    # A scalar parameter is a plain native-tilde draw (no submodel): NormalError's
    # σ has no `.θ` namespace path.
    draw = rand(as_turing_model(NormalError(), missing, fill(10.0, 5)))
    @test any(k -> occursin("σ", string(k)), keys(draw))
end
