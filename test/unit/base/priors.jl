@testitem "AbstractPriorModel role: wrapper and latent-as-prior" begin
    using ComposableTuringIDModels, Distributions
    @test AbstractLatentModel <: AbstractPriorModel
    @test BroadcastPrior(Normal()) isa AbstractPriorModel
    # Latent models satisfy the prior contract, so they are prior models.
    @test RandomWalk() isa AbstractPriorModel
    @test AR() isa AbstractPriorModel
    # Interface checker mirrors the other per-role checks.
    @test implements_prior_interface(BroadcastPrior(Normal()))
    @test implements_prior_interface(RandomWalk())
    @test !implements_prior_interface(Normal())   # a bare Distribution is not a prior model
end

@testitem "as_prior coerces distributions and passes prior models through" begin
    using ComposableTuringIDModels, Distributions
    @test as_prior(Normal()) isa BroadcastPrior
    @test as_prior([Normal(), Normal(2, 1)]) isa BroadcastPrior
    rw = RandomWalk()
    @test as_prior(rw) === rw            # a prior/latent model is accepted unchanged
    # Coercion is name-free: namespacing happens at the component's call site
    # (prefix-on `to_submodel`), not by carrying a name here (issue #80 is covered
    # by the linked-log-density test below).
    @test as_prior(Normal()) isa BroadcastPrior
    @test as_prior([Normal(), Normal()]) isa BroadcastPrior
end

@testitem "bare latent-model-as-prior threads under a linked log-density (#80)" begin
    using ComposableTuringIDModels, Distributions, Turing, Random
    using DynamicPPL: LogDensityFunction, VarInfo, link, getlogjoint
    import LogDensityProblems as LDP
    Random.seed!(180)
    # The exact #80 repro: a bare `AR(damp = RandomWalk())` used to sample via
    # `rand` but ERROR as a linked log-density, because the damping RandomWalk's
    # inner `std`/`ϵ_t`/`rw_init` collided with the AR innovation's. The prior slot
    # now namespaces the whole submodel (prefix-on `to_submodel`), so they cannot
    # collide.
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

@testitem "BroadcastPrior scalar (repeat-one) mode" begin
    using ComposableTuringIDModels, Distributions, Random
    Random.seed!(101)
    # length-1 scalar parameter: one draw, read back with `only`.
    v1 = as_turing_model(BroadcastPrior(Normal()), 1)()
    @test length(v1) == 1
    @test only(v1) isa Real
    # repeat-one: a single random variable repeated to length n, so all equal —
    # a global coefficient is not expanded into n i.i.d. draws.
    v = as_turing_model(BroadcastPrior(Normal()), 5)()
    @test length(v) == 5
    @test all(==(first(v)), v)
end

@testitem "BroadcastPrior vector mode gives one draw per element" begin
    using ComposableTuringIDModels, Distributions, Random
    Random.seed!(102)
    # Homogeneous vector: filldist (as _expand_dist did).
    vh = as_turing_model(BroadcastPrior([Normal(), Normal()]), 2)()
    @test length(vh) == 2
    # Heterogeneous vector: arraydist, element i drawn from distribution i.
    vt = as_turing_model(BroadcastPrior([Normal(0, 1), Normal(5, 0.1)]), 2)()
    @test length(vt) == 2
    @test vt[2] > vt[1]                  # second is tight around 5, first around 0
    # Length must match the vector length.
    @test_throws Exception as_turing_model(BroadcastPrior([Normal(), Normal()]), 3)()
end

@testitem "ergonomic constructor forms coerce priors" begin
    using ComposableTuringIDModels, Distributions
    using ComposableTuringIDModels: _prior_order
    # A single-distribution (repeat-one) or richer prior implies order 1, while a
    # vector wrapper fixes it to the vector length.
    @test _prior_order(BroadcastPrior(Normal())) == 1
    @test _prior_order(BroadcastPrior([Normal(), Normal()])) == 2
    # An already-wrapped prior model passes through unchanged.
    bp = BroadcastPrior(Normal())
    @test as_prior(bp) === bp
    # A single (non-vector) damping distribution coerces to an order-1 AR.
    @test AR(; damp = Normal()).p == 1
    # Positional bare-`Distribution` constructor forms coerce to the prior interface.
    for m in (AR(Normal(), Normal()),
        HierarchicalNormal(truncated(Normal(0, 1), 0, Inf)),
        HierarchicalNormal(0.5, truncated(Normal(0, 1), 0, Inf)),
        Intercept(; intercept = Normal()),
        MA(truncated(Normal(0.0, 0.05), -1, 1)),
        DiffLatentModel(RandomWalk(), Normal(); d = 2))
        @test m isa AbstractLatentModel
    end
    @test NegativeBinomialError(HalfNormal(0.1)).cluster_factor isa AbstractPriorModel
    @test NormalError(truncated(Normal(0, 1), 0, Inf)).std isa AbstractPriorModel
end

@testitem "priors compose as submodels (wrapper and latent)" begin
    using ComposableTuringIDModels, Distributions, Turing, Random
    Random.seed!(103)
    # A RandomWalk used directly as a (time-varying) prior returns a length-n path.
    @test length(as_turing_model(RandomWalk(), 8)()) == 8
    # Both the default wrapper and a latent model compose via to_submodel(..., false).
    @model function use_prior(prior, n)
        θ ~ to_submodel(as_turing_model(prior, n), false)
        return θ
    end
    @test length(use_prior(BroadcastPrior(Normal()), 4)()) == 4
    @test length(use_prior(RandomWalk(), 6)()) == 6
end
