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
