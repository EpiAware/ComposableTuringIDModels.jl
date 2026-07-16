# Role-type and interface-conformance tests. Two kinds, mirroring upstream's
# inline role assertions but driven through the reusable interface checkers:
#   1. role/type assertions — each concrete model, and each manipulator's OUTPUT,
#      `isa` its role (a wrapped latent stays latent, a wrapped observation stays
#      observation), plus an in-test extension proof.
#   2. behavioural conformance — build the model, draw, check output shape.

@testitem "role supertypes form the intended shallow hierarchy" begin
    using ComposableTuringIDModels
    for R in (AbstractLatentModel, AbstractInfectionModel, AbstractObservationModel)
        @test R <: AbstractComposableModel
    end
    # The error sub-role sits under the observation role.
    @test AbstractObservationErrorModel <: AbstractObservationModel
    @test AbstractObservationErrorModel <: AbstractComposableModel
    # Roles are distinct branches (no role is a subtype of another sibling).
    @test !(AbstractLatentModel <: AbstractObservationModel)
    @test !(AbstractInfectionModel <: AbstractLatentModel)
end

@testitem "latent models and latent manipulators are AbstractLatentModel" begin
    using ComposableTuringIDModels, Distributions
    # Concrete latent models.
    for m in (IID(Normal()), HierarchicalNormal(), RandomWalk(), AR(), MA(),
        Intercept(Normal()), FixedIntercept(0.1), Null())
        @test m isa AbstractLatentModel
    end
    # A wrapped/combined latent is still a latent (the key compositional contract).
    ar = AR()
    @test DiffLatentModel(; model = ar, init = [Normal(), Normal()]) isa
          AbstractLatentModel
    @test TransformLatentModel(ar, x -> exp.(x)) isa AbstractLatentModel
    @test PrefixLatentModel(; model = ar, prefix = "P") isa AbstractLatentModel
    @test RecordExpectedLatent(ar) isa AbstractLatentModel
    @test CombineLatentModels([Intercept(Normal()), AR()]) isa AbstractLatentModel
    @test ConcatLatentModels([Intercept(Normal()), AR()]) isa AbstractLatentModel
    @test BroadcastLatentModel(RandomWalk(), 7, RepeatEach()) isa AbstractLatentModel
    @test arma() isa AbstractLatentModel
    @test arima() isa AbstractLatentModel
end

@testitem "infection models are AbstractInfectionModel; ODE params are latent" begin
    using ComposableTuringIDModels, Distributions, OrdinaryDiffEq
    gen_int = [0.2, 0.3, 0.5]
    for m in (DirectInfections(; Z = RandomWalk()), ExpGrowthRate(; rt = RandomWalk()),
        Renewal(gen_int; rt = RandomWalk()))
        @test m isa AbstractInfectionModel
    end
    # ODE parameter structs play the latent role (they feed an ODEProcess slot).
    sir = SIRParams(tspan = (0.0, 30.0), infectiousness = LogNormal(log(0.3), 0.05),
        recovery_rate = LogNormal(log(0.1), 0.05),
        initial_prop_infected = Beta(1, 99))
    @test sir isa AbstractLatentModel
    proc = ODEProcess(; params = sir, sol2infs = sol -> sol[2, :])
    @test proc isa AbstractInfectionModel
end

@testitem "observation models and observation modifiers are AbstractObservationModel" begin
    using ComposableTuringIDModels, Distributions
    for m in (PoissonError(), NegativeBinomialError())
        @test m isa AbstractObservationModel
        @test m isa AbstractObservationErrorModel
    end
    # A wrapped observation model is still an observation model.
    @test LatentDelay(PoissonError(), [0.5, 0.5]) isa AbstractObservationModel
    @test Ascertainment(PoissonError(), FixedIntercept(0.1)) isa AbstractObservationModel
    @test Aggregate(PoissonError(), [0, 7]) isa AbstractObservationModel
    @test TransformObservationModel(PoissonError()) isa AbstractObservationModel
    @test PrefixObservationModel(; model = PoissonError(), prefix = "P") isa
          AbstractObservationModel
    @test RecordExpectedObs(PoissonError()) isa AbstractObservationModel
    @test Split((a = PoissonError(), b = PoissonError())) isa
          AbstractObservationModel
    @test Split(PoissonError()) isa AbstractObservationModel
end

@testitem "role slots reject wrong-role components at construction" begin
    using ComposableTuringIDModels, Distributions
    latent = RandomWalk()
    infection = DirectInfections(; Z = RandomWalk())
    obs = PoissonError()
    # Correct order constructs (infection then observation).
    @test IDModel(infection, obs) isa IDModel
    # Wrong order (the classic foot-gun) fails at construction, not at sampling.
    # Positional constructors fail dispatch (MethodError). The wrong-role
    # component is rejected before any sampling happens.
    @test_throws MethodError IDModel(obs, infection)
    # A latent model is not an infection model, so it cannot fill the infection
    # slot either.
    @test_throws MethodError IDModel(latent, obs)
    # IDProblem enforces the same role slots (keyword struct → TypeError).
    @test_throws Union{MethodError, TypeError} IDProblem(
        infection = obs, observation_model = obs, tspan = (1, 10))
    # A latent manipulator cannot wrap an observation model (slot is latent;
    # keyword constructor → TypeError).
    @test_throws Union{MethodError, TypeError} DiffLatentModel(; model = obs,
        init = [Normal(), Normal()])
    # An observation modifier cannot wrap a latent model (slot is observation;
    # positional constructor → MethodError).
    @test_throws MethodError LatentDelay(latent, [0.5, 0.5])
end

@testitem "reusable interface checkers confirm role conformance" begin
    using ComposableTuringIDModels, Distributions
    gen_int = [0.2, 0.3, 0.5]
    # Each checker is true for an in-role model implementing its as_turing_model.
    @test implements_latent_interface(RandomWalk())
    @test implements_latent_interface(AR(); n = 12)
    @test implements_infection_interface(DirectInfections(; Z = RandomWalk()))
    @test implements_infection_interface(Renewal(gen_int; rt = RandomWalk()); n = 20)
    @test implements_observation_interface(PoissonError())
    @test implements_observation_interface(NegativeBinomialError())
    # A model is NOT in a role it does not belong to.
    @test !implements_observation_interface(RandomWalk())
    @test !implements_latent_interface(PoissonError())
    @test !implements_infection_interface(RandomWalk())
end

@testitem "a user-defined struct in a role composes via its interface" begin
    using ComposableTuringIDModels, Distributions
    using DynamicPPL: @model

    # In-test extension proof: a tiny custom latent model. Subtyping the role and
    # implementing the role's as_turing_model is all that is required to compose.
    struct ConstantLatent <: AbstractLatentModel
        value::Float64
    end
    @model function ComposableTuringIDModels.as_turing_model(m::ConstantLatent, n)
        return fill(m.value, n)
    end

    custom = ConstantLatent(0.5)
    @test custom isa AbstractLatentModel
    @test implements_latent_interface(custom)

    # It slots into an infection model's latent position (the latent is now
    # folded into the infection model) and the composed model runs.
    infection = DirectInfections(; Z = custom, initialisation = Normal())
    model = IDModel(infection, PoissonError())
    @test model isa IDModel
    out = as_turing_model(model, missing, 10)()
    @test length(out.Z_t) == 10
    @test all(==(0.5), out.Z_t)
end

@testitem "behavioural conformance: each role's output has the expected shape" begin
    using ComposableTuringIDModels, Distributions, Random
    Random.seed!(101)
    n = 12
    # Latent: length-n path.
    z = as_turing_model(RandomWalk(), n)()
    @test length(z) == n
    @test eltype(z) <: Real
    # Infection: generates its own latent and maps it to an infection path,
    # returning (; I_t, Z_t).
    inf = as_turing_model(DirectInfections(; Z = RandomWalk()), n)()
    @test length(inf.I_t) == n
    @test length(inf.Z_t) == n
    @test all(>=(0), inf.I_t)
    # Observation: maps an expected series to observed counts, returning the
    # uniform `(; y_t, expected)` contract.
    y = as_turing_model(PoissonError(), missing, fill(10.0, n))()
    @test keys(y) == (:y_t, :expected)
    @test length(y.y_t) == n
    @test all(>=(0), y.y_t)
    @test y.expected == fill(10.0, n)
end
