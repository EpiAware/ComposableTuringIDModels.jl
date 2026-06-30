# Role-type and interface-conformance tests. Two kinds, mirroring upstream's
# inline role assertions but driven through the reusable interface checkers:
#   1. role/type assertions — each concrete model, and each manipulator's OUTPUT,
#      `isa` its role (a wrapped latent stays latent, a wrapped observation stays
#      observation), plus an in-test extension proof.
#   2. behavioural conformance — build the model, draw, check output shape.

@testitem "role supertypes form the intended shallow hierarchy" begin
    using EpiAwarePrototype
    for R in (AbstractLatentModel, AbstractInfectionModel, AbstractObservationModel)
        @test R <: AbstractEpiAwareModel
    end
    # The error sub-role sits under the observation role.
    @test AbstractObservationErrorModel <: AbstractObservationModel
    @test AbstractObservationErrorModel <: AbstractEpiAwareModel
    # Roles are distinct branches (no role is a subtype of another sibling).
    @test !(AbstractLatentModel <: AbstractObservationModel)
    @test !(AbstractInfectionModel <: AbstractLatentModel)
end

@testitem "latent models and latent manipulators are AbstractLatentModel" begin
    using EpiAwarePrototype, Distributions, Accessors
    # Concrete latent models.
    for m in (IID(Normal()), HierarchicalNormal(), RandomWalk(), AR(), MA(),
        Intercept(Normal()), FixedIntercept(0.1), Null())
        @test m isa AbstractLatentModel
    end
    # A wrapped/combined latent is still a latent (the key compositional contract).
    ar = AR()
    @test DiffLatentModel(; model = ar, init_priors = [Normal(), Normal()]) isa
          AbstractLatentModel
    @test TransformLatentModel(ar, x -> exp.(x)) isa AbstractLatentModel
    @test PrefixLatentModel(; model = ar, prefix = "P") isa AbstractLatentModel
    @test RecordExpectedLatent(ar) isa AbstractLatentModel
    @test CombineLatentModels([Intercept(Normal()), AR()]) isa AbstractLatentModel
    @test ConcatLatentModels([Intercept(Normal()), AR()]) isa AbstractLatentModel
    @test BroadcastLatentModel(RandomWalk(), 7, RepeatEach()) isa AbstractLatentModel
    @test arma() isa AbstractLatentModel
    @test arima() isa AbstractLatentModel
    # Hierarchy is a latent process over the grouping dimension.
    @test Hierarchy(FixedIntercept(0.0), (@optic _.intercept), AR()) isa
          AbstractLatentModel
end

@testitem "infection models are AbstractInfectionModel; ODE params are latent" begin
    using EpiAwarePrototype, Distributions, OrdinaryDiffEq
    data = EpiData([0.2, 0.3, 0.5], exp)
    for m in (DirectInfections(; Z = RandomWalk()), ExpGrowthRate(; rt = RandomWalk()),
        Renewal(; data = data, rt = RandomWalk()))
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
    using EpiAwarePrototype, Distributions
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
    @test StackObservationModels((a = PoissonError(), b = PoissonError())) isa
          AbstractObservationModel
end

@testitem "role slots reject wrong-role components at construction" begin
    using EpiAwarePrototype, Distributions
    latent = RandomWalk()
    infection = DirectInfections(; Z = RandomWalk())
    obs = PoissonError()
    # Correct order constructs (infection then observation).
    @test EpiAwareModel(infection, obs) isa EpiAwareModel
    # Wrong order (the classic foot-gun) fails at construction, not at sampling.
    # Positional constructors fail dispatch (MethodError). The wrong-role
    # component is rejected before any sampling happens.
    @test_throws MethodError EpiAwareModel(obs, infection)
    # A latent model is not an infection model, so it cannot fill the infection
    # slot either.
    @test_throws MethodError EpiAwareModel(latent, obs)
    # EpiProblem enforces the same role slots (keyword struct → TypeError).
    @test_throws Union{MethodError, TypeError} EpiProblem(
        epi_model = obs, observation_model = obs, tspan = (1, 10))
    # A latent manipulator cannot wrap an observation model (slot is latent;
    # keyword constructor → TypeError).
    @test_throws Union{MethodError, TypeError} DiffLatentModel(; model = obs,
        init_priors = [Normal(), Normal()])
    # An observation modifier cannot wrap a latent model (slot is observation;
    # positional constructor → MethodError).
    @test_throws MethodError LatentDelay(latent, [0.5, 0.5])
end

@testitem "reusable interface checkers confirm role conformance" begin
    using EpiAwarePrototype, Distributions
    data = EpiData([0.2, 0.3, 0.5], exp)
    # Each checker is true for an in-role model implementing its as_turing_model.
    @test implements_latent_interface(RandomWalk())
    @test implements_latent_interface(AR(); n = 12)
    @test implements_infection_interface(DirectInfections(; Z = RandomWalk()))
    @test implements_infection_interface(Renewal(; data = data, rt = RandomWalk()); n = 20)
    @test implements_observation_interface(PoissonError())
    @test implements_observation_interface(NegativeBinomialError())
    # A model is NOT in a role it does not belong to.
    @test !implements_observation_interface(RandomWalk())
    @test !implements_latent_interface(PoissonError())
    @test !implements_infection_interface(RandomWalk())
end

@testitem "a user-defined struct in a role composes via its interface" begin
    using EpiAwarePrototype, Distributions
    using DynamicPPL: @model

    # In-test extension proof: a tiny custom latent model. Subtyping the role and
    # implementing the role's as_turing_model is all that is required to compose.
    struct ConstantLatent <: AbstractLatentModel
        value::Float64
    end
    @model function EpiAwarePrototype.as_turing_model(m::ConstantLatent, n)
        return fill(m.value, n)
    end

    custom = ConstantLatent(0.5)
    @test custom isa AbstractLatentModel
    @test implements_latent_interface(custom)

    # It slots into an infection model's latent position (the latent is now
    # folded into the infection model) and the composed model runs.
    infection = DirectInfections(; Z = custom, initialisation_prior = Normal())
    model = EpiAwareModel(infection, PoissonError())
    @test model isa EpiAwareModel
    out = as_turing_model(model, missing, 10)()
    @test length(out.Z_t) == 10
    @test all(==(0.5), out.Z_t)
end

@testitem "behavioural conformance: each role's output has the expected shape" begin
    using EpiAwarePrototype, Distributions, Random
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
    # Observation: maps an expected series to observed counts.
    y = as_turing_model(PoissonError(), missing, fill(10.0, n))()
    @test length(y) == n
    @test all(>=(0), y)
end
