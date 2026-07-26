# Tests for the ImportedCases renewal modifier (issue #189).

@testitem "ImportedCases struct carries the importation prior" begin
    using ComposableTuringIDModels, Distributions
    # Constant unknown rate (Distribution)
    ic = ImportedCases(Normal(0.0, 0.1))
    @test ic.importation_rate isa Normal
    @test ic isa ComposableTuringIDModels.AbstractRenewalModifier

    # Time-varying rate (latent process)
    ic_tv = ImportedCases(RandomWalk())
    @test ic_tv.importation_rate isa RandomWalk
end

@testitem "ImportedCases composes into a Renewal step" begin
    using ComposableTuringIDModels, Distributions
    ic = ImportedCases(Normal(0.0, 0.1))
    gen_int = [0.2, 0.3, 0.5]
    # Via the positional modifier constructor
    r = Renewal(gen_int, ic; rt = RandomWalk(), initialisation = Normal())
    @test r.recurrent_step isa ComposableTuringIDModels.RenewalStep
    @test only(r.recurrent_step.modifiers) isa ImportedCases
    @test r.recurrent_step.modifiers[1].importation_rate isa Normal

    # Composed with another modifier
    ic2 = ImportedCases(Normal(0.0, 0.1))
    r2 = Renewal(gen_int, ic2,
        ComposableTuringIDModels.SusceptibleDepletion(1000.0);
        rt = RandomWalk(), initialisation = Normal())
    @test length(r2.recurrent_step.modifiers) == 2
    @test r2.recurrent_step.modifiers[1] isa ImportedCases
    @test r2.recurrent_step.modifiers[2] isa
          ComposableTuringIDModels.SusceptibleDepletion
end

@testitem "ImportedCases with constant rate adds to incidence" begin
    using ComposableTuringIDModels, Distributions, Random
    using DynamicPPL: fix
    Random.seed!(1891)
    gen_int = [0.2, 0.3, 0.5]
    logR = log(1.0)  # Rt = 1, no growth

    # No importation: flat incidence.
    plain = Renewal(gen_int; rt = FixedIntercept(logR),
        initialisation = Normal())
    fixinit = (init_incidence = log(1.0),)
    I_plain = fix(as_turing_model(plain, 30), fixinit)().I_t

    # With importation: incidence grows.
    ic = ImportedCases(Normal(0.5, 0.01))
    imported = Renewal(gen_int, ic; rt = FixedIntercept(logR),
        initialisation = Normal())
    I_imported = fix(as_turing_model(imported, 30), fixinit)().I_t

    @test all(isfinite, I_imported)
    @test all(>=(0), I_imported)
    @test I_imported[end] > I_plain[end]  # Importation adds to incidence
end

@testitem "ImportedCases with time-varying rate samples via as_turing_submodel" begin
    using ComposableTuringIDModels, Distributions, Random
    Random.seed!(1892)
    gen_int = [0.2, 0.3, 0.5]
    # Time-varying importation via RandomWalk
    ic_tv = ImportedCases(RandomWalk())
    r = Renewal(gen_int, ic_tv; rt = RandomWalk(),
        initialisation = Normal())
    out = as_turing_model(r, 20)()
    @test length(out.I_t) == 20
    @test all(isfinite, out.I_t)
    @test all(>=(0), out.I_t)

    # The importation rate is a sampled parameter.
    draw = rand(as_turing_model(r, 20))
    @test any(k -> startswith(string(k), "import_rates"), keys(draw))
end

@testitem "ImportedCases with fixed rate samples under NUTS" tags=[:sample] begin
    using ComposableTuringIDModels, Distributions, Turing, Random
    Random.seed!(1893)
    gen_int = [0.2, 0.3, 0.5]
    # An unconstrained prior is fine: `Renewal`'s transformation makes the rate
    # positive, so the sampler cannot drive incidence negative.
    ic = ImportedCases(Normal(0.0, 0.1))
    model = IDModel(
        Renewal(gen_int, ic; rt = RandomWalk(),
            initialisation = Normal()),
        PoissonError())
    y = as_turing_model(model, missing, 20)().generated_y_t
    chn = sample(as_turing_model(model, y, 20), NUTS(), 30;
        progress = false)
    @test chn !== nothing
end

@testitem "the modifier itself is a no-op; the rate is added by the step" begin
    using ComposableTuringIDModels, Distributions
    using ComposableTuringIDModels: modifier_init_state, apply_modifier
    ic = ImportedCases(Normal(0.0, 0.1))
    # The rate is a prior, not a constant on the step, so the modifier carries
    # no state and leaves the proposed incidence untouched.
    @test modifier_init_state(ic) == 0.0
    @test apply_modifier(ic, 3.5, modifier_init_state(ic)) == (3.5, 0.0)
end

@testitem "renewal steps without an ImportedCases modifier do not import" begin
    using ComposableTuringIDModels, Distributions
    using ComposableTuringIDModels: ConstantRenewalStep, RenewalStep,
                                    SusceptibleDepletion, _has_imported_cases,
                                    _import_model
    core = ConstantRenewalStep(reverse([0.2, 0.3, 0.5]))
    # The bare force-of-infection core, a modifier-free step, and a step with
    # only non-importing modifiers all take the plain (scalar-`Rt`) scan path.
    @test !_has_imported_cases(core)
    @test !_has_imported_cases(RenewalStep(core))
    @test !_has_imported_cases(RenewalStep(core, (SusceptibleDepletion(1000.0),)))
    # `nothing` is the inferred-generation-interval path, whose step is built
    # per draw rather than carried on the model.
    @test !_has_imported_cases(nothing)
    # Asking a non-importing step for its rate is a caller error, not a
    # silently-missing prior.
    @test_throws ArgumentError _import_model(RenewalStep(core))
end

@testitem "a negative unconstrained rate still imports rather than subtracting" begin
    using ComposableTuringIDModels, Distributions
    using DynamicPPL: fix
    gen_int = [0.2, 0.3, 0.5]
    logR = log(1.0)  # Rt = 1, no growth
    fixinit = (init_incidence = log(1.0),)

    plain = Renewal(gen_int; rt = FixedIntercept(logR),
        initialisation = Normal())
    I_plain = fix(as_turing_model(plain, 30), fixinit)().I_t

    # A strongly negative *unconstrained* rate: the transformation maps it to
    # exp(-2) ≈ 0.135 imports per step. Untransformed this subtracted 2.0 from
    # every incidence, driving the series onto the 1e-6 floor.
    imported = Renewal(gen_int, ImportedCases(FixedIntercept(-2.0));
        rt = FixedIntercept(logR), initialisation = Normal())
    I_imported = fix(as_turing_model(imported, 30), fixinit)().I_t

    @test all(>(0), I_imported)
    @test I_imported[end] > I_plain[end]
end

@testitem "a second ImportedCases modifier is rejected at construction" begin
    using ComposableTuringIDModels, Distributions
    gen_int = [0.2, 0.3, 0.5]
    ic1 = ImportedCases(Normal(0.0, 0.1))
    ic2 = ImportedCases(Normal(1.0, 0.1))
    # Only the first would ever be sampled, so a second is an error rather than
    # a silent drop.
    @test_throws AssertionError Renewal(gen_int, ic1, ic2;
        rt = RandomWalk(), initialisation = Normal())
    # One remains fine, alongside other modifiers.
    @test Renewal(gen_int, ic1,
        ComposableTuringIDModels.SusceptibleDepletion(1000.0);
        rt = RandomWalk(), initialisation = Normal()) isa Renewal
end

@testitem "ImportedCases bounds its slot to a prior" begin
    using ComposableTuringIDModels, Distributions
    # The slot is `PriorLike`, so a wrong-role value fails at construction
    # rather than deep inside sampling.
    @test_throws MethodError ImportedCases(1.0)
    @test_throws MethodError ImportedCases("constant")
    @test ImportedCases(Normal(0.0, 0.1)) isa ImportedCases
    @test ImportedCases([Normal(), Normal()]) isa ImportedCases
    @test ImportedCases(RandomWalk()) isa ImportedCases
end
