# Tests for the ImportedCases renewal modifier and the generic pre-scan
# modifier seam it uses (issue #189).

@testitem "ImportedCases struct carries the importation prior" begin
    using ComposableTuringIDModels, Distributions
    # Constant unknown rate (Distribution)
    ic = ImportedCases(Normal(0.0, 0.1))
    @test ic.importation_rate isa Normal
    @test ic.transformation === exp
    @test ic isa ComposableTuringIDModels.AbstractRenewalModifier

    # Time-varying rate (latent process)
    ic_tv = ImportedCases(RandomWalk())
    @test ic_tv.importation_rate isa RandomWalk

    # The map onto the positive rate is the modifier's own, not the host's.
    ic_id = ImportedCases(truncated(Normal(1.0, 0.1), 0, Inf);
        transformation = identity)
    @test ic_id.transformation === identity
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

    # Composed with another modifier, in the order given
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

@testitem "ImportedCases draws a time-varying rate through the seam" begin
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

    # The importation rate is a sampled parameter, namespaced by the modifier's
    # position on the step so modifiers cannot collide.
    draw = rand(as_turing_model(r, 20))
    @test any(k -> startswith(string(k), "modifier_1.import_rates"), keys(draw))
end

@testitem "ImportedCases samples under NUTS" tags=[:sample] begin
    using ComposableTuringIDModels, Distributions, Turing, Random
    Random.seed!(1893)
    gen_int = [0.2, 0.3, 0.5]
    # An unconstrained prior is fine: the modifier's transformation makes the
    # rate positive, so the sampler cannot drive incidence negative.
    ic = ImportedCases(Normal(0.0, 0.1))
    model = IDModel(
        Renewal(gen_int, ic; rt = RandomWalk(),
            initialisation = Normal()),
        PoissonError())
    y = as_turing_model(model, missing, 20)().generated_y_t
    chn = sample(as_turing_model(model, y, 20), NUTS(), 30;
        progress = false)
    # The importation rate is sampled under the positional name the seam gives
    # it, and its draws are usable.
    draws = vec(chn[@varname(modifier_1.import_rates)])
    @test length(draws) == 30
    @test all(isfinite, draws)
end

@testitem "ImportedCases resolves to an ImportedRate before the scan" begin
    using ComposableTuringIDModels, Distributions
    using ComposableTuringIDModels: ImportedRate, modifier_init_state,
                                    apply_modifier
    # The pre-scan seam draws the unconstrained slot and hands back the
    # modifier the scan actually uses.
    resolved = as_turing_model(ImportedCases(FixedIntercept(log(2.0))), 5)()
    @test resolved isa ImportedRate
    @test resolved.rate ≈ fill(2.0, 5)

    # A bare `Distribution` is one constant, so the drawn rate stays a scalar
    # rather than being materialised as `n` copies of itself.
    constant = as_turing_model(ImportedCases(Dirac(log(3.0))), 5)()
    @test constant.rate isa Number
    @test constant.rate ≈ 3.0
    @test first(apply_modifier(constant, 1.0, 7)) ≈ 4.0

    # The resolved modifier is a plain scan modifier: its substate is the step
    # counter, so step `t` adds the rate at `t`.
    inc, s0 = 3.5, modifier_init_state(resolved)
    @test s0 == 0
    inc1, s1 = apply_modifier(resolved, inc, s0)
    @test inc1 ≈ inc + 2.0
    @test s1 == 1
    inc2, s2 = apply_modifier(resolved, inc, s1)
    @test inc2 ≈ inc + 2.0
    @test s2 == 2
end

@testitem "a modifier with no pre-scan contribution resolves to itself" begin
    using ComposableTuringIDModels, Distributions, Random
    using ComposableTuringIDModels: ConstantRenewalStep, RenewalStep,
                                    SusceptibleDepletion
    using DynamicPPL: fix
    core = ConstantRenewalStep(reverse([0.2, 0.3, 0.5]))
    # `SusceptibleDepletion` implements nothing beyond the scan interface, so
    # the default seam returns it unchanged.
    dep = SusceptibleDepletion(1000.0)
    @test as_turing_model(dep, 10)() === dep
    # The bare force-of-infection core and a modifier-free step likewise.
    @test as_turing_model(core, 10)() === core
    resolved = as_turing_model(RenewalStep(core, (dep,)), 10)()
    @test resolved isa RenewalStep
    @test only(resolved.modifiers) === dep

    # And a depleting renewal is unchanged end to end by the seam.
    gen_int = [0.2, 0.3, 0.5]
    depleting = Renewal(gen_int, dep; rt = FixedIntercept(log(1.5)),
        initialisation = Normal())
    I_t = fix(as_turing_model(depleting, 30), (init_incidence = 0.0,))().I_t
    @test length(I_t) == 30
    @test all(>(0), I_t)
    @test sum(I_t) < 1000.0  # the pool bounds the epidemic
    # No importation parameters are sampled when no modifier declares any.
    draw = rand(as_turing_model(depleting, 10))
    @test !any(k -> occursin("import", string(k)), keys(draw))
end

@testitem "the pre-scan seam takes any modifier, not just ImportedCases" begin
    using ComposableTuringIDModels, Distributions, Random, Turing
    using ComposableTuringIDModels: AbstractRenewalModifier,
                                    SusceptibleDepletion, _at
    using DynamicPPL: fix
    Random.seed!(1895)

    # A modifier defined here rather than in the package: it carries a prior,
    # draws it through the seam, and resolves to a plain scan modifier. Nothing
    # in the renewal model knows what it is.
    struct ScaledIncidence{P} <: AbstractRenewalModifier
        log_scale::P
    end
    struct FixedScale{T} <: AbstractRenewalModifier
        scale::T
    end
    ComposableTuringIDModels.modifier_init_state(::FixedScale) = 0
    function ComposableTuringIDModels.apply_modifier(
            mod::FixedScale, incidence, t)
        return mod.scale * incidence, t + 1
    end
    Turing.@model function ComposableTuringIDModels.as_turing_model(
            mod::ScaledIncidence, n)
        log_scale ~ as_turing_submodel(mod.log_scale, n; prefix = true)
        return FixedScale(exp(_at(log_scale, 1)))
    end

    gen_int = [0.2, 0.3, 0.5]
    fixinit = (init_incidence = log(1.0),)
    args = (; rt = FixedIntercept(log(1.0)), initialisation = Normal())
    # The scan uses the resolved value: halving each step's incidence leaves
    # the series below the unscaled one.
    plain = Renewal(gen_int; args...)
    scaled = Renewal(gen_int, ScaledIncidence(FixedIntercept(log(0.5)));
        args...)
    I_plain = fix(as_turing_model(plain, 30), fixinit)().I_t
    I_scaled = fix(as_turing_model(scaled, 30), fixinit)().I_t
    @test all(>(0), I_scaled)
    @test I_scaled[end] < I_plain[end]

    # Composed alongside the package's own modifiers it is namespaced by its
    # position like any other, and only the modifiers carrying priors sample.
    r = Renewal(gen_int, SusceptibleDepletion(500.0),
        ScaledIncidence(Normal(-0.1, 0.05)),
        ImportedCases(Normal(-1.0, 0.1)); rt = RandomWalk(),
        initialisation = Normal())
    ks = string.(collect(keys(rand(as_turing_model(r, 15)))))
    @test any(k -> startswith(k, "modifier_2.log_scale"), ks)
    @test any(k -> startswith(k, "modifier_3.import_rates"), ks)
    @test !any(k -> startswith(k, "modifier_1."), ks)
    out = as_turing_model(r, 15)()
    @test length(out.I_t) == 15
    @test all(>(0), out.I_t)
end

@testitem "an unresolved modifier says so instead of a MethodError" begin
    using ComposableTuringIDModels, Distributions
    using ComposableTuringIDModels: modifier_init_state, apply_modifier
    # `ImportedCases` carries priors and has no scan interface of its own, so
    # scanning a hand-built step that was never resolved names the problem.
    ic = ImportedCases(Normal(0.0, 0.1))
    @test_throws "has no scan interface" modifier_init_state(ic)
    @test_throws "has no scan interface" apply_modifier(ic, 1.0, 0)
end

@testitem "several sampling modifiers compose without colliding" begin
    using ComposableTuringIDModels, Distributions, Random
    Random.seed!(1894)
    gen_int = [0.2, 0.3, 0.5]
    # Two importation streams, each with its own prior: the seam namespaces
    # each modifier by its position, so both are sampled.
    r = Renewal(gen_int, ImportedCases(Normal(-1.0, 0.1)),
        ImportedCases(RandomWalk()); rt = RandomWalk(),
        initialisation = Normal())
    draw = rand(as_turing_model(r, 15))
    ks = string.(collect(keys(draw)))
    @test any(k -> startswith(k, "modifier_1.import_rates"), ks)
    @test any(k -> startswith(k, "modifier_2.import_rates"), ks)
    out = as_turing_model(r, 15)()
    @test length(out.I_t) == 15
    @test all(>(0), out.I_t)
end

@testitem "modifier order sets how importation composes with depletion" begin
    using ComposableTuringIDModels, Distributions
    using ComposableTuringIDModels: SusceptibleDepletion
    using DynamicPPL: fix
    gen_int = [0.2, 0.3, 0.5]
    fixinit = (init_incidence = log(1.0),)
    ic() = ImportedCases(FixedIntercept(log(5.0)))
    # Depletion first: imports are added to the depleted incidence, so they are
    # not scaled by the susceptible fraction.
    after = Renewal(gen_int, SusceptibleDepletion(50.0), ic();
        rt = FixedIntercept(log(1.0)), initialisation = Normal())
    # Importation first: the imports are part of the incidence the pool
    # depletes, so they are scaled down as susceptibles run out.
    before = Renewal(gen_int, ic(), SusceptibleDepletion(50.0);
        rt = FixedIntercept(log(1.0)), initialisation = Normal())
    I_after = fix(as_turing_model(after, 40), fixinit)().I_t
    I_before = fix(as_turing_model(before, 40), fixinit)().I_t
    @test all(>(0), I_after)
    @test all(>(0), I_before)
    @test I_after[end] > I_before[end]
end

@testitem "a negative unconstrained rate still imports" begin
    using ComposableTuringIDModels, Distributions
    using DynamicPPL: fix
    gen_int = [0.2, 0.3, 0.5]
    logR = log(1.0)  # Rt = 1, no growth
    fixinit = (init_incidence = log(1.0),)

    plain = Renewal(gen_int; rt = FixedIntercept(logR),
        initialisation = Normal())
    I_plain = fix(as_turing_model(plain, 30), fixinit)().I_t

    # A strongly negative *unconstrained* rate: the transformation maps it to
    # exp(-2) ≈ 0.135 imports per step, so it still imports.
    imported = Renewal(gen_int, ImportedCases(FixedIntercept(-2.0));
        rt = FixedIntercept(logR), initialisation = Normal())
    I_imported = fix(as_turing_model(imported, 30), fixinit)().I_t

    @test all(>(0), I_imported)
    @test I_imported[end] > I_plain[end]

    # Nothing clamps the incidence: taking the same rate on its natural scale
    # subtracts 2.0 per step, which is a misspecified model and shows as one
    # rather than being silently floored at a small positive value.
    negative = Renewal(gen_int,
        ImportedCases(FixedIntercept(-2.0); transformation = identity);
        rt = FixedIntercept(logR), initialisation = Normal())
    I_negative = fix(as_turing_model(negative, 30), fixinit)().I_t
    @test any(<(0), I_negative)
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

@testitem "a renewal with importation is ForwardDiff-differentiable" begin
    using ComposableTuringIDModels, Distributions, ForwardDiff
    using ComposableTuringIDModels: ImportedRate, RenewalStep,
                                    ConstantRenewalStep, SusceptibleDepletion,
                                    accumulate_scan, _renewal_init_state
    gen_int = [0.2, 0.3, 0.5]
    core = ConstantRenewalStep(reverse(gen_int))
    n = 20
    # The scan runs over resolved modifiers, so the gradient flows through the
    # sampled importation rate exactly as it does through `Rt`.
    function total(θ)
        step = RenewalStep(core,
            (SusceptibleDepletion(1000.0), ImportedRate(fill(θ[2], n))))
        init = _renewal_init_state(step, 1.0, 0.0, length(gen_int))
        return sum(accumulate_scan(step, init, fill(θ[1], n)))
    end
    g = ForwardDiff.gradient(total, [1.1, 0.5])
    @test all(isfinite, g)
    @test all(>(0), g)
end
