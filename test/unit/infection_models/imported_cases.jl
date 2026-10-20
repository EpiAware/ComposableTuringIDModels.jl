# Tests for the ImportedCases renewal modifier and the generic pre-scan
# modifier seam it uses.

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
    ic_id = ImportedCases(
        truncated(Normal(1.0, 0.1), 0, Inf);
        transformation = identity
    )
    @test ic_id.transformation === identity
end

@testitem "ImportedCases composes into a Renewal" begin
    using ComposableTuringIDModels, Distributions
    ic = ImportedCases(Normal(0.0, 0.1))
    gen_int = [0.2, 0.3, 0.5]
    # Via the positional modifier constructor
    r = Renewal(gen_int, ic; rt = RandomWalk(), initialisation = Normal())
    @test only(r.modifiers) isa ImportedCases
    @test r.modifiers[1].importation_rate isa Normal

    # Composed with another modifier, in the order given
    ic2 = ImportedCases(Normal(0.0, 0.1))
    r2 = Renewal(
        gen_int, ic2,
        ComposableTuringIDModels.SusceptibleDepletion(1000.0);
        rt = RandomWalk(), initialisation = Normal()
    )
    @test length(r2.modifiers) == 2
    @test r2.modifiers[1] isa ImportedCases
    @test r2.modifiers[2] isa ComposableTuringIDModels.SusceptibleDepletion
end

@testitem "ImportedCases with constant rate adds to incidence" begin
    using ComposableTuringIDModels, Distributions, Random
    using DynamicPPL: fix
    Random.seed!(1891)
    gen_int = [0.2, 0.3, 0.5]
    logR = log(1.0)  # Rt = 1, no growth

    # No importation: flat incidence.
    plain = Renewal(
        gen_int; rt = FixedIntercept(logR),
        initialisation = Normal()
    )
    fixinit = (init_incidence = log(1.0),)
    I_plain = fix(as_turing_model(plain, 30), fixinit)().I_t

    # With importation: incidence grows.
    ic = ImportedCases(Normal(0.5, 0.01))
    imported = Renewal(
        gen_int, ic; rt = FixedIntercept(logR),
        initialisation = Normal()
    )
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
    r = Renewal(
        gen_int, ic_tv; rt = RandomWalk(),
        initialisation = Normal()
    )
    out = as_turing_model(r, 20)()
    @test length(out.I_t) == 20
    @test all(isfinite, out.I_t)
    @test all(>=(0), out.I_t)

    # The importation rate is a sampled parameter, namespaced by the modifier's
    # position on the step so modifiers cannot collide.
    draw = rand(as_turing_model(r, 20))
    @test any(k -> startswith(string(k), "modifier_1.import_rates"), keys(draw))
end

@testitem "ImportedCases samples under NUTS" tags = [:sample] begin
    using ComposableTuringIDModels, Distributions, Turing, Random
    Random.seed!(1893)
    gen_int = [0.2, 0.3, 0.5]
    # An unconstrained prior is fine: the modifier's transformation makes the
    # rate positive, so the sampler cannot drive incidence negative.
    ic = ImportedCases(Normal(0.0, 0.1))
    model = IDModel(
        Renewal(
            gen_int, ic; rt = RandomWalk(),
            initialisation = Normal()
        ),
        PoissonError()
    )
    y = as_turing_model(model, missing, 20)().generated_y_t
    chn = sample(
        as_turing_model(model, y, 20), NUTS(), 30;
        progress = false
    )
    # The importation rate is sampled under the positional name the seam gives
    # it, and its draws are usable.
    draws = vec(chn[@varname(modifier_1.import_rates)])
    @test length(draws) == 30
    @test all(isfinite, draws)
end

@testitem "ImportedCases resolves to an Add before the scan" begin
    using ComposableTuringIDModels, Distributions
    import ComposableRecurrences as CR
    # The pre-scan seam draws the unconstrained slot and hands back the
    # ComposableRecurrences modifier the scan actually uses.
    resolved = as_turing_model(ImportedCases(FixedIntercept(log(2.0))), 5)()
    @test resolved isa CR.Add
    @test resolved.b isa CR.TimeVarying
    @test resolved.b.x ≈ fill(2.0, 5)

    # A bare `Distribution` is one constant, so the drawn rate stays a scalar
    # rather than being materialised as `n` copies of itself.
    constant = as_turing_model(ImportedCases(Dirac(log(3.0))), 5)()
    @test constant.b isa Number
    @test constant.b ≈ 3.0

    # Step `t` adds the rate at `t`.
    rate = [1.0, 2.0, 3.0]
    path = as_turing_model(ImportedCases([Dirac(log(x)) for x in rate]), 3)()
    @test CR.Recurrence([0.0]; modifiers = (path,))(; stop = 3) ≈ rate
end

@testitem "SusceptibleDepletion resolves to a floored Depletion" begin
    using ComposableTuringIDModels, Distributions, Random
    using ComposableTuringIDModels: SusceptibleDepletion
    import ComposableRecurrences as CR
    using DynamicPPL: fix
    # `SusceptibleDepletion` samples nothing, and resolves to the
    # ComposableRecurrences depletion with the floored form.
    dep = SusceptibleDepletion(1000.0)
    resolved = as_turing_model(dep, 10)()
    @test resolved isa CR.Depletion
    @test resolved.form isa CR.Floor
    # One population per stratum is read per stratum.
    per_stratum = as_turing_model(SusceptibleDepletion([1.0, 2.0]), (2, 10))()
    @test per_stratum isa CR.Depletion

    # And a depleting renewal is unchanged end to end by the seam.
    gen_int = [0.2, 0.3, 0.5]
    depleting = Renewal(
        gen_int, dep; rt = FixedIntercept(log(1.5)),
        initialisation = Normal()
    )
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
    import ComposableRecurrences as CR
    using DynamicPPL: fix
    Random.seed!(1895)

    # A modifier defined here rather than in the package: it carries a prior,
    # draws it through the seam, and resolves to a ComposableRecurrences
    # modifier scaling each step's incidence. Nothing in the renewal model
    # knows what it is.
    struct ScaledIncidence{P} <: AbstractRenewalModifier
        log_scale::P
    end
    Turing.@model function ComposableTuringIDModels.as_turing_model(
            mod::ScaledIncidence, n
        )
        log_scale ~ as_turing_submodel(mod.log_scale, n; prefix = true)
        return CR.Transform((v, c) -> c * v, exp(_at(log_scale, 1)))
    end

    gen_int = [0.2, 0.3, 0.5]
    fixinit = (init_incidence = log(1.0),)
    args = (; rt = FixedIntercept(log(1.0)), initialisation = Normal())
    # The scan uses the resolved value: halving each step's incidence leaves
    # the series below the unscaled one.
    plain = Renewal(gen_int; args...)
    scaled = Renewal(
        gen_int, ScaledIncidence(FixedIntercept(log(0.5)));
        args...
    )
    I_plain = fix(as_turing_model(plain, 30), fixinit)().I_t
    I_scaled = fix(as_turing_model(scaled, 30), fixinit)().I_t
    @test all(>(0), I_scaled)
    @test I_scaled[end] < I_plain[end]

    # Composed alongside the package's own modifiers it is namespaced by its
    # position like any other, and only the modifiers carrying priors sample.
    r = Renewal(
        gen_int, SusceptibleDepletion(500.0),
        ScaledIncidence(Normal(-0.1, 0.05)),
        ImportedCases(Normal(-1.0, 0.1)); rt = RandomWalk(),
        initialisation = Normal()
    )
    ks = string.(collect(keys(rand(as_turing_model(r, 15)))))
    @test any(k -> startswith(k, "modifier_2.log_scale"), ks)
    @test any(k -> startswith(k, "modifier_3.import_rates"), ks)
    @test !any(k -> startswith(k, "modifier_1."), ks)
    out = as_turing_model(r, 15)()
    @test length(out.I_t) == 15
    @test all(>(0), out.I_t)
end

@testitem "several sampling modifiers compose without colliding" begin
    using ComposableTuringIDModels, Distributions, Random
    Random.seed!(1894)
    gen_int = [0.2, 0.3, 0.5]
    # Two importation streams, each with its own prior: the seam namespaces
    # each modifier by its position, so both are sampled.
    r = Renewal(
        gen_int, ImportedCases(Normal(-1.0, 0.1)),
        ImportedCases(RandomWalk()); rt = RandomWalk(),
        initialisation = Normal()
    )
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
    after = Renewal(
        gen_int, SusceptibleDepletion(50.0), ic();
        rt = FixedIntercept(log(1.0)), initialisation = Normal()
    )
    # Importation first: the imports are part of the incidence the pool
    # depletes, so they are scaled down as susceptibles run out.
    before = Renewal(
        gen_int, ic(), SusceptibleDepletion(50.0);
        rt = FixedIntercept(log(1.0)), initialisation = Normal()
    )
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

    plain = Renewal(
        gen_int; rt = FixedIntercept(logR),
        initialisation = Normal()
    )
    I_plain = fix(as_turing_model(plain, 30), fixinit)().I_t

    # A strongly negative *unconstrained* rate: the transformation maps it to
    # exp(-2) ≈ 0.135 imports per step, so it still imports.
    imported = Renewal(
        gen_int, ImportedCases(FixedIntercept(-2.0));
        rt = FixedIntercept(logR), initialisation = Normal()
    )
    I_imported = fix(as_turing_model(imported, 30), fixinit)().I_t

    @test all(>(0), I_imported)
    @test I_imported[end] > I_plain[end]

    # Nothing clamps the incidence: taking the same rate on its natural scale
    # subtracts 2.0 per step, which is a misspecified model and shows as one
    # rather than being silently floored at a small positive value.
    negative = Renewal(
        gen_int,
        ImportedCases(FixedIntercept(-2.0); transformation = identity);
        rt = FixedIntercept(logR), initialisation = Normal()
    )
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

@testitem "a per-element importation prior gives a length-n rate" begin
    using ComposableTuringIDModels, Distributions, Random
    using DynamicPPL: fix
    Random.seed!(1898)
    n = 8
    # A `Vector{<:Distribution}` is the third documented slot shape: one prior
    # per time, so the seam draws a length-`n` path rather than a constant.
    resolved = as_turing_model(
        ImportedCases([Dirac(log(float(t))) for t in 1:n]), n
    )()
    @test resolved.b.x ≈ collect(1.0:n)

    # It drives a renewal like any other slot shape, and the length is checked
    # against the series rather than silently recycled.
    gen_int = [0.2, 0.3, 0.5]
    fixinit = (init_incidence = log(1.0),)
    args = (; rt = FixedIntercept(log(1.0)), initialisation = Normal())
    plain = Renewal(gen_int; args...)
    seeded = Renewal(
        gen_int,
        ImportedCases([Normal(-1.0, 0.1) for _ in 1:n]); args...
    )
    I_plain = fix(as_turing_model(plain, n), fixinit)().I_t
    I_seeded = fix(as_turing_model(seeded, n), fixinit)().I_t
    @test all(>(0), I_seeded)
    @test I_seeded[end] > I_plain[end]

    mismatched = Renewal(
        gen_int,
        ImportedCases([Normal(-1.0, 0.1) for _ in 1:(n + 1)]); args...
    )
    @test_throws Exception as_turing_model(mismatched, n)()
end

@testitem "a renewal with importation is ForwardDiff-differentiable" begin
    using ComposableTuringIDModels, Distributions, ForwardDiff
    using ComposableTuringIDModels: _renewal_recurrence, SusceptibleDepletion
    import ComposableRecurrences as CR
    using LinearAlgebra: I
    gen_int = [0.2, 0.3, 0.5]
    n = 20
    window = ones(length(gen_int))
    depletion = as_turing_model(SusceptibleDepletion(1000.0), n)()
    # The recurrence runs over resolved modifiers, so the gradient flows
    # through the sampled importation rate exactly as it does through `Rt`.
    function total(θ)
        imports = CR.Add(CR.TimeVarying(fill(θ[2], n)))
        r = _renewal_recurrence(gen_int, I, (depletion, imports))
        return sum(r(fill(θ[1], n); history = window))
    end
    g = ForwardDiff.gradient(total, [1.1, 0.5])
    @test all(isfinite, g)
    @test all(>(0), g)
end
