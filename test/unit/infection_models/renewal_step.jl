# The renewal accumulation step.
#
# `RenewalStep` is the renewal step: an internal `ConstantRenewalStep` force-of-
# infection core with modifiers composing on top. With no modifiers it matches the
# plain renewal exactly; with a `SusceptibleDepletion` modifier it is a renewal
# with a fixed population. The `Renewal` helper composes modifiers onto the step.

@testitem "renewal_foi matches the ConstantRenewalStep recurrence" begin
    using ComposableTuringIDModels: ConstantRenewalStep, renewal_foi
    rev_gen = reverse([0.2, 0.3, 0.5])
    step = ConstantRenewalStep(rev_gen)
    window = [1.0, 2.0, 3.0]
    Rt = 1.4
    # The extracted FOI is exactly the new-incidence term the step commits.
    state = (; val = last(window), window = window)
    @test renewal_foi(step, window, Rt) ≈ step(state, Rt).val
end

@testitem "RenewalStep with no modifiers matches the plain renewal core" begin
    using ComposableTuringIDModels: ConstantRenewalStep, RenewalStep,
        accumulate_scan, renewal_init_state
    rev_gen = reverse([0.2, 0.3, 0.5])
    core = ConstantRenewalStep(rev_gen)
    step = RenewalStep(core)   # no modifiers -> plain renewal
    Rt = [1.5, 1.2, 1.0, 0.8]
    init = renewal_init_state(step, 5.0, 0.1, length(rev_gen))
    init_core = renewal_init_state(core, 5.0, 0.1, length(rev_gen))
    # Identical seed state and identical output: the default path is unchanged.
    @test init == init_core
    @test accumulate_scan(step, init, Rt) == accumulate_scan(core, init, Rt)
end

@testitem "RenewalStep with depletion matches the reference recurrence" begin
    using ComposableTuringIDModels: ConstantRenewalStep, RenewalStep,
        SusceptibleDepletion, accumulate_scan,
        renewal_init_state
    using LinearAlgebra: dot
    rev_gen = reverse([0.2, 0.3, 0.5])
    N = 1000.0
    Rt = [1.6, 1.4, 1.2, 1.0, 0.9, 0.8]

    core = ConstantRenewalStep(rev_gen)
    step = RenewalStep(core, (SusceptibleDepletion(N),))
    init = renewal_init_state(step, 5.0, 0.1, length(rev_gen))
    comp = accumulate_scan(step, init, Rt)

    # Explicit reference: the S/N-scaled renewal recurrence, hand-rolled. Guards
    # against drift in either the FOI core or the depletion modifier.
    function depletion_reference(Rt, window0, S0, rev_gen, N)
        window = copy(window0)
        S = S0
        ref = similar(Rt)
        for (k, r) in enumerate(Rt)
            inc = max(S / N, 1.0e-6) * r * dot(window, rev_gen)
            S -= inc
            window = vcat(window[2:end], inc)
            ref[k] = inc
        end
        return ref
    end
    ref = depletion_reference(Rt, init.window, only(init.substates), rev_gen, N)

    @test comp ≈ ref
    @test only(init.substates) ≈ N          # seeded with the full population
    # Depletion never raises incidence above the unconstrained renewal path.
    plain = accumulate_scan(
        core, (; val = last(init.window), window = init.window), Rt
    )
    @test all(comp .<= plain .+ 1.0e-8)
end

@testitem "RenewalStep with depletion is ForwardDiff-differentiable" begin
    using ComposableTuringIDModels: ConstantRenewalStep, RenewalStep,
        SusceptibleDepletion, accumulate_scan,
        renewal_init_state
    using ForwardDiff
    rev_gen = reverse([0.2, 0.3, 0.5])
    N = 1000.0
    step = RenewalStep(ConstantRenewalStep(rev_gen), (SusceptibleDepletion(N),))
    # Differentiate the summed infection path with respect to the Rt path; the
    # composed recurrence must stay AD-friendly (no mutation of tracked state).
    f = function (Rt)
        init = renewal_init_state(step, 5.0, 0.1, length(rev_gen))
        sum(accumulate_scan(step, init, Rt))
    end
    g = ForwardDiff.gradient(f, [1.6, 1.4, 1.2, 1.0])
    @test all(isfinite, g)
    @test length(g) == 4
end

@testitem "Renewal composes modifiers onto a RenewalStep" begin
    using ComposableTuringIDModels: RenewalStep, ConstantRenewalStep,
        SusceptibleDepletion, Renewal, RandomWalk
    gen_int = [0.2, 0.3, 0.5]
    plain = Renewal(gen_int; rt = RandomWalk())
    depleting = Renewal(
        gen_int, SusceptibleDepletion(1000.0);
        rt = RandomWalk()
    )
    # With no modifiers the step is the bare force-of-infection core, whichever
    # constructor built it.
    @test plain.recurrent_step isa ConstantRenewalStep
    @test isempty(plain.modifiers)
    @test Renewal(; generation_time = gen_int).recurrent_step isa
        ConstantRenewalStep
    @test depleting.recurrent_step isa RenewalStep
    @test depleting.recurrent_step.core isa ConstantRenewalStep
    @test only(depleting.recurrent_step.modifiers) isa SusceptibleDepletion
end

@testitem "Renewal susceptible depletion bends incidence below the no-depletion path" begin
    using ComposableTuringIDModels
    using DynamicPPL: fix
    using Random
    Random.seed!(48)
    gen_int = [0.2, 0.3, 0.5]
    # Pin R_t high and constant so, without depletion, incidence grows unbounded.
    logR = log(2.0)
    plain = Renewal(gen_int; rt = FixedIntercept(logR))
    depleting = Renewal(
        gen_int, SusceptibleDepletion(500.0);
        rt = FixedIntercept(logR)
    )
    fixinit = (init_incidence = log(1.0),)
    I_plain = fix(as_turing_model(plain, 40), fixinit)().I_t
    I_dep = fix(as_turing_model(depleting, 40), fixinit)().I_t
    @test all(isfinite, I_dep)
    @test all(>=(0), I_dep)
    # Depletion holds late incidence far below the unlimited-growth path.
    @test I_dep[end] < I_plain[end]
    # With a finite population the epidemic peaks and turns over, so late
    # incidence falls below the peak (impossible under unbounded growth).
    @test I_dep[end] < maximum(I_dep)
end

@testitem "Renewal with susceptible depletion samples under NUTS" tags = [:sample] begin
    using ComposableTuringIDModels, Distributions, Turing, Random
    Random.seed!(482)
    gen_int = [0.2, 0.3, 0.5]
    model = IDModel(
        Renewal(
            gen_int, SusceptibleDepletion(1000.0);
            rt = RandomWalk(), initialisation = Normal()
        ),
        PoissonError()
    )
    y = as_turing_model(model, missing, 20)().generated_y_t
    # A few NUTS steps exercise the composed-step gradient path (ForwardDiff).
    chn = sample(as_turing_model(model, y, 20), NUTS(), 30; progress = false)
    @test chn !== nothing
end

@testitem "Renewal combines a continuous generation time with a modifier" begin
    using ComposableTuringIDModels, Distributions
    using ComposableTuringIDModels: RenewalStep
    # A continuous generation time discretised by the constructor, carrying a
    # modifier: the two compose in either constructor form.
    renewal = Renewal(
        Gamma(2, 1.5), SusceptibleDepletion(1000.0);
        D_gen = 15.0, rt = RandomWalk(), initialisation = Normal(log(50), 0.2)
    )
    @test renewal.gen_int isa AbstractVector
    @test sum(renewal.gen_int) ≈ 1
    @test renewal.recurrent_step isa RenewalStep
    @test only(renewal.recurrent_step.modifiers) isa SusceptibleDepletion
    # One discretisation path: the interval matches the modifier-free model's,
    # so a modifier does not force a second discretisation outside the package.
    base = Renewal(; generation_time = Gamma(2, 1.5), D_gen = 15.0)
    @test renewal.gen_int == base.gen_int
    # The keyword form takes the same modifiers.
    kw = Renewal(;
        generation_time = Gamma(2, 1.5), D_gen = 15.0,
        modifiers = SusceptibleDepletion(1000.0), rt = RandomWalk(),
        initialisation = Normal(log(50), 0.2)
    )
    @test kw.gen_int == renewal.gen_int
    @test only(kw.modifiers) isa SusceptibleDepletion
    out = as_turing_model(renewal, 20)()
    @test length(out.I_t) == 20
    @test all(isfinite, out.I_t)
end

@testitem "Renewal composes a modifier onto an inferred generation interval" begin
    using ComposableTuringIDModels, Distributions, Random
    using DynamicPPL: fix
    gen = UncertainDelay(
        LogNormal,
        [Normal(1.9, 0.2), truncated(Normal(0.5, 0.2), 0, Inf)]; D = 14.0
    )
    # No interval is baked, so the modifiers are held on the model and the step
    # is rebuilt with them per draw.
    depleting = Renewal(
        gen, SusceptibleDepletion(50.0);
        rt = FixedIntercept(log(2.0)), initialisation = Normal(log(10), 0.1)
    )
    plain = Renewal(
        gen; rt = FixedIntercept(log(2.0)),
        initialisation = Normal(log(10), 0.1)
    )
    @test isnothing(depleting.recurrent_step)
    @test only(depleting.modifiers) isa SusceptibleDepletion
    @test isempty(plain.modifiers)
    fixinit = (init_incidence = log(10.0),)
    Random.seed!(51)
    depleted = fix(as_turing_model(depleting, 25), fixinit)().I_t
    Random.seed!(51)
    undepleted = fix(as_turing_model(plain, 25), fixinit)().I_t
    @test all(isfinite, depleted)
    # A small population bites: depletion holds incidence below the plain path.
    @test last(depleted) < last(undepleted)
end

@testitem "the modifiers keyword rejects a value that is not a modifier" begin
    using ComposableTuringIDModels, Distributions
    gen_int = [0.2, 0.3, 0.5]
    # The positional form is constrained by its signature; the keyword form has
    # to say so itself rather than failing inside the step's seam.
    @test_throws "AbstractRenewalModifier" Renewal(;
        generation_time = gen_int, modifiers = (Normal(),)
    )
    @test_throws "AbstractRenewalModifier" Renewal(;
        generation_time = gen_int, modifiers = [SusceptibleDepletion(1.0e3), 1]
    )
    # A single modifier, a tuple and a vector of them are all accepted.
    for mods in (
            SusceptibleDepletion(1.0e3), (SusceptibleDepletion(1.0e3),),
            [SusceptibleDepletion(1.0e3)],
        )
        r = Renewal(; generation_time = gen_int, modifiers = mods)
        @test only(r.modifiers) isa SusceptibleDepletion
    end
end

@testitem "Renewal widens its rt slot on both construction paths" begin
    using ComposableTuringIDModels, Distributions, Random
    using ComposableTuringIDModels: path_prior
    Random.seed!(344)

    # The positional constructor must widen a bare `Distribution` in the `rt`
    # PATH slot exactly as the keyword one does. Left to Julia's generated
    # constructor it does not: the slot draws a single scalar, which reaches
    # the scan as an `R_t` with no time axis.
    gi = [0.2, 0.3, 0.5]
    kw = Renewal(;
        generation_time = gi, rt = Normal(0, 0.1), initialisation = Normal()
    )
    pos = Renewal(
        gi, exp, Normal(0, 0.1), Normal(), kw.recurrent_step, kw.mixing, ()
    )

    @test typeof(pos.rt) == typeof(kw.rt)
    @test pos.rt isa Intercept
    n = 20
    @test length(as_turing_model(kw, n)().I_t) == n
    @test length(as_turing_model(pos, n)().I_t) == n
    # Rebuilding from stored fields is a no-op, so `Accessors.@set` is safe.
    @test path_prior(kw.rt) === kw.rt
    # Widening runs before `new`, so the slot's role guard still fires.
    @test_throws TypeError Renewal(
        gi, exp, PoissonError(), Normal(), kw.recurrent_step, kw.mixing, ()
    )
end

@testitem "the scan records the noise-free expectation alongside the draw" begin
    using ComposableTuringIDModels: ConstantRenewalStep, RenewalStep,
        SusceptibleDepletion, accumulate_scan, renewal_init_state,
        get_expected_state
    rev_gen = reverse([0.2, 0.3, 0.5])
    core = ConstantRenewalStep(rev_gen)
    Rt = [1.6, 1.4, 1.2, 1.0, 0.9, 0.8]

    # A step with no modifier transforms nothing, so its expectation series is
    # its committed series, on the bare core and on the step wrapping it.
    for plain in (core, RenewalStep(core, ()))
        init = renewal_init_state(plain, 5.0, 0.1, length(rev_gen))
        result = accumulate(plain, Rt; init = init)
        @test get_expected_state(plain, init, result) ==
            accumulate_scan(plain, init, Rt)
    end

    good = RenewalStep(core, (SusceptibleDepletion(1000.0),))
    init = renewal_init_state(good, 5.0, 0.1, length(rev_gen))
    result = accumulate(good, Rt; init = init)
    exp_val = get_expected_state(good, init, result)
    @test length(exp_val) == length(Rt)
    @test all(exp_val .>= 0)
    # The expectation is the pre-modifier incidence. Susceptible depletion
    # scales the committed draw down by the susceptible fraction, so the draw
    # stays at or below the expectation and equals it only while the pool is
    # full.
    @test all(accumulate_scan(good, init, Rt) .<= exp_val .+ 1.0e-8)
    @test maximum(accumulate_scan(good, init, Rt)) < maximum(exp_val)
end

@testitem "get_expected_state names a step that has no expectation" begin
    using ComposableTuringIDModels: ARStep, get_expected_state
    # Only the renewal steps carry an expectation. `get_expected_state` is
    # public with a generic name, so a caller reaching it on another
    # accumulation step gets the step named rather than a bare `MethodError`.
    step = ARStep([0.5])
    state = (; val = 1.0, window = [1.0])
    @test_throws "ARStep" get_expected_state(step, state, [state])
    @test_throws "has no noise-free expectation" get_expected_state(
        step, state, [state]
    )
end

@testitem "a plain renewal exposes exp_I_t identical to I_t" begin
    using ComposableTuringIDModels, Distributions, Random
    Random.seed!(481)
    gen_int = [0.2, 0.3, 0.5]
    model = Renewal(;
        generation_time = gen_int, rt = RandomWalk(),
        initialisation = Normal()
    )
    out = as_turing_model(model, 30)()
    @test haskey(out, :exp_I_t)
    @test out.exp_I_t ≈ out.I_t
    # The stratified shape reports one expectation per stratum.
    strat = Renewal(;
        generation_time = gen_int,
        rt = Stratify(RandomWalk(), FixedIntercept(0.0)),
        initialisation = Normal()
    )
    s = as_turing_model(strat, (2, 20))()
    @test size(s.exp_I_t) == (2, 20)
    @test s.exp_I_t ≈ s.I_t
end

@testitem "exp_I_t comes back from a chain in a composed model" tags = [:sample] begin
    using ComposableTuringIDModels, Distributions, Turing, Random
    Random.seed!(482)
    gen_int = [0.2, 0.3, 0.5]
    model = IDModel(
        Renewal(
            gen_int, SusceptibleDepletion(1000.0);
            rt = RandomWalk(), initialisation = Normal()
        ),
        PoissonError()
    )
    y = as_turing_model(model, missing, 20)().generated_y_t
    chn = sample(as_turing_model(model, y, 20), NUTS(), 15; progress = false)
    # `exp_I_t` is tracked via `:=` like `ρ`, so it comes out of the chain, one
    # series per draw.
    draws = vec(chn[:exp_I_t])
    @test length(draws) == 15
    @test all(d -> length(d) == 20, draws)
    @test all(isfinite, reduce(vcat, draws))
end

@testitem "exp_I_t is the incidence ImportedCases adds on top of" begin
    using ComposableTuringIDModels, Distributions, Random
    Random.seed!(483)
    # An importation modifier adds a spike on top of the force of infection, so
    # the committed draw exceeds the renewal expectation wherever it bites. The
    # portal is `Renewal`, whose `as_turing_model` resolves the modifier's
    # priors through the pre-scan seam before accumulating.
    gen_int = [0.2, 0.3, 0.5]
    model = Renewal(
        gen_int, ImportedCases(FixedIntercept(log(2.0)));
        rt = RandomWalk(), initialisation = Normal()
    )
    out = as_turing_model(model, 30)()
    @test haskey(out, :exp_I_t)
    @test all(out.I_t .>= out.exp_I_t)
    @test mean(out.I_t .- out.exp_I_t) > 0.0
end
