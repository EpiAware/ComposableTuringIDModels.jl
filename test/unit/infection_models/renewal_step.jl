# The renewal accumulation step (#48).
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
    @test renewal_foi(step, window, Rt) ≈ step(window, Rt)[end]
end

@testitem "RenewalStep with no modifiers matches the plain renewal core" begin
    using ComposableTuringIDModels: ConstantRenewalStep, RenewalStep,
                                    accumulate_scan, _renewal_init_state
    rev_gen = reverse([0.2, 0.3, 0.5])
    core = ConstantRenewalStep(rev_gen)
    step = RenewalStep(core)   # no modifiers -> plain renewal
    Rt = [1.5, 1.2, 1.0, 0.8]
    init = _renewal_init_state(step, 5.0, 0.1, length(rev_gen))
    init_core = _renewal_init_state(core, 5.0, 0.1, length(rev_gen))
    # Bare-window state and identical output: the default path is unchanged.
    @test init == init_core
    @test accumulate_scan(step, init, Rt) == accumulate_scan(core, init_core, Rt)
end

@testitem "RenewalStep with depletion matches the reference recurrence" begin
    using ComposableTuringIDModels: ConstantRenewalStep, RenewalStep,
                                    SusceptibleDepletion, accumulate_scan,
                                    _renewal_init_state
    using LinearAlgebra: dot
    rev_gen = reverse([0.2, 0.3, 0.5])
    N = 1000.0
    Rt = [1.6, 1.4, 1.2, 1.0, 0.9, 0.8]

    core = ConstantRenewalStep(rev_gen)
    step = RenewalStep(core, (SusceptibleDepletion(N),))
    init = _renewal_init_state(step, 5.0, 0.1, length(rev_gen))
    comp = accumulate_scan(step, init, Rt)

    # Explicit reference: the S/N-scaled renewal recurrence, hand-rolled. Guards
    # against drift in either the FOI core or the depletion modifier.
    function depletion_reference(Rt, window0, S0, rev_gen, N)
        window = copy(window0)
        S = S0
        ref = similar(Rt)
        for (k, r) in enumerate(Rt)
            inc = max(S / N, 1e-6) * r * dot(window, rev_gen)
            S -= inc
            window = vcat(window[2:end], inc)
            ref[k] = inc
        end
        return ref
    end
    ref = depletion_reference(Rt, init[1], init[2], rev_gen, N)

    @test comp ≈ ref
    @test init[2] ≈ N          # seeded with the full population
    # Depletion never raises incidence above the unconstrained renewal path.
    plain = accumulate_scan(core, init[1], Rt)
    @test all(comp .<= plain .+ 1e-8)
end

@testitem "RenewalStep with depletion is ForwardDiff-differentiable" begin
    using ComposableTuringIDModels: ConstantRenewalStep, RenewalStep,
                                    SusceptibleDepletion, accumulate_scan,
                                    _renewal_init_state
    using ForwardDiff
    rev_gen = reverse([0.2, 0.3, 0.5])
    N = 1000.0
    step = RenewalStep(ConstantRenewalStep(rev_gen), (SusceptibleDepletion(N),))
    # Differentiate the summed infection path with respect to the Rt path; the
    # composed recurrence must stay AD-friendly (no mutation of tracked state).
    f = function (Rt)
        init = _renewal_init_state(step, 5.0, 0.1, length(rev_gen))
        sum(accumulate_scan(step, init, Rt))
    end
    g = ForwardDiff.gradient(f, [1.6, 1.4, 1.2, 1.0])
    @test all(isfinite, g)
    @test length(g) == 4
end

@testitem "Renewal composes modifiers onto a RenewalStep" begin
    using ComposableTuringIDModels: RenewalStep, ConstantRenewalStep,
                                    SusceptibleDepletion
    data = IDData([0.2, 0.3, 0.5], exp)
    plain = Renewal(data; rt = RandomWalk())
    depleting = Renewal(data, SusceptibleDepletion(1000.0); rt = RandomWalk())
    @test plain.recurrent_step isa RenewalStep
    @test isempty(plain.recurrent_step.modifiers)
    @test depleting.recurrent_step isa RenewalStep
    @test depleting.recurrent_step.core isa ConstantRenewalStep
    @test only(depleting.recurrent_step.modifiers) isa SusceptibleDepletion
end

@testitem "Renewal susceptible depletion bends incidence below the no-depletion path" begin
    using ComposableTuringIDModels
    using DynamicPPL: fix
    using Random
    Random.seed!(48)
    data = IDData([0.2, 0.3, 0.5], exp)
    # Pin R_t high and constant so, without depletion, incidence grows unbounded.
    logR = log(2.0)
    plain = Renewal(data; rt = FixedIntercept(logR))
    depleting = Renewal(data, SusceptibleDepletion(500.0); rt = FixedIntercept(logR))
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

@testitem "Renewal with susceptible depletion samples under NUTS" tags=[:sample] begin
    using ComposableTuringIDModels, Distributions, Turing, Random
    Random.seed!(482)
    data = IDData([0.2, 0.3, 0.5], exp)
    model = IDModel(
        Renewal(data, SusceptibleDepletion(1000.0);
            rt = RandomWalk(), initialisation_prior = Normal()),
        PoissonError())
    y = as_turing_model(model, missing, 20)().generated_y_t
    # A few NUTS steps exercise the composed-step gradient path (ForwardDiff).
    chn = sample(as_turing_model(model, y, 20), NUTS(), 30; progress = false)
    @test chn !== nothing
end
