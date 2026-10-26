# The renewal recurrence.
#
# `Renewal` runs a ComposableRecurrences `Recurrence` whose kernel is the
# generation interval and whose gain is R_t, seeded on an exponential growth
# path. Its modifiers resolve to ComposableRecurrences modifiers before the
# scan, so a `SusceptibleDepletion` is a floored `Depletion`. These tests check
# the recurrence it builds against hand-written renewal loops.

@testsnippet RenewalReference begin
    using LinearAlgebra: dot
    # I_t = R_t Σ_i g_i I_{t-i}, from a window of past incidence, oldest
    # first, with optional floored susceptible depletion from a pool of `N`.
    function renewal_reference(g, window0, Rt; N = nothing)
        window = copy(window0)
        S = N
        out = similar(Rt, promote_type(eltype(Rt), eltype(window0)))
        for (k, r) in enumerate(Rt)
            inc = r * dot(window, reverse(g))
            if N !== nothing
                inc = max(S / N, 1.0e-6) * inc
                S -= inc
            end
            window = vcat(window[2:end], inc)
            out[k] = inc
        end
        return out
    end
end

@testitem "the renewal recurrence matches a hand-written renewal loop" setup = [RenewalReference] begin
    using ComposableTuringIDModels: _renewal_recurrence, _renewal_history,
        R_to_r
    using LinearAlgebra: I
    g = [0.2, 0.3, 0.5]
    Rt = [1.5, 1.2, 1.0, 0.8, 1.1, 1.3]
    window = _renewal_history(g, 5.0, first(Rt))
    # The window decays into the seed at the growth rate implied by R_1.
    r = R_to_r(first(Rt), g)
    @test window ≈ 5.0 .* exp.(-r .* [2, 1, 0])
    I_t = _renewal_recurrence(g, I, ())(Rt; history = window)
    @test I_t ≈ renewal_reference(g, window, Rt)
end

@testitem "the renewal recurrence with depletion matches a hand-written loop" setup = [RenewalReference] begin
    using ComposableTuringIDModels
    using ComposableTuringIDModels: _renewal_recurrence
    using LinearAlgebra: I
    g = [0.2, 0.3, 0.5]
    N = 1000.0
    Rt = [1.6, 1.4, 1.2, 1.0, 0.9, 0.8]
    window = 5.0 .* exp.(-0.1 .* [2, 1, 0])
    depletion = as_turing_model(SusceptibleDepletion(N), 6)()
    comp = _renewal_recurrence(g, I, (depletion,))(Rt; history = window)
    @test comp ≈ renewal_reference(g, window, Rt; N = N)
    # Depletion never raises incidence above the unconstrained renewal path.
    plain = _renewal_recurrence(g, I, ())(Rt; history = window)
    @test all(comp .<= plain .+ 1.0e-8)
    # The floor binds once the pool is exhausted.
    small = as_turing_model(SusceptibleDepletion(20.0), 6)()
    @test _renewal_recurrence(g, I, (small,))(10 .* Rt; history = window) ≈
        renewal_reference(g, window, 10 .* Rt; N = 20.0)
end

@testitem "the renewal recurrence with depletion is ForwardDiff-differentiable" setup = [RenewalReference] begin
    using ComposableTuringIDModels
    using ComposableTuringIDModels: _renewal_recurrence
    using ForwardDiff
    using LinearAlgebra: I
    g = [0.2, 0.3, 0.5]
    window = 5.0 .* exp.(-0.1 .* [2, 1, 0])
    depletion = as_turing_model(SusceptibleDepletion(1000.0), 4)()
    f(Rt) = sum(_renewal_recurrence(g, I, (depletion,))(Rt; history = window))
    ref(Rt) = sum(renewal_reference(g, window, Rt; N = 1000.0))
    Rt = [1.6, 1.4, 1.2, 1.0]
    @test ForwardDiff.gradient(f, Rt) ≈ ForwardDiff.gradient(ref, Rt)
end

@testitem "Renewal holds its modifiers in order" begin
    using ComposableTuringIDModels
    gen_int = [0.2, 0.3, 0.5]
    plain = Renewal(gen_int; rt = RandomWalk())
    depleting = Renewal(
        gen_int, SusceptibleDepletion(1000.0);
        rt = RandomWalk()
    )
    @test isempty(plain.modifiers)
    @test isempty(Renewal(; generation_time = gen_int).modifiers)
    @test only(depleting.modifiers) isa SusceptibleDepletion
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
    # A continuous generation time discretised by the constructor, carrying a
    # modifier: the two compose in either constructor form.
    renewal = Renewal(
        Gamma(2, 1.5), SusceptibleDepletion(1000.0);
        D_gen = 15.0, rt = RandomWalk(), initialisation = Normal(log(50), 0.2)
    )
    @test renewal.gen_int isa AbstractVector
    @test sum(renewal.gen_int) ≈ 1
    @test only(renewal.modifiers) isa SusceptibleDepletion
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
    # The interval is drawn per draw and the recurrence built from it and the
    # modifiers.
    depleting = Renewal(
        gen, SusceptibleDepletion(50.0);
        rt = FixedIntercept(log(2.0)), initialisation = Normal(log(10), 0.1)
    )
    plain = Renewal(
        gen; rt = FixedIntercept(log(2.0)),
        initialisation = Normal(log(10), 0.1)
    )
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
    # to say so itself rather than failing inside the modifier seam.
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
    # the recurrence as an `R_t` with no time axis.
    gi = [0.2, 0.3, 0.5]
    kw = Renewal(;
        generation_time = gi, rt = Normal(0, 0.1), initialisation = Normal()
    )
    pos = Renewal(gi, exp, Normal(0, 0.1), Normal(), kw.mixing, ())

    @test typeof(pos.rt) == typeof(kw.rt)
    @test pos.rt isa Intercept
    n = 20
    @test length(as_turing_model(kw, n)().I_t) == n
    @test length(as_turing_model(pos, n)().I_t) == n
    # Rebuilding from stored fields is a no-op, so `Accessors.@set` is safe.
    @test path_prior(kw.rt) === kw.rt
    # Widening runs before `new`, so the slot's role guard still fires.
    @test_throws TypeError Renewal(
        gi, exp, PoissonError(), Normal(), kw.mixing, ()
    )
end
