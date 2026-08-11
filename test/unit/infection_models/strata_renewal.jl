# Tests for a stratified `Renewal`: one recursion per stratum sharing a
# window contract with the single-series case, per-stratum generation
# intervals, and the two modifiers that widen to a strata axis
# (`SusceptibleDepletion`'s pool, `ImportedCases`' exogenous rate).

@testitem "the uncoupled path is unchanged by adding a strata axis" begin
    using ComposableTuringIDModels, Distributions, Random
    # The single most important guarantee in this file: a one-stratum
    # `Dims{2}` shape must reproduce a plain `Int`-shaped series exactly, so a
    # `Stratify`-driven model with no strata gain never changes existing
    # behaviour.
    gen_int = [0.2, 0.3, 0.5]
    single = Renewal(;
        generation_time = gen_int, rt = RandomWalk(),
        initialisation = Normal()
    )
    strat = Renewal(;
        generation_time = gen_int,
        rt = Stratify(RandomWalk(), FixedIntercept(0.0)),
        initialisation = Normal()
    )
    n = 20

    Random.seed!(42)
    I_single = as_turing_model(single, n)().I_t
    Random.seed!(42)
    I_strat = as_turing_model(strat, (1, n))().I_t

    @test size(I_strat) == (1, n)
    @test vec(I_strat) ≈ I_single
end

@testitem "a stratified renewal runs one recursion per stratum" begin
    using ComposableTuringIDModels, Distributions, Random, Statistics
    Random.seed!(702)
    gen_int = [0.2, 0.3, 0.5]
    n_strata, n_time = 3, 30
    # A vector-valued (IID) initialisation gives each stratum its own seed;
    # a bare Distribution would instead broadcast one shared seed.
    model = Renewal(;
        generation_time = gen_int,
        rt = Stratify(RandomWalk(), Hierarchy()),
        initialisation = IID(Normal(log(20.0), 1.0))
    )
    out = as_turing_model(model, (n_strata, n_time))()

    @test size(out.I_t) == (n_strata, n_time)
    # Per-stratum seeds differ.
    @test length(unique(round.(out.I_t[:, 1]; digits = 8))) == n_strata
    # Strata diverge dynamically rather than staying proportional: the ratio
    # between two strata's rows is not constant over time.
    ratio = out.I_t[1, :] ./ out.I_t[2, :]
    @test std(ratio) / abs(mean(ratio)) > 0.01
end

@testitem "a strata x lags interval convolves each stratum on its own" begin
    using ComposableTuringIDModels: ConstantRenewalStep, accumulate_scan
    G = [0.6 0.3 0.1; 0.1 0.3 0.6]   # 2 strata, already-reversed lags
    core = ConstantRenewalStep(G)
    # Identical seed values in both rows: any divergence is down to the
    # per-stratum generation interval alone.
    window0 = [1.0 2.0 3.0; 1.0 2.0 3.0]
    Rt = [fill(1.5, 2) for _ in 1:10]
    I_t = accumulate_scan(
        core, (; val = window0[:, end], window = window0), Rt
    )
    @test size(I_t) == (2, 10)
    @test I_t[1, :] != I_t[2, :]
end

@testitem "a Hierarchy-parameterised UncertainDelay pools across strata" begin
    using ComposableTuringIDModels, Distributions, Random
    Random.seed!(708)
    gen = UncertainDelay(
        LogNormal,
        [Hierarchy(), truncated(Normal(0.3, 0.05), 0, Inf)]; D = 10.0
    )
    model = Renewal(;
        generation_time = gen,
        rt = Stratify(RandomWalk(), FixedIntercept(0.0)),
        initialisation = IID(Normal(log(20.0), 0.3))
    )
    n_strata, n_time = 3, 15
    out = as_turing_model(model, (n_strata, n_time))()
    @test size(out.I_t) == (n_strata, n_time)
    @test all(isfinite, out.I_t)
    @test all(>=(0), out.I_t)
end

@testitem "a constant-parameter UncertainDelay is shared by every stratum" begin
    using ComposableTuringIDModels, Distributions, Random
    Random.seed!(709)
    gen = UncertainDelay(
        LogNormal,
        [Normal(1.5, 0.2), truncated(Normal(0.3, 0.05), 0, Inf)]; D = 10.0
    )
    model = Renewal(;
        generation_time = gen,
        rt = Stratify(RandomWalk(), FixedIntercept(0.0)),
        initialisation = IID(Normal(log(20.0), 0.3))
    )
    n_strata, n_time = 3, 15
    out = as_turing_model(model, (n_strata, n_time))()
    @test size(out.I_t) == (n_strata, n_time)
    @test all(isfinite, out.I_t)
end

@testitem "SusceptibleDepletion, scalar pop_size: one pool, all strata" begin
    using ComposableTuringIDModels: ConstantRenewalStep, RenewalStep,
        SusceptibleDepletion, accumulate_scan,
        renewal_init_state
    rev_gen = reverse([0.2, 0.3, 0.5])
    core = ConstantRenewalStep(rev_gen)
    step = RenewalStep(core, (SusceptibleDepletion(1000.0),))
    n_strata = 2
    I0, r = fill(5.0, n_strata), fill(0.1, n_strata)
    init = renewal_init_state(step, I0, r, length(rev_gen))
    # The scalar pool is spread to a separate, equal-sized copy per stratum.
    @test only(init.substates) == fill(1000.0, n_strata)

    Rt = [fill(1.6, n_strata) for _ in 1:10]
    I_t = accumulate_scan(step, init, Rt)
    plain = accumulate_scan(
        core, (; val = init.window[:, end], window = init.window), Rt
    )
    @test size(I_t) == (n_strata, 10)
    @test all(I_t .<= plain .+ 1.0e-8)
end

@testitem "SusceptibleDepletion, per-stratum pop_size: separate pools" begin
    using ComposableTuringIDModels: ConstantRenewalStep, RenewalStep,
        SusceptibleDepletion, accumulate_scan,
        renewal_init_state
    rev_gen = reverse([0.2, 0.3, 0.5])
    core = ConstantRenewalStep(rev_gen)
    pop = [200.0, 20_000.0]
    step = RenewalStep(core, (SusceptibleDepletion(pop),))
    I0, r = fill(5.0, 2), fill(0.2, 2)
    init = renewal_init_state(step, I0, r, length(rev_gen))
    @test only(init.substates) == pop

    Rt = [fill(1.6, 2) for _ in 1:20]
    I_t = accumulate_scan(step, init, Rt)
    plain = accumulate_scan(
        core, (; val = init.window[:, end], window = init.window), Rt
    )
    @test all(I_t .<= plain .+ 1.0e-8)
    # The small pool depletes harder: it ends further below the unconstrained
    # path than the (barely touched) large pool.
    rel_small = I_t[1, end] / plain[1, end]
    rel_large = I_t[2, end] / plain[2, end]
    @test rel_small < rel_large

    mismatched = RenewalStep(core, (SusceptibleDepletion([1.0, 2.0, 3.0]),))
    @test_throws AssertionError renewal_init_state(
        mismatched, I0, r, length(rev_gen)
    )
end

@testitem "_seed rejects a mismatched per-stratum initialisation" begin
    using ComposableTuringIDModels: _seed
    # A vector `initialisation` must draw one seed per stratum; a length that
    # disagrees with the stratum count is a build-time mistake, not silently
    # broadcast or truncated.
    @test_throws AssertionError _seed([1.0, 2.0, 3.0], (2, 10))
end

@testitem "SusceptibleDepletion and ImportedCases together: order matters" begin
    using ComposableTuringIDModels, Distributions, Random
    Random.seed!(712)
    # Each modifier is tested alone at `Dims{2}`, and the two are tested
    # together at single-series (imported_cases.jl), but never both together
    # on a stratified renewal. Their order decides whether the imports are
    # scaled by the susceptible fraction (depletion first) or added on top
    # of it untouched (imports first), so the two orderings must differ.
    gen_int = [0.2, 0.3, 0.5]
    n_strata, n_time = 2, 20
    rt = Stratify(FixedIntercept(log(1.2)), FixedIntercept(0.0))
    seeds = [Dirac(log(5.0)), Dirac(log(5.0))]
    pop = [200.0, 200.0]
    rate = Stratify(FixedIntercept(log(0.5)), FixedIntercept(0.0))

    deplete_then_import = Renewal(
        gen_int, SusceptibleDepletion(pop),
        ImportedCases(rate); rt = rt, initialisation = seeds
    )
    import_then_deplete = Renewal(
        gen_int, ImportedCases(rate),
        SusceptibleDepletion(pop); rt = rt, initialisation = seeds
    )

    Random.seed!(99)
    I_deplete_first = as_turing_model(
        deplete_then_import, (n_strata, n_time)
    )().I_t
    Random.seed!(99)
    I_import_first = as_turing_model(
        import_then_deplete, (n_strata, n_time)
    )().I_t

    @test size(I_deplete_first) == (n_strata, n_time)
    @test all(isfinite, I_deplete_first)
    @test all(isfinite, I_import_first)
    @test I_deplete_first != I_import_first
end

@testitem "ImportedCases with a Stratify rate: one stream per stratum" begin
    using ComposableTuringIDModels, Distributions, Random
    Random.seed!(711)
    gen_int = [0.2, 0.3, 0.5]
    n_strata, n_time = 2, 25
    rt = Stratify(FixedIntercept(log(1.0)), FixedIntercept(0.0))
    seeds = [Dirac(log(1.0)), Dirac(log(1.0))]

    plain = Renewal(gen_int; rt = rt, initialisation = seeds)
    rate = Stratify(FixedIntercept(log(0.5)), Hierarchy(; mean = Dirac(0.0)))
    imported = Renewal(
        gen_int, ImportedCases(rate); rt = rt,
        initialisation = seeds
    )

    I_plain = as_turing_model(plain, (n_strata, n_time))().I_t
    I_imported = as_turing_model(imported, (n_strata, n_time))().I_t

    @test size(I_imported) == (n_strata, n_time)
    # The importation rate is strictly positive (its `exp` transformation),
    # so it lifts every stratum's incidence at every step.
    @test all(I_imported .>= I_plain .- 1.0e-8)
    @test all(I_imported[:, end] .> I_plain[:, end])
end
