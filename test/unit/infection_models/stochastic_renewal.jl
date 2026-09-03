# Tests for the centred stochastic renewal: the per-step moments, infections
# as a sampled parameter, and the propose/commit split the loop shares with
# the deterministic scan.

@testitem "InfectionNoise matches the negative binomial's first two moments" begin
    using ComposableTuringIDModels, Distributions
    import ComposableTuringIDModels: _noise_dist, _noise_sd

    # Mean iota, variance iota(1 + iota xi^2), which is what the family is
    # matched to. Uncapped, so the exact variance is kept.
    noise = InfectionNoise(; cv_cap = Inf)
    for ι in (5.0, 250.0, 4000.0), ξ in (0.0, 0.1, 0.4)
        d = _noise_dist(InfectionNoise(; overdispersion = ξ, cv_cap = Inf), ι, ξ)
        @test mean(d) ≈ ι rtol = 1.0e-10
        @test var(d) ≈ ι * (1 + ι * ξ^2) rtol = 1.0e-10
        @test _noise_sd(noise, ι, ξ) ≈ sqrt(ι * (1 + ι * ξ^2)) rtol = 1.0e-10
    end
    # The default family is the positive one, so no truncation is needed.
    @test InfectionNoise().dist === LogNormal
    @test minimum(_noise_dist(InfectionNoise(), 100.0, 0.1)) == 0.0
end

@testitem "InfectionNoise caps the coefficient of variation" begin
    using ComposableTuringIDModels, Distributions
    import ComposableTuringIDModels: _noise_sd, _soft_upper

    # `sqrt(1/iota + xi^2)` diverges as the expectation approaches zero: an
    # arbitrarily small expectation would get an arbitrarily wide draw.
    ξ, ι = 0.1, 1.0e-4
    capped = InfectionNoise(; overdispersion = ξ)
    uncapped = InfectionNoise(; overdispersion = ξ, cv_cap = Inf)
    @test _noise_sd(capped, ι, ξ) / ι ≈ _soft_upper(sqrt(1 / ι + ξ^2), 0.5, 10.0)
    # Saturated at the cap there, two orders of magnitude below the exact
    # coefficient of variation.
    @test _noise_sd(capped, ι, ξ) / ι ≈ 0.5 atol = 1.0e-6
    @test _noise_sd(capped, ι, ξ) < 0.01 * _noise_sd(uncapped, ι, ξ)
    # At an expectation where the cap does not bite the two are close, but
    # the smooth limit still costs an absolute offset of
    # `log1p(exp(-k(u - c)))/k` — the price of it being smooth. That offset
    # approaches `log1p(exp(-k u))/k` as the coefficient of variation falls,
    # so with `overdispersion = 0` and a large enough expectation it takes it
    # to zero. `cv_cap = Inf` is the way out.
    u, k = 0.5, 10.0
    c = sqrt(1 / 1.0e4 + ξ^2)
    @test _noise_sd(capped, 1.0e4, ξ) ≈ (c - log1p(exp(-k * (u - c))) / k) * 1.0e4
    @test _noise_sd(capped, 1.0e4, ξ) > 0.97 * _noise_sd(uncapped, 1.0e4, ξ)
    # Past that point the cap would return a negative coefficient of
    # variation, so the exact one is used instead and the draw stays valid.
    @test _soft_upper(sqrt(1 / 1.0e8), 0.5, 10.0) < 0
    @test _noise_sd(InfectionNoise(; overdispersion = 0.0), 1.0e8, 0.0) > 0

    # An expectation the recursion cannot come back from has no moment pair,
    # and must not raise from `sqrt` mid-gradient.
    @test isnan(_noise_sd(capped, -1.0, ξ))
    @test isnan(_noise_sd(capped, 0.0, ξ))
end

@testitem "StochasticRenewal makes infections a sampled parameter" begin
    using ComposableTuringIDModels, Distributions, Random
    using DynamicPPL: DebugUtils

    sr = StochasticRenewal(
        [0.2, 0.3, 0.5]; rt = FixedIntercept(0.0),
        initialisation = Normal(log(500.0), 0.0)
    )
    mdl = as_turing_model(sr, 12)
    drawn = rand(Xoshiro(1), mdl)
    names = string.(collect(keys(drawn)))
    # One variable per time point: this is the centred form, so infections are
    # parameters the likelihood informs rather than a transform of standard
    # normals.
    for t in 1:12
        @test any(contains("I_t[$t]"), names)
    end
    @test !any(contains("I_raw"), names)
    @test DebugUtils.check_model(mdl; error_on_failure = false)

    # Stochastic: with R_t and the seeding both fixed the deterministic
    # process is constant across draws and this one is not.
    a = mdl(Xoshiro(1)).I_t
    b = mdl(Xoshiro(2)).I_t
    @test !isapprox(a, b)
    @test all(>(0), a)
    @test all(isfinite, a)
end

@testitem "StochasticRenewal centres on the renewal expectation" begin
    using ComposableTuringIDModels, Distributions, Random
    using DynamicPPL: fix
    using Statistics: mean

    gen = [0.2, 0.3, 0.5]
    fixinit = (init_incidence = log(1.0e5),)
    det = fix(
        as_turing_model(
            Renewal(gen; rt = FixedIntercept(log(1.2)), initialisation = Normal()),
            10
        ), fixinit
    )().I_t
    # At a large expectation the coefficient of variation is small, so the
    # centred draw tracks the deterministic recursion it is centred on.
    sr = StochasticRenewal(
        gen; rt = FixedIntercept(log(1.2)), initialisation = Normal(),
        noise = InfectionNoise(; overdispersion = 0.0)
    )
    paths = [fix(as_turing_model(sr, 10), fixinit)(Xoshiro(s)).I_t for s in 1:60]
    got = mean(paths)
    @test all(isapprox.(got, det; rtol = 0.05))
end

@testitem "StochasticRenewal takes the modifiers Renewal takes" begin
    using ComposableTuringIDModels, Distributions, Random
    using DynamicPPL: fix

    gen = [0.2, 0.3, 0.5]
    fixinit = (init_incidence = log(50.0),)
    # A modifier transforms the expectation the draw is centred on, so
    # susceptible depletion still bends the epidemic over.
    growing = StochasticRenewal(
        gen; rt = FixedIntercept(log(2.0)), initialisation = Normal(),
        noise = InfectionNoise(; overdispersion = 0.0)
    )
    depleting = StochasticRenewal(
        gen, SusceptibleDepletion(2000.0);
        rt = FixedIntercept(log(2.0)), initialisation = Normal(),
        noise = InfectionNoise(; overdispersion = 0.0)
    )
    free = fix(as_turing_model(growing, 20), fixinit)(Xoshiro(4)).I_t
    bounded = fix(as_turing_model(depleting, 20), fixinit)(Xoshiro(4)).I_t
    @test bounded[end] < free[end]
    @test all(>(0), bounded)

    # A prior-carrying modifier resolves through the same pre-scan seam.
    imported = StochasticRenewal(
        gen, ImportedCases(Normal(0.5, 0.01));
        rt = FixedIntercept(log(1.0)), initialisation = Normal()
    )
    names = string.(collect(keys(rand(Xoshiro(5), as_turing_model(imported, 12)))))
    @test any(contains("modifier_1.import_rates"), names)
end

@testitem "StochasticRenewal resolves an inferred generation interval" begin
    using ComposableTuringIDModels, Distributions, Random

    gen = UncertainDelay(
        LogNormal, [Normal(1.9, 0.2), truncated(Normal(0.5, 0.2), 0, Inf)];
        D = 14.0
    )
    sr = StochasticRenewal(;
        generation_time = gen, rt = RandomWalk(), initialisation = Normal()
    )
    @test isnothing(sr.recurrent_step)
    drawn = rand(Xoshiro(6), as_turing_model(sr, 10))
    names = string.(collect(keys(drawn)))
    @test any(contains("gen"), names)
    @test any(contains("I_t[10]"), names)
end

@testitem "StochasticRenewal keeps a fixed overdispersion out of the chain" begin
    using ComposableTuringIDModels, Distributions, Random

    gen = [0.2, 0.3, 0.5]
    slots(noise) = string.(
        collect(
            keys(
                rand(
                    Xoshiro(7),
                    as_turing_model(
                        StochasticRenewal(
                            gen; rt = FixedIntercept(0.0),
                            initialisation = Normal(), noise = noise
                        ), 6
                    )
                )
            )
        )
    )
    # A fixed scalar is a constant, not a point-mass parameter.
    @test !any(contains("ξ"), slots(InfectionNoise()))
    # A prior on it adds one slot.
    @test any(contains("ξ"), slots(InfectionNoise(; overdispersion = HalfNormal(0.1))))
end

@testitem "StochasticRenewal runs one series" begin
    using ComposableTuringIDModels, Distributions

    sr = StochasticRenewal([0.2, 0.3, 0.5]; rt = Stratify(RandomWalk(), Hierarchy()))
    @test_throws ErrorException as_turing_model(sr, (3, 10))
end

@testitem "StochasticRenewal scores and differentiates as a posterior" begin
    using ComposableTuringIDModels, Distributions, Random, ForwardDiff
    using DynamicPPL: LogDensityFunction, VarInfo, link, getlogjoint
    import LogDensityProblems as LDP

    gen = [0.2, 0.3, 0.5]
    idm = IDModel(
        StochasticRenewal(gen; rt = RandomWalk(), initialisation = Normal()),
        NegativeBinomialError()
    )
    n = 10
    y = as_turing_model(idm, missing, n)(Xoshiro(8)).generated_y_t
    mdl = as_turing_model(idm, y, n)
    ldf = LogDensityFunction(mdl, getlogjoint, link(VarInfo(mdl), mdl))
    x = zeros(LDP.dimension(ldf))
    @test isfinite(LDP.logdensity(ldf, x))
    g = ForwardDiff.gradient(v -> LDP.logdensity(ldf, v), x)
    @test all(isfinite, g)
    # Infections carry their own gradient dimension, which is what a centred
    # parameterisation means: `n` of them, plus the latent path and seeding.
    @test LDP.dimension(ldf) >= n + 1
    @test !all(iszero, g)
end

@testitem "the propose/commit split leaves the deterministic scan unchanged" begin
    using ComposableTuringIDModels, Distributions, Random
    import ComposableTuringIDModels: ConstantRenewalStep, RenewalStep,
        SusceptibleDepletion, renewal_init_state, _propose, _commit

    # The stochastic loop and the scan run the same arithmetic, so a step must
    # be exactly its own proposal committed.
    rev = reverse([0.2, 0.3, 0.5])
    for step in (
            ConstantRenewalStep(rev),
            RenewalStep(ConstantRenewalStep(rev)),
            RenewalStep(ConstantRenewalStep(rev), (SusceptibleDepletion(1000.0),)),
        )
        state = renewal_init_state(step, 10.0, 0.1, 3)
        inc, subs = _propose(step, state, 1.4)
        got = _commit(step, state, inc, subs)
        want = step(state, 1.4)
        @test got.val ≈ want.val
        @test got.window ≈ want.window
    end
end

@testitem "the renewal summaries serve both renewal models" begin
    using ComposableTuringIDModels, Distributions

    # `R_to_r` and `expected_Rt` read the generation interval off the model,
    # which the stochastic renewal carries exactly as the deterministic one
    # does.
    gen = [0.2, 0.3, 0.5]
    det = Renewal(gen; rt = RandomWalk(), initialisation = Normal())
    sto = StochasticRenewal(gen; rt = RandomWalk(), initialisation = Normal())
    @test R_to_r(1.5, sto) ≈ R_to_r(1.5, det)
    I_t = [100.0, 200.0, 300.0, 400.0, 500.0]
    @test expected_Rt(sto, I_t) ≈ expected_Rt(det, I_t)
    # An inferred interval has no single value for either to use.
    inferred = StochasticRenewal(;
        generation_time = UncertainDelay(
            LogNormal, [Normal(1.9, 0.2), truncated(Normal(0.5, 0.2), 0, Inf)];
            D = 14.0
        )
    )
    @test_throws ArgumentError expected_Rt(inferred, I_t)
end

@testitem "the shared setup keeps both renewal models' variable names" begin
    using ComposableTuringIDModels, Distributions, DynamicPPL, Random

    # Both models draw their setup through one submodel with prefixing off, so
    # every name is the one it carries when the setup is written out inline.
    # Pinned here because a prefix slipping back on renames every chain and
    # nothing else would fail.
    varnames(m, n) = string.(collect(keys(VarInfo(Xoshiro(1), as_turing_model(m, n)))))
    gen = [0.2, 0.3, 0.5]
    uncertain = UncertainDelay(
        LogNormal, [Normal(1.9, 0.2), truncated(Normal(0.5, 0.2), 0, Inf)];
        D = 14.0
    )
    latent = ["init", "std", "ϵ_t"]
    seed_path = SeedingPath(RandomWalk(; init = Normal(log(50), 0.5)))
    seed_names = ["init_incidence.init", "init_incidence.std", "init_incidence.ϵ_t"]

    det(; kwargs...) = Renewal(; rt = RandomWalk(), kwargs...)
    @test varnames(det(; generation_time = gen, initialisation = Normal()), 20) ==
        [latent; "init_incidence"]
    @test varnames(det(; generation_time = gen, initialisation = seed_path), 20) ==
        [latent; seed_names]
    @test varnames(det(; generation_time = uncertain, initialisation = Normal()), 20) ==
        [latent; "gen.θ"; "init_incidence"]
    @test varnames(
        Renewal(gen, ImportedCases(Normal(-1, 0.5)); rt = RandomWalk()), 20
    ) == [latent; "init_incidence"; "modifier_1.import_rates"]
    @test varnames(
        det(;
            generation_time = gen, initialisation = Normal(),
            rt = Stratify(RandomWalk(), Hierarchy())
        ), (3, 20)
    ) == [
        latent; "across.mean"; "across.group_effects.ϵ_t"; "init_incidence"
    ]

    n = 12
    I_names = ["I_t[$t]" for t in 1:n]
    sto(; kwargs...) = StochasticRenewal(; rt = RandomWalk(), kwargs...)
    @test varnames(sto(; generation_time = gen, initialisation = Normal()), n) ==
        [latent; "init_incidence"; I_names]
    @test varnames(sto(; generation_time = gen, initialisation = seed_path), n) ==
        [latent; seed_names; I_names]
    @test varnames(sto(; generation_time = uncertain, initialisation = Normal()), n) ==
        [latent; "gen.θ"; "init_incidence"; I_names]
    # A drawn overdispersion sits between the setup and the recursion, as it
    # does when the setup is inline.
    @test varnames(
        sto(;
            generation_time = gen, initialisation = Normal(),
            noise = InfectionNoise(; overdispersion = Exponential(0.1))
        ), n
    ) == [latent; "init_incidence"; "ξ"; I_names]
end
