# Tests for the non-centred `InfectionNoise` renewal modifier: the draw it
# applies inside the scan, the seam it resolves through, and the `is_noise`
# trait that marks where the renewal expectation stops being one.

@testitem "InfectionNoise draws the log-scale non-centred transform" begin
    using ComposableTuringIDModels, Distributions
    import ComposableTuringIDModels: InfectionNoiseDraws, apply_modifier,
        modifier_init_state

    # sigma2 = log1p(cv^2); I = exp(log(iota) - sigma2/2 + raw*sigma), with the
    # coefficient of variation uncapped so it is the exact one.
    raw = [0.5, -1.25, 2.0]
    ξ = 0.1
    noise = InfectionNoise(; overdispersion = ξ, cv_cap = Inf)
    mod = InfectionNoiseDraws(raw, noise, ξ)
    ιs = (100.0, 2500.0, 7.0)
    # Thread the substate by folding, so each step reads its own draw.
    got = accumulate(
        (acc, i) -> apply_modifier(mod, ιs[i], acc[2]), eachindex(ιs);
        init = (0.0, modifier_init_state(mod, nothing))
    )
    for (i, ι) in enumerate(ιs)
        σ² = log1p(1 / ι + ξ^2)
        @test got[i][1] ≈ exp(log(ι) - σ² / 2 + raw[i] * sqrt(σ²))
    end
    # The substate is the step counter, so it advances once per step.
    @test last(got)[2] == length(ιs)
end

@testitem "InfectionNoise takes the noise family as an argument" begin
    using ComposableTuringIDModels, Distributions
    import ComposableTuringIDModels: InfectionNoiseDraws, apply_modifier

    ξ, ι, z = 0.1, 400.0, 1.25
    cv = 0.5 - log1p(exp(10.0 * (0.5 - sqrt(1 / ι + ξ^2)))) / 10.0

    linear = InfectionNoise(; dist = Normal, overdispersion = ξ)
    got_n, _ = apply_modifier(InfectionNoiseDraws([z], linear, ξ), ι, 0)
    # The linear form: I = iota + z cv iota.
    @test got_n ≈ ι + z * cv * ι

    positive = InfectionNoise(; overdispersion = ξ)
    got_l, _ = apply_modifier(InfectionNoiseDraws([z], positive, ξ), ι, 0)
    σ² = log1p(cv^2)
    @test got_l ≈ exp(log(ι) - σ² / 2 + z * sqrt(σ²))

    # The positive family is the default, because the linear one can take the
    # incidence below zero and the modifier cannot bound it.
    below, _ = apply_modifier(
        InfectionNoiseDraws([-6.0], InfectionNoise(; dist = Normal, cv_cap = Inf), ξ),
        1.0, 0
    )
    @test below < 0
    for z in (-50.0, -8.0, 0.0, 8.0), ι in (1.0e-3, 1.0, 1.0e4)
        above, _ = apply_modifier(
            InfectionNoiseDraws([z], InfectionNoise(; cv_cap = Inf), ξ), ι, 0
        )
        @test above > 0
        @test isfinite(above)
    end
end

@testitem "InfectionNoise samples one standard normal per time point" begin
    using ComposableTuringIDModels, Distributions, Random

    r = Renewal(
        [0.2, 0.3, 0.5], InfectionNoise();
        rt = RandomWalk(), initialisation = Normal()
    )
    names = string.(collect(keys(rand(Xoshiro(1), as_turing_model(r, 10)))))
    # Namespaced under the modifier's position, so it is recoverable.
    @test any(contains("modifier_1.I_raw"), names)
    # A fixed overdispersion is a constant, not a sampled parameter.
    @test !any(contains("ξ"), names)

    withprior = Renewal(
        [0.2, 0.3, 0.5], InfectionNoise(; overdispersion = HalfNormal(0.1));
        rt = RandomWalk(), initialisation = Normal()
    )
    @test any(
        contains("ξ"),
        string.(collect(keys(rand(Xoshiro(1), as_turing_model(withprior, 10)))))
    )
end

@testitem "InfectionNoise makes a renewal's infections stochastic" begin
    using ComposableTuringIDModels, Distributions, Random

    gen = [0.2, 0.3, 0.5]
    fixed = Renewal(
        gen; rt = FixedIntercept(0.0), initialisation = Normal(log(500.0), 0.0)
    )
    noisy = Renewal(
        gen, InfectionNoise(); rt = FixedIntercept(0.0),
        initialisation = Normal(log(500.0), 0.0)
    )
    # With R_t and the seeding fixed the deterministic process is constant
    # across draws and the noisy one is not.
    @test as_turing_model(fixed, 12)(Xoshiro(1)).I_t ≈
        as_turing_model(fixed, 12)(Xoshiro(2)).I_t
    c = as_turing_model(noisy, 12)(Xoshiro(1)).I_t
    e = as_turing_model(noisy, 12)(Xoshiro(2)).I_t
    @test !isapprox(c, e)
    @test all(>(0), c)
end

@testitem "the two parameterisations draw the same noise" begin
    using ComposableTuringIDModels, Distributions, Random
    using DynamicPPL: fix
    using Statistics: mean, std

    # Centred and non-centred differ in geometry, not in the model they
    # express: the first step's marginal is the same distribution either way.
    gen = [0.2, 0.3, 0.5]
    noise = InfectionNoise(; overdispersion = 0.2, cv_cap = Inf)
    args = (; rt = FixedIntercept(0.0), initialisation = Normal())
    fixinit = (init_incidence = log(200.0),)
    centred = StochasticRenewal(gen; noise = noise, args...)
    non_centred = Renewal(gen, noise; args...)
    firsts(m) = [
        fix(as_turing_model(m, 6), fixinit)(Xoshiro(s)).I_t[1] for s in 1:4000
    ]
    a, b = firsts(centred), firsts(non_centred)
    # The renewal expectation the first step is centred on, and the
    # negative-binomial standard deviation at it.
    ι = 200.0
    @test mean(a) ≈ ι rtol = 0.05
    @test mean(b) ≈ ι rtol = 0.05
    @test std(a) ≈ sqrt(ι * (1 + ι * 0.2^2)) rtol = 0.1
    @test std(b) ≈ sqrt(ι * (1 + ι * 0.2^2)) rtol = 0.1
end

@testitem "is_noise marks the modifier that draws infections" begin
    using ComposableTuringIDModels, Distributions
    import ComposableTuringIDModels: is_noise, InfectionNoiseDraws,
        SusceptibleDepletion, ImportedRate

    # A trait beside `apply_modifier`, false for every modifier that only
    # transforms the expectation.
    @test !is_noise(SusceptibleDepletion(1000.0))
    @test !is_noise(ImportedCases(Normal()))
    @test !is_noise(ImportedRate([1.0]))
    # True for the one that replaces it with a draw, resolved or not, since
    # the scan sees the resolved modifier.
    @test is_noise(InfectionNoise())
    @test is_noise(InfectionNoiseDraws([0.0], InfectionNoise(), 0.1))
end

@testitem "InfectionNoise composes with the other modifiers" begin
    using ComposableTuringIDModels, Distributions, Random
    using DynamicPPL: fix

    gen = [0.2, 0.3, 0.5]
    fixinit = (init_incidence = log(50.0),)
    # Order is the tuple's: depletion scales the expectation, then the noise
    # draws around the depleted value.
    r = Renewal(
        gen, SusceptibleDepletion(2000.0), InfectionNoise();
        rt = FixedIntercept(log(2.0)), initialisation = Normal()
    )
    free = Renewal(
        gen, InfectionNoise(); rt = FixedIntercept(log(2.0)),
        initialisation = Normal()
    )
    bounded_I = fix(as_turing_model(r, 20), fixinit)(Xoshiro(4)).I_t
    free_I = fix(as_turing_model(free, 20), fixinit)(Xoshiro(4)).I_t
    @test bounded_I[end] < free_I[end]
    @test all(>(0), bounded_I)
    @test length(r.recurrent_step.modifiers) == 2
end

@testitem "the noisy renewal scan step infers concretely" begin
    using ComposableTuringIDModels, Distributions
    import ComposableTuringIDModels: InfectionNoiseDraws, apply_modifier

    # The scan runs this once per timestep inside the differentiated log
    # density, so an abstract return type here is paid on every gradient.
    for family in (LogNormal, Normal)
        noise = InfectionNoise(; dist = family)
        # `isconcretetype` is the wrong predicate for a field holding a
        # `Type{...}`; `isdispatchelem` asks whether a method can be selected
        # on it at compile time, which is what matters.
        ft = fieldtype(typeof(noise), :dist)
        @test Base.isdispatchelem(ft)
        @test !(ft isa UnionAll)

        draws = InfectionNoiseDraws(randn(8), noise, noise.overdispersion)
        rt = only(
            Base.return_types(apply_modifier, (typeof(draws), Float64, Int))
        )
        @test rt === Tuple{Float64, Int}
    end
end

@testitem "the noisy renewal scores and differentiates as a posterior" begin
    using ComposableTuringIDModels, Distributions, Random, ForwardDiff
    using DynamicPPL: LogDensityFunction, VarInfo, link, getlogjoint
    import LogDensityProblems as LDP

    idm = IDModel(
        Renewal(
            [0.2, 0.3, 0.5], InfectionNoise();
            rt = RandomWalk(), initialisation = Normal()
        ),
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
    @test !all(iszero, g)
end

@testitem "an epidemic-scale expectation draws, rather than going quiet" begin
    using ComposableTuringIDModels, Distributions, Random
    import ComposableTuringIDModels: _noise_sd, _soft_upper

    # The soft cap is smooth, and smoothness costs it an absolute offset.
    # Below that offset it returns zero or less, which is no coefficient of
    # variation at all, reachable at a few million infections with a small
    # overdispersion, which is an ordinary national-scale model. The cap is a
    # limit, so it is not allowed to annihilate the value it limits.
    noise = InfectionNoise(; overdispersion = 0.0)
    for ι in (1.0e5, 1.0e6, 3.0e6, 1.0e8)
        sd = _noise_sd(noise, ι, 0.0)
        @test isfinite(sd)
        @test sd > 0
    end
    # Where the smooth cap is still positive it is used unchanged, so R's
    # constants are matched wherever they mean anything.
    raw = sqrt(1 / 1.0e5)
    @test _noise_sd(noise, 1.0e5, 0.0) ≈ _soft_upper(raw, 0.5, 10.0) * 1.0e5
    @test _soft_upper(raw, 0.5, 10.0) < raw
    # Past the crossover the exact coefficient of variation is used instead.
    @test _noise_sd(noise, 3.0e6, 0.0) ≈ sqrt(1 / 3.0e6) * 3.0e6
    @test _soft_upper(sqrt(1 / 3.0e6), 0.5, 10.0) < 0

    # The whole model, at national-scale infections.
    gi = [0.2, 0.3, 0.5]
    args = (;
        rt = FixedIntercept(log(1.5)),
        initialisation = FixedIntercept(log(3.0e6)),
    )
    non_centred = Renewal(gi, noise; args...)
    centred = StochasticRenewal(gi; noise = noise, args...)
    for m in (non_centred, centred)
        I_t = as_turing_model(m, 12)(Xoshiro(3)).I_t
        @test all(isfinite, I_t)
        @test all(>(0), I_t)
    end
end

@testitem "a model that cannot produce a value says so when asked for one" begin
    using ComposableTuringIDModels, Distributions, Random
    using DynamicPPL: @varname
    import ComposableTuringIDModels: _moment_dist

    # Scoring and sampling want opposite things from an invalid moment pair.
    # A rejected proposal must stay silent and cheap; a simulated value has no
    # proposal to reject, so `Inf` there would be a wrong answer told to
    # nobody.
    d = _moment_dist(LogNormal, 100.0, -1.0)
    @test logpdf(d, 50.0) == -Inf
    @test_throws ArgumentError rand(d)
    @test occursin(
        "cannot produce a value", sprint(
            showerror, try
                rand(d)
            catch e
                e
            end
        )
    )

    # A `missing` observation on a model that is fine draws as before.
    y = rand(Xoshiro(4), as_turing_model(NormalError(), missing, fill(10.0, 5)))
    @test all(isfinite, y[@varname(y_t)])
end
