# Per-step parameters wired through the single seam + `_at`: a bare `Distribution`
# draws a SCALAR (constant, no length-n allocation) while a process draws a
# length-n path (time-varying / hierarchical), with no per-component special-casing.

@testitem "wired params: Distribution ⇒ scalar, process ⇒ length-n path" begin
    using ComposableTuringIDModels, Distributions, Random
    Random.seed!(1)
    # count the RVs a slot draws under its own namespace prefix
    nkeys(m, pre) = count(k -> startswith(string(k), pre), keys(rand(m)))

    # NegativeBinomialError.cluster_factor
    for n in (8, 40)
        @test nkeys(
            as_turing_model(NegativeBinomialError(), missing,
                fill(10.0, n)), "cluster_factor") == 1        # scalar, flat in n
    end
    @test nkeys(
        as_turing_model(
            NegativeBinomialError(; cluster_factor =
            RandomWalk()), missing, fill(10.0, 20)),
        "cluster_factor") > 1

    # NormalError.std
    @test nkeys(as_turing_model(NormalError(), missing, fill(10.0, 30)), "σ") == 1
    @test nkeys(
        as_turing_model(NormalError(; std = RandomWalk()), missing,
            fill(10.0, 20)), "σ") > 1

    # HierarchicalNormal.std
    @test nkeys(as_turing_model(HierarchicalNormal(), 30), "std") == 1
    @test nkeys(
        as_turing_model(
            HierarchicalNormal(; std = RandomWalk()), 20), "std") > 1

    # MA.θ (order-1 coefficient)
    @test nkeys(as_turing_model(MA(), 30), "θ") == 1
    @test nkeys(as_turing_model(MA(; θ = RandomWalk()), 20), "θ") > 1
end

@testitem "wired params: process prior yields a length-n path" begin
    using ComposableTuringIDModels, Distributions, Random
    Random.seed!(2)
    # HierarchicalNormal with a (positive) time-varying scale returns length n
    tv_scale = TransformLatentModel(RandomWalk(), x -> exp.(x))
    @test length(as_turing_model(HierarchicalNormal(; std = tv_scale), 15)()) == 15
    # MA(1) with a time-varying coefficient returns length n
    @test length(as_turing_model(MA(; θ = RandomWalk()), 15)()) == 15
end

@testitem "time-varying observation-error params differentiate (ForwardDiff)" begin
    using ComposableTuringIDModels, Distributions, Turing, Random
    using DynamicPPL: LogDensityFunction, VarInfo, link, getlogjoint
    import LogDensityProblems as LDP
    import DifferentiationInterface as DI
    Random.seed!(3)
    grad_finite(model, n) = begin
        y = as_turing_model(model, missing, n)().generated_y_t
        m = as_turing_model(model, y, n)
        vi = link(VarInfo(m), m)
        ldf = LogDensityFunction(m, getlogjoint, vi)
        # A representative point, as the package AD harness uses (ADFixtures),
        # not the all-zeros origin: a squared time-varying `cluster_factor`
        # is exactly 0 there (the negative binomial's degenerate Poisson
        # limit), a measure-zero singularity the sampler never visits.
        θ = 0.3 .* randn(MersenneTwister(1), LDP.dimension(ldf))
        grad = DI.gradient(
            x -> LDP.logdensity(ldf, x), DI.AutoForwardDiff(), θ)
        all(isfinite, grad) && length(grad) == length(θ)
    end
    # time-varying negative-binomial overdispersion
    @test grad_finite(
        IDModel(DirectInfections(; Z = RandomWalk(), initialisation = Normal()),
            NegativeBinomialError(; cluster_factor = RandomWalk())), 10)
    # time-varying Gaussian observation noise (kept positive by a log transform)
    tv_sd = TransformLatentModel(RandomWalk(), x -> exp.(x))
    @test grad_finite(
        IDModel(DirectInfections(; Z = RandomWalk(), initialisation = Normal()),
            NormalError(; std = tv_sd)), 10)
end

@testitem "time-varying latent params differentiate (ForwardDiff)" begin
    using ComposableTuringIDModels, Distributions, Turing, Random
    using DynamicPPL: LogDensityFunction, VarInfo, link, getlogjoint
    import LogDensityProblems as LDP
    import DifferentiationInterface as DI
    Random.seed!(4)
    grad_finite(m) = begin
        vi = link(VarInfo(m), m)
        ldf = LogDensityFunction(m, getlogjoint, vi)
        grad = DI.gradient(x -> LDP.logdensity(ldf, x),
            DI.AutoForwardDiff(), zeros(LDP.dimension(ldf)))
        all(isfinite, grad)
    end
    # time-varying innovation scale (stochastic volatility)
    tv_scale = TransformLatentModel(RandomWalk(), x -> exp.(x))
    @test grad_finite(as_turing_model(HierarchicalNormal(; std = tv_scale), 12))
    # time-varying MA(1) coefficient
    @test grad_finite(as_turing_model(MA(; θ = RandomWalk()), 12))
end
