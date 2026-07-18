@testitem "ExpGrowthRate generates a growth-rate path and maps it to infections" begin
    using ComposableTuringIDModels, Distributions, Random
    Random.seed!(41)
    egr = ExpGrowthRate(; rt = RandomWalk(), initialisation = Normal())
    out = as_turing_model(egr, 20)()
    @test length(out.I_t) == 20
    @test length(out.Z_t) == 20
    @test all(>(0), out.I_t)
end

@testitem "Renewal generates an Rt path and maps it to infections" begin
    using ComposableTuringIDModels, Distributions, Random
    Random.seed!(42)
    gen_int = [0.2, 0.3, 0.5]
    renewal = Renewal(; generation_time = gen_int, rt = RandomWalk(), initialisation = Normal())
    out = as_turing_model(renewal, 20)()
    @test length(out.I_t) == 20
    @test length(out.Z_t) == 20
    @test all(isfinite, out.I_t)
    @test all(>=(0), out.I_t)
end

@testitem "fixed generation_time (vector / distribution) bakes an interval" begin
    using ComposableTuringIDModels, Distributions, Random
    Random.seed!(420)
    # A pmf vector and a continuous distribution both bake a fixed interval and
    # renewal step (the inferred path stores a prior model and no step instead).
    r_vec = Renewal(; generation_time = [0.2, 0.3, 0.5], rt = RandomWalk(),
        initialisation = Normal())
    @test r_vec.gen_int isa AbstractVector
    @test !isnothing(r_vec.recurrent_step)
    r_dist = Renewal(; generation_time = Gamma(2.0, 1.0), D_gen = 10.0,
        rt = RandomWalk(), initialisation = Normal())
    @test r_dist.gen_int isa AbstractVector
    @test isapprox(sum(r_dist.gen_int), 1.0)
    @test !isnothing(r_dist.recurrent_step)
    for r in (r_vec, r_dist)
        out = as_turing_model(r, 20)()
        @test length(out.I_t) == 20
        @test all(isfinite, out.I_t)
    end
end

@testitem "uncertain generation_time infers the interval per draw" begin
    using ComposableTuringIDModels, Distributions, Random
    Random.seed!(44)
    gen = UncertainDelay(LogNormal,
        [Normal(1.9, 0.2), truncated(Normal(0.5, 0.2), 0, Inf)]; D = 14.0)
    renewal = Renewal(; generation_time = gen, rt = RandomWalk(),
        initialisation = Normal())
    # The inferred path holds the prior model and bakes no fixed interval/step.
    @test renewal.gen_int isa UncertainDelay
    @test isnothing(renewal.recurrent_step)
    out = as_turing_model(renewal, 20)()
    @test length(out.I_t) == 20
    @test length(out.Z_t) == 20
    @test all(isfinite, out.I_t)
    @test all(>=(0), out.I_t)
    # The generation interval's distribution parameters are inferred RVs
    # (namespaced under the `gen` slot).
    draw = rand(as_turing_model(renewal, 20))
    @test any(k -> startswith(string(k), "gen"), keys(draw))
end

@testitem "uncertain generation_time differentiates (ForwardDiff)" begin
    using ComposableTuringIDModels, Distributions, Random
    using DynamicPPL: LogDensityFunction, VarInfo, link, getlogjoint
    import LogDensityProblems as LDP
    import DifferentiationInterface as DI
    Random.seed!(45)
    gen = UncertainDelay(LogNormal,
        [Normal(1.0, 0.3), truncated(Normal(0.4, 0.2), 0, Inf)]; D = 6.0)
    model = IDModel(
        Renewal(; generation_time = gen, rt = RandomWalk(),
            initialisation = Normal()),
        NegativeBinomialError())
    y = as_turing_model(model, missing, 12)().generated_y_t
    m = as_turing_model(model, y, 12)
    vi = link(VarInfo(m), m)
    ldf = LogDensityFunction(m, getlogjoint, vi)
    dim = LDP.dimension(ldf)
    # A representative point, not the all-zeros origin (see the note in
    # test/unit/base/timevarying_params.jl).
    θ = 0.3 .* randn(MersenneTwister(1), dim)
    grad = DI.gradient(x -> LDP.logdensity(ldf, x), DI.AutoForwardDiff(), θ)
    @test all(isfinite, grad)
    @test length(grad) == dim
end

@testitem "infection models fix their latent to a deterministic path" begin
    using ComposableTuringIDModels, Distributions, Random
    using DynamicPPL: fix
    Random.seed!(421)
    gen_int = [0.2, 0.3, 0.5]
    # Pinning the latent to a known (log) Rt trajectory via a FixedIntercept
    # latent makes the renewal infection path deterministic given I₀ — the
    # standalone-style illustration under the folded interface.
    logR = log(1.5)
    renewal = Renewal(; generation_time = gen_int, rt = FixedIntercept(logR),
        initialisation = Normal())
    mdl = fix(as_turing_model(renewal, 30), (init_incidence = 0.0,))
    out = mdl()
    @test all(≈(logR), out.Z_t)
    @test all(>(0), out.I_t)
    # A constant Rt > 1 grows incidence.
    @test out.I_t[end] > out.I_t[1]
end

@testitem "growth-rate / reproduction-number conversions round-trip" begin
    using ComposableTuringIDModels
    w = [0.2, 0.3, 0.5]
    r = R_to_r(1.5, w)
    @test r_to_R(r, w) ≈ 1.5 rtol=1e-3
    # r and R move in the same direction.
    @test R_to_r(2.0, w) > R_to_r(1.2, w)
end

@testitem "composed Renewal model runs a short NUTS sample" tags=[:sample] begin
    using ComposableTuringIDModels, Distributions, Turing, Random
    Random.seed!(43)
    gen_int = [0.2, 0.3, 0.5]
    model = IDModel(
        Renewal(; generation_time = gen_int, rt = RandomWalk(), initialisation = Normal()),
        PoissonError())
    y = as_turing_model(model, missing, 20)().generated_y_t
    chn = sample(as_turing_model(model, y, 20), NUTS(), 30; progress = false)
    @test chn !== nothing
end
