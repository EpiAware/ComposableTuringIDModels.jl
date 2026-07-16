@testitem "TimeVaryingAR returns a length-n numeric path" begin
    using ComposableTuringIDModels, Random
    Random.seed!(1)
    tv = TimeVaryingAR()
    @test tv isa AbstractLatentModel
    z = as_turing_model(tv, 12)()
    @test z isa AbstractVector
    @test length(z) == 12
    # n must exceed 1 (there is at least one transition)
    @test_throws Exception as_turing_model(TimeVaryingAR(), 1)()
end

@testitem "TVARStep threads a per-step coefficient through the recursion" begin
    using ComposableTuringIDModels
    using ComposableTuringIDModels: TVARStep, accumulate_scan
    z1 = 0.5
    ρ = [0.6, -0.3, 0.9, 0.1]
    ϵ = [0.1, -0.2, 0.05, 0.3]
    got = accumulate_scan(TVARStep(), z1, collect(zip(ρ, ϵ)))
    # explicit time-varying AR(1) recursion z_t = ρ_t z_{t-1} + ϵ_t
    want = [z1]
    for t in 1:4
        push!(want, ρ[t] * want[end] + ϵ[t])
    end
    @test got ≈ want
    @test length(got) == 5
end

@testitem "TimeVaryingAR leaves AR constant-coefficient semantics unchanged" begin
    using ComposableTuringIDModels, Distributions
    # A process-valued damping prior still gives an order-1 constant-coefficient
    # AR (a structured prior over the constant coefficient), NOT a time path.
    @test AR(; damp = RandomWalk()).p == 1
    # A length-2 vector damping prior still gives an order-2 AR (init length
    # must match the order).
    @test AR(;
        damp = [truncated(Normal(0.5, 0.1), 0, 1),
            truncated(Normal(0.2, 0.1), 0, 1)],
        init = [Normal(), Normal()]).p == 2
end

@testitem "TimeVaryingAR composes as a latent in the stack" begin
    using ComposableTuringIDModels, Distributions, Random
    Random.seed!(2)
    # It returns a numeric path, so it drops into a bare-vector latent slot such as
    # a Renewal's rt inside a composed IDModel.
    gen_int = [0.2, 0.3, 0.5]
    idmodel = IDModel(
        Renewal(gen_int; rt = TimeVaryingAR(), initialisation = Normal()),
        PoissonError())
    y = as_turing_model(idmodel, missing, 12)().generated_y_t
    @test length(y) == 12
end

@testitem "TimeVaryingAR differentiates under ForwardDiff" begin
    using ComposableTuringIDModels, Turing
    using DynamicPPL: LogDensityFunction, VarInfo, link, getlogjoint
    import LogDensityProblems as LDP
    import DifferentiationInterface as DI
    m = as_turing_model(TimeVaryingAR(), 15)
    vi = link(VarInfo(m), m)
    ldf = LogDensityFunction(m, getlogjoint, vi)
    θ = zeros(LDP.dimension(ldf))
    grad = DI.gradient(x -> LDP.logdensity(ldf, x), DI.AutoForwardDiff(), θ)
    @test length(grad) == length(θ)
    @test all(isfinite, grad)
end

@testitem "TimeVaryingAR recovers a time-varying damping path" tags=[:sample] begin
    using ComposableTuringIDModels, Distributions, Turing, Random, Statistics
    using Turing: to_submodel
    Random.seed!(80)
    # a genuine ρ_t ramp crossing zero, driving one observed series
    n = 50
    tgrid = range(0, 1; length = n - 1)
    ρ_true = 0.85 .- 1.15 .* tgrid
    z = zeros(n)
    z[1] = randn()
    for t in 2:n
        z[t] = ρ_true[t - 1] * z[t - 1] + 0.3 * randn()
    end
    @model function observe_path(y, n)
        latent ~ as_turing_submodel(TimeVaryingAR(), n)
        for t in 1:n
            y[t] ~ Normal(latent[t], 0.01)
        end
    end
    model = observe_path(z, n)
    chain = sample(model, NUTS(0.8; adtype = Turing.AutoForwardDiff()), 300;
        progress = false)
    # ρ is tracked (via `:=`) so it is recovered from the chain: `chain[:ρ]` is a
    # per-draw vector of the coefficient path.
    ρ_draws = reduce(hcat, vec(chain[:ρ]))   # (n-1) × draws
    ρ_mean = vec(mean(ρ_draws; dims = 2))
    @test length(ρ_mean) == n - 1
    @test cor(ρ_mean, ρ_true) > 0.6
end
