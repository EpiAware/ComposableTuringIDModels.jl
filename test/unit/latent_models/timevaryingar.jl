@testitem "TimeVaryingAR is a thin alias over AR(damp = process)" begin
    using ComposableTuringIDModels, Random
    Random.seed!(1)
    tv = TimeVaryingAR()
    # No separate type: TimeVaryingAR builds an AR with a process damp slot.
    @test tv isa AR
    @test tv isa AbstractLatentModel
    @test tv.damp isa RandomWalk
    @test tv.p == 1
    @test tv.transform === tanh
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
    # A vector coefficient is read per step via `_at`; the driving sequence is
    # `(t, ϵ_t)` pairs.
    got = accumulate_scan(TVARStep(ρ), z1, collect(zip(eachindex(ϵ), ϵ)))
    # explicit time-varying AR(1) recursion z_t = ρ_t z_{t-1} + ϵ_t
    want = [z1]
    for t in 1:4
        push!(want, ρ[t] * want[end] + ϵ[t])
    end
    @test got ≈ want
    @test length(got) == 5
    # A scalar coefficient is the constant AR(1): `_at` ignores `t`.
    gotc = accumulate_scan(TVARStep(0.6), z1, collect(zip(eachindex(ϵ), ϵ)))
    wantc = [z1]
    for t in 1:4
        push!(wantc, 0.6 * wantc[end] + ϵ[t])
    end
    @test gotc ≈ wantc
end

@testitem "AR damp slot: Distribution ⇒ constant, process ⇒ time-varying" begin
    using ComposableTuringIDModels, Distributions
    # A process-valued damping prior makes an order-1 AR genuinely time-varying
    # (a per-step coefficient path), while staying order 1.
    tvar = AR(; damp = RandomWalk())
    @test tvar.p == 1
    @test tvar.transform === tanh          # unbounded process ⇒ tanh band
    # A bare Distribution damping prior is a constant order-1 coefficient.
    cvar = AR(; damp = Normal())
    @test cvar.p == 1
    @test cvar.transform === identity      # bounded prior ⇒ used as-is
    # A length-2 vector damping prior still gives an order-2 (constant) AR (init
    # length must match the order).
    @test AR(;
        damp = [truncated(Normal(0.5, 0.1), 0, 1),
            truncated(Normal(0.2, 0.1), 0, 1)],
        init = [Normal(), Normal()]).p == 2
end

@testitem "AR constant damp stays a scalar RV (no length-n allocation)" begin
    using ComposableTuringIDModels, Distributions, Random
    Random.seed!(3)
    # The whole point of the scalar path: a Distribution damp is ONE scalar RV,
    # so the number of damping parameters does NOT grow with the series length.
    for nn in (10, 100)
        draw = rand(as_turing_model(AR(; damp = Normal()), nn))
        dk = [k for k in keys(draw) if occursin("damp_AR", string(k))]
        @test length(dk) == 1
        @test draw[only(dk)] isa Real          # scalar, not a length-n vector
    end
    # By contrast a process damp draws a length-(n-1) path, so its parameter count
    # DOES scale with n — the flexibility is there when asked for.
    Random.seed!(4)
    proc = rand(as_turing_model(AR(; damp = RandomWalk()), 30))
    pk = [k for k in keys(proc) if occursin("damp_AR", string(k))]
    @test length(pk) > 1
end

@testitem "TimeVaryingAR composes as a latent in the stack" begin
    using ComposableTuringIDModels, Distributions, Random
    Random.seed!(2)
    # It returns a numeric path, so it drops into a bare-vector latent slot such as
    # a Renewal's rt inside a composed IDModel.
    idmodel = IDModel(
        Renewal([0.2, 0.3, 0.5]; rt = TimeVaryingAR(),
            initialisation = Normal()),
        PoissonError())
    y = as_turing_model(idmodel, missing, 12)().generated_y_t
    @test length(y) == 12
end

@testitem "TimeVaryingAR differentiates under ForwardDiff" begin
    using ComposableTuringIDModels, Turing
    using DynamicPPL: LogDensityFunction, VarInfo, link, getlogjoint
    import LogDensityProblems as LDP
    import DifferentiationInterface as DI
    m = as_turing_model(AR(; damp = RandomWalk()), 15)
    vi = link(VarInfo(m), m)
    ldf = LogDensityFunction(m, getlogjoint, vi)
    θ = zeros(LDP.dimension(ldf))
    grad = DI.gradient(x -> LDP.logdensity(ldf, x), DI.AutoForwardDiff(), θ)
    @test length(grad) == length(θ)
    @test all(isfinite, grad)
end

@testitem "AR(damp = RandomWalk()) recovers a time-varying damping path" tags=[:sample] begin
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
        latent ~ as_turing_submodel(AR(; damp = RandomWalk()), n)
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
