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
    data = IDData([0.2, 0.3, 0.5], exp)
    renewal = Renewal(data; rt = RandomWalk(), initialisation = Normal())
    out = as_turing_model(renewal, 20)()
    @test length(out.I_t) == 20
    @test length(out.Z_t) == 20
    @test all(isfinite, out.I_t)
    @test all(>=(0), out.I_t)
end

@testitem "infection models fix their latent to a deterministic path" begin
    using ComposableTuringIDModels, Distributions, Random
    using DynamicPPL: fix
    Random.seed!(421)
    data = IDData([0.2, 0.3, 0.5], exp)
    # Pinning the latent to a known (log) Rt trajectory via a FixedIntercept
    # latent makes the renewal infection path deterministic given I₀ — the
    # standalone-style illustration under the folded interface.
    logR = log(1.5)
    renewal = Renewal(data; rt = FixedIntercept(logR),
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
    data = IDData([0.2, 0.3, 0.5], exp)
    model = IDModel(
        Renewal(data; rt = RandomWalk(), initialisation = Normal()),
        PoissonError())
    y = as_turing_model(model, missing, 20)().generated_y_t
    chn = sample(as_turing_model(model, y, 20), NUTS(), 30; progress = false)
    @test chn !== nothing
end
