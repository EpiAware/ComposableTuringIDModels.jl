@testitem "continuous distributions discretise into valid PMFs (CensoredDistributions)" begin
    using ComposableTuringIDModels, Distributions
    # Renewal and LatentDelay discretise a continuous distribution into a PMF via
    # CensoredDistributions.double_interval_censored (internal `_discretised_pmf`).
    pmf = ComposableTuringIDModels._discretised_pmf(Gamma(2.0, 1.0); D = 10.0)
    @test isapprox(sum(pmf), 1.0)
    @test all(>=(0), pmf)

    renewal = Renewal(; gen_distribution = Gamma(2.0, 1.0), D_gen = 10.0)
    @test isapprox(sum(renewal.gen_int), 1.0)
    @test all(>=(0), renewal.gen_int)

    obs = LatentDelay(PoissonError(), truncated(Normal(5.0, 2.0), 0.0, Inf))
    @test isapprox(sum(obs.delay), 1.0)
    @test all(>=(0), obs.delay)
end

@testitem "expected_Rt inverts the renewal relationship" begin
    using ComposableTuringIDModels
    gen_int = [0.2, 0.3, 0.5]
    rt = expected_Rt(gen_int, [100.0, 200, 300, 400, 500])
    @test length(rt) == 2
    @test all(>(0), rt)
end

@testitem "DirectSample draws from the prior" begin
    using ComposableTuringIDModels, Distributions, Turing, Random
    Random.seed!(81)
    @model g() = (x ~ Normal())
    # apply_method wraps the solution in IDObservables; `.samples` is the
    # raw inference result.
    chain = apply_method(g(), DirectSample(; n_samples = 10))
    @test chain isa IDObservables
    @test chain.samples !== nothing
    single = apply_method(g(), DirectSample())
    @test haskey(single.samples, @varname(x))
end

@testitem "get_param_array reshapes a Chains into (draws, chains)" begin
    using ComposableTuringIDModels, Distributions, Turing, MCMCChains, Random
    Random.seed!(82)
    @model g() = (x ~ Normal())
    chn = MCMCChains.Chains(sample(g(), Prior(), MCMCSerial(), 3, 2; progress = false))
    A = get_param_array(chn)
    @test size(A) == (3, 2)
end
