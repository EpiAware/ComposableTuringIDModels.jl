@testitem "continuous distributions discretise into valid PMFs (CensoredDistributions)" begin
    using ComposableTuringIDModels, Distributions
    # Renewal and LatentDelay discretise a continuous distribution into a PMF via
    # CensoredDistributions.double_interval_censored (internal `_discretised_pmf`).
    pmf = ComposableTuringIDModels._discretised_pmf(Gamma(2.0, 1.0); D = 10.0)
    @test isapprox(sum(pmf), 1.0)
    @test all(>=(0), pmf)

    renewal = Renewal(; generation_time = Gamma(2.0, 1.0), D_gen = 10.0)
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

@testitem "_safe_int_floor rejects non-finite input with a clear DomainError" begin
    using ComposableTuringIDModels: _safe_int_floor
    # A non-finite value must raise a `DomainError` naming the offending value,
    # not an `InexactError` from the eventual `BigInt` conversion.
    for x in (Inf, -Inf, NaN)
        @test_throws DomainError _safe_int_floor(x)
    end
    # Finite values, including ones outside the `Int` range, are unaffected.
    @test _safe_int_floor(3.7) == 3
    @test _safe_int_floor(exp(48.0)) == floor(BigInt, exp(48.0))
end

@testitem "SafePoisson and SafeNegativeBinomial reject non-finite rates" begin
    using ComposableTuringIDModels, Distributions, Random
    Random.seed!(140)
    # A non-finite `SafePoisson` rate must raise a `DomainError`.
    @test_throws DomainError rand(SafePoisson(Inf))
    @test_throws DomainError rand(SafePoisson(NaN))
    # `SafeNegativeBinomial` samples through an inner `SafePoisson` whose rate
    # is a Gamma mixing draw. Degenerate `r`/`p` make that draw non-finite, so
    # it reaches the same guard.
    @test_throws DomainError rand(SafeNegativeBinomial(Inf, 0.5))
    @test_throws DomainError rand(SafeNegativeBinomial(1.0, 0.0))
end

@testitem "quantile and mode reject non-finite parameters, not just rand" begin
    using ComposableTuringIDModels, Distributions
    # `quantile` (both types) and `SafeNegativeBinomial`'s `mode` delegate to
    # the wrapped `Distributions.Poisson` / `NegativeBinomial`, which floor
    # internally without passing through `_safe_int_floor`. They must raise the
    # same `DomainError` naming the offending value on a non-finite parameter.
    for λ in (Inf, NaN)
        @test_throws DomainError quantile(SafePoisson(λ), 0.5)
    end
    for (r, p) in ((Inf, 0.5), (NaN, 0.5), (1.0, Inf), (1.0, NaN))
        @test_throws DomainError quantile(SafeNegativeBinomial(r, p), 0.5)
        @test_throws DomainError mode(SafeNegativeBinomial(r, p))
    end
    # The guard must be transparent for finite parameters: each guarded entry
    # point still delegates to the wrapped distribution and returns its value.
    for q in (0.1, 0.5, 0.9)
        @test quantile(SafePoisson(4.0), q) == quantile(Poisson(4.0), q)
        @test quantile(SafeNegativeBinomial(3.0, 0.4), q) ==
            quantile(NegativeBinomial(3.0, 0.4), q)
    end
    @test mode(SafeNegativeBinomial(3.0, 0.4)) ==
        mode(NegativeBinomial(3.0, 0.4))
end
