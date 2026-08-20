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

@testitem "_discretised_pmf takes D from a truncated distribution's support" begin
    using ComposableTuringIDModels, Distributions
    # When `D` is not supplied, a finite `maximum(dist)` (e.g. a
    # caller-truncated distribution) should set the horizon directly rather
    # than being re-derived from a quantile, which silently shortens the PMF.
    gen = Gamma(1.5625, 1.92) + 1.0 # mean 3, sd 2.4, shifted
    truncated_gen = truncated(gen, nothing, 15.0)

    explicit = ComposableTuringIDModels._discretised_pmf(gen; Δd = 1.0, D = 15.0)
    from_support = ComposableTuringIDModels._discretised_pmf(truncated_gen; Δd = 1.0)
    explicit_upper = ComposableTuringIDModels._discretised_pmf(
        truncated_gen; Δd = 1.0, upper = 1.0
    )

    @test length(explicit) == 15
    @test length(from_support) == 15
    @test length(explicit_upper) == 15
    @test isapprox(explicit, from_support; atol = 1.0e-12)
    @test isapprox(explicit, explicit_upper; atol = 1.0e-12)
end

@testitem "_discretised_pmf rounds a misaligned support bound up to Δd" begin
    using ComposableTuringIDModels, Distributions
    # A bounded distribution's `maximum` need not land on a Δd grid point
    # (e.g. truncating at 15.5 with Δd = 1.0). The PMF horizon must still be
    # a multiple of Δd, or the bins `0, Δd, …, D-Δd` misalign with the
    # right-truncation applied inside `double_interval_censored` and the
    # final bin covers a partial, non-Δd-wide interval.
    gen = Gamma(1.5625, 1.92) + 1.0
    truncated_gen = truncated(gen, nothing, 15.5)

    pmf = ComposableTuringIDModels._discretised_pmf(truncated_gen; Δd = 1.0)

    # Horizon rounds 15.5 up to 16.0, giving 16 unit-width bins, so no mass
    # between the true bound (15.5) and the rounded horizon (16.0) is lost:
    # the distribution has no density past its own maximum.
    @test length(pmf) == 16
    @test isapprox(sum(pmf), 1.0)
    @test all(>=(0), pmf)

    # A non-unit Δd that also does not divide the bound exactly.
    truncated_small = truncated(gen, nothing, 15.5)
    pmf_small = ComposableTuringIDModels._discretised_pmf(
        truncated_small; Δd = 0.3
    )
    # ceil(15.5 / 0.3) = 52, so the horizon is 15.6 and there are 52 bins.
    @test length(pmf_small) == 52
    @test isapprox(sum(pmf_small), 1.0)
    @test all(>=(0), pmf_small)
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
