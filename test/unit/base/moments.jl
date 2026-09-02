# Tests for the moment-matching core: the closed-form fast paths, the guards
# that keep an invalid moment pair off `rand` and `quantile`, and the soft cap.

@testitem "the fast draw agrees with the generic moment-matched path" begin
    using ComposableTuringIDModels, Distributions, ReparameterisedDistributions
    import ComposableTuringIDModels: _moment_draw

    # `_moment_draw` short-circuits the generic construct-then-quantile path
    # for the two closed-form families. It exists for speed, so it must agree.
    for m in (1.0e-3, 1.0, 500.0, 1.0e5), s in (0.1, 1.0, 50.0),
            z in (-3.0, -1.25, 0.0, 0.5, 2.5)
        s >= m && continue
        p = cdf(Normal(), z)
        @test _moment_draw(LogNormal, m, s, z) ≈ quantile(
            reparameterise(LogNormal; mean = m, sd = s, check_args = false), p
        ) rtol = 1.0e-12
        # `Normal` is already native in `(mean, sd)`, so it has no registered
        # reparameterisation; the closed form is checked against the native
        # distribution's own inverse instead.
        @test _moment_draw(Normal, m, s, z) ≈ quantile(Normal(m, s), p) rtol =
            1.0e-12
    end
end

@testitem "the closed-form draw keeps the tails the generic path loses" begin
    using ComposableTuringIDModels, Distributions, ReparameterisedDistributions
    import ComposableTuringIDModels: _moment_draw

    # The generic path round-trips through `cdf` then `quantile`, and the round
    # trip has no resolution left in the tails: `cdf` saturates at exactly 1
    # by `z = 8.3`, so the draw comes back as the support's endpoint, and it
    # underflows to zero on the way down. The closed forms are exact there.
    m, s = 100.0, 20.0
    generic(z) = quantile(
        reparameterise(LogNormal; mean = m, sd = s, check_args = false),
        cdf(Normal(), z)
    )
    for z in (8.3, 9.0, 40.0)
        @test generic(z) == Inf
        @test isfinite(_moment_draw(LogNormal, m, s, z))
    end
    @test generic(-40.0) == 0.0
    @test _moment_draw(LogNormal, m, s, -40.0) > 0.0
    # The loss sets in well before the endpoint is reached: three digits are
    # already gone at eight standard deviations.
    @test !isapprox(generic(8.0), _moment_draw(LogNormal, m, s, 8.0); rtol = 1.0e-3)
    # A family with no closed form has only the generic path, so it collapses.
    @test _moment_draw(Gamma, m, s, 9.0) == Inf
end

@testitem "invalid moments score -Inf rather than throwing" begin
    using ComposableTuringIDModels, Distributions
    import ComposableTuringIDModels: _moment_dist

    # A diverging proposal reaches arguments that throw. `logpdf` must be
    # `-Inf` so the point is rejected, and `rand` must not raise, because
    # imputing a `missing` observation samples from whatever this returns.
    invalid = (
        (LogNormal, -1.0, 1.0), (LogNormal, 1.0, -1.0), (LogNormal, NaN, 1.0),
        (LogNormal, Inf, 1.0), (LogNormal, 1.0, Inf), (Normal, Inf, 1.0),
        (Normal, 1.0, -1.0), (Gamma, -1.0, 1.0), (Gamma, 1.0, Inf),
    )
    for (family, m, s) in invalid
        d = _moment_dist(family, m, s)
        @test logpdf(d, 3.0) == -Inf
        @test isinf(rand(d))
    end

    # The log-scale conversion squares `sd / mean`, which overflows to `Inf`
    # once the ratio passes `sqrt(floatmax)` — and `Inf` then reaches the
    # constructor looking like a valid argument.
    big = sqrt(floatmax(Float64))
    @test isinf((10 * big)^2)
    d = _moment_dist(LogNormal, 1.0, 10 * big)
    @test logpdf(d, 3.0) == -Inf
    @test isinf(rand(d))
    # Just inside the overflow the pair is honoured rather than rejected.
    @test isfinite(logpdf(_moment_dist(LogNormal, 1.0, 0.1 * big), 3.0))
end

@testitem "an invalid draw is NaN rather than a DomainError" begin
    using ComposableTuringIDModels, Distributions
    import ComposableTuringIDModels: _moment_draw

    # `quantile` raises on an invalid pair exactly as `rand` does, so the
    # generic path is guarded; the closed forms are guarded with it so a
    # negative mean cannot reach `log`.
    for family in (LogNormal, Gamma)
        @test isnan(_moment_draw(family, -1.0, 1.0, 0.5))
        @test isnan(_moment_draw(family, 1.0, NaN, 0.5))
    end
    @test isnan(_moment_draw(Normal, NaN, 1.0, 0.5))
    # A valid pair still draws.
    @test _moment_draw(LogNormal, 10.0, 1.0, 0.0) ≈ 10.0 / sqrt(1.01)
end

@testitem "the soft cap holds and does not overflow" begin
    using ComposableTuringIDModels
    import ComposableTuringIDModels: _softplus, _soft_upper

    # `soft_upper(x, u, k) = u - softplus(u - x, k)`, tracking `x` well below
    # the cap and saturating at it above.
    @test _soft_upper(0.01, 0.5, 10.0) ≈ 0.01 atol = 2.0e-3
    @test _soft_upper(0.1, 0.5, 10.0) ≈ 0.1 atol = 2.0e-2
    @test _soft_upper(50.0, 0.5, 10.0) ≈ 0.5 atol = 1.0e-6
    # An infinite cap is no cap at all.
    @test _soft_upper(50.0, Inf, 10.0) == 50.0

    # The naive `log1p(exp(k x)) / k` overflows to `Inf` once `k x` passes
    # about 709, which a diverging proposal reaches. The stable form is
    # asymptotically `x` there.
    @test isinf(log1p(exp(800.0)))
    @test _softplus(80.0, 10.0) ≈ 80.0
    @test _softplus(1.0e6, 10.0) ≈ 1.0e6
    @test isfinite(_soft_upper(-1.0e6, 0.5, 10.0))
    @test _soft_upper(-1.0e6, 0.5, 10.0) ≈ -1.0e6
    # Below the overflow it is the naive form, exactly.
    @test _softplus(0.3, 10.0) ≈ log1p(exp(3.0)) / 10.0
end

@testitem "the moment core is type stable and gradient safe" begin
    using ComposableTuringIDModels, Distributions, ForwardDiff
    import ComposableTuringIDModels: _moment_draw, _moment_dist

    for family in (Normal, LogNormal)
        rt = only(
            Base.return_types(_moment_draw, (Type{family}, Float64, Float64, Float64))
        )
        @test rt === Float64
    end
    # The draw differentiates through both closed forms.
    for family in (Normal, LogNormal)
        g = ForwardDiff.derivative(m -> _moment_draw(family, m, 0.2 * m, 1.5), 40.0)
        @test isfinite(g)
        @test g > 0
    end
    # And the guarded distribution scores under AD without raising.
    g = ForwardDiff.derivative(
        σ -> logpdf(_moment_dist(LogNormal, 10.0, σ), 9.0), 2.0
    )
    @test isfinite(g)
end

@testitem "a family is a type, checked where it is passed" begin
    using ComposableTuringIDModels, Distributions

    # An instance handed in place of the family would otherwise be stored and
    # fail deep inside the moment solve, so it is rejected at the constructor.
    @test_throws ArgumentError InfectionNoise(; dist = LogNormal(0.0, 1.0))
    @test_throws ArgumentError ObservationError(Normal(0.0, 1.0))
    # The family itself is stored in the type domain, which is what keeps the
    # draw's return type concrete.
    @test fieldtype(typeof(InfectionNoise(; dist = Normal)), :dist) === Type{Normal}
end
