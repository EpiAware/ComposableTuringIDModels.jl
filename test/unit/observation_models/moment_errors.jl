# Tests for the moment-parameterised observation-error models: the generic
# `ObservationError`, the two named location-scale families, and the guards
# they inherit from the moment core.

@testitem "ObservationError parameterises a family by its moments" begin
    using ComposableTuringIDModels, Distributions

    # Absolute spread: the sampled `σ` is the standard deviation.
    abs_err = ObservationError(Normal; sd = HalfNormal(0.1))
    @test !abs_err.relative
    d = observation_error(abs_err, 100.0, 2.0)
    @test mean(d) ≈ 100.0
    @test std(d) ≈ 2.0

    # Relative spread: it is a coefficient of variation instead.
    rel_err = ObservationError(LogNormal; sd = HalfNormal(0.1), relative = true)
    @test rel_err.relative
    d = observation_error(rel_err, 100.0, 0.2)
    @test mean(d) ≈ 100.0 rtol = 1.0e-10
    @test std(d) ≈ 20.0 rtol = 1.0e-10

    # Any family registered for `(:mean, :sd)`, which is the point: a new
    # family is a call, not a new type.
    gamma_err = ObservationError(Gamma; sd = HalfNormal(0.2), relative = true)
    d = observation_error(gamma_err, 50.0, 0.2)
    @test mean(d) ≈ 50.0 rtol = 1.0e-10
    @test std(d) ≈ 10.0 rtol = 1.0e-10
    @test minimum(d) ≈ 0.0
end

@testitem "ObservationError samples its spread through the one seam" begin
    using ComposableTuringIDModels, Distributions, Random
    using DynamicPPL: @varname

    oe = ObservationError(LogNormal; sd = HalfNormal(0.2), relative = true)
    mdl = as_turing_model(oe, missing, fill(100.0, 6))
    drawn = rand(Xoshiro(1), mdl)
    @test any(contains("σ"), string.(collect(keys(drawn))))
    @test all(>(0), mdl(Xoshiro(1)).y_t)

    # A process-valued prior makes the spread time-varying, with no other
    # change: one `σ` per time point rather than one constant. The process
    # has to keep the spread positive, as any spread prior does — an
    # unconstrained one draws a negative standard deviation, which is no
    # distribution and simulating from it says so.
    tv = ObservationError(
        LogNormal;
        sd = TransformLatentModel(HierarchicalNormal(), x -> exp.(x)),
        relative = true
    )
    drawn_tv = rand(Xoshiro(1), as_turing_model(tv, missing, fill(100.0, 6)))
    @test length(drawn_tv[@varname(σ.ϵ_t)]) == 6
    unconstrained = ObservationError(
        LogNormal; sd = HierarchicalNormal(), relative = true
    )
    @test_throws ArgumentError rand(
        Xoshiro(1), as_turing_model(unconstrained, missing, fill(100.0, 6))
    )
end

@testitem "LogNormalError is relative noise about the expected value" begin
    using ComposableTuringIDModels, Distributions

    lne = LogNormalError()
    @test lne.cv isa HalfNormal
    # The two moments the component states, whatever the log-scale parameters
    # they imply.
    for Y_t in (1.0, 100.0, 1.0e4), σ in (0.05, 0.2, 1.5)
        d = observation_error(lne, Y_t, σ)
        @test mean(d) ≈ Y_t rtol = 1.0e-10
        @test std(d) ≈ σ * Y_t rtol = 1.0e-10
        # The log-scale spread does not depend on the expected value, which is
        # what makes the noise relative.
        @test d.σ ≈ sqrt(log1p(σ^2)) rtol = 1.0e-12
    end
    # It is the same component the generic form builds.
    generic = ObservationError(LogNormal; sd = lne.cv, relative = true)
    @test observation_error(generic, 100.0, 0.2) == observation_error(lne, 100.0, 0.2)
end

@testitem "LogNormalError names its prior cv and simulates positive data" begin
    using ComposableTuringIDModels, Distributions, Random

    mdl = as_turing_model(LogNormalError(), missing, fill(50.0, 8))
    @test any(contains("cv"), string.(collect(keys(rand(Xoshiro(2), mdl)))))
    y = mdl(Xoshiro(2)).y_t
    @test all(>(0), y)
    @test length(y) == 8
end

@testitem "NormalError keeps its interface and gains the guard" begin
    using ComposableTuringIDModels, Distributions, Random

    # Unchanged to a caller: same field, same prior name, same distribution.
    ne = NormalError()
    @test ne.std isa HalfNormal
    @test observation_error(ne, 10.0, 2.0) == Normal(10.0, 2.0)
    drawn = rand(Xoshiro(3), as_turing_model(ne, missing, fill(10.0, 5)))
    @test any(contains("σ"), string.(collect(keys(drawn))))

    # A spread a diverging proposal drives non-finite now scores `-Inf`
    # instead of raising, and refuses to be drawn from.
    d = observation_error(ne, 10.0, Inf)
    @test logpdf(d, 9.0) == -Inf
    @test_throws ArgumentError rand(d)
end

@testitem "a moment error model rejects rather than raises on a bad proposal" begin
    using ComposableTuringIDModels, Distributions

    # Every reachable failure of the log-scale conversion, from a sampler that
    # has wandered: the ratio overflow, a non-finite expectation, a
    # non-positive spread. All score `-Inf`, so the proposal is rejected with
    # nothing raised on the differentiated path.
    lne = LogNormalError()
    for (Y_t, σ) in (
            (100.0, 10 * sqrt(floatmax(Float64))), (Inf, 0.2), (NaN, 0.2),
            (100.0, -0.2), (-100.0, 0.2),
        )
        d = observation_error(lne, Y_t, σ)
        @test logpdf(d, 50.0) == -Inf
        # Asking the same point for a value is a different question, and it
        # has no answer.
        @test_throws ArgumentError rand(d)
    end
end

@testitem "a moment error model scores and differentiates a series" begin
    using ComposableTuringIDModels, Distributions, ForwardDiff
    using DynamicPPL: LogDensityFunction, VarInfo, link, getlogjoint
    import LogDensityProblems as LDP

    y = [95.0, 104.0, 99.0, 110.0, 88.0]
    for m in (LogNormalError(), ObservationError(Gamma; relative = true))
        mdl = as_turing_model(m, y, fill(100.0, 5))
        ldf = LogDensityFunction(mdl, getlogjoint, link(VarInfo(mdl), mdl))
        x = zeros(LDP.dimension(ldf))
        @test isfinite(LDP.logdensity(ldf, x))
        g = ForwardDiff.gradient(v -> LDP.logdensity(ldf, v), x)
        @test all(isfinite, g)
        @test !all(iszero, g)
    end
end

@testitem "the moment error models infer concretely" begin
    using ComposableTuringIDModels, Distributions

    # `observation_error` runs once per scored time point inside the
    # differentiated log density, so the family has to be resolvable at
    # compile time rather than stored as a `UnionAll`.
    for m in (
            ObservationError(LogNormal), ObservationError(Normal),
            LogNormalError(), NormalError(),
        )
        rt = only(Base.return_types(observation_error, (typeof(m), Float64, Float64)))
        @test isconcretetype(rt) || rt isa Union
        @test !(rt === Any)
    end
    ft = fieldtype(typeof(ObservationError(LogNormal)), :dist)
    @test Base.isdispatchelem(ft)
    @test !(ft isa UnionAll)
end
