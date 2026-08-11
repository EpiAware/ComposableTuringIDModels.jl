# Tests for `Stratify`: adding a strata axis to a shared latent path by
# broadcasting a once-over-time draw against a once-over-strata (or
# once-over-both) draw, plus the `across_shape` extension point and the
# path-model guard it and every other strata-shaped call relies on.

@testitem "Stratify constructs and wraps a bare Distribution slot" begin
    using ComposableTuringIDModels, Distributions
    s = Stratify(RandomWalk(), Hierarchy())
    @test s isa Stratify
    @test s isa AbstractLatentModel
    @test s.combine === (+)
    # A bare `Distribution` in either slot is auto-wrapped in an `Intercept`,
    # the same PATH convention used elsewhere.
    wrapped = Stratify(Normal(0, 1), Normal(0, 2))
    @test wrapped.shared isa Intercept
    @test wrapped.across isa Intercept
end

@testitem "Stratify broadcasts shared and across via combine" begin
    using ComposableTuringIDModels, Distributions
    # Fully deterministic (no random draws), so the broadcast arithmetic is
    # checked exactly rather than statistically.
    default = Stratify(FixedIntercept(1.0), FixedIntercept(2.0))
    out = as_turing_model(default, (3, 4))()
    @test out == fill(3.0, 3, 4)   # default combine is `+`

    mult = Stratify(FixedIntercept(2.0), FixedIntercept(3.0); combine = *)
    outm = as_turing_model(mult, (2, 5))()
    @test outm == fill(6.0, 2, 5)
end

@testitem "across_shape: stratum count, or the full shape for Replicate" begin
    using ComposableTuringIDModels
    # `across_shape` is public but not exported, so it is imported by name.
    using ComposableTuringIDModels: across_shape
    @test across_shape(Hierarchy(), (4, 10)) == 4
    @test across_shape(FixedIntercept(0.0), (4, 10)) == 4
    @test across_shape(Replicate(RandomWalk()), (4, 10)) == (4, 10)
end

@testitem "a vector across is time-constant; Replicate varies over time" begin
    using ComposableTuringIDModels, Distributions, Random
    Random.seed!(601)
    # `Hierarchy` draws once per stratum (across_shape ⇒ n_strata), so the gap
    # between two strata's rows is the same at every time point.
    constant_offset = Stratify(RandomWalk(), Hierarchy())
    out1 = as_turing_model(constant_offset, (3, 15))()
    diffs1 = out1[1, :] .- out1[2, :]
    @test all(x -> isapprox(x, diffs1[1]; atol = 1.0e-8), diffs1)

    Random.seed!(602)
    # `Replicate` spans both axes (across_shape ⇒ full shape), so the gap
    # varies over time instead of staying fixed.
    time_varying_offset = Stratify(RandomWalk(), Replicate(RandomWalk()))
    out2 = as_turing_model(time_varying_offset, (3, 15))()
    diffs2 = out2[1, :] .- out2[2, :]
    @test !all(x -> isapprox(x, diffs2[1]; atol = 1.0e-6), diffs2)
end

@testitem "across_shape is a third-party extension point" begin
    using ComposableTuringIDModels, Distributions, Random, Turing
    using ComposableTuringIDModels: across_shape, AbstractPriorModel
    Random.seed!(603)
    # A user-defined `across` model, added here rather than in the package:
    # it needs the full `(n_strata, n_time)` shape, exactly like `Replicate`,
    # and adds its own `across_shape` method to say so. `Stratify` picks it
    # up unchanged — mirrors the "the pre-scan seam takes any modifier"
    # extension test for `AbstractRenewalModifier` in imported_cases.jl.
    struct TimeVaryingAcross{M} <: AbstractPriorModel
        model::M
    end
    ComposableTuringIDModels.across_shape(::TimeVaryingAcross, n::Dims{2}) = n
    Turing.@model function ComposableTuringIDModels.as_turing_model(
            m::TimeVaryingAcross, n::Dims{2}
        )
        out ~ as_turing_submodel(m.model, n)
        return out
    end

    @test across_shape(TimeVaryingAcross(Replicate(RandomWalk())), (4, 10)) ==
        (4, 10)

    strat = Stratify(RandomWalk(), TimeVaryingAcross(Replicate(RandomWalk())))
    out = as_turing_model(strat, (3, 15))()
    @test size(out) == (3, 15)
    diffs = out[1, :] .- out[2, :]
    @test !all(x -> isapprox(x, diffs[1]; atol = 1.0e-6), diffs)
end

@testitem "a path model at a strata shape errors, naming Stratify" begin
    using ComposableTuringIDModels, Distributions
    @test_throws "Stratify" as_turing_model(RandomWalk(), (3, 10))
    @test_throws "Stratify" as_turing_model(IID(Normal()), (2, 5))
end
