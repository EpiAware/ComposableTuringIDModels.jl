# Tests for `Replicate`: `n_strata` fully independent draws of the same path
# model, each prefixed by its stratum index, stacked into a strata x time
# matrix. Used on its own, and (see stratify.jl) in `Stratify`'s `across` slot.

@testitem "Replicate constructs and wraps a bare Distribution" begin
    using ComposableTuringIDModels, Distributions
    rep = Replicate(RandomWalk())
    @test rep isa Replicate
    @test rep isa AbstractLatentModel
    @test rep.model isa RandomWalk
    wrapped = Replicate(Normal(0, 1))
    @test wrapped.model isa Intercept
end

@testitem "Replicate returns a strata x time matrix" begin
    using ComposableTuringIDModels, Distributions, Random
    Random.seed!(603)
    rep = Replicate(RandomWalk())
    out = as_turing_model(rep, (4, 12))()
    @test size(out) == (4, 12)
    @test all(x -> x isa Real, out)
end

@testitem "each stratum's variables are prefixed by its index" begin
    using ComposableTuringIDModels, Distributions
    rep = Replicate(RandomWalk())
    mdl = as_turing_model(rep, (3, 8))
    names = string.(collect(keys(rand(mdl))))
    @test any(startswith(n, "stratum1.") for n in names)
    @test any(startswith(n, "stratum2.") for n in names)
    @test any(startswith(n, "stratum3.") for n in names)
end

@testitem "Replicate draws every stratum independently" begin
    using ComposableTuringIDModels, Distributions, Random, Statistics
    Random.seed!(604)
    rep = Replicate(RandomWalk())
    reps, n_time = 400, 8
    draws = [as_turing_model(rep, (2, n_time))() for _ in 1:reps]
    row1 = reduce(hcat, [d[1, :] for d in draws])   # n_time x reps
    row2 = reduce(hcat, [d[2, :] for d in draws])
    per_time_cor = [cor(row1[t, :], row2[t, :]) for t in 1:n_time]
    # No shared parameters between strata, so the two rows are uncorrelated.
    @test mean(abs.(per_time_cor)) < 0.15
end

@testitem "Replicate in Stratify's across slot overrides across_shape" begin
    using ComposableTuringIDModels, Distributions, Random
    using ComposableTuringIDModels: across_shape
    Random.seed!(605)
    strat = Stratify(RandomWalk(), Replicate(RandomWalk()))
    @test across_shape(strat.across, (3, 10)) == (3, 10)
    out = as_turing_model(strat, (3, 10))()
    @test size(out) == (3, 10)
end
