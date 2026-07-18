@testitem "latent modifiers generate length-n paths" begin
    using ComposableTuringIDModels, Distributions, Random
    Random.seed!(31)
    n = 10
    trans = TransformLatentModel(Intercept(Normal(2, 0.2)), x -> exp.(x))
    @test length(as_turing_model(trans, n)()) == n
    @test all(>(0), as_turing_model(trans, n)())  # exp transform is positive

    rec = RecordExpectedLatent(FixedIntercept(0.1))
    @test length(as_turing_model(rec, 5)()) == 5
end

@testitem "latent modifiers accept a raw member (Distribution / IID / process)" begin
    using ComposableTuringIDModels, Distributions, Random
    Random.seed!(310)
    n = 5

    # Each modifier struct stores a bare `Distribution` member raw (it draws a
    # single scalar through the single seam).
    @test TransformLatentModel(Normal(), x -> exp.(x)).model isa Distribution
    @test RecordExpectedLatent(Normal()).model isa Distribution
    @test PrefixLatentModel(Normal(), "Test").model isa Distribution
    @test PrefixLatentModel(; model = Normal(), prefix = "Test").model isa
          Distribution
    @test BroadcastLatentModel(Normal(), 7, RepeatEach()).model isa Distribution
    @test BroadcastLatentModel(Normal(); period = 7,
        broadcast_rule = RepeatEach()).model isa Distribution

    # For a length-n member use IID() (n i.i.d.) or a process; it composes like a
    # model and produces a length-n path.
    trans = TransformLatentModel(IID(Normal()), x -> exp.(x))
    @test length(as_turing_model(trans, n)()) == n
    rec = RecordExpectedLatent(IID(Normal()))
    @test length(as_turing_model(rec, n)()) == n
    bcast = BroadcastLatentModel(IID(Normal()), 7, RepeatEach())
    @test length(as_turing_model(bcast, 10)()) == 10

    # A vector of `Distribution`s is stored as the raw vector.
    vec_member = TransformLatentModel([Normal(), Normal()], identity)
    @test vec_member.model isa AbstractVector{<:Distribution}

    # A richer prior model is stored unchanged (no double-wrapping).
    proc = TransformLatentModel(RandomWalk(), identity)
    @test proc.model isa RandomWalk
end

@testitem "PrefixLatentModel prefixes inner variable names" begin
    using ComposableTuringIDModels, Distributions, Random
    Random.seed!(32)
    pm = PrefixLatentModel(; model = HierarchicalNormal(), prefix = "Test")
    names = string.(collect(keys(rand(as_turing_model(pm, 10)))))
    @test all(startswith("Test."), names)
end

@testitem "CombineLatentModels sums components" begin
    using ComposableTuringIDModels, Distributions, Random
    Random.seed!(33)
    combined = CombineLatentModels([Intercept(Normal(2, 0.2)), AR()])
    out = as_turing_model(combined, 10)()
    @test length(out) == 10
    @test all(isfinite, out)
end

@testitem "ConcatLatentModels concatenates segments" begin
    using ComposableTuringIDModels, Distributions, Random
    Random.seed!(34)
    combined = ConcatLatentModels([Intercept(Normal(2, 0.2)), AR()])
    out = as_turing_model(combined, 10)()
    @test length(out) == 10
    @test all(isfinite, out)
    # equal_dimensions splits 10 across 2 models as [5, 5].
    @test ComposableTuringIDModels.equal_dimensions(10, 2) == [5, 5]
end

@testitem "equal_dimensions always sums to n (>=3 segments)" begin
    using ComposableTuringIDModels
    # Regression: the remainder must be distributed across the leading segments,
    # not dumped onto the first, so the lengths always sum to exactly n. The old
    # `vcat(ceil(n/m), fill(floor(n/m), m-1))` gave [3,2,2,2]=9 for (10, 4) and
    # tripped ConcatLatentModels' `@assert sum(dims) == n`.
    ed = ComposableTuringIDModels.equal_dimensions
    @test ed(10, 4) == [3, 3, 2, 2]
    @test ed(11, 3) == [4, 4, 3]
    for (n, m) in ((10, 4), (11, 3), (13, 5), (7, 7), (100, 6), (5, 2))
        dims = ed(n, m)
        @test length(dims) == m
        @test sum(dims) == n
        @test all(>(0), dims)
        @test maximum(dims) - minimum(dims) <= 1
    end
end

@testitem "ConcatLatentModels builds and runs with >=3 segments" begin
    using ComposableTuringIDModels, Distributions, Random
    Random.seed!(340)
    # Exercises the default equal_dimensions on a partition with a non-trivial
    # remainder (10 across 3 -> [4, 3, 3]); previously threw at model build.
    combined = ConcatLatentModels([Intercept(Normal(2, 0.2)), AR(), RandomWalk()])
    out = as_turing_model(combined, 10)()
    @test length(out) == 10
    @test all(isfinite, out)
    @test ComposableTuringIDModels.equal_dimensions(10, 3) == [4, 3, 3]
end

@testitem "MAStep convolves innovations in correct lag order" begin
    using ComposableTuringIDModels: MAStep, accumulate_scan
    # MA(2): Z_t = ϵ_t + θ1 ϵ_{t-1} + θ2 ϵ_{t-2}; the first q outputs are the raw
    # warm-up innovations in natural order. Drives the same seeding MA.jl uses
    # (newest-first buffer via `reverse(ϵ[1:q])`).
    θ = [0.5, 0.2]
    q = 2
    ϵ = [1.0, 2.0, 3.0, 4.0, 5.0]
    ma = accumulate_scan(
        MAStep(θ), (; val = 0.0, state = reverse(ϵ[1:q])), ϵ[(q + 1):end])
    expected = [1.0, 2.0,
        3.0 + 0.5 * 2.0 + 0.2 * 1.0,
        4.0 + 0.5 * 3.0 + 0.2 * 2.0,
        5.0 + 0.5 * 4.0 + 0.2 * 3.0]
    @test ma ≈ expected
end

@testitem "ARStep applies damping in correct lag order" begin
    using ComposableTuringIDModels: ARStep, accumulate_scan
    # AR reverses the damping coefficients before building the step (see AR.jl)
    # so ρ[i] multiplies the lag-i term. With ρ = [ρ1, ρ2], init [Z1, Z2] and
    # zero innovations: Z3 = ρ1 Z2 + ρ2 Z1, Z4 = ρ1 Z3 + ρ2 Z2, …
    ρ = [0.6, 0.3]
    init = [0.0, 1.0]        # oldest→newest: Z1 = 0, Z2 = 1
    ϵ = [0.0, 0.0, 0.0]
    ar = accumulate_scan(ARStep(reverse(ρ)), init, ϵ)
    z3 = 0.6 * 1.0 + 0.3 * 0.0
    z4 = 0.6 * z3 + 0.3 * 1.0
    z5 = 0.6 * z4 + 0.3 * z3
    @test ar ≈ [0.0, 1.0, z3, z4, z5]
end

@testitem "broadcasting expands a shorter process to length n" begin
    using ComposableTuringIDModels, Random
    Random.seed!(35)
    each = BroadcastLatentModel(RandomWalk(), 7, RepeatEach())
    @test length(as_turing_model(each, 10)()) == 10
    block = BroadcastLatentModel(RandomWalk(), 3, RepeatBlock())
    @test length(as_turing_model(block, 10)()) == 10
    # Rule behaviour.
    @test broadcast_rule(RepeatEach(), [1, 2], 5, 2) == [1, 2, 1, 2, 1]
    @test broadcast_rule(RepeatBlock(), [1, 2], 4, 2) == [1, 1, 2, 2]
    @test length(as_turing_model(broadcast_dayofweek(RandomWalk()), 14)()) == 14
    @test length(as_turing_model(broadcast_weekly(RandomWalk()), 14)()) == 14
end

@testitem "arma and arima build composable latent processes" begin
    using ComposableTuringIDModels, Distributions, Random
    Random.seed!(36)
    @test length(as_turing_model(arma(), 10)()) == 10
    @test length(as_turing_model(arima(), 10)()) == 10
    # arima is a differenced arma.
    @test arima() isa DiffLatentModel
end
