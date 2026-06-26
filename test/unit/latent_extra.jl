@testitem "latent modifiers generate length-n paths" begin
    using EpiAwarePrototype, Distributions, Random
    Random.seed!(31)
    n = 10
    trans = TransformLatentModel(Intercept(Normal(2, 0.2)), x -> exp.(x))
    @test length(as_turing_model(trans, n)()) == n
    @test all(>(0), as_turing_model(trans, n)())  # exp transform is positive

    rec = RecordExpectedLatent(FixedIntercept(0.1))
    @test length(as_turing_model(rec, 5)()) == 5
end

@testitem "PrefixLatentModel prefixes inner variable names" begin
    using EpiAwarePrototype, Distributions, Random
    Random.seed!(32)
    pm = PrefixLatentModel(; model = HierarchicalNormal(), prefix = "Test")
    names = string.(collect(keys(rand(as_turing_model(pm, 10)))))
    @test all(startswith("Test."), names)
end

@testitem "CombineLatentModels sums components" begin
    using EpiAwarePrototype, Distributions, Random
    Random.seed!(33)
    combined = CombineLatentModels([Intercept(Normal(2, 0.2)), AR()])
    out = as_turing_model(combined, 10)()
    @test length(out) == 10
    @test all(isfinite, out)
end

@testitem "ConcatLatentModels concatenates segments" begin
    using EpiAwarePrototype, Distributions, Random
    Random.seed!(34)
    combined = ConcatLatentModels([Intercept(Normal(2, 0.2)), AR()])
    out = as_turing_model(combined, 10)()
    @test length(out) == 10
    @test all(isfinite, out)
    # equal_dimensions splits 10 across 2 models as [5, 5].
    @test EpiAwarePrototype.equal_dimensions(10, 2) == [5, 5]
end

@testitem "broadcasting expands a shorter process to length n" begin
    using EpiAwarePrototype, Random
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
    using EpiAwarePrototype, Distributions, Random
    Random.seed!(36)
    @test length(as_turing_model(arma(), 10)()) == 10
    @test length(as_turing_model(arima(), 10)()) == 10
    # arima is a differenced arma.
    @test arima() isa DiffLatentModel
end
