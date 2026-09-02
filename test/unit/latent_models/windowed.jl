@testitem "InWindow places a series on its window and zeroes the rest" begin
    using ComposableTuringIDModels
    @test broadcast_n(InWindow(20:30), 40, 1) == 11
    @test broadcast_rule(InWindow(2:3), [1.0, 2.0], 5, 1) == [0.0, 1.0, 2.0, 0.0, 0.0]
    # A window touching either end needs no special case.
    @test broadcast_rule(InWindow(1:2), [1.0, 2.0], 4, 1) == [1.0, 2.0, 0.0, 0.0]
    @test broadcast_rule(InWindow(3:4), [1.0, 2.0], 4, 1) == [0.0, 0.0, 1.0, 2.0]
    # The whole series is a window too.
    @test broadcast_rule(InWindow(1:3), [1.0, 2.0, 3.0], 3, 1) == [1.0, 2.0, 3.0]
    @test_throws AssertionError InWindow(3:2)
    @test_throws AssertionError InWindow(0:2)
    @test_throws AssertionError broadcast_rule(InWindow(3:4), [1.0, 2.0], 3, 1)
end

@testitem "broadcast_window confines an effect to its window" begin
    using ComposableTuringIDModels, Distributions, Random
    base = RandomWalk(; init = Normal(0, 0.1), ϵ_t = IID(Normal(0, 0.05)))
    win = 20:30
    outside = setdiff(1:40, win)

    # A known effect size lands exactly on the window.
    fixed = CombineLatentModels(
        [base, broadcast_window(FixedIntercept(2.0), win)], ["", "Window"]
    )
    Random.seed!(101)
    plain = as_turing_model(base, 40)()
    Random.seed!(101)
    shifted = as_turing_model(fixed, 40)()
    @test shifted .- plain ≈ 2.0 .* in.(1:40, Ref(win))

    # An estimated effect size contributes exactly zero outside the window.
    est = CombineLatentModels(
        [base, broadcast_window(Normal(1, 0.2), win)], ["", "Window"]
    )
    Random.seed!(102)
    combined = as_turing_model(est, 40)()
    Random.seed!(102)
    inner_only = as_turing_model(base, 40)()
    @test length(combined) == 40
    @test combined[outside] ≈ inner_only[outside]
    @test all(combined[win] .!= inner_only[win])
    # The effect draws one parameter under the prefix it was given.
    @test string.(keys(rand(as_turing_model(est, 40)))) ==
        ["init", "ϵ_t", "Window.intercept"]
end

@testitem "broadcast_window accepts a window touching either end" begin
    using ComposableTuringIDModels, Distributions, Random
    base = RandomWalk(; init = Normal(0, 0.1), ϵ_t = IID(Normal(0, 0.05)))
    # This is where a segment-partitioning composition cannot go, because the
    # segment on the untouched side has length zero.
    for win in (1:10, 31:40, 1:40)
        model = CombineLatentModels(
            [base, broadcast_window(FixedIntercept(3.0), win)], ["", "Window"]
        )
        Random.seed!(103)
        plain = as_turing_model(base, 40)()
        Random.seed!(103)
        shifted = as_turing_model(model, 40)()
        @test shifted .- plain ≈ 3.0 .* in.(1:40, Ref(win))
    end
end

@testitem "a windowed process is generated over the window only" begin
    using ComposableTuringIDModels, Distributions, Random
    # A process member varies within the window and draws its own parameters
    # under the prefix, sized by the window rather than by the series.
    inner = RandomWalk(; init = Normal(), ϵ_t = IID(Normal(0, 0.1)))
    model = CombineLatentModels(
        [FixedIntercept(0.0), broadcast_window(inner, 20:30)], ["", "Window"]
    )
    Random.seed!(104)
    path = as_turing_model(model, 40)()
    @test length(path) == 40
    @test all(iszero, path[setdiff(1:40, 20:30)])
    @test length(unique(path[20:30])) > 1
    names = string.(keys(rand(as_turing_model(model, 40))))
    @test "Window.init" in names
    @test "Window.ϵ_t" in names
end

@testitem "broadcast_window shows its window in the model tree" begin
    using ComposableTuringIDModels, Distributions
    tree = sprint(
        show, MIME"text/plain"(), broadcast_window(Normal(1, 0.2), 20:30)
    )
    @test occursin("BroadcastLatentModel(InWindow(20:30))", tree)
    # The label carries the rule for every broadcast model, so two rules that
    # behave differently no longer print the same.
    weekly = sprint(show, MIME"text/plain"(), broadcast_weekly(RandomWalk()))
    @test occursin("BroadcastLatentModel(RepeatBlock)", weekly)
end

@testitem "the from-scratch composition matches broadcast_window" begin
    using ComposableTuringIDModels, Distributions, Random
    # An index-reading transform summed on by `CombineLatentModels` is the
    # route the helper replaces, and the docs still show it. For a constant
    # effect the two give the same path.
    base = RandomWalk(; init = Normal(0, 0.1), ϵ_t = IID(Normal(0, 0.05)))
    win = 20:30
    helper = CombineLatentModels(
        [base, broadcast_window(FixedIntercept(2.0), win)], ["", "Window"]
    )
    manual = CombineLatentModels(
        [
            base,
            TransformLatentModel(
                FixedIntercept(2.0), x -> x .* in.(eachindex(x), Ref(win))
            ),
        ], ["", "Window"]
    )
    Random.seed!(105)
    from_helper = as_turing_model(helper, 40)()
    Random.seed!(105)
    from_manual = as_turing_model(manual, 40)()
    @test from_helper ≈ from_manual
    # With an estimated effect the manual route draws its parameter under the
    # same prefixed name, so a model can move to the helper without renaming.
    est_manual = CombineLatentModels(
        [
            base,
            TransformLatentModel(
                Intercept(Normal(1, 0.2)),
                x -> x .* in.(eachindex(x), Ref(win))
            ),
        ], ["", "Window"]
    )
    est_helper = CombineLatentModels(
        [base, broadcast_window(Normal(1, 0.2), win)], ["", "Window"]
    )
    @test string.(keys(rand(as_turing_model(est_manual, 40)))) ==
        string.(keys(rand(as_turing_model(est_helper, 40))))
end

@testitem "a windowed effect multiplies under an exp transform" begin
    using ComposableTuringIDModels, Distributions, Random
    # `CombineLatentModels` sums, so a multiplicative window is the same
    # composition in log space. The documented form is asserted here rather
    # than left as a claim on the page.
    base = RandomWalk(; init = Normal(0, 0.1), ϵ_t = IID(Normal(0, 0.05)))
    win = 20:30
    multiplied = TransformLatentModel(
        CombineLatentModels(
            [base, broadcast_window(FixedIntercept(log(2)), win)], ["", "Window"]
        ), x -> exp.(x)
    )
    plain = TransformLatentModel(base, x -> exp.(x))
    Random.seed!(106)
    with_window = as_turing_model(multiplied, 40)()
    Random.seed!(106)
    without = as_turing_model(plain, 40)()
    ratio = with_window ./ without
    @test ratio[win] ≈ fill(2.0, length(win))
    @test ratio[setdiff(1:40, win)] ≈ ones(40 - length(win))
end
