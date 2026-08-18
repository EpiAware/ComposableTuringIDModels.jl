@testitem "latent components generate length-n paths" begin
    using ComposableTuringIDModels, Distributions, Random
    Random.seed!(1)
    n = 12
    for m in (
            IID(Normal()), HierarchicalNormal(), RandomWalk(), AR(), MA(),
            Intercept(Normal()), FixedIntercept(2.0), HilbertSpaceGP(),
        )
        path = as_turing_model(m, n)()
        @test length(path) == n
    end
    # Null generates nothing.
    @test as_turing_model(Null(), n)() === nothing
end

@testitem "AR and MA respect their order via priors" begin
    using ComposableTuringIDModels, Distributions, Random
    Random.seed!(2)
    ar2 = AR(;
        damp = [
            truncated(Normal(0, 0.05), 0, 1),
            truncated(Normal(0, 0.05), 0, 1),
        ],
        init = [Normal(), Normal()]
    )
    @test ar2.p == 2
    @test length(as_turing_model(ar2, 10)()) == 10

    ma2 = MA(;
        θ = [
            truncated(Normal(0, 0.05), -1, 1),
            truncated(Normal(0, 0.05), -1, 1),
        ]
    )
    @test ma2.q == 2
    @test length(as_turing_model(ma2, 10)()) == 10
end

@testitem "DiffLatentModel composes an ARIMA-style latent process" begin
    using ComposableTuringIDModels, Distributions, Random
    Random.seed!(3)
    arima = DiffLatentModel(; model = AR(), init = [Normal(), Normal()])
    @test arima.d == 2
    path = as_turing_model(arima, 20)()
    @test length(path) == 20
    @test all(isfinite, path)
end

@testitem "rand from a latent model namespaces prior variables" begin
    using ComposableTuringIDModels, Distributions, Random
    using DynamicPPL: VarName
    Random.seed!(4)
    draw = rand(as_turing_model(RandomWalk(), 10))
    names = string.(collect(keys(draw)))
    # The init prior slot is prefixed at the call site (prefix-on
    # `as_turing_submodel`), so a RandomWalk exposes its init under a namespace
    # path (e.g. `init.θ`); the inner HierarchicalNormal's `std` is a flat
    # native-tilde scalar draw.
    @test any(startswith("init"), names)
    @test any(contains("std"), names)
end

@testitem "latent components name their parameters generically" begin
    using ComposableTuringIDModels, Distributions
    using DynamicPPL: VarInfo, DebugUtils
    # A parameter is named for what it is, not for the component that draws
    # it, so the same quantity reads the same way whichever process supplies
    # it and a component that appears twice is told apart by a prefix rather
    # than by a longer name.
    varnames(model, n) = Symbol.(
        string.(
            keys(
                VarInfo(
                    as_turing_model(model, n)
                )
            )
        )
    )

    ar = varnames(AR(), 20)
    @test :init ∈ ar
    @test :damp ∈ ar
    @test :ar_init ∉ ar
    @test :damp_AR ∉ ar

    @test :init ∈ varnames(RandomWalk(), 20)
    @test :mean ∈ varnames(Hierarchy(), 3)

    # Two processes that now share `init` are told apart by the prefix their
    # composition already applies.
    combined = varnames(CombineLatentModels([AR(), RandomWalk()]), 20)
    @test Symbol("Combine.1.init") ∈ combined
    @test Symbol("Combine.2.init") ∈ combined

    # `DiffLatentModel` keeps `latent_init` rather than taking `init` too: it
    # wraps a single inner model prefix-off, so a generic name there would
    # land on the inner process's own `init` and the shipped
    # `DiffLatentModel(model = AR())` would stop working out of the box.
    diffed = varnames(DiffLatentModel(model = AR()), 20)
    @test :latent_init ∈ diffed
    @test :init ∈ diffed
    @test DebugUtils.check_model(
        as_turing_model(DiffLatentModel(model = AR()), 20);
        error_on_failure = false
    )
end
