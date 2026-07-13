@testitem "latent components generate length-n paths" begin
    using ComposableTuringIDModels, Distributions, Random
    Random.seed!(1)
    n = 12
    for m in (IID(Normal()), HierarchicalNormal(), RandomWalk(), AR(), MA(),
        Intercept(Normal()), FixedIntercept(2.0))
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
        damp = [truncated(Normal(0, 0.05), 0, 1),
            truncated(Normal(0, 0.05), 0, 1)],
        init = [Normal(), Normal()])
    @test ar2.p == 2
    @test length(as_turing_model(ar2, 10)()) == 10

    ma2 = MA(;
        θ = [truncated(Normal(0, 0.05), -1, 1),
        truncated(Normal(0, 0.05), -1, 1)])
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
    # Prior slots are prefixed at the call site (prefix-on `to_submodel`), so a
    # RandomWalk exposes its init and the inner HierarchicalNormal's std under a
    # namespace path (e.g. `rw_init.θ`, `ϵ_t.std.θ`).
    @test any(startswith("rw_init"), names)
    @test any(contains("std"), names)
end
