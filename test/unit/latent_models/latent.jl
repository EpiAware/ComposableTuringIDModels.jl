@testitem "latent components generate length-n paths" begin
    using EpiAwarePrototype, Distributions, Random
    Random.seed!(1)
    n = 12
    for m in (IID(Normal()), HierarchicalNormal(), RandomWalk(), AR(), MA(),
        Intercept(Normal()), FixedIntercept(2.0), HilbertSpaceGP())
        path = as_turing_model(m, n)()
        @test length(path) == n
    end
    # Null generates nothing.
    @test as_turing_model(Null(), n)() === nothing
end

@testitem "AR and MA respect their order via priors" begin
    using EpiAwarePrototype, Distributions, Random
    Random.seed!(2)
    ar2 = AR(;
        damp_priors = [truncated(Normal(0, 0.05), 0, 1),
            truncated(Normal(0, 0.05), 0, 1)],
        init_priors = [Normal(), Normal()])
    @test ar2.p == 2
    @test length(as_turing_model(ar2, 10)()) == 10

    ma2 = MA(;
        θ_priors = [truncated(Normal(0, 0.05), -1, 1),
        truncated(Normal(0, 0.05), -1, 1)])
    @test ma2.q == 2
    @test length(as_turing_model(ma2, 10)()) == 10
end

@testitem "DiffLatentModel composes an ARIMA-style latent process" begin
    using EpiAwarePrototype, Distributions, Random
    Random.seed!(3)
    arima = DiffLatentModel(; model = AR(), init_priors = [Normal(), Normal()])
    @test arima.d == 2
    path = as_turing_model(arima, 20)()
    @test length(path) == 20
    @test all(isfinite, path)
end

@testitem "HilbertSpaceGP draws a length-n path with named GP parameters" begin
    using EpiAwarePrototype, Distributions, Random
    Random.seed!(5)
    gp = HilbertSpaceGP(; m = 8)
    @test gp isa AbstractLatentModel
    @test implements_latent_interface(gp; n = 25)
    n = 25
    path = as_turing_model(gp, n)()
    @test length(path) == n
    @test all(isfinite, path)
    # The sampled parameters are the length scale, the marginal sd, and the
    # m basis weights — no inner error-model variables leak through.
    draw = rand(as_turing_model(gp, n))
    pairs_dict = Dict(string(k) => v for (k, v) in pairs(draw))
    @test haskey(pairs_dict, "ℓ")
    @test haskey(pairs_dict, "σ")
    @test haskey(pairs_dict, "β")
    @test length(pairs_dict["β"]) == 8
end

@testitem "HilbertSpaceGP basis approximates the squared-exponential kernel" begin
    using EpiAwarePrototype, Distributions, LinearAlgebra
    using EpiAwarePrototype: hsgp_basis, se_spectral_density
    # With enough basis functions the reconstructed covariance
    # Φ diag(S(√λ)) Φ' converges to the exact SE-kernel Gram matrix.
    n, σ, ℓ, c = 20, 1.0, 1.0, 2.0
    x = collect(1:n) .- (n + 1) / 2
    K_exact = [σ^2 * exp(-(xi - xj)^2 / (2ℓ^2)) for xi in x, xj in x]
    Φ, sqrt_λ = hsgp_basis(n, 40, c)
    sd = sqrt.(se_spectral_density(sqrt_λ, σ, ℓ))
    K_approx = Φ * Diagonal(sd .^ 2) * Φ'
    @test size(Φ) == (n, 40)
    @test norm(K_approx - K_exact) / norm(K_exact) < 0.05
end

@testitem "HilbertSpaceGP rejects invalid m and c" begin
    using EpiAwarePrototype, Distributions
    @test_throws AssertionError HilbertSpaceGP(; m = 0)
    @test_throws AssertionError HilbertSpaceGP(; c = 1.0)
    # n must exceed 1 for a meaningful basis.
    @test_throws Exception as_turing_model(HilbertSpaceGP(), 1)()
end

@testitem "rand from a latent model uses flat (unprefixed) names" begin
    using EpiAwarePrototype, Distributions, Random
    using DynamicPPL: VarName
    Random.seed!(4)
    draw = rand(as_turing_model(RandomWalk(), 10))
    names = string.(collect(keys(draw)))
    # to_submodel(..., false) keeps inner variable names flat: a RandomWalk
    # exposes its init and the inner HierarchicalNormal's std/ϵ_t without a
    # path prefix.
    @test "rw_init" in names
    @test any(startswith("std"), names) || "std" in names
end
