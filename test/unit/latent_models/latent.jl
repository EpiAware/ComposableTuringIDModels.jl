@testitem "latent components generate length-n paths" begin
    using ComposableTuringIDModels, Distributions, Random
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

@testitem "HilbertSpaceGP draws a length-n path with named GP parameters" begin
    using ComposableTuringIDModels, Distributions, Random
    Random.seed!(5)
    gp = HilbertSpaceGP(; m = 8)
    @test gp isa AbstractLatentModel
    @test implements_latent_interface(gp; n = 25)
    n = 25
    path = as_turing_model(gp, n)()
    @test length(path) == n
    @test all(isfinite, path)
    draw = rand(as_turing_model(gp, n))
    pairs_dict = Dict(string(k) => v for (k, v) in pairs(draw))
    @test haskey(pairs_dict, "ℓ")
    @test haskey(pairs_dict, "σ")
    @test haskey(pairs_dict, "β")
    @test length(pairs_dict["β"]) == 8
end

@testitem "HilbertSpaceGP basis approximates the squared-exponential kernel" begin
    using ComposableTuringIDModels, Distributions, LinearAlgebra
    using ComposableTuringIDModels: hsgp_basis, se_spectral_density,
                                    _hsgp_standardised_index
    using KernelFunctions: with_lengthscale, kernelmatrix
    n, σ, ℓ, c = 20, 1.0, 0.5, 2.0
    x = _hsgp_standardised_index(n)
    K_exact = kernelmatrix(σ^2 * with_lengthscale(SqExponentialKernel(), ℓ), x)
    Φ, sqrt_λ = hsgp_basis(n, 40, c)
    sd = sqrt.(se_spectral_density(sqrt_λ, σ, ℓ))
    K_approx = Φ * Diagonal(sd .^ 2) * Φ'
    @test size(Φ) == (n, 40)
    @test norm(K_approx - K_exact) / norm(K_exact) < 0.05
end

@testitem "HilbertSpaceGP rejects invalid m and c" begin
    using ComposableTuringIDModels, Distributions
    @test_throws AssertionError HilbertSpaceGP(; m = 0)
    @test_throws AssertionError HilbertSpaceGP(; c = 1.0)
    @test_throws Exception as_turing_model(HilbertSpaceGP(), 1)()
end

@testitem "HilbertSpaceGP supports squared-exponential and Matern kernels" begin
    using ComposableTuringIDModels, Distributions, Random
    using KernelFunctions: Kernel
    Random.seed!(6)
    n = 25
    @test HilbertSpaceGP().kernel isa SqExponentialKernel
    for K in (SqExponentialKernel(), Matern32Kernel(), Matern52Kernel())
        @test K isa Kernel
        gp = HilbertSpaceGP(; m = 10, kernel = K)
        @test gp.kernel === K
        @test implements_latent_interface(gp; n = n)
        path = as_turing_model(gp, n)()
        @test length(path) == n
        @test all(isfinite, path)
    end
end

@testitem "HilbertSpaceGP spectral densities are positive, finite, and kernel-specific" begin
    using ComposableTuringIDModels
    using ComposableTuringIDModels: spectral_density, se_spectral_density
    ω = collect(range(0, 5; length = 12))
    σ, ℓ = 1.0, 1.0
    for K in (SqExponentialKernel(), Matern32Kernel(), Matern52Kernel())
        S = spectral_density(K, ω, σ, ℓ)
        @test length(S) == length(ω)
        @test all(>(0), S)
        @test all(isfinite, S)
        @test S[1] > S[end]
    end
    @test se_spectral_density(ω, σ, ℓ) ≈
          spectral_density(SqExponentialKernel(), ω, σ, ℓ)
    @test spectral_density(Matern32Kernel(), ω, σ, ℓ) !=
          spectral_density(SqExponentialKernel(), ω, σ, ℓ)
end

@testitem "HilbertSpaceGP Matern bases approximate their kernel covariance" begin
    using ComposableTuringIDModels, LinearAlgebra
    using ComposableTuringIDModels: hsgp_basis, spectral_density,
                                    _hsgp_standardised_index
    using KernelFunctions: with_lengthscale, kernelmatrix
    n, σ, ℓ, c = 20, 1.0, 0.8, 3.0
    x = _hsgp_standardised_index(n)
    K_exact = kernelmatrix(σ^2 * with_lengthscale(Matern52Kernel(), ℓ), x)
    Φ, sqrt_λ = hsgp_basis(n, 60, c)
    sd = sqrt.(spectral_density(Matern52Kernel(), sqrt_λ, σ, ℓ))
    K_approx = Φ * Diagonal(sd .^ 2) * Φ'
    @test norm(K_approx - K_exact) / norm(K_exact) < 0.15
end

@testitem "HilbertSpaceGP builds its basis once, outside the model body" begin
    using ComposableTuringIDModels, Distributions, Random
    using DynamicPPL: DynamicPPL
    Random.seed!(7)
    gp = HilbertSpaceGP(; m = 12)
    mdl = as_turing_model(gp, 30)
    @test mdl isa DynamicPPL.Model
    @test length(mdl()) == 30
    @test length(mdl()) == 30
end

@testitem "HilbertSpaceGP samples in the DEFAULT ℓ/m regime" tags=[:sample] begin
    using ComposableTuringIDModels, Distributions, Turing, Random
    Random.seed!(8)
    gp = HilbertSpaceGP()
    @test gp.m == 20
    @test minimum(gp.length_scale_prior) > 0
    n = 30
    chn = sample(as_turing_model(gp, n), NUTS(), 40; progress = false)
    @test size(chn, 1) == 40
    @test all(isfinite, Array(chn))
    ℓ = vec(chn[:ℓ])
    @test all(>(0), ℓ)
    path = as_turing_model(gp, n)()
    @test length(path) == n
    @test all(isfinite, path)
end

@testitem "ExactGP draws a length-n path with named GP parameters" begin
    using ComposableTuringIDModels, Distributions, Random
    Random.seed!(9)
    gp = ExactGP()
    @test gp isa AbstractLatentModel
    @test implements_latent_interface(gp; n = 25)
    n = 25
    path = as_turing_model(gp, n)()
    @test length(path) == n
    @test all(isfinite, path)
    draw = rand(as_turing_model(gp, n))
    pairs_dict = Dict(string(k) => v for (k, v) in pairs(draw))
    @test haskey(pairs_dict, "ℓ")
    @test haskey(pairs_dict, "σ")
    @test haskey(pairs_dict, "z")
    @test length(pairs_dict["z"]) == n
end

@testitem "ExactGP prior covariance is the exact kernel Gram matrix" begin
    using ComposableTuringIDModels, Distributions, Random, LinearAlgebra
    using ComposableTuringIDModels: _hsgp_standardised_index
    using DynamicPPL: fix
    using KernelFunctions: with_lengthscale, kernelmatrix
    using Statistics: cov, mean
    Random.seed!(10)
    n, σ, ℓ = 10, 1.0, 0.6
    x = _hsgp_standardised_index(n)
    K_exact = kernelmatrix(σ^2 * with_lengthscale(SqExponentialKernel(), ℓ), x)
    mdl = fix(as_turing_model(ExactGP(), n), (ℓ = ℓ, σ = σ))
    draws = reduce(hcat, (mdl() for _ in 1:4000))
    K_emp = cov(draws; dims = 2)
    @test norm(K_emp - K_exact) / norm(K_exact) < 0.1
end

@testitem "ExactGP rejects invalid jitter and n" begin
    using ComposableTuringIDModels, Distributions
    @test_throws AssertionError ExactGP(; jitter = 0.0)
    @test_throws Exception as_turing_model(ExactGP(), 1)()
end

@testitem "ExactGP supports squared-exponential and Matern kernels" begin
    using ComposableTuringIDModels, Distributions, Random
    using KernelFunctions: Kernel
    Random.seed!(11)
    n = 25
    @test ExactGP().kernel isa SqExponentialKernel
    for K in (SqExponentialKernel(), Matern32Kernel(), Matern52Kernel())
        gp = ExactGP(; kernel = K)
        @test gp.kernel === K
        @test implements_latent_interface(gp; n = n)
        path = as_turing_model(gp, n)()
        @test length(path) == n
        @test all(isfinite, path)
    end
end

@testitem "ExactGP samples in the DEFAULT regime" tags=[:sample] begin
    using ComposableTuringIDModels, Distributions, Turing, Random
    Random.seed!(12)
    gp = ExactGP()
    @test minimum(gp.length_scale_prior) > 0
    n = 20
    chn = sample(as_turing_model(gp, n), NUTS(), 40; progress = false)
    @test size(chn, 1) == 40
    @test all(isfinite, Array(chn))
    @test all(>(0), vec(chn[:ℓ]))
    path = as_turing_model(gp, n)()
    @test length(path) == n
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
    # path (e.g. `rw_init.θ`); the inner HierarchicalNormal's `std` is a flat
    # native-tilde scalar draw.
    @test any(startswith("rw_init"), names)
    @test any(contains("std"), names)
end
