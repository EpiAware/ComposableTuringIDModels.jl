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
    @test implements_prior_interface(gp; n = 25)
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
    using ComposableTuringIDModels: hsgp_basis, spectral_density
    using KernelFunctions: with_lengthscale, kernelmatrix
    n, σ, ℓ, c = 20, 1.0, 0.5, 2.0
    x = standardised_index(n)
    K_exact = kernelmatrix(σ^2 * with_lengthscale(SqExponentialKernel(), ℓ), x)
    Φ, sqrt_λ = hsgp_basis(n, 40, c)
    sd = sqrt.(spectral_density(SqExponentialKernel(), sqrt_λ, σ, ℓ))
    K_approx = Φ * Diagonal(sd .^ 2) * Φ'
    @test size(Φ) == (n, 40)
    @test norm(K_approx - K_exact) / norm(K_exact) < 0.05
end

@testitem "HilbertSpaceGP rejects invalid m and c" begin
    using ComposableTuringIDModels, Distributions
    @test_throws AssertionError HilbertSpaceGP(; m = 0)
    @test_throws AssertionError HilbertSpaceGP(; c = 1.0)
    @test_throws AssertionError as_turing_model(HilbertSpaceGP(), 1)()
end

@testitem "GP latents reject a hyperprior that strays off its support" begin
    using ComposableTuringIDModels, Distributions
    # An unbounded prior is caught at construction. Left to the sampler the
    # first negative proposal would abort the chain from inside the spectral
    # density (HSGP) or the Cholesky factorisation (exact), far from the cause.
    for GP in (HilbertSpaceGP, ExactGP)
        @test_throws AssertionError GP(; length_scale = Normal())
        @test_throws AssertionError GP(;
            length_scale = truncated(Normal(), -1, Inf))
        @test_throws AssertionError GP(; marginal_std = Normal())
        # the shipped defaults, and a hand-rolled prior whose support is open
        # at zero, both pass
        @test GP() isa GP
        @test GP(; length_scale = Gamma(2, 0.2),
            marginal_std = Exponential(1.0)) isa GP
    end
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
        @test implements_prior_interface(gp; n = n)
        path = as_turing_model(gp, n)()
        @test length(path) == n
        @test all(isfinite, path)
    end
end

@testitem "HilbertSpaceGP spectral densities are positive, finite, and kernel-specific" begin
    using ComposableTuringIDModels
    using ComposableTuringIDModels: spectral_density
    ω = collect(range(0, 5; length = 12))
    σ, ℓ = 1.0, 1.0
    for K in (SqExponentialKernel(), Matern32Kernel(), Matern52Kernel())
        S = spectral_density(K, ω, σ, ℓ)
        @test length(S) == length(ω)
        @test all(>(0), S)
        @test all(isfinite, S)
        @test S[1] > S[end]
    end
    @test spectral_density(Matern32Kernel(), ω, σ, ℓ) !=
          spectral_density(SqExponentialKernel(), ω, σ, ℓ)
end

@testitem "spectral_density rejects a zero length scale rather than returning NaN" begin
    using ComposableTuringIDModels
    using ComposableTuringIDModels: spectral_density
    # ν = √(2p+1)/ℓ is Inf at ℓ = 0 and the Matern densities become Inf/Inf,
    # so an unguarded call would hand a NaN basis weight to the sampler.
    for K in (SqExponentialKernel(), Matern32Kernel(), Matern52Kernel())
        @test_throws AssertionError spectral_density(K, [1.0], 1.0, 0.0)
        @test_throws AssertionError spectral_density(K, [1.0], 1.0, -0.5)
        @test_throws AssertionError spectral_density(K, [1.0], -1.0, 1.0)
    end
end

@testitem "HilbertSpaceGP Matern bases approximate their kernel covariance" begin
    using ComposableTuringIDModels, LinearAlgebra
    using ComposableTuringIDModels: hsgp_basis, spectral_density
    using KernelFunctions: with_lengthscale, kernelmatrix
    n, σ, ℓ, c = 20, 1.0, 0.8, 3.0
    x = standardised_index(n)
    K_exact = kernelmatrix(σ^2 * with_lengthscale(Matern52Kernel(), ℓ), x)
    Φ, sqrt_λ = hsgp_basis(n, 60, c)
    sd = sqrt.(spectral_density(Matern52Kernel(), sqrt_λ, σ, ℓ))
    K_approx = Φ * Diagonal(sd .^ 2) * Φ'
    @test norm(K_approx - K_exact) / norm(K_exact) < 0.15
end

@testitem "HilbertSpaceGP captures its basis as a model argument" begin
    using ComposableTuringIDModels, Distributions, Random
    using ComposableTuringIDModels: hsgp_basis
    using DynamicPPL: DynamicPPL
    Random.seed!(7)
    n, m, c = 30, 12, 1.5
    gp = HilbertSpaceGP(; m = m, c = c)
    mdl = as_turing_model(gp, n)
    @test mdl isa DynamicPPL.Model
    # The basis is an ARGUMENT of the returned model, not something the model
    # body builds: `as_turing_model` computes it once and captures it, so no
    # log-density evaluation differentiates through the sinusoids.
    Φ, sqrt_λ = hsgp_basis(n, m, c)
    @test mdl.args.Φ == Φ
    @test mdl.args.sqrt_λ == sqrt_λ
    @test length(mdl()) == n
end

@testitem "HilbertSpaceGP samples in the DEFAULT ℓ/m regime" tags=[:sample] begin
    using ComposableTuringIDModels, Distributions, Turing, Random
    Random.seed!(8)
    gp = HilbertSpaceGP()
    @test gp.m == 20
    @test minimum(gp.length_scale) > 0
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

@testitem "HilbertSpaceGP composes into a renewal and recovers the latent" tags=[:sample] begin
    using ComposableTuringIDModels, Distributions, Turing, Random
    using DynamicPPL: InitFromPrior
    using Statistics: cor, mean
    Random.seed!(13)
    # An epidemic large enough for the counts to identify the latent path: with
    # only a handful of cases a day the posterior collapses onto the prior and
    # nothing about recovery is testable.
    n, ndraws = 40, 150
    model = IDModel(
        Renewal(; generation_time = Gamma(6.5, 0.62),
            rt = HilbertSpaceGP(m = 12),
            initialisation = Normal(log(50), 0.1)),
        NegativeBinomialError(cluster_factor = HalfNormal(0.1)))
    prior = as_turing_model(model, fill(missing, n), n)
    sim = fix(prior, (ℓ = 0.5, σ = 0.5))()
    y_obs = sim.generated_y_t
    posterior = as_turing_model(model, y_obs, n)
    chain = sample(posterior, NUTS(0.9), ndraws;
        initial_params = InitFromPrior(), progress = false)
    gen = vec(generated_observables(posterior, y_obs, chain).generated)
    @test length(gen) == ndraws
    @test all(g -> length(g.Z_t) == n && all(isfinite, g.Z_t), gen)
    # A degenerate fit (σ initialised far out in the prior tail, step size
    # collapsed, chain stuck) still returns finite paths of the right length,
    # just wrong ones, so also check the posterior mean tracks the simulated
    # latent and stays on a plausible scale. Recovery holds well away from these
    # thresholds; a sweep over seeds gives correlations of 0.86 to 0.98 and a
    # largest posterior-mean log Rt of about 1.1.
    Z_mean = vec(mean(reduce(hcat, (g.Z_t for g in gen)); dims = 2))
    @test cor(Z_mean, sim.Z_t) > 0.7
    @test maximum(abs, Z_mean) < 3
    # Every reported-count index must come back from `predict`: the credible
    # bands in the Gaussian-process case study read `y_t[i]` for every `i`, and
    # a missing index would otherwise show up as a silently shortened band.
    pred = predict(prior, chain)
    @test all(i -> length(vec(pred[@varname(y_t[i])])) == ndraws, 1:n)
end

@testitem "ExactGP draws a length-n path with named GP parameters" begin
    using ComposableTuringIDModels, Distributions, Random
    Random.seed!(9)
    gp = ExactGP()
    @test gp isa AbstractLatentModel
    @test implements_prior_interface(gp; n = 25)
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
    using DynamicPPL: fix
    using KernelFunctions: with_lengthscale, kernelmatrix
    using Statistics: cov, mean
    Random.seed!(10)
    n, σ, ℓ = 10, 1.0, 0.6
    x = standardised_index(n)
    K_exact = kernelmatrix(σ^2 * with_lengthscale(SqExponentialKernel(), ℓ), x)
    mdl = fix(as_turing_model(ExactGP(), n), (ℓ = ℓ, σ = σ))
    draws = reduce(hcat, (mdl() for _ in 1:4000))
    K_emp = cov(draws; dims = 2)
    @test norm(K_emp - K_exact) / norm(K_exact) < 0.1
end

@testitem "ExactGP rejects invalid jitter and n" begin
    using ComposableTuringIDModels, Distributions
    @test_throws AssertionError ExactGP(; jitter = 0.0)
    @test_throws AssertionError as_turing_model(ExactGP(), 1)()
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
        @test implements_prior_interface(gp; n = n)
        path = as_turing_model(gp, n)()
        @test length(path) == n
        @test all(isfinite, path)
    end
end

@testitem "ExactGP samples in the DEFAULT regime" tags=[:sample] begin
    using ComposableTuringIDModels, Distributions, Turing, Random
    Random.seed!(12)
    gp = ExactGP()
    @test minimum(gp.length_scale) > 0
    n = 20
    chn = sample(as_turing_model(gp, n), NUTS(), 40; progress = false)
    @test size(chn, 1) == 40
    @test all(isfinite, Array(chn))
    @test all(>(0), vec(chn[:ℓ]))
    path = as_turing_model(gp, n)()
    @test length(path) == n
    @test all(isfinite, path)
end

@testitem "ExactGP factorises at extreme marginal standard deviations" begin
    using ComposableTuringIDModels, Distributions, Random
    using DynamicPPL: fix
    Random.seed!(16)
    # The nugget is relative to σ². With a fixed absolute nugget the Cholesky
    # factorisation throws a `PosDefException` once σ is large enough that the
    # covariance is only numerically indefinite — mid-chain, which kills the
    # sampler rather than returning a bad draw.
    n = 40
    for σ in (0.0, 1e-8, 1.0, 1e3, 1e5), ℓ in (0.05, 0.5, 5.0)

        path = fix(as_turing_model(ExactGP(), n), (ℓ = ℓ, σ = σ))()
        @test length(path) == n
        @test all(isfinite, path)
    end
end

@testitem "ExactGP captures its input grid as a model argument" begin
    using ComposableTuringIDModels, Distributions, Random
    using DynamicPPL: DynamicPPL
    Random.seed!(14)
    n = 30
    mdl = as_turing_model(ExactGP(), n)
    @test mdl isa DynamicPPL.Model
    # The standardised grid depends only on `n`, so `as_turing_model` builds it
    # once and captures it rather than rebuilding it on every evaluation.
    @test mdl.args.x == standardised_index(n)
    @test length(mdl()) == n
end

@testitem "standardised_index is zero mean, unit sd, and guards n = 1" begin
    using ComposableTuringIDModels
    using Statistics: mean, std
    x = standardised_index(40)
    @test length(x) == 40
    @test mean(x)≈0 atol=1e-12
    @test std(x)≈1 atol=1e-12
    # The half-range approaches √3 from below as n grows, which is what makes a
    # fixed length-scale prior meaningful across series lengths.
    @test maximum(abs, x) < sqrt(3)
    @test maximum(abs, standardised_index(400)) > maximum(abs, x)
    # n = 1 has zero standard deviation and would give a NaN grid.
    @test_throws AssertionError standardised_index(1)
end

@testitem "ExactGP composes into a renewal and recovers the latent" tags=[:sample] begin
    using ComposableTuringIDModels, Distributions, Turing, Random
    using DynamicPPL: InitFromPrior
    using Statistics: cor, mean
    Random.seed!(15)
    # The exact GP counterpart of the HilbertSpaceGP composed fit above, kept
    # short because its covariance factorisation is O(n^3) per evaluation. The
    # Gaussian-process case study drives this same composition under Mooncake.
    n, ndraws = 30, 150
    model = IDModel(
        Renewal(; generation_time = Gamma(6.5, 0.62), rt = ExactGP(),
            initialisation = Normal(log(50), 0.1)),
        NegativeBinomialError(cluster_factor = HalfNormal(0.1)))
    prior = as_turing_model(model, fill(missing, n), n)
    sim = fix(prior, (ℓ = 0.5, σ = 0.5))()
    y_obs = sim.generated_y_t
    posterior = as_turing_model(model, y_obs, n)
    chain = sample(posterior, NUTS(0.9), ndraws;
        initial_params = InitFromPrior(), progress = false)
    gen = vec(generated_observables(posterior, y_obs, chain).generated)
    @test length(gen) == ndraws
    @test all(g -> length(g.Z_t) == n && all(isfinite, g.Z_t), gen)
    Z_mean = vec(mean(reduce(hcat, (g.Z_t for g in gen)); dims = 2))
    @test cor(Z_mean, sim.Z_t) > 0.7
    @test maximum(abs, Z_mean) < 3
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
