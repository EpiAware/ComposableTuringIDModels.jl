# Gaussian-process latent models: the Hilbert-space approximation
# (`HilbertSpaceGP`), the exact GP it approximates (`ExactGP`), and the
# shared pieces they are built from (`standardised_index`, `hsgp_basis`,
# `spectral_density`).

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

@testitem "the Hilbert-space basis converges to the exact kernel Gram matrix" begin
    using ComposableTuringIDModels, LinearAlgebra
    using ComposableTuringIDModels: hsgp_basis, spectral_density
    using KernelFunctions: with_lengthscale, kernelmatrix
    # The approximation is Φ diag(S(√λ)) Φ'. Comparing it against the Gram
    # matrix KernelFunctions builds from the same kernel pins the eigenvalues,
    # the eigenfunction normalisation and the spectral-density constants at
    # once: a wrong constant anywhere leaves an error that does not shrink.
    # σ ≠ 1 deliberately: at σ = 1 a σ-versus-σ² mix-up in the spectral density
    # is exactly invisible, so a reconstruction test fixed at σ = 1 says nothing
    # about the marginal-variance convention.
    n, σ, ℓ = 60, 2.5, 0.5
    x = standardised_index(n)
    ladder = ((20, 1.5), (40, 2.0), (80, 3.0))
    function relerr(m, c, kernel)
        K_exact = kernelmatrix(σ^2 * with_lengthscale(kernel, ℓ), x)
        Φ, sqrt_λ = hsgp_basis(n, m, c)
        @test size(Φ) == (n, m)
        K_approx = Φ * Diagonal(spectral_density(kernel, sqrt_λ, σ, ℓ)) * Φ'
        return norm(K_approx - K_exact) / norm(K_exact)
    end
    # Tolerances are the measured error at the top of each ladder, rounded up an
    # order of magnitude. They are deterministic (no RNG anywhere here), and the
    # squared-exponential floor is set by the basis, not by round-off.
    for (kernel, tol) in (
            (SqExponentialKernel(), 1.0e-10),
            (Matern32Kernel(), 1.0e-3), (Matern52Kernel(), 1.0e-4),
        )
        errs = [relerr(m, c, kernel) for (m, c) in ladder]
        @test all(<(0), diff(errs))   # strictly better with more m and wider L
        @test errs[end] < tol
    end
end

@testitem "the Hilbert-space basis reproduces the marginal variance σ²" begin
    using ComposableTuringIDModels, LinearAlgebra
    using ComposableTuringIDModels: hsgp_basis, spectral_density
    # k(x, x) = σ², so the diagonal of Φ diag(S) Φ' must be σ². This is the
    # sharpest single check on the eigenfunction normalisation: a factor α in
    # φ_j rescales the whole diagonal by α².
    Φ, sqrt_λ = hsgp_basis(50, 100, 3.0)
    K(σ) = Φ * Diagonal(spectral_density(SqExponentialKernel(), sqrt_λ, σ, 0.5)) *
        Φ'
    for σ in (0.5, 1.0, 3.0)
        @test all(d -> isapprox(d, σ^2; rtol = 1.0e-6), diag(K(σ)))
    end
    # The spectral density carries σ² as a bare factor, so the covariance scales
    # exactly quadratically in σ.
    @test K(2.0) ≈ 4 .* K(1.0) rtol = 1.0e-12
end

@testitem "HilbertSpaceGP rejects invalid m and c" begin
    using ComposableTuringIDModels, Distributions
    @test_throws AssertionError HilbertSpaceGP(; m = 0)
    @test_throws AssertionError as_turing_model(HilbertSpaceGP(), 1)()
    # c below 1.2 leaves a boundary error floor no m can clear, so it is
    # rejected rather than accepted as an unusable configuration.
    @test_throws AssertionError HilbertSpaceGP(; c = 1.0)
    @test_throws AssertionError HilbertSpaceGP(; c = 1.1)
    @test HilbertSpaceGP(; c = 1.2) isa HilbertSpaceGP
end

@testitem "GP hyperpriors must be univariate distributions" begin
    using ComposableTuringIDModels, Distributions, LinearAlgebra
    # The support guard calls `minimum` and the model body draws with `~`, both
    # of which need a univariate distribution. Bounding the slot rejects a wrong
    # kind of prior at construction, rather than letting `minimum` fail with a
    # `MethodError` from inside the guard whose message never reaches the user.
    for GP in (HilbertSpaceGP, ExactGP)
        @test_throws TypeError GP(; length_scale = MvNormal(zeros(2), I))
        @test_throws TypeError GP(; marginal_std = MvNormal(zeros(2), I))
    end
end

@testitem "ExactGP's nugget stays negligible across the σ range" begin
    using ComposableTuringIDModels, LinearAlgebra
    using DynamicPPL: fix
    # The nugget is there to keep the Cholesky factorisation defined, not to
    # change the model. Scaled by σ² it stays a fixed relative perturbation; a
    # nugget with an absolute floor would dominate the covariance at small σ.
    n, ℓ, jitter = 15, 0.6, 1.0e-6
    mdl = as_turing_model(ExactGP(; jitter = jitter), n)
    for σ in (1.0e-3, 1.0e-2, 1.0, 1.0e2)
        L = reduce(
            hcat,
            begin
                    e = zeros(n)
                    e[j] = 1.0
                    fix(mdl, (ℓ = ℓ, σ = σ, z = e))()
                end for j in 1:n
        )
        inflation = maximum(diag(L * L')) / σ^2 - 1
        @test 0 <= inflation < 10 * jitter
    end
    # σ = 0 is degenerate but must still factorise rather than throw, leaving
    # only the absolute floor (√(jitter·eps) ≈ 1.5e-11) in the path.
    @test all(<(1.0e-8) ∘ abs, fix(mdl, (ℓ = ℓ, σ = 0.0))())
end

@testitem "GP latents reject a hyperprior that strays off its support" begin
    using ComposableTuringIDModels, Distributions
    # An unbounded prior is caught at construction. Left to the sampler the
    # first negative proposal would abort the chain from inside the spectral
    # density (HSGP) or the Cholesky factorisation (exact), far from the cause.
    for GP in (HilbertSpaceGP, ExactGP)
        @test_throws AssertionError GP(; length_scale = Normal())
        @test_throws AssertionError GP(;
            length_scale = truncated(Normal(), -1, Inf)
        )
        @test_throws AssertionError GP(; marginal_std = Normal())
        # the shipped defaults, and a hand-rolled prior whose support is open
        # at zero, both pass
        @test GP() isa GP
        @test GP(;
            length_scale = Gamma(2, 0.2),
            marginal_std = Exponential(1.0)
        ) isa GP
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

@testitem "spectral densities match the Fourier transform of their kernel" begin
    using ComposableTuringIDModels
    using ComposableTuringIDModels: spectral_density
    using KernelFunctions: with_lengthscale
    # The closed forms are only right if they are the Fourier transform of the
    # kernel KernelFunctions actually builds, under the convention the
    # Hilbert-space weights assume, S(ω) = ∫ k(τ) e^{-iωτ} dτ. Integrating that
    # numerically checks every constant (the √(2π)ℓ, the 4ν³, the 16ν⁵/3, and
    # the σ² and ℓ powers) against the kernel rather than against itself.
    function fourier_transform(kernel, ω, σ, ℓ)
        k = σ^2 * with_lengthscale(kernel, ℓ)
        τ = range(0.0, 25.0 * ℓ; length = 200_001)   # k(τ) is negligible beyond
        f = [k(t, 0.0) * cos(ω * t) for t in τ]      # k and cos are even, so
        return 2 * step(τ) * (sum(f) - (f[1] + f[end]) / 2)   # 2∫₀^∞ suffices
    end
    for kernel in (SqExponentialKernel(), Matern32Kernel(), Matern52Kernel()),
            ℓ in (0.4, 1.3), σ in (1.0, 2.5), ω in (0.0, 0.7, 3.1)
        @test spectral_density(kernel, ω, σ, ℓ) ≈ fourier_transform(
            kernel, ω, σ, ℓ
        ) rtol = 1.0e-6
    end
end

@testitem "HilbertSpaceGP converges to ExactGP as m grows and L widens" begin
    using ComposableTuringIDModels, LinearAlgebra
    using DynamicPPL: fix
    # Both models are linear in their standard-normal weights, so feeding unit
    # vectors recovers the loading matrix A (the Cholesky factor for `ExactGP`,
    # Φ diag(√S) for `HilbertSpaceGP`) and A A' is the implied covariance —
    # deterministic, with no Monte Carlo noise to loosen the tolerance for.
    n, ℓ, σ = 30, 0.5, 1.0
    function implied_cov(mdl, weight::Symbol, k)
        A = reduce(
            hcat,
            begin
                    e = zeros(k)
                    e[j] = 1.0
                    fix(mdl, (ℓ = ℓ, σ = σ, weight => e))()
                end for j in 1:k
        )
        return A * A'
    end
    ladder = ((10, 1.5), (20, 1.5), (40, 2.0), (120, 4.0))
    # `ExactGP` adds a relative nugget τ(σ² + 1) to its diagonal, so the two
    # cannot agree below about 2τ/σ² ≈ 2e-6 however large m gets.
    for (kernel, tol) in ((SqExponentialKernel(), 1.0e-5), (Matern52Kernel(), 1.0e-4))
        exact = implied_cov(as_turing_model(ExactGP(; kernel = kernel), n), :z, n)
        errs = map(ladder) do (m, c)
            hs = implied_cov(
                as_turing_model(HilbertSpaceGP(; m = m, c = c, kernel = kernel), n),
                :β, m
            )
            return norm(hs - exact) / norm(exact)
        end
        @test all(<(0), diff(collect(errs)))
        @test errs[end] < tol
    end
end

@testitem "HilbertSpaceGP degrades as documented when m or c is too small" begin
    using ComposableTuringIDModels, LinearAlgebra
    using ComposableTuringIDModels: hsgp_basis, spectral_density
    using KernelFunctions: with_lengthscale, kernelmatrix
    # `m` and `c` set the floor on accuracy at opposite ends of the length-scale
    # range, and neither fixes the other's failure. The model stays finite and
    # samples happily in both regimes, so a test has to look at the error.
    n = 60
    function relerr(m, c, ℓ)
        x = standardised_index(n)
        K_exact = kernelmatrix(with_lengthscale(SqExponentialKernel(), ℓ), x)
        Φ, sqrt_λ = hsgp_basis(n, m, c)
        K_approx = Φ *
            Diagonal(
            spectral_density(
                SqExponentialKernel(), sqrt_λ, 1.0,
                ℓ
            )
        ) * Φ'
        return norm(K_approx - K_exact) / norm(K_exact)
    end
    # A single basis function cannot represent the covariance at any ℓ.
    @test relerr(1, 1.5, 0.5) > 0.5
    # A short ℓ needs a larger m: the default m = 20 is 25% out at ℓ = 0.1.
    @test relerr(20, 1.5, 0.1) > 0.1
    @test relerr(60, 1.5, 0.1) < 1.0e-3
    # A long ℓ needs a larger c instead — adding basis functions does not help,
    # because the approximation is periodic on [-L, L]. Even at the smallest c
    # the constructor allows, a long ℓ has an error floor m cannot clear.
    @test relerr(400, 1.2, 1.5) > 0.05
    @test relerr(80, 4.0, 1.5) < 1.0e-8
    # The smallest well-defined basis (n = 2, m = 1) still evaluates.
    @test length(as_turing_model(HilbertSpaceGP(; m = 1), 2)()) == 2
end

@testitem "HilbertSpaceGP is under-dispersed at the default m for a short ℓ" begin
    using ComposableTuringIDModels, LinearAlgebra
    using ComposableTuringIDModels: hsgp_basis, spectral_density
    # A truncated basis can only lose variance, never add it, so where m is too
    # small for ℓ the implied marginal standard deviation falls below σ. This is
    # invisible in a plot of the fit but biases the inferred σ upward, so the
    # docstring quotes it and this pins the numbers it quotes.
    gp = HilbertSpaceGP()
    @test (gp.m, gp.c) == (20, 1.5)
    function implied_sd(m, ℓ, kernel)
        Φ, sqrt_λ = hsgp_basis(60, m, gp.c)
        K = Φ * Diagonal(spectral_density(kernel, sqrt_λ, 1.0, ℓ)) * Φ'
        return sqrt(sum(diag(K)) / size(K, 1))
    end
    for (kernel, at_02, at_01) in (
            (SqExponentialKernel(), 0.99, 0.89),
            (Matern32Kernel(), 0.95, 0.84), (Matern52Kernel(), 0.97, 0.86),
        )
        @test implied_sd(gp.m, 0.2, kernel) ≈ at_02 atol = 0.02
        @test implied_sd(gp.m, 0.1, kernel) ≈ at_01 atol = 0.02
        # It is m, not the kernel, that is the problem: a larger basis recovers
        # the marginal variance at the same ℓ.
        @test implied_sd(60, 0.1, kernel) > implied_sd(gp.m, 0.1, kernel)
        @test implied_sd(60, 0.2, kernel) > 0.99
    end
end

@testitem "a composed GP keeps its parameters unprefixed" begin
    using ComposableTuringIDModels, Distributions
    using DynamicPPL: VarInfo
    # `Renewal` threads its `rt` slot through the prefix-off submodel seam, so a
    # GP latent contributes flat `ℓ`, `σ` and weight names to the composed
    # model. The composed-fit test items below `fix` on those flat names, and
    # the case study reads them out of the chain.
    for (gp, weight) in ((HilbertSpaceGP(; m = 5), :β), (ExactGP(), :z))
        model = IDModel(
            Renewal(;
                generation_time = Gamma(6.5, 0.62), rt = gp,
                initialisation = Normal(log(50), 0.1)
            ),
            NegativeBinomialError()
        )
        names = Symbol.(
            string.(
                keys(
                    VarInfo(
                        as_turing_model(
                            model, fill(missing, 8), 8
                        )
                    )
                )
            )
        )
        @test :ℓ ∈ names
        @test :σ ∈ names
        @test weight ∈ names
    end
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

@testitem "HilbertSpaceGP samples in the DEFAULT ℓ/m regime" tags = [:sample] begin
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

@testitem "HilbertSpaceGP composes into a renewal and recovers the latent" tags = [:sample] begin
    using ComposableTuringIDModels, Distributions, Turing, Random
    using DynamicPPL: InitFromPrior
    using Statistics: cor, mean
    Random.seed!(13)
    # An epidemic large enough for the counts to identify the latent path: with
    # only a handful of cases a day the posterior collapses onto the prior and
    # nothing about recovery is testable.
    n, ndraws = 40, 150
    model = IDModel(
        Renewal(;
            generation_time = Gamma(6.5, 0.62),
            rt = HilbertSpaceGP(m = 12),
            initialisation = Normal(log(50), 0.1)
        ),
        NegativeBinomialError(cluster_factor = HalfNormal(0.1))
    )
    prior = as_turing_model(model, fill(missing, n), n)
    sim = fix(prior, (ℓ = 0.5, σ = 0.5))()
    y_obs = sim.generated_y_t
    posterior = as_turing_model(model, y_obs, n)
    chain = sample(
        posterior, NUTS(0.9), ndraws;
        initial_params = InitFromPrior(), progress = false
    )
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
    using ComposableTuringIDModels, LinearAlgebra
    using DynamicPPL: fix
    using KernelFunctions: with_lengthscale, kernelmatrix
    # The path is L z with L the Cholesky factor, so it is linear in z: feeding
    # unit vectors recovers L exactly and L L' is the prior covariance. That is
    # a deterministic reading of the covariance, so the tolerance can be set by
    # the nugget rather than by Monte Carlo error.
    n, ℓ = 20, 0.6
    x = standardised_index(n)
    # σ ≠ 1 deliberately: at σ = 1 the model would agree with the Gram matrix
    # just as well if it scaled the kernel by σ rather than σ².
    for (σ, kernel) in (
            (2.5, SqExponentialKernel()), (0.4, Matern32Kernel()),
            (1.7, Matern52Kernel()),
        )
        K_exact = kernelmatrix(σ^2 * with_lengthscale(kernel, ℓ), x)
        mdl = as_turing_model(ExactGP(; kernel = kernel), n)
        L = reduce(
            hcat,
            begin
                    e = zeros(n)
                    e[j] = 1.0
                    fix(mdl, (ℓ = ℓ, σ = σ, z = e))()
                end for j in 1:n
        )
        # `ExactGP` adds a relative nugget to the diagonal, so it reproduces the
        # Gram matrix only up to that; everything else must be exact.
        @test norm(L * L' - K_exact) / norm(K_exact) < 1.0e-5
    end
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

@testitem "ExactGP samples in the DEFAULT regime" tags = [:sample] begin
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
    for σ in (0.0, 1.0e-8, 1.0, 1.0e3, 1.0e5), ℓ in (0.05, 0.5, 5.0)

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
    @test mean(x) ≈ 0 atol = 1.0e-12
    @test std(x) ≈ 1 atol = 1.0e-12
    # The half-range approaches √3 from below as n grows, which is what makes a
    # fixed length-scale prior meaningful across series lengths.
    @test maximum(abs, x) < sqrt(3)
    @test maximum(abs, standardised_index(400)) > maximum(abs, x)
    # n = 1 has zero standard deviation and would give a NaN grid.
    @test_throws AssertionError standardised_index(1)
end

@testitem "ExactGP composes into a renewal and recovers the latent" tags = [:sample] begin
    using ComposableTuringIDModels, Distributions, Turing, Random
    using DynamicPPL: InitFromPrior
    using Statistics: cor, mean
    Random.seed!(15)
    # The exact GP counterpart of the HilbertSpaceGP composed fit above, kept
    # short because its covariance factorisation is O(n^3) per evaluation. The
    # Gaussian-process case study drives this same composition under Mooncake.
    n, ndraws = 30, 150
    model = IDModel(
        Renewal(;
            generation_time = Gamma(6.5, 0.62), rt = ExactGP(),
            initialisation = Normal(log(50), 0.1)
        ),
        NegativeBinomialError(cluster_factor = HalfNormal(0.1))
    )
    prior = as_turing_model(model, fill(missing, n), n)
    sim = fix(prior, (ℓ = 0.5, σ = 0.5))()
    y_obs = sim.generated_y_t
    posterior = as_turing_model(model, y_obs, n)
    chain = sample(
        posterior, NUTS(0.9), ndraws;
        initial_params = InitFromPrior(), progress = false
    )
    gen = vec(generated_observables(posterior, y_obs, chain).generated)
    @test length(gen) == ndraws
    @test all(g -> length(g.Z_t) == n && all(isfinite, g.Z_t), gen)
    Z_mean = vec(mean(reduce(hcat, (g.Z_t for g in gen)); dims = 2))
    @test cor(Z_mean, sim.Z_t) > 0.7
    @test maximum(abs, Z_mean) < 3
end
