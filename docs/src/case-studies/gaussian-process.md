# [A Gaussian-process latent process](@id case-study-gp)

One of the design claims of the prototype is that *any* latent process — anything
implementing `as_turing_model(model, n)` and returning a length-`n` path — can
drive an infection model, without the latent process knowing anything about the
rest of the package. This case study makes good on that claim with a **Gaussian
process** (GP). A GP is a prior over functions, so it is a natural, flexible
choice for a smoothly varying quantity such as ``\log R_t``: instead of assuming
a parametric form (a random walk, an autoregression) we let the data choose the
shape, subject only to a smoothness assumption.

The obstacle is cost. An exact GP over ``n`` time points needs an ``n \times n``
covariance matrix and an ``O(n^3)`` Cholesky factorisation *at every leapfrog
step* of the sampler, which makes it impractical inside NUTS at realistic ``n``.
The package ships an **approximate** GP — the [`HilbertSpaceGP`](@ref) latent
model — that sidesteps the factorisation entirely and is fast and stable under
gradient-based sampling.

## The Hilbert-space approximation

The Hilbert-space approximate GP (HSGP) of [riutortmayol2023practical](@citet),
building on [solin2020hilbert](@citep), replaces the GP with a short, fixed sum of
basis functions. On an interval ``[-L, L]`` the eigenfunctions and eigenvalues of
the Laplacian are

```math
\phi_j(x) = \sqrt{\tfrac{1}{L}}\, \sin\!\Big(\sqrt{\lambda_j}\,(x + L)\Big),
\qquad \sqrt{\lambda_j} = \frac{\pi j}{2 L}, \qquad j = 1, \ldots, m,
```

and the GP draw is approximated as

```math
f(x) \approx \sum_{j=1}^{m} \phi_j(x)\, \sqrt{S(\sqrt{\lambda_j})}\; \beta_j,
\qquad \beta_j \sim \mathrm{Normal}(0, 1),
```

where ``S`` is the spectral density of the chosen kernel. For the
squared-exponential kernel with marginal standard deviation ``\sigma`` and length
scale ``\ell``,

```math
S(\omega) = \sigma^2 \sqrt{2\pi}\, \ell \, \exp\!\Big(-\tfrac{1}{2}\ell^2\omega^2\Big).
```

Two properties make this fast and sampler-friendly:

  - **The basis is fixed.** ``\phi_j`` and ``\lambda_j`` depend only on the series
    length ``n``, the number of basis functions ``m``, and the boundary factor
    ``c`` (the domain is extended to ``L = c\,S`` beyond the half-range ``S`` of
    the inputs). They depend on none of the sampled parameters, so the basis lies
    entirely outside the differentiated path. Each evaluation of the latent path is
    then a single ``n \times m`` matrix–vector product.
  - **It is non-centred.** The only sampled quantities are ``\ell``, ``\sigma``,
    and the ``m`` standard-normal weights ``\beta``. A non-centred parameterisation
    like this is exactly what NUTS handles well, and the gradient is cheap, so the
    model differentiates cleanly under reverse-mode AD — including
    [Mooncake](https://chalk-lab.github.io/Mooncake.jl/), which we use below.

The accuracy/speed trade-off is set by two numbers: more basis functions ``m``
resolve shorter length scales (at linear cost), and the boundary factor ``c``
pushes the artificial boundary far enough away that it does not distort the fit
near the data [riutortmayol2023practical](@citep).

## A GP latent process on its own

The latent process is a model in its own right. Constructing a
[`HilbertSpaceGP`](@ref) and giving it a series length produces a length-`n`
draw, just like any other latent component.

```@example gp
using EpiAwarePrototype, Distributions, Random
Random.seed!(202)

gp = HilbertSpaceGP(
    length_scale_prior = truncated(Normal(0.0, 0.5), 0, Inf),
    marginal_std_prior = truncated(Normal(0.0, 0.5), 0, Inf),
    m = 20, c = 1.5)
```

```@example gp
draw = as_turing_model(gp, 60)()
(length = length(draw), extrema = round.(extrema(draw), digits = 2))
```

Because the basis is fixed, the reconstructed covariance of the approximation
converges to the exact squared-exponential kernel as the number of basis
functions grows. The package exposes the basis builder
[`EpiAwarePrototype.hsgp_basis`](@ref) and spectral density
[`EpiAwarePrototype.se_spectral_density`](@ref) used internally, so we can check
the approximation directly against the kernel it targets:

```@example gp
using LinearAlgebra
n = 40
σ, ℓ, c = 1.0, 1.0, 2.0
x = collect(1:n) .- (n + 1) / 2
K_exact = [σ^2 * exp(-(xi - xj)^2 / (2ℓ^2)) for xi in x, xj in x]

Φ, sqrt_λ = EpiAwarePrototype.hsgp_basis(n, 40, c)
sd = sqrt.(EpiAwarePrototype.se_spectral_density(sqrt_λ, σ, ℓ))
K_approx = Φ * Diagonal(sd .^ 2) * Φ'

round(norm(K_approx - K_exact) / norm(K_exact), digits = 4)
```

The relative error is a fraction of a percent: with enough basis functions the
approximation reproduces the kernel it stands in for.

## Composing it into an infection model

To use the GP as the time-varying reproduction number we hand it to a
[`Renewal`](@ref) infection model as its `rt` latent process. Nothing about the
renewal process changes — it asks its latent slot for a length-`n` ``\log R_t``
path and gets one. The generation interval is a ``\mathrm{Gamma}(6.5, 0.62)``
serial interval discretised by [`EpiData`](@ref), and reported cases are
overdispersed counts via [`NegativeBinomialError`](@ref).

```@example gp
data = EpiData(gen_distribution = Gamma(6.5, 0.62))
renewal = Renewal(data; rt = gp, initialisation_prior = Normal(log(2.0), 0.1))
obs = NegativeBinomialError(cluster_factor_prior = HalfNormal(0.1))
model = EpiAwareModel(renewal, obs)
```

## Simulate

Passing `missing` observations turns the composed model into a prior simulator.
We draw one trajectory with a clear epidemic wave to treat as our observed data;
the returned generated quantities include the reported cases `generated_y_t`, the
latent infections `I_t`, and the GP latent process `Z_t = \log R_t`.

```@example gp
n = 70
local sim
for _ in 1:60
    global sim = as_turing_model(model, fill(missing, n), n)()
    (sum(sim.generated_y_t) > 500 && maximum(exp.(sim.Z_t)) > 1.3) && break
end
y_obs = sim.generated_y_t
Z_true = sim.Z_t
(total_cases = sum(y_obs), peak = maximum(y_obs),
    Rt_range = round.(extrema(exp.(Z_true)), digits = 2))
```

## Fit with Mooncake

Conditioning on the observed counts and sampling with NUTS recovers the
posterior. We choose [Mooncake](https://chalk-lab.github.io/Mooncake.jl/) as the
automatic-differentiation backend through `NUTS(; adtype = AutoMooncake(...))`.
The non-centred basis-weight parameterisation is well suited to reverse-mode AD,
so the gradient is cheap and the chain mixes well. A short run keeps the page
quick to build.

```@example gp
using Turing, Mooncake
using ADTypes: AutoMooncake

posterior = as_turing_model(model, y_obs, n)
chain = sample(posterior,
    NUTS(0.9; adtype = AutoMooncake(; config = nothing)), 300;
    progress = false)
nothing # hide
```

The sampled parameters are the GP hyperparameters — the length scale ``\ell`` and
marginal standard deviation ``\sigma`` — the ``m`` basis weights ``\beta``, and
the observation cluster factor. The hyperparameters and the cluster factor are
identified from the single simulated series:

```@example gp
using MCMCChains, Statistics
mc = MCMCChains.Chains(chain)
summarystats(mc[[:ℓ, :σ, :cluster_factor]])
```

## Recover the latent process

The reproduction number ``R_t = \exp(Z_t)`` is a generated quantity rather than a
sampled parameter, but because the GP latent is a deterministic function of the
sampled weights we can reconstruct it directly from the posterior draws of
``\ell``, ``\sigma`` and ``\beta`` using the same fixed basis. This is the whole
point of the basis-function form: the posterior over functions is just the
posterior over a handful of weights, pushed back through ``\Phi``.

```@example gp
Φ_post, sqrt_λ_post = EpiAwarePrototype.hsgp_basis(n, gp.m, gp.c)
ℓ_draws = vec(mc[:ℓ].data)
σ_draws = vec(mc[:σ].data)
β_draws = reduce(hcat, [vec(mc[Symbol("β[$j]")].data) for j in 1:gp.m])

Z_draws = map(1:length(ℓ_draws)) do d
    w = sqrt.(EpiAwarePrototype.se_spectral_density(sqrt_λ_post, σ_draws[d], ℓ_draws[d]))
    Φ_post * (w .* β_draws[d, :])
end
Z_post_mean = mean(Z_draws)

(correlation = round(cor(Z_post_mean, Z_true), digits = 2),
    rmse = round(sqrt(mean((Z_post_mean .- Z_true) .^ 2)), digits = 3))
```

The posterior-mean log-``R_t`` tracks the truth closely, and the implied
reproduction number recovers the simulated range:

```@example gp
(Rt_true = round.(extrema(exp.(Z_true)), digits = 2),
    Rt_posterior_mean = round.(extrema(exp.(Z_post_mean)), digits = 2))
```

The GP never had to know it was modelling a reproduction number, and the renewal
and observation models never had to know their latent process was a GP. The two
sides met only through the length-`n` latent contract — which is exactly the
composability the package is built around. Swapping the GP for an
[`AR`](@ref) or a [`RandomWalk`](@ref) latent is a one-line change to the
`rt` argument; the rest of the model is untouched.

!!! note "Prototype"
    This example uses a short sampler run and simulated data to stay fast to
    build. For a real analysis you would use more iterations, check convergence
    diagnostics, tune the number of basis functions `m` to the expected
    smoothness, and supply observed data.

## References

```@bibliography
Pages = ["gaussian-process.md"]
Canonical = false
```
