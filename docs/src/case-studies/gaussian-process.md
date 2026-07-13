# [A Gaussian-process latent process](@id case-study-gp)

Any latent process that implements `as_turing_model(model, n)` and returns a
length-`n` path can drive an infection model. This case study uses a **Gaussian
process** (GP) — a prior over functions, and so a flexible choice for a smoothly
varying quantity such as ``\log R_t``: rather than assume a parametric form we
let the data choose the shape, subject only to a smoothness assumption from the
kernel.

An exact GP over ``n`` points needs an ``n \times n`` covariance factorisation
that costs ``O(n^3)`` *at every leapfrog step* of the sampler. The package ships
two GP latent models built on the same ecosystem kernels:

  - [`ExactGP`](@ref) — the exact GP. Accurate, but ``O(n^3)`` per evaluation, so
    best for short series. It is the accuracy reference.
  - [`HilbertSpaceGP`](@ref) — the Hilbert-space basis-function approximation of
    [riutortmayol2023practical](@citet), building on [solin2020hilbert](@citep).
    A fixed basis makes each evaluation an ``n \times m`` matrix–vector product,
    which is fast and stable under gradient-based sampling.

We build both, check them against the GP ecosystem, then fit both to the same
simulated data and compare accuracy and speed.

## Two GP latent models

Both are latent models in their own right: give one a series length and it
returns a length-`n` draw, like any other component.

```@example gp
using ComposableTuringIDModels, Distributions, Random
Random.seed!(202)

hsgp = HilbertSpaceGP(m = 20, c = 1.5)
exact = ExactGP()
(hsgp = length(as_turing_model(hsgp, 60)()),
    exact = length(as_turing_model(exact, 60)()))
```

Both sample a length scale ``\ell`` and marginal standard deviation ``\sigma``.
[`HilbertSpaceGP`](@ref) then draws ``m`` basis weights ``\beta`` and forms the
path from the fixed basis; [`ExactGP`](@ref) draws ``n`` non-centred weights
``z`` and pushes them through the Cholesky factor of the full covariance. Both
parameterisations are non-centred, which NUTS handles well.

## The kernel and the GP ecosystem

The kernels are the standard
[KernelFunctions.jl](https://juliagaussianprocesses.github.io/KernelFunctions.jl/)
types — the same ones
[AbstractGPs.jl](https://juliagaussianprocesses.github.io/AbstractGPs.jl/) builds
exact GPs from — so the model reuses the ecosystem's kernels rather than defining
its own. The default `SqExponentialKernel` gives very smooth paths;
`Matern32Kernel` and `Matern52Kernel` give progressively rougher ones. A kernel
enters [`HilbertSpaceGP`](@ref) only through its spectral density, so adding one
is a single `ComposableTuringIDModels.spectral_density(::MyKernel, ω, σ, ℓ)`
method; [`ExactGP`](@ref) uses the kernel's Gram matrix directly.

Because both models are built on the same kernel, we can check them against
AbstractGPs. On standardised inputs, [`ExactGP`](@ref)'s prior covariance *is* the
AbstractGPs Gram matrix by construction, and the [`HilbertSpaceGP`](@ref) basis
reconstructs it to a fraction of a percent with enough basis functions:

```@example gp
using LinearAlgebra
using KernelFunctions: with_lengthscale, kernelmatrix
using AbstractGPs: GP, cov as gp_cov
using ComposableTuringIDModels: hsgp_basis, se_spectral_density,
    _hsgp_standardised_index

n0, σ0, ℓ0, c0 = 40, 1.0, 1.0, 2.0
x = _hsgp_standardised_index(n0)
k = σ0^2 * with_lengthscale(SqExponentialKernel(), ℓ0)
K_ecosystem = gp_cov(GP(k)(x))          # AbstractGPs' own Gram matrix

Φ, sqrt_λ = hsgp_basis(n0, 40, c0)
sd = sqrt.(se_spectral_density(sqrt_λ, σ0, ℓ0))
K_hsgp = Φ * Diagonal(sd .^ 2) * Φ'

(exact = norm(kernelmatrix(k, x) - K_ecosystem) / norm(K_ecosystem),
    hsgp = round(norm(K_hsgp - K_ecosystem) / norm(K_ecosystem), digits = 4))
```

The exact GP matches AbstractGPs to machine precision, and the Hilbert-space
weights reproduce the same kernel: the spectral density the HSGP applies is the
Fourier transform of exactly this KernelFunctions kernel.

## Simulate from an exact GP

To use a GP as the reproduction number we hand it to a [`Renewal`](@ref)
infection model as its `rt` latent process. The generation interval is a
``\mathrm{Gamma}(6.5, 0.62)`` serial interval discretised by [`IDData`](@ref),
and reported cases are overdispersed counts via [`NegativeBinomialError`](@ref).

We simulate a ground truth by driving the renewal model with an [`ExactGP`](@ref)
whose hyperparameters are *fixed* to a known length scale, then drawing one
trajectory. Passing `missing` observations makes the composed model a prior
simulator; `fix` pins the GP hyperparameters and a single seeded draw gives the
observed data — no rejection loop. The returned quantities include the reported
cases `generated_y_t`, the latent infections `I_t`, and the GP path
`Z_t = \log R_t`.

```@example gp
using Turing: fix

data = IDData(gen_distribution = Gamma(6.5, 0.62))
obs = NegativeBinomialError(cluster_factor_prior = HalfNormal(0.1))
n = 70

truth = IDModel(Renewal(data; rt = ExactGP(),
        initialisation_prior = Normal(log(2.0), 0.1)), obs)
Random.seed!(10)
sim = fix(as_turing_model(truth, fill(missing, n), n), (ℓ = 0.55, σ = 0.55))()
y_obs = sim.generated_y_t
Z_true = sim.Z_t
(total_cases = sum(y_obs), peak = maximum(y_obs),
    Rt_range = round.(extrema(exp.(Z_true)), digits = 2))
```

## Fit both GPs and compare

Conditioning on the observed counts and sampling with NUTS recovers the
posterior. We differentiate with
[Mooncake](https://chalk-lab.github.io/Mooncake.jl/), the recommended backend for
this package. A small helper fits a chosen GP latent, times the run, recovers the
posterior-mean ``\log R_t`` with [`generated_observables`](@ref), and scores it
against the truth. `sample` returns a
[FlexiChains](https://github.com/penelopeysm/FlexiChains.jl) chain, read directly
— no conversion.

```@example gp
using Turing, Mooncake, Statistics
using ADTypes: AutoMooncake

function fit_gp(latent)
    model = IDModel(Renewal(data; rt = latent,
            initialisation_prior = Normal(log(2.0), 0.1)), obs)
    posterior = as_turing_model(model, y_obs, n)
    time = @elapsed chain = sample(posterior,
        NUTS(0.9; adtype = AutoMooncake(; config = nothing)), 300;
        progress = false)
    gen = vec(generated_observables(posterior, y_obs, chain).generated)
    Z_mean = mean([g.Z_t for g in gen])
    (; model, chain, Z_mean, time,
        cor = cor(Z_mean, Z_true),
        rmse = sqrt(mean((Z_mean .- Z_true) .^ 2)))
end

Random.seed!(1)
hs = fit_gp(HilbertSpaceGP(m = 20))
Random.seed!(1)
ex = fit_gp(ExactGP())

(hsgp = (cor = round(hs.cor, digits = 2), rmse = round(hs.rmse, digits = 3),
        seconds = round(hs.time, digits = 1)),
    exact = (cor = round(ex.cor, digits = 2), rmse = round(ex.rmse, digits = 3),
        seconds = round(ex.time, digits = 1)))
```

Both recover the latent reproduction number closely. The approximate model is the
faster of the two — the Hilbert-space basis is fixed, so each evaluation is a
matrix–vector product, whereas the exact GP rebuilds and factorises the full
covariance at ``O(n^3)`` on every evaluation. That gap is the reason the
approximation exists, and it widens with the series length; at short ``n`` like
this the exact GP is still affordable and gives the reference the approximation is
judged against.

## Posterior trajectories

Following the plotting convention of the other case studies, two small helpers
reduce the per-draw trajectories to credible bands.

```@setup gp
using Statistics

const CI_QS = [0.025, 0.25, 0.5, 0.75, 0.975]

function credible_bands(mat; qs = CI_QS)
    reduce(hcat, (map(eachrow(mat)) do row
        vals = collect(skipmissing(row))
        isempty(vals) ? missing : quantile(vals, q)
    end for q in qs))
end

function ci_ribbon!(ax, ts, bands; color, label)
    keep = findall(!ismissing, view(bands, :, 3))
    x, b = ts[keep], Float64.(bands[keep, :])
    band!(ax, x, b[:, 1], b[:, 5]; color = (color, 0.15))
    band!(ax, x, b[:, 2], b[:, 4]; color = (color, 0.3))
    lines!(ax, x, b[:, 3]; color = color, linewidth = 2, label = label)
end

function predictive_bands(pred, n)
    ndraws = length(vec(pred[@varname(y_t[n])]))
    rows = map(1:n) do i
        try
            permutedims(vec(pred[@varname(y_t[i])]))
        catch
            fill(missing, 1, ndraws)
        end
    end
    credible_bands(reduce(vcat, rows))
end
```

The reproduction number ``R_t = \exp(Z_t)`` comes from the returned `Z_t` draws;
the posterior-predictive case counts come from `predict` on the model with the
observations set to `missing`. We overlay both GP fits on the simulated truth.

```@example gp
using CairoMakie
using Turing: @varname

ts = 1:n
fig = Figure(; size = (760, 620))
ax1 = Axis(fig[1, 1]; ylabel = "Reproduction number R(t)")
for (fit, color, label) in
    ((hs, :teal, "HSGP"), (ex, :purple, "exact"))
    gen = vec(generated_observables(fit.model, y_obs, fit.chain).generated)
    Rt = credible_bands(reduce(hcat, (exp.(g.Z_t) for g in gen)))
    ci_ribbon!(ax1, ts, Rt; color = color, label = label)
end
lines!(ax1, ts, exp.(Z_true); color = :black, linewidth = 2,
    linestyle = :dash, label = "truth")
axislegend(ax1; position = :rt)

ax2 = Axis(fig[2, 1]; xlabel = "Day", ylabel = "Reported cases")
pred = predict(as_turing_model(hs.model, fill(missing, n), n), hs.chain)
ci_ribbon!(ax2, ts, predictive_bands(pred, n); color = :teal,
    label = "HSGP posterior predictive")
scatter!(ax2, ts, y_obs; color = :black, markersize = 7, label = "observed")
axislegend(ax2; position = :rt)
fig
```

Both posterior ``R_t`` bands bracket the simulated truth, and the
posterior-predictive counts cover the observed epidemic curve. Neither GP ever
had to know it was modelling a reproduction number, and the renewal and
observation models never had to know their latent process was a GP: the two sides
met only through the length-`n` latent contract. Swapping the GP for an
[`AR`](@ref) or a [`RandomWalk`](@ref) is a one-line change to the `rt` argument.

Use [`HilbertSpaceGP`](@ref) by default — it stays cheap under NUTS as the series
grows. Reach for [`ExactGP`](@ref) on short series when you want the exact GP as a
reference, or a check on the approximation.

!!! note "Illustrative run"
    This example uses a short sampler run and simulated data to stay fast to
    build. For a real analysis you would use more iterations, check convergence
    diagnostics, tune the number of basis functions `m` to the expected
    smoothness, and supply observed data.

## References

```@bibliography
Pages = ["gaussian-process.md"]
Canonical = false
```
