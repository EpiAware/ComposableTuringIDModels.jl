# [A Gaussian-process latent process](@id case-study-gp)

One of the design claims of the package is that *any* latent process — anything
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
    entirely outside the differentiated path. [`as_turing_model`](@ref) therefore
    builds the basis **once** when the model is constructed and captures it, rather
    than rebuilding it on every gradient evaluation; each log-density evaluation is
    then just a single ``n \times m`` matrix–vector product. Caching the basis this
    way cuts the per-gradient Mooncake cost by an order of magnitude relative to
    rebuilding it inside the model.
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
using ComposableTuringIDModels, Distributions, Random
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
functions grows. The kernel itself is the ecosystem-standard
[KernelFunctions.jl](https://juliagaussianprocesses.github.io/KernelFunctions.jl/)
`SqExponentialKernel`, so we can check the approximation directly against that
kernel's own Gram matrix. The package exposes the basis builder
[`ComposableTuringIDModels.hsgp_basis`](@ref) and spectral density
[`ComposableTuringIDModels.se_spectral_density`](@ref) used internally; the inputs
are standardised, so we build the target Gram matrix on the same standardised
coordinates:

```@example gp
using LinearAlgebra
using KernelFunctions: with_lengthscale, kernelmatrix
n = 40
σ, ℓ, c = 1.0, 1.0, 2.0
x = ComposableTuringIDModels._hsgp_standardised_index(n)
K_exact = kernelmatrix(σ^2 * with_lengthscale(SqExponentialKernel(), ℓ), x)

Φ, sqrt_λ = ComposableTuringIDModels.hsgp_basis(n, 40, c)
sd = sqrt.(ComposableTuringIDModels.se_spectral_density(sqrt_λ, σ, ℓ))
K_approx = Φ * Diagonal(sd .^ 2) * Φ'

round(norm(K_approx - K_exact) / norm(K_exact), digits = 4)
```

The relative error is a fraction of a percent: with enough basis functions the
Hilbert-space weights reproduce the KernelFunctions kernel they stand in for. This
is the concrete link to the GP ecosystem — the spectral density the HSGP applies
is the Fourier transform of exactly this KernelFunctions kernel.

## Choosing a kernel

The kernels are the standard [KernelFunctions.jl](https://juliagaussianprocesses.github.io/KernelFunctions.jl/)
types, so the model reuses the ecosystem's kernels rather than defining its own.
Only the spectral density changes between kernels, so the kernel is a one-field
choice on the model. The default `SqExponentialKernel` gives very smooth paths;
`Matern32Kernel` and `Matern52Kernel` give progressively rougher ones, which can
suit a less smooth latent process. These three match the kernels offered by the
EpiNow2 Gaussian-process implementation. The basis is shared, so swapping the
kernel reuses everything else:

```@example gp
using KernelFunctions: SqExponentialKernel, Matern52Kernel
gp_se = HilbertSpaceGP(m = 20, kernel = SqExponentialKernel())
gp_matern = HilbertSpaceGP(m = 20, kernel = Matern52Kernel())
(se = length(as_turing_model(gp_se, 60)()),
    matern = length(as_turing_model(gp_matern, 60)()))
```

A kernel enters the HSGP only through its spectral density, so adding a new kernel
is a single `ComposableTuringIDModels.spectral_density(::MyKernel, ω, σ, ℓ)`
method — any KernelFunctions `Kernel` for which that method exists can drive the
GP. This case study uses the squared-exponential kernel; the renewal example below
is identical for a Matérn kernel bar that one argument.

### Where this sits in the Julia GP ecosystem

The HSGP is deliberately a *NUTS-friendly approximation*, not a hand-rolled GP in
isolation. It connects to the wider ecosystem at two points, and has two natural
siblings:

  - **[KernelFunctions.jl](https://juliagaussianprocesses.github.io/KernelFunctions.jl/)**
    supplies the covariance kernels used here. KernelFunctions defines the kernels
    and their Gram matrices (which we validated against above) but not their
    spectral densities; the HSGP adds the one-dimensional spectral density each
    kernel needs for the basis-function approximation.
  - **[AbstractGPs.jl](https://juliagaussianprocesses.github.io/AbstractGPs.jl/)**
    is the ecosystem's interface for *exact* GPs built from those same kernels. It
    gives an exact posterior in closed form for Gaussian likelihoods, but the
    ``O(n^3)`` factorisation it performs at every evaluation is what makes an exact
    GP impractical inside NUTS at realistic ``n`` — the obstacle the HSGP exists to
    sidestep. It is the right tool when an exact GP is affordable; the HSGP is the
    approximation that stays cheap under gradient-based sampling.
  - **[TemporalGPs.jl](https://github.com/JuliaGaussianProcesses/TemporalGPs.jl)**
    is the closest sibling: for a one-dimensional input (like time here) it
    reformulates an AbstractGPs GP as a linear-Gaussian state-space model and gives
    an *exact* GP at ``O(n)`` cost via Kalman filtering. That makes it an appealing
    alternative latent process for ``\log R_t`` — exact rather than approximate,
    and still linear in the series length — though the HSGP's non-centred,
    fixed-basis form is the more direct fit for a Turing model differentiated under
    Mooncake. Both would plug into the same length-`n` latent contract.

## Composing it into an infection model

To use the GP as the time-varying reproduction number we hand it to a
[`Renewal`](@ref) infection model as its `rt` latent process. Nothing about the
renewal process changes — it asks its latent slot for a length-`n` ``\log R_t``
path and gets one. The generation interval is a ``\mathrm{Gamma}(6.5, 0.62)``
serial interval discretised by [`IDData`](@ref), and reported cases are
overdispersed counts via [`NegativeBinomialError`](@ref).

```@example gp
data = IDData(gen_distribution = Gamma(6.5, 0.62))
renewal = Renewal(data; rt = gp, initialisation_prior = Normal(log(2.0), 0.1))
obs = NegativeBinomialError(cluster_factor_prior = HalfNormal(0.1))
model = IDModel(renewal, obs)
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
the observation cluster factor. `sample` returns a
[FlexiChains](https://github.com/penelopeysm/FlexiChains.jl) chain, so we read the
posterior from it directly — no conversion to another chain type. `summarystats`
works natively on the chain and gives the usual per-parameter summary — point
estimates and their uncertainty, alongside the effective sample size and
``\hat{R}`` convergence diagnostic:

```@example gp
using MCMCChains, Statistics
summarystats(chain)
```

Individual parameters are read by name straight off the FlexiChains chain with
`vec(chain[@varname(...)])` — no conversion step — from which posterior summaries
of the GP hyperparameters follow directly:

```@example gp
using Turing: @varname
ℓ_post = vec(chain[@varname(ℓ)])
σ_post = vec(chain[@varname(σ)])
(ℓ = round(mean(ℓ_post), digits = 2), σ = round(mean(σ_post), digits = 2))
```

The length scale ``\ell``, marginal standard deviation ``\sigma``, and the
observation cluster factor are all identified from the single simulated series.

## Recover the latent process

The reproduction number ``R_t = \exp(Z_t)`` is a generated quantity rather than a
sampled parameter, so we do not rebuild it by hand from the basis weights. The
composed model already *returns* `Z_t` (alongside `I_t` and `generated_y_t`), and
[`generated_observables`](@ref) re-runs the fitted model over the chain and hands
back those returned quantities for every posterior draw. It is the counterpart of
`predict`: `predict` samples the *observed* variables (the reported counts `y_t`)
forward, whereas `generated_observables` collects the model's *returned* latent
quantities — which is what we want for ``Z_t``. Both push the fitted model forward
over the posterior; neither reaches into the basis. This is the basis-function
form paying off — the posterior over the latent function is just the model pushed
forward over the posterior over its weights.

```@example gp
gen = vec(generated_observables(posterior, y_obs, chain).generated)
Z_draws = [g.Z_t for g in gen]
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

## Posterior trajectories

Point summaries only go so far. Following the same plotting convention as the
other case studies, two small helpers reduce the per-draw trajectories to
credible bands, which we overlay on the simulated truth.

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
observations set to `missing`.

```@example gp
using CairoMakie

ts = 1:n
Rt = credible_bands(reduce(hcat, (exp.(z) for z in Z_draws)))

pred = predict(as_turing_model(model, fill(missing, n), n), chain)
yt = predictive_bands(pred, n)

fig = Figure(; size = (760, 620))
ax1 = Axis(fig[1, 1]; ylabel = "Reproduction number R(t)")
ci_ribbon!(ax1, ts, Rt; color = :purple, label = "posterior")
lines!(ax1, ts, exp.(Z_true); color = :black, linewidth = 2,
    linestyle = :dash, label = "truth")
axislegend(ax1; position = :rt)
ax2 = Axis(fig[2, 1]; xlabel = "Day", ylabel = "Reported cases")
ci_ribbon!(ax2, ts, yt; color = :teal, label = "posterior predictive")
scatter!(ax2, ts, y_obs; color = :black, markersize = 7, label = "observed")
axislegend(ax2; position = :rt)
fig
```

The posterior credible band for ``R_t`` brackets the simulated truth, and the
posterior-predictive case counts cover the observed epidemic curve — the GP has
recovered the latent reproduction number through the renewal and observation
models without ever seeing them directly.

The GP never had to know it was modelling a reproduction number, and the renewal
and observation models never had to know their latent process was a GP. The two
sides met only through the length-`n` latent contract — which is exactly the
composability the package is built around. Swapping the GP for an
[`AR`](@ref) or a [`RandomWalk`](@ref) latent is a one-line change to the
`rt` argument; the rest of the model is untouched.

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
