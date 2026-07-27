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

The kernels are
[KernelFunctions.jl](https://juliagaussianprocesses.github.io/KernelFunctions.jl/)
types — the ones
[AbstractGPs.jl](https://juliagaussianprocesses.github.io/AbstractGPs.jl/) builds
exact GPs from. A kernel enters [`ExactGP`](@ref) through its Gram matrix and
[`HilbertSpaceGP`](@ref) only through its spectral density, so either model takes
`SqExponentialKernel` (the default), `Matern32Kernel` or `Matern52Kernel`.

Both models are linear in a vector of standard-normal weights (`z` for the exact
GP, ``\beta`` for the Hilbert-space one), so fixing the hyperparameters and
feeding in unit vectors traces out that linear map; the Gram matrix of the
result is the model's implied prior covariance, which we can compare against the
covariance AbstractGPs builds from the same kernel.

```@example gp
using LinearAlgebra, Statistics
using AbstractGPs: GP, cov as gp_cov
using KernelFunctions: with_lengthscale
using Turing: fix

n0, ℓ0, σ0, m0 = 40, 1.0, 1.0, 40
x = (collect(1:n0) .- mean(1:n0)) ./ std(1:n0)   # models standardise the index
K_ref = gp_cov(GP(σ0^2 * with_lengthscale(SqExponentialKernel(), ℓ0))(x))

unit(k, j) = [i == j ? 1.0 : 0.0 for i in 1:k]
relerr(K) = norm(K - K_ref) / norm(K_ref)

function gram(cols)
    M = reduce(hcat, cols)
    return M * M'
end

K_exact = gram([fix(as_turing_model(ExactGP(), n0),
                    (ℓ = ℓ0, σ = σ0, z = unit(n0, j)))() for j in 1:n0])
K_hsgp = gram([fix(as_turing_model(HilbertSpaceGP(m = m0, c = 2.0), n0),
                   (ℓ = ℓ0, σ = σ0, β = unit(m0, j)))() for j in 1:m0])
(exact = round(relerr(K_exact), sigdigits = 2),
    hsgp = round(relerr(K_hsgp), sigdigits = 2))
```

The exact GP reproduces the ecosystem covariance to the jitter it adds for a
stable Cholesky factor, and the Hilbert-space basis reproduces it to a fraction
of a percent.

## Simulate from an exact GP

To use a GP as the reproduction number we hand it to a [`Renewal`](@ref)
infection model as its `rt` latent process. The generation interval is a
``\mathrm{Gamma}(6.5, 0.62)`` serial interval discretised by [`Renewal`](@ref),
and reported cases are overdispersed counts via [`NegativeBinomialError`](@ref).

We simulate a ground truth by driving the renewal model with an [`ExactGP`](@ref)
whose hyperparameters are *fixed* to a known length scale, then drawing one
trajectory. Passing `missing` observations makes the composed model a prior
simulator; `fix` pins the GP hyperparameters and a single seeded draw gives the
observed data — no rejection loop. The returned quantities include the reported
cases `generated_y_t`, the latent infections `I_t`, and the GP path
`Z_t = \log R_t`.

```@example gp
si = Gamma(6.5, 0.62)
obs = NegativeBinomialError(cluster_factor = HalfNormal(0.1))
n = 70

truth = IDModel(Renewal(; generation_time = si, rt = ExactGP(),
        initialisation = Normal(log(2.0), 0.1)), obs)
Random.seed!(10)
sim = fix(as_turing_model(truth, fill(missing, n), n), (ℓ = 0.55, σ = 0.55))()
y_obs = sim.generated_y_t
Z_true = sim.Z_t
(total_cases = sum(y_obs), peak = maximum(y_obs),
    Rt_range = round.(extrema(exp.(Z_true)), digits = 2))
```

The simulated ``R_t`` path starts above one, dips, peaks again around day 40 and
falls through one near day 50, so the epidemic grows, turns over and declines —
the shape the fits have to recover.

```@example gp
using CairoMakie

fig_sim = Figure(; size = (760, 480))
ax_rt = Axis(fig_sim[1, 1]; ylabel = "Reproduction number Rₜ")
lines!(ax_rt, 1:n, exp.(Z_true); color = :black, linewidth = 2,
    label = "simulated truth")
hlines!(ax_rt, [1.0]; color = :grey, linestyle = :dash)
axislegend(ax_rt; position = :rt)
ax_y = Axis(fig_sim[2, 1]; xlabel = "Day", ylabel = "Reported cases")
scatter!(ax_y, 1:n, y_obs; color = :black, markersize = 7,
    label = "simulated observations")
axislegend(ax_y; position = :lt)
fig_sim
```

## Fit both GPs and compare

Conditioning on the observed counts and sampling with NUTS recovers the
posterior. We differentiate with
[Mooncake](https://chalk-lab.github.io/Mooncake.jl/), the recommended backend for
this package.

The chain is started from a prior draw (`InitFromPrior()`). Turing's default
initialisation instead draws each parameter uniformly on ``[-2, 2]`` in
unconstrained space, which for a positive scale parameter means a starting value
up to ``e^2 \approx 7.4``: seven prior standard deviations out for the marginal
standard deviation ``\sigma``. A GP with ``\sigma \approx 7`` puts ``\log R_t``
in the tens, where the renewal likelihood is astronomically bad and its curvature
enormous, so NUTS shrinks the step size to ``10^{-8}`` and the chain never leaves
its starting point. Any GP hyperprior with a long tail has this failure mode, and
starting from the prior avoids it.

A small helper fits a chosen GP latent, times the run, recovers the per-draw
``\log R_t`` with [`generated_observables`](@ref), and scores the posterior mean
against the truth. `sample` returns a
[FlexiChains](https://github.com/penelopeysm/FlexiChains.jl) chain, read directly
— no conversion.

```@example gp
using Turing, Mooncake, Statistics
using ADTypes: AutoMooncake

function fit_gp(latent)
    model = IDModel(Renewal(; generation_time = si, rt = latent,
            initialisation = Normal(log(2.0), 0.1)), obs)
    posterior = as_turing_model(model, y_obs, n)
    time = @elapsed chain = sample(posterior,
        NUTS(0.9; adtype = AutoMooncake(; config = nothing)), 300;
        initial_params = InitFromPrior(), progress = false)
    gen = vec(generated_observables(posterior, y_obs, chain).generated)
    Z = reduce(hcat, (g.Z_t for g in gen))          # time × draw
    Z_mean = vec(mean(Z; dims = 2))
    (; model, chain, Z, Z_mean, time, cor = cor(Z_mean, Z_true),
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
reduce the per-draw trajectories to credible bands and draw a median line with
50% and 95% ribbons.

```@example gp
CI_QS = [0.025, 0.25, 0.5, 0.75, 0.975]

# time × 5 credible bands from a time × draws matrix
credible_bands(mat; qs = CI_QS) = reduce(hcat,
    (map(row -> quantile(row, q), eachrow(mat)) for q in qs))

# median line with 50% and 95% ribbons
function ci_ribbon!(ax, ts, b; color, label)
    band!(ax, ts, b[:, 1], b[:, 5]; color = (color, 0.15))
    band!(ax, ts, b[:, 2], b[:, 4]; color = (color, 0.3))
    lines!(ax, ts, b[:, 3]; color = color, linewidth = 2, label = label)
end
```

The reproduction number ``R_t = \exp(Z_t)`` comes from the returned `Z_t` draws.
The reported counts are scored element-wise, so their posterior *predictive*
distribution — fresh counts under each posterior parameter set — comes from
`predict` on the same model with the observations set to `missing`, stacking the
`y_t[i]` draws into a time × draws matrix. We overlay both GP fits on the
simulated truth.

```@example gp
ts = 1:n
fig = Figure(; size = (760, 620))
ax1 = Axis(fig[1, 1]; ylabel = "Reproduction number Rₜ")
for (fit, color, label) in ((hs, :teal, "HSGP"), (ex, :purple, "exact"))
    ci_ribbon!(ax1, ts, credible_bands(exp.(fit.Z)); color = color,
        label = label)
end
lines!(ax1, ts, exp.(Z_true); color = :black, linewidth = 2,
    linestyle = :dash, label = "truth")
hlines!(ax1, [1.0]; color = :grey, linestyle = :dash)
axislegend(ax1; position = :rt)

ax2 = Axis(fig[2, 1]; xlabel = "Day", ylabel = "Reported cases")
pred = predict(as_turing_model(hs.model, fill(missing, n), n), hs.chain)
yt = reduce(vcat, (permutedims(vec(pred[@varname(y_t[i])])) for i in ts))
ci_ribbon!(ax2, ts, credible_bands(yt); color = :teal,
    label = "HSGP posterior predictive")
scatter!(ax2, ts, y_obs; color = :black, markersize = 7, label = "observed")
axislegend(ax2; position = :lt)
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
