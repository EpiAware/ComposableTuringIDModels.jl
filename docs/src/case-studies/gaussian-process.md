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

They differ in how much of their work depends on the sampled parameters. The
Hilbert-space basis is set by the series length, the number of basis functions
and the boundary factor alone, so [`as_turing_model`](@ref) builds it before
returning the model and captures it: the `@model` body never sees it, and one
``n \times m`` matrix–vector product is all that is differentiated. The exact
GP's covariance and its Cholesky factor depend on ``\ell`` and ``\sigma``, so
both are rebuilt *and* differentiated on every evaluation. That is where the
timing gap below comes from.

Composed into a [`Renewal`](@ref), the enclosing `@model` reconstructs its
submodels on every evaluation, so the basis is rebuilt inside the traced call.
It still depends on no sampled parameter, but it is executed and taped along
with everything else, so it is not free. That cost is measured below.

## The kernel and the GP ecosystem

The kernels are
[KernelFunctions.jl](https://juliagaussianprocesses.github.io/KernelFunctions.jl/)
types — the ones
[AbstractGPs.jl](https://juliagaussianprocesses.github.io/AbstractGPs.jl/) builds
exact GPs from. A kernel enters [`ExactGP`](@ref) through its Gram matrix and
[`HilbertSpaceGP`](@ref) only through its spectral density, so either model takes
`SqExponentialKernel` (the default), `Matern32Kernel` or `Matern52Kernel`.
[`ExactGP`](@ref) forms that Gram matrix itself rather than wrapping an
AbstractGPs `GP`, because it needs the factor inside a Turing model where the
hyperparameters are sampled; the kernels are shared either way.

Both models are linear in a vector of standard-normal weights, so unit vectors at
fixed hyperparameters trace out that map and the Gram matrix of the result is the
model's implied prior covariance. Against AbstractGPs' own:

```@example gp
using LinearAlgebra
using AbstractGPs: GP, cov as gp_cov
using KernelFunctions: with_lengthscale
using Turing: fix

n0, σ0, ℓ0 = 40, 1.0, 0.3          # ℓ0 is the median of the default prior
unit(k, j) = [i == j ? 1.0 : 0.0 for i in 1:k]
gram(cols) = (M = reduce(hcat, cols); M * M')
relerr(K) = norm(K - K_ref) / norm(K_ref)

K_ref = gp_cov(GP(σ0^2 * with_lengthscale(SqExponentialKernel(), ℓ0))(
    standardised_index(n0)))
K_exact = gram([fix(as_turing_model(ExactGP(), n0),
                   (ℓ = ℓ0, σ = σ0, z = unit(n0, j)))() for j in 1:n0])
K_hsgp = gram([fix(as_turing_model(HilbertSpaceGP(m = 20), n0),
                  (ℓ = ℓ0, σ = σ0, β = unit(20, j)))() for j in 1:20])

err = (exact = relerr(K_exact), hsgp = relerr(K_hsgp))
@assert err.exact < 1e-5 && err.hsgp < 1e-3   # drift fails the build
map(e -> round(e, sigdigits = 2), err)
```

The exact GP matches to the jitter it adds for a stable Cholesky factor, and the
basis to two parts in ten thousand. That accuracy is not uniform in ``\ell``: the
basis cannot resolve wiggles finer than ``2L/m``, so a short length scale needs a
larger `m`, and a long one feels the ``[-L, L]`` boundary and needs a wider `c`.

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

The epidemic is seeded at around fifty infections rather than a handful, which
matters for what the page can claim. With only a few cases a day the first
fortnight of counts would carry almost no information about ``R_t``, both fits
would fall back on the prior there, and a comparison against the truth would
mostly measure the seed size rather than either model.

```@example gp
si = Gamma(6.5, 0.62)
obs = NegativeBinomialError(cluster_factor = HalfNormal(0.1))
n = 70

id_model(latent) = IDModel(Renewal(; generation_time = si, rt = latent,
        initialisation = Normal(log(50), 0.1)), obs)

Random.seed!(10)
sim = fix(as_turing_model(id_model(ExactGP()), fill(missing, n), n),
    (ℓ = 0.55, σ = 0.55))()
y_obs = sim.generated_y_t
Z_true = sim.Z_t
(total_cases = sum(y_obs), peak = maximum(y_obs),
    Rt_range = round.(extrema(exp.(Z_true)), digits = 2))
```

The simulated ``R_t`` path starts around two, dips, peaks again around day 40 and
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

A small helper fits a chosen GP latent to the first `n_fit` days, times the run,
recovers the per-draw ``\log R_t`` with [`generated_observables`](@ref), and
scores the posterior mean against the truth. `sample` returns a
[FlexiChains](https://github.com/penelopeysm/FlexiChains.jl) chain, read directly
— no conversion.

```@example gp
using Turing, Mooncake, Statistics
using ADTypes: AutoMooncake

function fit_gp(latent, n_fit)
    model = id_model(latent)
    y = y_obs[1:n_fit]
    posterior = as_turing_model(model, y, n_fit)
    time = @elapsed chain = sample(posterior,
        NUTS(0.9; adtype = AutoMooncake(; config = nothing)), 300;
        initial_params = InitFromPrior(), progress = false)
    gen = vec(generated_observables(posterior, y, chain).generated)
    Z = reduce(hcat, (g.Z_t for g in gen))          # time × draw
    Z_mean = vec(mean(Z; dims = 2))
    Z_ref = Z_true[1:n_fit]
    (; model, chain, n = n_fit, Z, Z_mean, time, cor = cor(Z_mean, Z_ref),
        rmse = sqrt(mean((Z_mean .- Z_ref) .^ 2)))
end

Random.seed!(1)
n_ex = 60
hs = fit_gp(HilbertSpaceGP(m = 20), n)
Random.seed!(1)
ex = fit_gp(ExactGP(), n_ex)

score(f) = (days = f.n, cor = round(f.cor, digits = 2),
    rmse = round(f.rmse, digits = 3), seconds = round(f.time, digits = 1))
scores() = "$(score(hs)), $(score(ex))"

@assert hs.cor>0.9 && ex.cor>0.85 "posterior mean lost the latent: $(scores())"
@assert hs.rmse<0.15 && ex.rmse<0.2 "posterior mean off the truth: $(scores())"
(hsgp = score(hs), exact = score(ex))
```

Both posterior means correlate with the simulated ``\log R_t`` above 0.9 on the
days they cover, and the trajectory figure below shows the two agreeing with
each other over the days they share — the check the exact GP is here to provide.

The exact GP is fit to the first `n_ex` days rather than all `n`: its ``O(n^3)``
factorisation is rebuilt and differentiated at every leapfrog step, and at the
full length it would dominate the cost of building this page. The wall-clock
times above are single un-warmed runs that include the model's compilation and
cover different series lengths, so read them as indicative. Timing one gradient
of each composed log-density on the *same* `n_ex` days is the like-for-like
comparison, and it is what the sampler actually pays per leapfrog step:

```@example gp
using Chairmarks: @b
using DynamicPPL: LogDensityFunction, VarInfo
using LogDensityProblems: logdensity_and_gradient

composed(latent, n_fit) = as_turing_model(
    id_model(latent), y_obs[1:n_fit], n_fit)

function gradient_time(model)
    ldf = LogDensityFunction(model; adtype = AutoMooncake(; config = nothing))
    θ = VarInfo(model)[:]
    return (@b logdensity_and_gradient(ldf, θ)).time
end

grad = (hsgp = gradient_time(composed(HilbertSpaceGP(m = 20), n_ex)),
    exact = gradient_time(composed(ExactGP(), n_ex)))
round(grad.exact / grad.hsgp; sigdigits = 2)
```

The exact GP costs more per gradient even on the shorter series, and the gap
widens with the series length: that is the reason the approximation exists.

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
for (fit, color, label) in ((hs, :teal, "HSGP"),
    (ex, :purple, "exact (first $(n_ex) days)"))
    ci_ribbon!(ax1, 1:(fit.n), credible_bands(exp.(fit.Z)); color = color,
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
