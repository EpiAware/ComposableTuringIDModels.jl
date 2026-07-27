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

Both models are linear in a vector of standard-normal weights (`z` for the exact
GP, ``\beta`` for the Hilbert-space one), so fixing the hyperparameters and
feeding in unit vectors traces out that linear map. The Gram matrix of the result
is the model's implied prior covariance, and we compare it against the covariance
AbstractGPs builds from the same kernel on the same inputs — both models hand the
kernel [`standardised_index`](@ref).

```@example gp
using LinearAlgebra
using AbstractGPs: GP, cov as gp_cov
using KernelFunctions: with_lengthscale
using Turing: fix

n0, σ0 = 40, 1.0
x = standardised_index(n0)
unit(k, j) = [i == j ? 1.0 : 0.0 for i in 1:k]
gram(cols) = (M = reduce(hcat, cols); M * M')

exact_cov(ℓ) = gram([fix(as_turing_model(ExactGP(), n0),
                        (ℓ = ℓ, σ = σ0, z = unit(n0, j)))() for j in 1:n0])
hsgp_cov(ℓ, m, c) = gram([fix(as_turing_model(HilbertSpaceGP(m = m, c = c), n0),
                             (ℓ = ℓ, σ = σ0, β = unit(m, j)))() for j in 1:m])

function relerr(K, ℓ)
    K_ref = gp_cov(GP(σ0^2 * with_lengthscale(SqExponentialKernel(), ℓ))(x))
    return norm(K - K_ref) / norm(K_ref)
end

# ℓ = 0.3 sits in the bulk of the default prior, ℓ = 1.0 far out in its tail
err = (exact = relerr(exact_cov(0.3), 0.3),
    default_short = relerr(hsgp_cov(0.3, 20, 1.5), 0.3),
    default_long = relerr(hsgp_cov(1.0, 20, 1.5), 1.0),
    more_basis_long = relerr(hsgp_cov(1.0, 40, 1.5), 1.0),
    wider_domain_long = relerr(hsgp_cov(1.0, 20, 2.0), 1.0))

# asserted, so drift fails the build rather than printing a worse number
@assert err.exact < 1e-5 && err.default_short < 1e-3
@assert err.default_long < 0.06 && err.wider_domain_long < 1e-3
map(e -> round(e, sigdigits = 2), err)
```

The exact GP matches to the jitter it adds for a stable Cholesky factor.
At its defaults the Hilbert-space basis matches to two parts in ten thousand at
the shorter length scale, the regime this page fits in.
At the longer one it is out by a few percent and more basis functions do not
help: there the boundary factor `c` sets the floor, because the approximation is
periodic on ``[-L, L]`` and a slowly varying path feels that boundary. Widening
to `c = 2` recovers two orders of magnitude.

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

The epidemic is seeded at around fifty infections. That matters for what the
page can claim: with a handful of cases a day the first fortnight of counts
carries almost no information about ``R_t``, both fits fall back on the prior
there, and a comparison against the truth would mostly be measuring the seed
size rather than either model.

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

The chain is started from a prior draw (`InitFromPrior()`). Turing's default
initialisation instead draws each parameter uniformly on ``[-2, 2]`` in
unconstrained space, which for a positive scale parameter means a starting value
up to ``e^2 \approx 7.4``. The prior on the marginal standard deviation
``\sigma`` is a half-normal with unit scale, which puts almost all of its mass
below 2, so that is far out in the tail. A GP with ``\sigma \approx 7`` puts
``\log R_t`` in the tens, where the renewal likelihood is astronomically bad and
its curvature enormous, so NUTS shrinks the step size by six or seven orders of
magnitude and the chain never leaves its starting point. Any GP hyperprior with a
long tail has this failure mode, and starting from the prior avoids it.

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

n_ex = 40   # the exact GP is fit to the first n_ex days only
Random.seed!(1)
hs = fit_gp(HilbertSpaceGP(m = 20), n)
Random.seed!(1)
ex = fit_gp(ExactGP(), n_ex)

# the page claims recovery, so it asserts recovery: a fit that had collapsed
# onto the prior would fail the build rather than render a wrong figure
@assert hs.cor > 0.9 && ex.cor > 0.9
@assert hs.rmse < 0.15 && ex.rmse < 0.15
score(f) = (days = f.n, cor = round(f.cor, digits = 2),
    rmse = round(f.rmse, digits = 3), seconds = round(f.time, digits = 1))
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
@assert grad.exact > grad.hsgp
round(grad.exact / grad.hsgp; sigdigits = 2)
```

The exact GP costs more per gradient even on the shorter series, and the gap
widens with the series length: that is the reason the approximation exists. At
short ``n`` the exact GP is still affordable and gives the reference the
approximation is judged against.

That leaves the basis rebuild the composed model does on every evaluation. It is
not differentiated with respect to any sampled parameter, but it does run inside
the traced call, so it is taped like the rest of the body. The honest measurement
is therefore a gradient with the basis rebuilt inside against a gradient with it
captured outside — not a bare call to the constructor:

```@example gp
using DynamicPPL: to_submodel

# the same GP, but with the basis rebuilt inside the traced body, exactly as the
# enclosing Renewal does when it reconstructs its submodels
@model function rebuilt_basis(gp, n)
    z ~ to_submodel(as_turing_model(gp, n), false)
    return z
end

gp = HilbertSpaceGP(m = 20)
overhead = gradient_time(rebuilt_basis(gp, n)) -
           gradient_time(as_turing_model(gp, n))
share = overhead / gradient_time(composed(gp, n))
round(share; sigdigits = 2)
```

That is the fraction of each gradient a cached basis would recover. Caching
would tie a latent model to one series length, which is the composability the
rest of the package is built on, so the repeat stands; the number above is what
it costs, measured rather than assumed.

## Posterior trajectories

Following the plotting convention of the other case studies, two small helpers
reduce the per-draw trajectories to credible bands and draw a median line with
50% and 95% ribbons. Neither skips a missing or failed value: nothing here is
scored on a delayed scale, so a band shorter than the days it is drawn against
is a bug and raises rather than plotting a truncated ribbon.

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

The share of days on which the 95% band contains the quantity it is estimating
is a calibration check on the bands the figure draws. It is not on its own a
check on fit quality — a fit that had collapsed onto a wide, near-prior band
would score well on it — which is why the accuracy of the posterior means is
asserted separately above.

```@example gp
covered(b, x) = round(mean(b[:, 1] .<= x .<= b[:, 5]); digits = 2)
coverage = (Rt_hsgp = covered(credible_bands(exp.(hs.Z)), exp.(Z_true)),
    Rt_exact = covered(credible_bands(exp.(ex.Z)), exp.(Z_true[1:(ex.n)])),
    cases = covered(credible_bands(yt), y_obs))
@assert all(>=(0.8), coverage) "a 95% band covers under 80% of days: $coverage"
coverage
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
