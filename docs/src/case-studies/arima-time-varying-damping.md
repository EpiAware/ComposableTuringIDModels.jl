# [Time-varying damping in an AR process](@id case-study-arima-tvdamp)

Every parameter slot of a component is an [`AbstractPriorModel`](@ref): as well as
a bare `Distribution`, it accepts a *latent process*, so a coefficient can carry a
structured prior instead of a single global draw.
This page uses that on the damping of an autoregressive process, the workhorse
inside [`AR`](@ref), [`arma`](@ref) and [`arima`](@ref).
It does two things: it checks that a latent process dropped straight into the
[`AR`](@ref) damping slot threads correctly under a gradient-based sampler, and it
builds a genuinely *time-varying* damping ``\rho_t`` from a latent process,
simulates from it, fits it under NUTS and shows the ``\rho_t`` trajectory is
recovered.

This is early-stage, actively developed software; the API may change.

## A damping coefficient with a latent-process prior

Passing a latent model to the damping slot — `AR(damp = RandomWalk())` — makes the
damping coefficient carry a random-walk prior rather than a single distribution.
The same works through [`arima`](@ref).
The latent model is auto-prefixed inside the [`as_prior`](@ref) coercion seam, so
its internal variables (`std`, `ϵ_t`, `rw_init`) cannot collide with the AR
innovation's own latent under the prefix-off submodel convention.
Without that prefixing the bare form sampled via `rand` but errored when evaluated
as a *linked* log-density (the target a gradient sampler differentiates).

```@example tvdamp
using ComposableTuringIDModels, Distributions, Turing, Random, Statistics
using Turing: to_submodel, returned
using DynamicPPL: LogDensityFunction, VarInfo, link, getlogjoint
import LogDensityProblems as LDP
Random.seed!(80)

for (name, mdl) in ("AR" => AR(; damp = RandomWalk()),
        "arima" => arima(; damp = RandomWalk()))
    m = as_turing_model(mdl, 40)
    vi = link(VarInfo(m), m)
    ldf = LogDensityFunction(m, getlogjoint, vi)
    finite = isfinite(LDP.logdensity(ldf, zeros(LDP.dimension(ldf))))
    println("$name: rand ok = ", rand(m) !== nothing,
        ", linked log-density finite = ", finite)
end
```

Both models draw a value *and* evaluate as a linked log-density, so they sample
under NUTS end-to-end:

```@example tvdamp
m = as_turing_model(AR(; damp = RandomWalk()), 40)
chain = sample(m, NUTS(0.8; adtype = Turing.AutoForwardDiff()), 300;
    progress = false)
size(chain, 1)
```

With an order-one [`AR`](@ref) the damping is a single coefficient drawn from that
process.
The rest of this page makes the damping vary over the series.

## A genuinely time-varying damping

An AR(1) recursion with a *time-varying* damping is

```math
z_t = \rho_t\, z_{t-1} + \sigma\, \epsilon_t, \qquad \epsilon_t \sim
\mathrm{Normal}(0, 1),
```

where ``\rho_t`` is itself a latent path.
We take that path from a package latent model — a [`RandomWalk`](@ref) — mapped
through `tanh` into the stationary band ``(-1, 1)``.
A single series barely identifies a whole ``\rho_t`` path (each ``\rho_t`` sees one
transition), so we simulate a small panel of independent series that share the
same damping path, which sharpens the recovery — the standard identification
route for a time-varying coefficient.

The true damping ramps from strong positive persistence through zero to mild
anti-persistence:

```@example tvdamp
n = 30
tgrid = range(0, 1; length = n)
ρ_true = 0.85 .- 1.15 .* tgrid          # 0.85 → −0.30, crossing zero
σ_true = 0.3
K = 20                                  # independent series sharing ρ_t
Y = zeros(n, K)
for k in 1:K
    Y[1, k] = randn()
    for t in 2:n
        Y[t, k] = ρ_true[t] * Y[t - 1, k] + σ_true * randn()
    end
end
(n = n, series = K, ρ_start = round(ρ_true[1], digits = 2),
    ρ_end = round(ρ_true[end], digits = 2))
```

The model draws the damping path from the latent process, squashes it into
``(-1, 1)``, and applies the AR(1) recursion to every series.
`n_groups`-style dimensions — here the series count `K` and the length `n` — are
read from the data matrix, not stored on the model:

```@example tvdamp
@model function time_varying_ar(Y, damping_process)
    n, K = size(Y)
    ρ_latent ~ to_submodel(as_turing_model(damping_process, n), false)
    ρ = tanh.(ρ_latent)                 # stationary band (−1, 1)
    σ ~ truncated(Normal(0, 0.5), 0, Inf)
    for k in 1:K, t in 2:n
        Y[t, k] ~ Normal(ρ[t] * Y[t - 1, k], σ)
    end
    return (; ρ, σ)
end
nothing # hide
```

The damping process is supplied as a component, so swapping the prior over
``\rho_t`` (a smoother walk, a different innovation scale, an [`AR`](@ref) instead
of a [`RandomWalk`](@ref)) is a one-line change.
We fit with a [`RandomWalk`](@ref) damping path:

```@example tvdamp
model = time_varying_ar(Y, RandomWalk())
fit = sample(model, NUTS(0.8; adtype = Turing.AutoForwardDiff()), 300;
    progress = false)
size(fit, 1)
```

The ``\rho_t`` path is a generated quantity of the model, recovered per draw with
`returned`:

```@example tvdamp
ρ_draws = reduce(hcat, [g.ρ for g in vec(returned(model, fit))])   # n × draws
ρ_mean = vec(mean(ρ_draws; dims = 2))
σ_post = mean(fit[:σ])
(cor_with_truth = round(cor(ρ_mean, ρ_true), digits = 3),
    σ_posterior_mean = round(σ_post, digits = 3), σ_true = σ_true)
```

The posterior mean damping path tracks the true ramp closely and the innovation
scale ``\sigma`` is recovered.
Plotting the path against the truth with an 80% credible band makes the recovery
visible:

```@example tvdamp
using CairoMakie
qs = [quantile(ρ_draws[t, :], [0.1, 0.5, 0.9]) for t in 1:n]
lo = getindex.(qs, 1)
md = getindex.(qs, 2)
hi = getindex.(qs, 3)

fig = Figure(; size = (760, 420))
ax = Axis(fig[1, 1]; xlabel = "Time step", ylabel = "Damping ρₜ")
band!(ax, 1:n, lo, hi; color = (:seagreen, 0.25))
lines!(ax, 1:n, md; color = :seagreen, linewidth = 2, label = "posterior mean")
lines!(ax, 1:n, ρ_true; color = :black, linestyle = :dash, linewidth = 2,
    label = "true ρₜ")
hlines!(ax, [0.0]; color = :grey, linestyle = :dot)
axislegend(ax; position = :lb)
fig
```

The posterior band covers the true trajectory across the whole series, including
the sign change, so the time-varying damping is recovered from data.
The same construction takes any latent process in the damping slot: the AR
recursion is fixed, the prior over how the damping evolves is a component you
choose.
