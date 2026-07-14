# [Time-varying damping in an AR process](@id case-study-tvdamp)

Every parameter slot of a component is an [`AbstractPriorModel`](@ref): as well as
a bare `Distribution`, a slot accepts a *prior model*, including a latent process.
On the damping of an autoregressive process this buys two things, both through the
struct's own constructor with no hand-written recursion:

  - a **structured prior** over the (time-constant) damping coefficient, via
    [`AR`](@ref); and
  - a **genuinely time-varying** coefficient path ``\rho_t``, via
    [`TimeVaryingAR`](@ref).

This is early-stage, actively developed software; the API may change.

## A structured prior on the constant damping

[`AR`](@ref) applies its damping as a constant length-`p` coefficient, so a prior
model in its `damp` slot enriches the *prior* over that constant rather than making
it vary. A [`HierarchicalNormal`](@ref) there gives the coefficient an adaptive
(inferred-scale) prior:

```@example tvdamp
using ComposableTuringIDModels, Distributions, Turing, Random, Statistics
using Turing: to_submodel
Random.seed!(80)

ar = AR(; damp = HierarchicalNormal())
(order = ar.p, damping_is_constant = true)
```

The coefficient is one value shared across the series. To let it *change* over the
series, reach for [`TimeVaryingAR`](@ref).

## A genuinely time-varying damping

[`TimeVaryingAR`](@ref) is a first-order AR whose coefficient is a path:

```math
z_t = \rho_t\, z_{t-1} + \epsilon_t, \qquad \rho_t = \tanh(u_t),
```

where the unconstrained path ``u_t`` is drawn from a latent process (a
[`RandomWalk`](@ref) by default) and `tanh` maps it into the stationary band
``(-1, 1)``. The coefficient path is a component you choose: swapping the prior
over how the damping evolves is a one-line change to the `damp` argument, and the
AR recursion is untouched.

Built with `as_turing_model(m, n)` it returns the numeric length-`n` path (like
every other latent model), and tracks the coefficient path ``\rho_t`` as a
generated quantity `ρ` (recovered from the chain, below):

```@example tvdamp
tv = TimeVaryingAR()
length(as_turing_model(tv, 8)())
```

Because it returns a plain path it drops straight into any latent slot — here as
the ``r_t`` process of a [`Renewal`](@ref) inside a composed [`IDModel`](@ref),
with no glue code:

```@example tvdamp
data = IDData([0.2, 0.3, 0.5], exp)
nested = IDModel(
    Renewal(data; rt = TimeVaryingAR(), initialisation = Normal()),
    PoissonError())
length(as_turing_model(nested, missing, 12)().generated_y_t)
```

We simulate one series whose true damping ramps from strong positive persistence
through zero to mild anti-persistence, then recover it. A single series informs
each ``\rho_t`` through one transition, so recovery leans on the smoothness of the
[`RandomWalk`](@ref) damping prior (a panel of series sharing one ``\rho_t`` draw
would sharpen it further):

```@example tvdamp
n = 50
tgrid = range(0, 1; length = n - 1)
ρ_true = 0.85 .- 1.15 .* tgrid          # 0.85 → −0.30, crossing zero
z = zeros(n)
z[1] = randn()
for t in 2:n
    z[t] = ρ_true[t - 1] * z[t - 1] + 0.3 * randn()
end
(n = n, ρ_start = round(ρ_true[1], digits = 2),
    ρ_end = round(ρ_true[end], digits = 2))
```

The model wraps [`TimeVaryingAR`](@ref) in a thin observation of the path and fits
under NUTS. The coefficient path is tracked as the generated quantity `ρ`, so it is
recovered straight from the chain:

```@example tvdamp
@model function observe_path(y, n)
    latent ~ to_submodel(as_turing_model(TimeVaryingAR(), n), false)
    for t in 1:n
        y[t] ~ Normal(latent[t], 0.01)
    end
end

model = observe_path(z, n)
fit = sample(model, NUTS(0.8; adtype = Turing.AutoForwardDiff()), 100;
    progress = false)
# ρ is tracked as a generated quantity: `fit[:ρ]` is a per-draw coefficient path
ρ_draws = reduce(hcat, vec(fit[:ρ]))     # (n-1) × draws
ρ_mean = vec(mean(ρ_draws; dims = 2))
(correlation_with_truth = round(cor(ρ_mean, ρ_true), digits = 2),)
```

The posterior mean damping path tracks the true ramp, including the sign change.
Plotting it against the truth with an 80% credible band makes the recovery
visible:

```@example tvdamp
using CairoMakie
qs = [quantile(ρ_draws[t, :], [0.1, 0.5, 0.9]) for t in 1:(n - 1)]
lo = getindex.(qs, 1)
md = getindex.(qs, 2)
hi = getindex.(qs, 3)

fig = Figure(; size = (760, 420))
ax = Axis(fig[1, 1]; xlabel = "Time step", ylabel = "Damping ρₜ")
band!(ax, 1:(n - 1), lo, hi; color = (:seagreen, 0.25))
lines!(ax, 1:(n - 1), md; color = :seagreen, linewidth = 2,
    label = "posterior mean")
lines!(ax, 1:(n - 1), ρ_true; color = :black, linestyle = :dash, linewidth = 2,
    label = "true ρₜ")
hlines!(ax, [0.0]; color = :grey, linestyle = :dot)
axislegend(ax; position = :lb)
fig
```

The band covers the true trajectory across the series, so the time-varying damping
is recovered from data — through the [`TimeVaryingAR`](@ref) struct, with the
coefficient's evolution supplied as a component rather than coded by hand. Only the
order-1 case is built so far; time-varying coefficients for higher-order AR(`p`)
are tracked in [#113](https://github.com/EpiAware/ComposableTuringIDModels.jl/issues/113).
