# [Partial pooling across groups in a composed model](@id case-study-hierarchy)

A multi-group epidemic is a panel: one shared infection process seen by several
groups, each reporting it at its own level.
[`GroupedIDModel`](@ref) expresses this panel as a composed model rather than a
hand-written `@model`.
It holds a shared infection process, a `group_effect` prior over the grouping
axis, an observation model, and a combiner mapping a group's effect to its
expected series.

This page drives a per-group reporting level with a [`Hierarchy`](@ref) inside a
[`GroupedIDModel`](@ref), simulates from it and fits it end-to-end under NUTS,
recovering the per-group levels.
The group dimension threads from the data and the group prior is namespaced by the
component, so the panel composes with no hand-orchestration.

## The composed model

The shared epidemic is a [`DirectInfections`](@ref) process carrying a
[`RandomWalk`](@ref) latent, observed with a [`PoissonError`](@ref).
Every group sees the *same* infection curve ``I_t`` but reports it at its own
level: a per-group log-reporting level ``\ell_g`` scales the expected counts
before the observation model.
Those per-group levels are partially pooled with a [`Hierarchy`](@ref), supplied
as the `group_effect`.

```@example hier
using ComposableTuringIDModels, Distributions, Turing, Random, Statistics
using Turing: returned
Random.seed!(77)

idmodel = IDModel(
    DirectInfections(; Z = RandomWalk(), initialisation = Normal(log(50.0), 0.2)),
    PoissonError())
hierarchy = Hierarchy(; mean = Normal(0.0, 0.5), across = IID(Normal(0.0, 0.5)))
model = GroupedIDModel(idmodel, hierarchy)
```

The grouping dimension is **not** a field of any component.
[`GroupedIDModel`](@ref) reads `n_time` and `n_groups` from the shape of the data
matrix `Y` (rows are time, columns are groups) and passes `n_groups` to the
[`Hierarchy`](@ref) and `n_time` to the infection process, the same way a series
length is passed to `as_turing_model(latent, n)`.

The group prior carries its own innovations (the [`IID`](@ref) `across` process
samples an `ϵ_t`), which under the prefix-off submodel convention would collide
with the infection [`RandomWalk`](@ref)'s own `ϵ_t`.
[`GroupedIDModel`](@ref) namespaces the group prior through the prior-slot prefix
convention so the two never collide, and prefixes each group's observation with
`group<g>`.
The mapping from a group's effect to its expected series is the `combiner` field.
Its default is a multiplicative effect on the exponential scale,
``\text{expected}_g = e^{\ell_g}\, I_t``, so `exp(ℓ_g)` scales the shared curve.
Swap `combiner` for a different mapping the way [`Ascertainment`](@ref) swaps its
`transform`.

## Simulate

Passing an all-`missing` matrix makes the model a prior simulator; the group
dimension threads from its column count.
We simulate eight groups over 24 time steps:

```@example hier
n_time, n_groups = 24, 8
Ymiss = Matrix{Union{Missing, Float64}}(missing, n_time, n_groups)
sim = as_turing_model(model, Ymiss)()
Ydata = reduce(hcat, [Int.(sim.y[g].y_t) for g in 1:n_groups])
true_levels = sim.group_levels
(n_time = n_time, n_groups = n_groups, data_size = size(Ydata),
    true_levels = round.(true_levels, digits = 2))
```

The shared infection curve is common to all groups; the columns of `Ydata` differ
only through the per-group level and Poisson noise.

## Fit

Conditioning on the simulated counts and sampling with NUTS recovers the
posterior end-to-end.
`n_groups` again threads from the data matrix, so nothing about the group
dimension is hard-coded in the components:

```@example hier
posterior = as_turing_model(model, Float64.(Ydata))
chain = sample(posterior, NUTS(0.8; adtype = Turing.AutoForwardDiff()), 300;
    progress = false)
size(chain, 1)
```

The per-group levels are a generated quantity, recovered per draw with `returned`
and compared with the simulated truth:

```@example hier
level_draws = reduce(hcat, [g.group_levels for g in vec(returned(posterior, chain))])
post_mean = vec(mean(level_draws; dims = 2))
(true_levels = round.(true_levels, digits = 2),
    posterior_means = round.(post_mean, digits = 2),
    correlation = round(cor(true_levels, post_mean), digits = 3))
```

The posterior per-group levels line up with the simulated truth.
A plot with 80% credible intervals against the ``y = x`` line makes the recovery
visible:

```@example hier
using CairoMakie
qs = [quantile(level_draws[g, :], [0.1, 0.5, 0.9]) for g in 1:n_groups]
lo = getindex.(qs, 1)
md = getindex.(qs, 2)
hi = getindex.(qs, 3)

fig = Figure(; size = (620, 460))
ax = Axis(fig[1, 1]; xlabel = "True group level ℓ_g",
    ylabel = "Posterior level ℓ_g")
lims = (minimum(true_levels) - 0.4, maximum(true_levels) + 0.4)
lines!(ax, [lims...], [lims...]; color = :grey, linestyle = :dash)
rangebars!(ax, true_levels, lo, hi; color = :seagreen, whiskerwidth = 10)
scatter!(ax, true_levels, md; color = :seagreen, markersize = 12)
fig
```

Each group's credible interval covers the ``y = x`` line, so the partially pooled
per-group levels are recovered inside a full composed panel.
[`GroupedIDModel`](@ref) supplied the panel structure, the [`Hierarchy`](@ref)
supplied the per-group levels, and the group dimension threaded from the data with
the group prior namespaced by the component.
Swapping `across = RandomWalk()` in the [`Hierarchy`](@ref) would instead pool
*neighbouring* groups (correlated ordered strata), and swapping the `group_effect`
for a bare [`IID`](@ref) or a `Distribution` would drop the shared level for
independent per-group levels, each with no other change.
