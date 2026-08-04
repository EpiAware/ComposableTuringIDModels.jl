# [Partial pooling across groups in a composed model](@id case-study-hierarchy)

A multi-group epidemic is a panel: one shared infection process seen by several
groups, each reporting it at its own level.
[`IDModel`](@ref) expresses this panel with the same constructor as a
single-group model, just with a `group_effect` prior threaded between the
infection and observation arguments:
`IDModel(infection_model, group_effect, observation_model)`.
Internally this builds a [`GroupedInfections`](@ref) (the shared curve
replicated per group) observed through a [`Split`](@ref) (every group through
the same observation model, namespaced by group) — see the [Composable
design](@ref) page for how the two compose.

This page drives a per-group reporting level with a [`Hierarchy`](@ref) inside
a grouped [`IDModel`](@ref), simulates from it and fits it end-to-end under
NUTS, recovering the per-group levels.
The group dimension threads from the data and the group prior is namespaced by
the component, so the panel composes with no hand-orchestration.

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

hierarchy = Hierarchy(; mean = Normal(0.0, 0.5), across = IID(Normal(0.0, 0.5)))
model = IDModel(
    DirectInfections(; Z = RandomWalk(), initialisation = Normal(log(50.0), 0.2)),
    hierarchy,
    PoissonError())
```

The grouping dimension is **not** a field of any component.
Passing `as_turing_model(model, Y)` reads `n_groups` and `n_time` from the
shape of the data matrix `Y` (rows are groups, columns are time) and passes
`n_groups` to the [`Hierarchy`](@ref) and `n_time` to the infection process,
the same way a series length is passed to `as_turing_model(latent, n)`.

The group prior carries its own innovations (the [`IID`](@ref) `across` process
samples an `ϵ_t`), which under the prefix-off submodel convention would collide
with the infection [`RandomWalk`](@ref)'s own `ϵ_t`.
[`GroupedInfections`](@ref) namespaces the group prior through the prior-slot
prefix convention so the two never collide, and [`Split`](@ref) prefixes each
group's observation with `group<g>`.
The mapping from a group's effect to its row of the infection matrix is the
`combiner` keyword argument.
Its default is a multiplicative effect on the exponential scale,
``\text{row}_g = e^{\ell_g}\, I_t``, so `exp(ℓ_g)` scales the shared curve.
Swap `combiner` for a different mapping the way [`Ascertainment`](@ref) swaps its
`transform`.

## Simulate

Passing an all-`missing` matrix makes the model a prior simulator; the group
dimension threads from its row count.
We simulate eight groups over 24 time steps:

```@example hier
n_time, n_groups = 24, 8
Ymiss = Matrix{Union{Missing, Float64}}(missing, n_groups, n_time)
sim = as_turing_model(model, Ymiss)()
Ydata = Float64.(reduce(vcat, [permutedims(sim.generated_y_t[g]) for g in 1:n_groups]))
true_levels = sim.Z_t.group_levels
(n_time = n_time, n_groups = n_groups, data_size = size(Ydata),
    true_levels = round.(true_levels, digits = 2))
```

The shared infection curve is common to all groups; the rows of `Ydata` differ
only through the per-group level and Poisson noise.

## Fit

Conditioning on the simulated counts and sampling with NUTS recovers the
posterior end-to-end.
`n_groups` again threads from the data matrix, so nothing about the group
dimension is hard-coded in the components:

```@example hier
posterior = as_turing_model(model, Ydata)
chain = sample(posterior, NUTS(0.85; adtype = Turing.AutoForwardDiff()), 300;
    progress = false)
size(chain, 1)
```

The per-group levels are a generated quantity, recovered per draw with `returned`
and compared with the simulated truth:

```@example hier
level_draws = reduce(hcat, [g.Z_t.group_levels for g in vec(returned(posterior, chain))])
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
`IDModel`'s grouped constructor supplied the panel structure, the
[`Hierarchy`](@ref) supplied the per-group levels, and the group dimension
threaded from the data with the group prior namespaced by the component.
Swapping `across = RandomWalk()` in the [`Hierarchy`](@ref) would instead pool
*neighbouring* groups (correlated ordered strata), and swapping the `group_effect`
for a bare [`IID`](@ref) or a `Distribution` would drop the shared level for
independent per-group levels, each with no other change.

When the groups are genuinely **separate** infection processes rather than one
shared curve — several distinct regions, say, each with its own latent — see
[`CombineInfections`](@ref) instead, described on the [Multiple observation
streams](@ref case-study-split) page.
