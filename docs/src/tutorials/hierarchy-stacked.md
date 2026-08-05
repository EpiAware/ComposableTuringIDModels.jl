# [Partial pooling across groups in a composed model](@id tutorial-hierarchy)

A multi-group epidemic is a panel: one shared infection process seen by several
groups, each reporting it at its own level.
[`Stratify`](@ref) expresses this panel directly: a shared process drawn once
over time, and a per-group deviation drawn once over the group axis, combined
into one `groups × time` latent matrix.
[`Split`](@ref) then observes every group through the same observation model,
namespaced by group.
See the [Composable design](@ref) page for how the two compose.

This page drives a per-group reporting level with a [`Hierarchy`](@ref) inside
a [`Stratify`](@ref), simulates from it and fits it end-to-end under NUTS,
recovering the per-group levels.
The group dimension threads from the data and the group prior is namespaced by
the component, so the panel composes with no hand-orchestration.

## The composed model

The shared epidemic is a [`DirectInfections`](@ref) process carrying a
[`RandomWalk`](@ref) latent, observed with a [`PoissonError`](@ref).
[`Stratify`](@ref) puts a per-group level on that latent: every group shares
the same random walk, offset by its own log-level ``\ell_g``, so each group's
infection curve is the shared curve scaled by ``e^{\ell_g}``.
Those per-group levels are partially pooled with a [`Hierarchy`](@ref),
supplied as `Stratify`'s `across` slot.

```@example hier
using ComposableTuringIDModels, Distributions, Turing, Random, Statistics
using Turing: returned
Random.seed!(77)

hierarchy = Hierarchy(; mean = Normal(0.0, 0.5), across = IID(Normal(0.0, 0.5)))
model = IDModel(
    DirectInfections(; Z = Stratify(RandomWalk(), hierarchy),
        initialisation = Normal(log(50.0), 0.2)),
    Split(PoissonError()))
```

The grouping dimension is **not** a field of any component.
Passing `as_turing_model(model, Y)` reads `n_groups` and `n_time` from the
shape of the data matrix `Y` (rows are groups, columns are time) and passes
`n_groups` to the [`Hierarchy`](@ref) (through `Stratify`'s `across` slot) and
`n_time` to the shared random walk, the same way a series length is passed to
`as_turing_model(latent, n)`.

The group prior carries its own innovations (the [`IID`](@ref) `across` process
samples an `ϵ_t`), which under the prefix-off submodel convention would collide
with the infection [`RandomWalk`](@ref)'s own `ϵ_t`.
[`Stratify`](@ref) prefixes its `across` slot automatically so the two never
collide, and [`Split`](@ref) prefixes each group's observation with
`group<g>`.
`Stratify`'s `combine` argument maps a group's level and the shared path onto
that group's row of the latent matrix, additively by default:
``Z_{g,t} = \text{shared}_t + \ell_g``, so after the model's `exp`
transformation each group's curve is the shared curve scaled by ``e^{\ell_g}``.
Swap `combine` for a different mapping the way [`Ascertainment`](@ref) swaps its
`transform`.

## Simulate

Passing an all-`missing` matrix makes the model a prior simulator; the group
dimension threads from its row count.
We simulate eight groups over 24 time steps:

```@example hier
n_time, n_groups = 24, 8
Ymiss = Matrix{Union{Missing, Float64}}(missing, n_groups, n_time)
sim = as_turing_model(model, Ymiss)()
Ydata = Float64.(reduce(vcat,
    [permutedims(sim.generated_y_t[Symbol(:group, g)]) for g in 1:n_groups]))
# Z_t[g, :] is the shared random walk plus group g's pooled level. Averaging
# over time and subtracting the grand mean cancels the shared component and
# leaves each group's level relative to the others. Only that relative level
# is identified: a constant added to every group's level and subtracted from
# the shared path gives the same Z_t, so the two are confounded.
group_level(Z) = vec(mean(Z; dims = 2)) .- mean(Z)
true_levels = group_level(sim.Z_t)
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

The relative per-group levels are recovered per draw from the generated `Z_t`
with
`returned`, and compared with the simulated truth:

```@example hier
level_draws = reduce(hcat,
    [group_level(g.Z_t) for g in vec(returned(posterior, chain))])
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

Each group's credible interval covers the ``y = x`` line, so the partially
pooled per-group levels are recovered inside a full composed panel.
[`Stratify`](@ref) supplied the panel structure, the [`Hierarchy`](@ref)
supplied the per-group levels, and the group dimension threaded from the data
with the group prior namespaced by the component.
Swapping `across = RandomWalk()` in the [`Hierarchy`](@ref) would instead pool
*neighbouring* groups (correlated ordered strata), and swapping `Stratify`'s
`across` slot for a bare [`IID`](@ref) or a `Distribution` would drop the
shared level for independent per-group levels, each with no other change.

When the groups are genuinely **separate** infection processes rather than one
shared curve — several distinct regions, say, each with its own latent — see
[`CombineInfections`](@ref) instead, described on the [Multiple observation
streams](@ref tutorial-split) page.
