# [Partial pooling across groups in a composed model](@id case-study-hierarchy)

[`Hierarchy`](@ref) is a partial-pooling latent process over a grouping dimension.
Built with `as_turing_model(h, n_groups)` it draws one shared level ``\mu`` and
`n_groups` deviations and returns the **numeric** per-group vector
``\ell_g = \mu + \delta_g``.
Because it returns numbers like any other latent model, it threads into a latent
slot of a full composed model with no group-axis contract change.

This page drives a per-group quantity with a [`Hierarchy`](@ref) inside a complete
infection + observation model, simulates from it and fits it end-to-end under
NUTS, recovering the per-group levels.
It makes explicit how the group dimension threads from the data and how the group
prior must be namespaced so it composes.

## The composed model

The shared epidemic is an [`IDModel`](@ref): a [`DirectInfections`](@ref) process
carrying a [`RandomWalk`](@ref) latent, observed with a [`PoissonError`](@ref).
Every group sees the *same* infection curve ``I_t`` but reports it at its own
level: a per-group log-ascertainment ``\ell_g`` scales the expected counts before
the observation model.
Those per-group levels are partially pooled with a [`Hierarchy`](@ref).

```@example hier
using ComposableTuringIDModels, Distributions, Turing, Random, Statistics
using Turing: to_submodel, returned
Random.seed!(77)

idmodel = IDModel(
    DirectInfections(; Z = RandomWalk(), initialisation = Normal(log(50.0), 0.2)),
    PoissonError())
hierarchy = Hierarchy(; mean = Normal(0.0, 0.5), across = IID(Normal(0.0, 0.5)))
nothing # hide
```

The grouping dimension is **not** a field of any component.
The top-level model reads `n_time` and `n_groups` from the shape of the data
matrix `Y` (rows are time, columns are groups) and passes `n_groups` to the
[`Hierarchy`](@ref) and `n_time` to the infection process — the same way a series
length is passed to `as_turing_model(latent, n)`.

The group prior carries its own innovations (the [`IID`](@ref) `across` process
samples an `ϵ_t`), which under the prefix-off submodel convention would collide
with the infection [`RandomWalk`](@ref)'s own `ϵ_t`.
Namespacing the group prior with a [`PrefixLatentModel`](@ref) keeps the two
apart — this is the group-requirement threading made explicit:

```@example hier
@model function grouped_epidemic(idmodel, hierarchy, Y)
    n_time, n_groups = size(Y)
    # Per-group pooled log-ascertainment levels; namespaced so the group prior's
    # innovations do not collide with the infection process's own latent.
    group_levels ~ to_submodel(
        as_turing_model(PrefixLatentModel(hierarchy, "groups"), n_groups), false)
    # One shared infection process (the IDModel's infection component).
    infections ~ to_submodel(
        as_turing_model(idmodel.infection_model, n_time), false)
    I_t = infections.I_t
    ys = Vector{Any}(undef, n_groups)
    for g in 1:n_groups
        expected_g = exp(group_levels[g]) .* I_t            # group-specific level
        og = PrefixObservationModel(idmodel.observation_model, "group$g")
        y_g ~ to_submodel(as_turing_model(og, Y[:, g], expected_g), false)
        ys[g] = y_g
    end
    return (; I_t, group_levels, y = ys)
end
nothing # hide
```

## Simulate

Passing an all-`missing` matrix makes the model a prior simulator; the group
dimension threads from its column count.
We simulate four groups over 24 time steps:

```@example hier
n_time, n_groups = 24, 4
Ymiss = Matrix{Union{Missing, Float64}}(missing, n_time, n_groups)
sim = grouped_epidemic(idmodel, hierarchy, Ymiss)()
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
`n_groups` again threads from the data matrix — nothing about the group dimension
is hard-coded in the components:

```@example hier
posterior = grouped_epidemic(idmodel, hierarchy, Float64.(Ydata))
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
per-group levels are recovered inside a full composed model.
The [`Hierarchy`](@ref) supplied the per-group quantity numerically, the group
dimension threaded from the data, and namespacing the group prior let it compose
with the infection process's own latent.
Swapping `across = RandomWalk()` would instead pool *neighbouring* groups
(correlated ordered strata) with no other change.
