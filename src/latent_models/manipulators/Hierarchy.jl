# Partial-pooling / hierarchical structure across groups, returning a numeric
# per-group latent path.

@doc raw"
Partially pool a per-group level across groups, where the cross-group
relationship is itself a prior model.

`Hierarchy` is a non-centred partial-pooling latent process over a grouping
dimension. Built with `as_turing_model(h, n_groups)` it draws a single shared
level ``\mu`` from `mean` and `n_groups` group deviations from `across`, and
returns the **numeric** length-`n_groups` vector

```math
\ell_g = \mu + \delta_g, \qquad g = 1, \ldots, n_{\text{groups}},
```

with ``\delta`` the group deviations. It returns numeric values (like every other
latent model), so it threads straight into any latent slot — e.g. an infection
model's `Z` — with no group-axis contract change.

The number of groups is **not** a field of the struct: it is supplied at build
time (read from the grouping dimension of the data), exactly the way a series
length `n` is passed to `as_turing_model(latent, n)`.

Both slots take a raw prior (a bare `Distribution`, or a latent/prior model),
sampled through [`as_turing_submodel`](@ref), so the pooling behaviour is
parameterised rather than hard-coded through `across`:

  - an i.i.d.-Normal ([`IID`](@ref)`(Normal())`, the default) gives **classic
    exchangeable partial pooling** — each group's deviation is an independent
    draw shrunk toward the shared level;
  - a [`RandomWalk`](@ref) relates **neighbouring** groups (adjacent age-bands /
    ordered strata), so the group effects are correlated along the grouping
    dimension;
  - any other [`AbstractLatentModel`](@ref) works — the hierarchy is a prior
    process *over the grouping dimension*.

This is the numeric, contract-compliant partial-pooling construct: it returns
length-`n_groups` values rather than per-group model variants, and takes its
cross-group relationship through the prior interface.

## Fields

  - `mean`: prior for the shared level ``\mu`` (a `Distribution` or prior model).
  - `across`: the cross-group relationship generating the group deviations
    (default `IID(Normal())`). A length-`n_groups` PATH slot: a bare
    `Distribution` here is auto-wrapped in an [`Intercept`](@ref) — a constant
    deviation shared across groups — so use [`IID`](@ref) for exchangeable
    (independent) group deviations.

# Examples
```@example Hierarchy
using ComposableTuringIDModels, Distributions
# Partially pool a per-group level across 3 groups with classic (exchangeable)
# pooling; n_groups is supplied at build time.
h = Hierarchy(; across = IID(Normal(0.0, 1.0)))
length(as_turing_model(h, 3)())
```
"
struct Hierarchy{M <: PriorLike, A <: PriorLike} <: AbstractLatentModel
    "Prior for the shared level ``\\mu``."
    mean::M
    "Cross-group relationship generating the group deviations."
    across::A

    function Hierarchy(mean, across)
        # `across` is a length-`n_groups` PATH slot (one deviation per group), so
        # a bare `Distribution` is wrapped in an `Intercept` (a constant, shared
        # deviation broadcast across groups), never left as a scalar. `mean` is a
        # single shared level, so it stays a bare scalar prior.
        wrapped = _path_prior(across)
        new{typeof(mean), typeof(wrapped)}(mean, wrapped)
    end
end

function Hierarchy(; mean = Normal(), across = IID(Normal()))
    return Hierarchy(mean, across)
end

@model function as_turing_model(h::Hierarchy, n_groups::Int)
    @assert n_groups>0 "n_groups must be greater than 0"
    # One shared level μ (a single global draw), plus n_groups deviations from the
    # cross-group prior: the latent parameterises the pooling (iid ⇒ exchangeable,
    # RandomWalk ⇒ correlated neighbours). The result is a numeric per-group path.
    hierarchy_mean ~ as_turing_submodel(h.mean, 1; prefix = true)
    group_effects ~ as_turing_submodel(h.across, n_groups; prefix = true)
    return only(hierarchy_mean) .+ group_effects
end
