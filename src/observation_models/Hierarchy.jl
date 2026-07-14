# Partial pooling / hierarchical structure across observation streams, built from
# a single base model with Accessors.jl lenses (port of upstream #255; see #38).
#
# The construction mechanism: keep one fully-specified `base` observation model,
# draw a per-group vector from a `latent` relationship model, then overwrite a
# deeply-nested field of `base` with each group's value using an Accessors optic
# (`set(base, lens, value)`). This yields one variant per group with no shadow
# structs, and the variants are stacked as named streams with `Split`. Because
# lenses compose, the pooled field can sit at any nesting depth (e.g. inside a
# `LatentDelay`- or `Ascertainment`-wrapped stream).

@doc raw"
Partially pool a field of an observation model across several data streams.

`Hierarchy` builds one model variant per data stream from a single **base**
observation model by overwriting a chosen deeply-nested field with a per-group
value drawn from a **latent relationship model**, then stacks the variants as
named streams (via [`Split`](@ref)). This shares structure across streams while
letting the pooled field vary — the middle ground between full pooling
(identical streams) and no pooling (independent streams).

The cross-group relationship is *any* [`AbstractLatentModel`](@ref) over the
grouping dimension, so partial pooling is one case among many:

  - `IID(Normal())` → independent group values (no pooling);
  - [`HierarchicalNormal`](@ref) → non-centred partial pooling (a shared spread
    hyperprior with per-group offsets);
  - [`RandomWalk`](@ref) → neighbouring groups/age-bands related.

The number of groups is **data-derived**: it is read from the `y_t` NamedTuple at
model-build time (one group per stream), not fixed at construction.

## The pooled field must be consumed as a value (injection-point constraint)

The lens target is overwritten with a plain sampled *value*, so the base model
must **use** that field as a value rather than sampling it with `~`. Fields that
feed a computation directly — a [`FixedIntercept`](@ref) ascertainment effect, a
scale, a transform argument — pool cleanly and are gradient-safe under HMC/NUTS.

A field that the base samples with `~` as a bare `Distribution` prior (e.g.
`NegativeBinomialError`'s `cluster_factor_prior`) is **not** poolable this way:
writing a value onto it would degenerate the `~` into a point mass, which
DynamicPPL rejects under gradient-based inference as non-differentiable. Pooling
such prior fields (the headline pooled-overdispersion case in #255) needs the
priors-as-submodels seam (#37/#372) so the field becomes a value-pass-through
submodel; `Hierarchy` then pools it with the same lens. See issue #38 for the
sequencing and the open API questions (lens ergonomics, placeholder priors).

## Fields

  - `base`: the base observation model whose field is pooled across streams.
  - `lens`: an `Accessors` optic onto the field to pool (e.g.
    `Accessors.@optic _.latent_model.intercept`). Lenses compose, so the field
    may be nested arbitrarily deep.
  - `latent`: the [`AbstractLatentModel`](@ref) describing the cross-group
    relationship; sampled at length `n_groups` to produce one value per stream.

## Data contract

`y_t` is a NamedTuple of observed series keyed by stream name (or a NamedTuple of
`missing` to simulate); its keys name and count the groups. `Y_t` is the expected
series passed on to [`Split`](@ref) (a shared vector broadcast to every stream, or
a per-stream NamedTuple). The return value is `Split`'s `(; y_t, expected)`.

# Examples
```@example Hierarchy
using ComposableTuringIDModels, Distributions, Accessors
# Pool a reporting-fraction (ascertainment) effect across cases and deaths, with
# a non-centred normal relationship between the two streams' log-ascertainment.
base = Ascertainment(PoissonError(), FixedIntercept(0.0); latent_prefix = \"\")
pooled = Hierarchy(base, (@optic _.latent_model.intercept), HierarchicalNormal())
mdl = as_turing_model(pooled, (cases = missing, deaths = missing), fill(100.0, 12))
rand(mdl)
```
"
struct Hierarchy{B <: AbstractObservationModel, L, M <: AbstractLatentModel} <:
       AbstractObservationModel
    "The base observation model whose field is pooled across streams."
    base::B
    "An `Accessors` optic onto the field to pool (may be nested arbitrarily deep)."
    lens::L
    "The latent model describing the cross-group relationship (sampled at `n_groups`)."
    latent::M
end

# Ordered group names from the data: a Hierarchy fans out one stream per `y_t`
# entry, so the grouping dimension is data-derived (never stored on the struct).
function _hierarchy_names(y_t)
    y_t isa NamedTuple ||
        error("A Hierarchy needs a NamedTuple `y_t` to derive its groups")
    return collect(string.(keys(y_t)))
end

@model function as_turing_model(m::Hierarchy, y_t, Y_t)
    snames = _hierarchy_names(y_t)
    n_groups = length(snames)

    # Draw one value per group from the cross-group relationship. Its internal
    # variables stay flat (prefix off) — they are shared hyperparameters sitting
    # above the per-stream Split.
    group_values ~ to_submodel(as_turing_model(m.latent, n_groups), false)

    # Build one variant per group by overwriting the pooled field with its value,
    # then stack the variants as named streams. `set` returns a modified copy, so
    # `base` is untouched and the field can sit at any nesting depth.
    variants = map(v -> Accessors.set(m.base, m.lens, v), group_values)
    streams = NamedTuple{Tuple(Symbol.(snames))}(Tuple(variants))

    obs ~ to_submodel(as_turing_model(Split(streams), y_t, Y_t), false)
    return obs
end
