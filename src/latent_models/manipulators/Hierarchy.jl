# Partial-pooling / hierarchical structure across groups via an Accessors.jl lens.

@doc raw"
Vary one field of a base model across groups, where the cross-group relationship
is itself a **latent model**.

`Hierarchy` takes a per-unit `base` model, an `Accessors.jl` `lens` identifying
the single field that varies by group, and a `latent` model describing how the
group-level values relate to one another. The number of groups is **not** a field
of the struct: it is supplied at build time (read from the grouping dimension of
the data), exactly the way a series length `n` is passed to
`as_turing_model(latent, n)` rather than stored in the latent struct.

When built with `as_turing_model(h, n_groups)` it

  1. draws `n_groups` group-level values as a length-`n_groups` latent path from
     `latent` (a submodel), then
  2. writes each value onto its own copy of `base` with the lens
     (`Accessors.set(base, lens, value_g)`), returning the vector of `n_groups`
     per-group base-model variants.

Because the cross-group structure is *any* latent model, the pooling behaviour is
parameterised rather than hard-coded:

  - an i.i.d.-Normal latent ([`IID`](@ref)`(Normal())`) or a
    [`HierarchicalNormal`](@ref) gives **classic partial pooling** — each group's
    value is an exchangeable draw shrunk toward a shared mean;
  - a [`RandomWalk`](@ref) latent relates **neighbouring** groups (adjacent
    age-bands / ordered strata), so the group effects are correlated along the
    grouping dimension;
  - any other [`AbstractLatentModel`](@ref) works — the hierarchy is simply a
    latent process *over the grouping dimension*.

The lens may point at any field whose type accepts the drawn value. Pointing it
at a plain **value** field (e.g. `FixedIntercept`'s `intercept`) writes the group
value directly and is the cleanest, AD-friendly case. Pointing it at a **prior**
field (e.g. `NegativeBinomialError`'s `cluster_factor_prior::S <: Sampleable`) is
the case @seabbs flagged on the design: the base declares a concrete placeholder
prior there, and because the field is typed `<: Sampleable` the written value must
itself be a `Sampleable` (e.g. a `Dirac` point mass at the drawn value), so the
group's parameter is sampled deterministically equal to it. That placeholder
prior is never sampled — the documented \"unused prior\" trade-off. The struct
keeps `lens` untyped so either target works.

`Hierarchy` composes with a future prior-role abstraction (`AbstractPriorModel`):
the `latent` slot is the cross-group prior, so when a latent model can play the
prior role directly, a `Hierarchy` is unchanged — it already takes the
cross-group relationship as a pluggable submodel.

# Arguments

  - `h`: the [`Hierarchy`](@ref) model.
  - `n_groups`: the number of groups, read from the grouping dimension of the data
    at build time.

# Examples
```@example Hierarchy
using EpiAwarePrototype, Distributions, Accessors
# Partially pool a per-group intercept (level) across groups: an iid-Normal
# cross-group latent gives classic partial pooling. The lens points at the
# `intercept` value field, so each group's drawn value is written straight onto
# its own variant. n_groups (3 here) is supplied at build time.
base = FixedIntercept(0.0)
h = Hierarchy(base, (@optic _.intercept), IID(Normal(0.0, 1.0)))
variants = as_turing_model(h, 3)()
length(variants)
```

## Fields

  - `base`: the per-unit base model whose `lens` field varies by group.
  - `lens`: an `Accessors.jl` optic identifying the field that varies by group.
  - `latent`: the latent model describing the cross-group relationship.
"
struct Hierarchy{M <: AbstractEpiAwareModel, L, T <: AbstractLatentModel} <:
       AbstractLatentModel
    "The per-unit base model whose `lens` field varies by group."
    base::M
    "An `Accessors.jl` optic identifying the field that varies by group."
    lens::L
    "The latent model describing the cross-group relationship."
    latent::T

    function Hierarchy(base::AbstractEpiAwareModel, lens, latent::AbstractLatentModel)
        # Fail fast if the lens does not point at a field of `base`: `get` it once
        # at construction so a bad optic errors here, not at sample time.
        Accessors.getall(base, lens)
        return new{typeof(base), typeof(lens), typeof(latent)}(base, lens, latent)
    end
end

function Hierarchy(; base::AbstractEpiAwareModel, lens,
        latent::AbstractLatentModel = IID(Normal()))
    return Hierarchy(base, lens, latent)
end

@doc raw"
Build the `n_groups` per-group base-model variants of a [`Hierarchy`](@ref) by
writing each group-level value onto its own copy of `base` with `lens`.

`group_values` is the length-`n_groups` cross-group latent path; `set_per_group`
applies `Accessors.set(base, lens, group_values[g])` for each group and returns
the vector of variants. This is the construction mechanism — separated from the
sampling of `group_values` so it can be reused and tested directly.

# Arguments

  - `base`: the per-unit base model.
  - `lens`: the `Accessors.jl` optic for the field that varies by group.
  - `group_values`: the per-group values to write onto `lens`.

# Examples
```@example
using EpiAwarePrototype, Distributions, Accessors
base = FixedIntercept(0.0)
EpiAwarePrototype.set_per_group(base, (@optic _.intercept), [0.1, 0.2])
```
"
function set_per_group(base::AbstractEpiAwareModel, lens, group_values)
    return map(v -> Accessors.set(base, lens, v), group_values)
end

@model function as_turing_model(h::Hierarchy, n_groups::Int)
    @assert n_groups>0 "n_groups must be greater than 0"
    # Draw the n_groups group-level values as a latent path over the grouping
    # dimension; the latent model parameterises the cross-group relationship
    # (iid-Normal ⇒ partial pooling, RandomWalk ⇒ correlated neighbours, ...).
    group_values ~ to_submodel(as_turing_model(h.latent, n_groups), false)
    # Write each value onto its own variant of the base model via the lens.
    variants = set_per_group(h.base, h.lens, group_values)
    return variants
end
