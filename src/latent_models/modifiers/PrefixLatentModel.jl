# Prefix-the-variables latent modifier, via `DynamicPPL.prefix`.

@doc raw"
Wrap an inner latent model so its sampled variables are prefixed with `prefix`.

The inner model is prefixed with `DynamicPPL.prefix` before being sampled as a
submodel, so its variables appear as `prefix.varname`.

# Arguments

  - `model`: the inner latent model.
  - `n`: the shape to generate — a length or an `(n_strata, n_time)` shape,
    whatever the inner model accepts.

# Examples
```@example PrefixLatentModel
using ComposableTuringIDModels
pm = PrefixLatentModel(; model = HierarchicalNormal(), prefix = \"Test\")
rand(as_turing_model(pm, 10))
```

The `model` slot takes a raw component: a latent model, or a `Distribution` (or a
vector of them).

## Fields

  - `model`: the latent model to prefix.
  - `prefix`: the string prefix applied to the inner model's variables.
"
@kwdef struct PrefixLatentModel{M <: PriorLike, P <: String} <:
    AbstractLatentModel
    "The latent model."
    model::M
    "The prefix for the latent model."
    prefix::P
end

# `PrefixLatentModel` wraps whatever its inner model returns, so it serves both a
# length-`n` path and an `(n_strata, n_time)` shape.
# The two methods below delegate to one shared `@model`, which avoids a dispatch
# ambiguity against the `AbstractPriorModel`/`Dims{2}` guard.
@model function _prefix_latent(model::PrefixLatentModel, n)
    submodel ~ to_submodel(
        prefix(as_turing_model(model.model, n), Symbol(model.prefix)), false
    )
    return submodel
end
as_turing_model(model::PrefixLatentModel, n::Int) = _prefix_latent(model, n)
as_turing_model(model::PrefixLatentModel, n::Dims{2}) = _prefix_latent(model, n)

# Widen a raw prior into a path and namespace it under `prefix`, for the
# components that own a prefix slot and store the prior already wrapped.
# Applying it twice gives the same object, which is what lets the same call sit
# in the public constructor and in the `ConstructionBase.constructorof` rebuild:
# a raw prior set into the slot is wrapped, and the stored wrapped prior is left
# alone.
# The slot's prefix is the single place its name is set, so a
# `PrefixLatentModel` put in the slot is re-prefixed to it rather than nested
# under it.
_prefixed_path(model, prefix::AbstractString) =
    prefix == "" ? path_prior(model) : PrefixLatentModel(path_prior(model), prefix)

function _prefixed_path(model::PrefixLatentModel, prefix::AbstractString)
    prefix == "" && return model.model
    prefix == model.prefix && return model
    return PrefixLatentModel(model.model, prefix)
end
