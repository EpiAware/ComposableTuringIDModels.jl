# Prefix-the-variables latent modifier (replaces the upstream `prefix_submodel`
# helper via `DynamicPPL.prefix`).

@doc raw"
Wrap an inner latent model so its sampled variables are prefixed with `prefix`.

This replaces the original `prefix_submodel` helper: the inner model is prefixed
with `DynamicPPL.prefix` before being sampled as a submodel, so its variables
appear as `prefix.varname`.

# Arguments

  - `model`: the inner latent model.
  - `n`: the length of the latent series to generate.

# Examples
```@example PrefixLatentModel
using ComposableTuringIDModels
pm = PrefixLatentModel(; model = HierarchicalNormal(), prefix = \"Test\")
rand(as_turing_model(pm, 10))
```

The `model` slot is an [`AbstractPriorModel`](@ref): a bare `Distribution` (or a
vector of them) is coerced via [`as_prior`](@ref), as at the top-level slots.

## Fields

  - `model`: the latent model to prefix.
  - `prefix`: the string prefix applied to the inner model's variables.
"
@kwdef struct PrefixLatentModel{M <: AbstractLatentModel, P <: String} <:
              AbstractLatentModel
    "The latent model."
    model::M
    "The prefix for the latent model."
    prefix::P
end

# Coerce a bare `Distribution` (or vector) member to the prior interface so it is
# accepted alongside a process, matching the top-level slots and Combine/Concat.
function PrefixLatentModel(model, prefix)
    return PrefixLatentModel(as_prior(model), prefix)
end

@model function as_turing_model(model::PrefixLatentModel, n)
    submodel ~ to_submodel(
        prefix(as_turing_model(model.model, n), Symbol(model.prefix)), false)
    return submodel
end
