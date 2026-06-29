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
using EpiAwarePrototype
pm = PrefixLatentModel(; model = HierarchicalNormal(), prefix = \"Test\")
rand(as_turing_model(pm, 10))
```

## Fields

  - `model`: the latent model to prefix.
  - `prefix`: the string prefix applied to the inner model's variables.
"
@kwdef struct PrefixLatentModel{M <: AbstractEpiAwareModel, P <: String} <:
              AbstractEpiAwareModel
    "The latent model."
    model::M
    "The prefix for the latent model."
    prefix::P
end

@model function as_turing_model(model::PrefixLatentModel, n)
    submodel ~ to_submodel(
        prefix(as_turing_model(model.model, n), Symbol(model.prefix)), false)
    return submodel
end
