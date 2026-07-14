# Transform-the-output latent modifier.

@doc raw"
Apply a transformation function to the output of an inner latent model.

# Arguments

  - `model`: the inner latent model whose output is transformed.
  - `n`: the length of the latent series to generate.

# Examples
```@example TransformLatentModel
using ComposableTuringIDModels, Distributions
trans = TransformLatentModel(Intercept(Normal(2, 0.2)), x -> exp.(x))
rand(as_turing_model(trans, 5))
```

The `model` slot is an [`AbstractPriorModel`](@ref): a bare `Distribution` (or a
vector of them) is coerced via [`as_prior`](@ref), as at the top-level slots.

## Fields

  - `model`: the latent model to transform.
  - `transform`: the transformation function applied to the latent vector.
"
@kwdef struct TransformLatentModel{M <: AbstractLatentModel, F <: Function} <:
              AbstractLatentModel
    "The latent model to transform."
    model::M
    "The transformation function."
    transform::F
end

# Coerce a bare `Distribution` (or vector) member to the prior interface so it is
# accepted alongside a process, matching the top-level slots and Combine/Concat.
function TransformLatentModel(model, transform)
    return TransformLatentModel(as_prior(model), transform)
end

@model function as_turing_model(model::TransformLatentModel, n)
    untransformed ~ to_submodel(as_turing_model(model.model, n))
    return model.transform(untransformed)
end
