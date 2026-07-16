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

The `model` slot takes a raw component: a latent model, or a `Distribution` (or a
vector of them), composed through [`as_turing_submodel`](@ref).

## Fields

  - `model`: the latent model to transform.
  - `transform`: the transformation function applied to the latent vector.
"
@kwdef struct TransformLatentModel{M <: PriorLike, F <: Function} <:
              AbstractLatentModel
    "The latent model to transform."
    model::M
    "The transformation function."
    transform::F
end

@model function as_turing_model(model::TransformLatentModel, n)
    untransformed ~ as_turing_submodel(model.model, n)
    return model.transform(untransformed)
end
