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

The `model` slot is a length-`n` PATH slot: a bare `Distribution` there is
auto-wrapped in an [`Intercept`](@ref), giving a constant inner path; a process,
an [`IID`](@ref), or a vector passes through. Use [`IID`](@ref) for `n`
independent draws. It is composed through [`as_turing_submodel`](@ref).

## Fields

  - `model`: the latent model to transform.
  - `transform`: the transformation function applied to the latent vector.
"
struct TransformLatentModel{M <: PriorLike, F <: Function} <:
       AbstractLatentModel
    "The latent model to transform."
    model::M
    "The transformation function."
    transform::F

    function TransformLatentModel(model, transform::F) where {F <: Function}
        # `model` is a length-`n` PATH slot: a bare `Distribution` is wrapped in
        # an `Intercept` (a constant inner path), never left as a scalar.
        wrapped = _path_prior(model)
        new{typeof(wrapped), F}(wrapped, transform)
    end
end

TransformLatentModel(; model, transform) = TransformLatentModel(model, transform)

@model function as_turing_model(model::TransformLatentModel, n)
    untransformed ~ as_turing_submodel(model.model, n)
    return model.transform(untransformed)
end
