# Transform-the-output latent modifier.

@doc raw"
Apply a transformation function to the output of an inner latent model.

# Arguments

  - `model`: the inner latent model whose output is transformed.
  - `n`: the shape to generate — a length or an `(n_strata, n_time)` shape,
    whatever the inner model accepts.

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
        return new{typeof(wrapped), F}(wrapped, transform)
    end
end

TransformLatentModel(; model, transform) = TransformLatentModel(model, transform)

# `TransformLatentModel` wraps whatever its inner model returns, so it serves
# both a length-`n` path and an `(n_strata, n_time)` shape. The two methods
# below delegate to one shared `@model`, which avoids both the duplication and
# a dispatch ambiguity against the `AbstractPriorModel`/`Dims{2}` guard.
@model function _transform_latent(model::TransformLatentModel, n)
    untransformed ~ as_turing_submodel(model.model, n)
    return model.transform(untransformed)
end
function as_turing_model(model::TransformLatentModel, n::Int)
    return _transform_latent(model, n)
end
function as_turing_model(model::TransformLatentModel, n::Dims{2})
    return _transform_latent(model, n)
end
