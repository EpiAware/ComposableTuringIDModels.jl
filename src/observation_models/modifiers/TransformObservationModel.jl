# Transform-the-expected-observations modifier.

@doc raw"
Apply a transformation function to the expected observations before passing them
to an inner observation model.

The expected observations `Y_t` are mapped through `transform` and the result is
passed to the inner `model`. The default `transform` applies a softplus
(`x -> log1pexp.(x)`), keeping the transformed expected observations positive.

# Arguments

  - `obs`: the [`TransformObservationModel`](@ref).
  - `y_t`: the observed series (or `missing` when simulating predictively).
  - `Y_t`: the expected-observation series.

# Examples
```@example TransformObservationModel
using EpiAwarePrototype
obs = TransformObservationModel(PoissonError(), x -> x .* 2)
mdl = as_turing_model(obs, missing, fill(10.0, 5))
rand(mdl)
```

## Fields

  - `model`: the inner observation model the transformed expected observations are
    passed to.
  - `transform`: the transformation applied to the expected observations.
"
@kwdef struct TransformObservationModel{M <: AbstractObservationModel, F <: Function} <:
              AbstractObservationModel
    "The inner observation model."
    model::M
    "The transformation applied to the expected observations."
    transform::F = x -> log1pexp.(x)
end

function TransformObservationModel(
        model::M; transform = x -> log1pexp.(x)) where {
        M <: AbstractObservationModel}
    return TransformObservationModel(model, transform)
end

@model function as_turing_model(obs::TransformObservationModel, y_t, Y_t)
    transformed_Y_t = obs.transform(Y_t)
    y_t ~ to_submodel(as_turing_model(obs.model, y_t, transformed_Y_t), false)
    return y_t
end
