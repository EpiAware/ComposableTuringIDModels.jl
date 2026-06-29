# Prefix-the-variables observation modifier (replaces the upstream
# `prefix_submodel` helper via `DynamicPPL.prefix`).

@doc raw"
Wrap an inner observation model so its sampled variables are prefixed with
`prefix`.

This replaces the original `prefix_submodel` helper for observation models: the
inner model is prefixed with `DynamicPPL.prefix` before being sampled as a
submodel, so its variables appear as `prefix.varname`.

# Arguments

  - `observation_model`: the [`PrefixObservationModel`](@ref).
  - `y_t`: the observed series (or `missing` when simulating predictively).
  - `Y_t`: the expected-observation series.

# Examples
```@example PrefixObservationModel
using EpiAwarePrototype
pm = PrefixObservationModel(; model = PoissonError(), prefix = \"Test\")
mdl = as_turing_model(pm, missing, fill(10.0, 5))
rand(mdl)
```

## Fields

  - `model`: the inner observation model to prefix.
  - `prefix`: the string prefix applied to the inner model's variables.
"
@kwdef struct PrefixObservationModel{M <: AbstractEpiAwareModel, P <: String} <:
              AbstractEpiAwareModel
    "The observation model."
    model::M
    "The prefix for the observation model."
    prefix::P
end

@model function as_turing_model(observation_model::PrefixObservationModel, y_t, Y_t)
    submodel ~ to_submodel(
        prefix(as_turing_model(observation_model.model, y_t, Y_t),
            Symbol(observation_model.prefix)), false)
    return submodel
end
