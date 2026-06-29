# Record-the-expected-observations modifier.

@doc raw"
Record the expected observations as a tracked generated quantity (`exp_y_t`).

The expected observations `Y_t` are tracked via the `:=` syntax before the inner
`model` is applied unchanged, so the expected observations are available in the
returned chain alongside the inner model's variables.

# Arguments

  - `model`: the [`RecordExpectedObs`](@ref) model.
  - `y_t`: the observed series (or `missing` when simulating predictively).
  - `Y_t`: the expected-observation series.

# Examples
```@example RecordExpectedObs
using EpiAwarePrototype
obs = RecordExpectedObs(PoissonError())
mdl = as_turing_model(obs, missing, fill(10.0, 5))
rand(mdl)
```

## Fields

  - `model`: the inner observation model whose expected observations are recorded.
"
struct RecordExpectedObs{M <: AbstractEpiAwareModel} <: AbstractEpiAwareModel
    "The inner observation model whose expected observations are recorded."
    model::M
end

@model function as_turing_model(model::RecordExpectedObs, y_t, Y_t)
    exp_y_t := Y_t
    y_t ~ to_submodel(as_turing_model(model.model, y_t, Y_t), false)
    return y_t
end
