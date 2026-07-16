# `IDProblem`: an infection + observation model over a time span.

@doc raw"
A full epidemiological inference problem: an infection process, an observation
model, and a time span. The latent (parameter) process is owned by the infection
model, so it is not a separate slot here.

`as_turing_model(problem, data)` assembles the corresponding [`IDModel`](@ref)
over `tspan` and conditions it on `data.y_t`.

# Arguments

  - `idproblem`: the [`IDProblem`](@ref).
  - `data`: a value with a `y_t` field holding the observations (or `missing`).

# Examples
```@example IDProblem
using ComposableTuringIDModels, Distributions
data = IDData([0.2, 0.3, 0.5], exp)
problem = IDProblem(
    infection = DirectInfections(; Z = RandomWalk(), initialisation = Normal()),
    observation_model = PoissonError(),
    tspan = (1, 20))
rand(as_turing_model(problem, (; y_t = missing)))
```

## Fields

  - `infection`: the infection process model.
  - `observation_model`: the observation model.
  - `tspan`: the `(first, last)` time span of the series.
"
@kwdef struct IDProblem{I <: AbstractInfectionModel, O <: AbstractObservationModel}
    "The infection process model."
    infection::I
    "The observation model."
    observation_model::O
    "The `(first, last)` time span of the series."
    tspan::Tuple{Int, Int}
end

@model function as_turing_model(idproblem::IDProblem, data)
    y_t = data.y_t
    time_steps = idproblem.tspan[end] - idproblem.tspan[1] + 1
    model = IDModel(idproblem.infection, idproblem.observation_model)
    out ~ as_turing_submodel(model, y_t, time_steps)
    return out
end
