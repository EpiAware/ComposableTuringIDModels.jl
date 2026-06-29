# `EpiProblem`: an infection + observation model over a time span.

@doc raw"
A full epidemiological inference problem: an infection process, an observation
model, and a time span. The latent (parameter) process is owned by the infection
model, so it is not a separate slot here.

`as_turing_model(problem, data)` assembles the corresponding [`EpiAwareModel`](@ref)
over `tspan` and conditions it on `data.y_t`.

# Arguments

  - `epiproblem`: the [`EpiProblem`](@ref).
  - `data`: a value with a `y_t` field holding the observations (or `missing`).

# Examples
```@example EpiProblem
using EpiAwarePrototype, Distributions
data = EpiData([0.2, 0.3, 0.5], exp)
problem = EpiProblem(
    epi_model = DirectInfections(; Z = RandomWalk(), initialisation_prior = Normal()),
    observation_model = PoissonError(),
    tspan = (1, 20))
rand(as_turing_model(problem, (; y_t = missing)))
```

## Fields

  - `epi_model`: the infection process model.
  - `observation_model`: the observation model.
  - `tspan`: the `(first, last)` time span of the series.
"
@kwdef struct EpiProblem{I <: AbstractInfectionModel, O <: AbstractObservationModel}
    "The infection process model."
    epi_model::I
    "The observation model."
    observation_model::O
    "The `(first, last)` time span of the series."
    tspan::Tuple{Int, Int}
end

@model function as_turing_model(epiproblem::EpiProblem, data)
    y_t = data.y_t
    time_steps = epiproblem.tspan[end] - epiproblem.tspan[1] + 1
    model = EpiAwareModel(epiproblem.epi_model, epiproblem.observation_model)
    out ~ to_submodel(as_turing_model(model, y_t, time_steps), false)
    return out
end
