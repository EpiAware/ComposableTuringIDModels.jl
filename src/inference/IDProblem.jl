# `IDProblem`: an infection + observation model over a time span.

@doc raw"
A full epidemiological inference problem: an infection process, an observation
model, and a time span. The latent (parameter) process is owned by the infection
model, so it is not a separate slot here.

`as_turing_model(problem, data)` assembles the corresponding [`IDModel`](@ref)
over `tspan` and conditions it on `data.y_t`. The infection process's shape is
fixed at build time rather than stored on the problem. The observation model
states what it can of its own stratum count (see [`infection_strata`](@ref)), and
`data.y_t` is read only where the model leaves that count open, along the axis
the component consuming the data calls a stream axis. A plain vector (or
`missing`) gives a single-series infection process, while a data matrix or a
`NamedTuple` of streams gives one infection stratum per stream. Data whose array
dimensions mean something else does not stratify anything: a
[`ReportTriangle`](@ref)'s rows are reference days, and a [`BinomialError`](@ref)'s
`N` field is a trials covariate rather than a second stream.

# Arguments

  - `idproblem`: the [`IDProblem`](@ref).
  - `data`: a value with a `y_t` field holding the observations (or `missing`).

# Examples
```@example IDProblem
using ComposableTuringIDModels, Distributions
problem = IDProblem(
    infection = DirectInfections(; Z = RandomWalk(), initialisation = Normal()),
    observation_model = PoissonError(),
    tspan = (1, 20))
rand(as_turing_model(problem, (; y_t = missing)))
```

## Fields

  - `infection`: the infection process model.
  - `observation_model`: the observation model.
  - `tspan`: the `(first, last)` time span of the observations. An observation
    chain with delays consumes the head of the series it is handed, so the
    infection process runs over a longer span than this; see
    [`observation_lead_in`](@ref).
"
@kwdef struct IDProblem{I <: AbstractInfectionModel, O <: AbstractObservationModel}
    "The infection process model."
    infection::I
    "The observation model."
    observation_model::O
    "The `(first, last)` time span of the observations."
    tspan::Tuple{Int, Int}
end

# The problem's lead-in is its observation model's, and `tspan` is the span of
# the OBSERVATIONS, so a requirements report needs no separate `n`.
function observation_lead_in(idproblem::IDProblem)
    return observation_lead_in(idproblem.observation_model)
end

_observation_chain(idproblem::IDProblem) = idproblem.observation_model

data_requirements(idproblem::IDProblem, data) = data_requirements(
    idproblem, data.y_t, _tspan_length(idproblem)
)

data_requirements(idproblem::IDProblem) = data_requirements(
    idproblem, missing, _tspan_length(idproblem)
)

# Disambiguation: a bare `ModelShape` second argument is a length, not a data
# carrier, so it overrides the problem's own `tspan`.
data_requirements(idproblem::IDProblem, n::ModelShape) = data_requirements(
    idproblem, missing, n
)

_tspan_length(idproblem::IDProblem) =
    idproblem.tspan[end] - idproblem.tspan[1] + 1

@model function as_turing_model(idproblem::IDProblem, data)
    y_t = data.y_t
    time_steps = _tspan_length(idproblem)
    model = IDModel(idproblem.infection, idproblem.observation_model)
    shape = _obs_data_shape(idproblem.observation_model, y_t, time_steps)
    out ~ as_turing_submodel(model, y_t, shape)
    return out
end
