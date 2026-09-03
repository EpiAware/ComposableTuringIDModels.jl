# `IDProblem`: an infection + observation model over a time span.

@doc raw"
A full epidemiological inference problem: an infection process, an observation
model, and a time span. The latent (parameter) process is owned by the infection
model, so it is not a separate slot here.

`as_turing_model(problem, data)` assembles the corresponding [`IDModel`](@ref)
over `tspan` and conditions it on `data.y_t`. The infection process's shape is
read from the observation model and `data.y_t` at build time (see
[`infection_strata`](@ref)), not stored on the problem: a plain vector (or
`missing`) gives a single-series infection process, exactly as today, while a
data matrix or a `NamedTuple` of streams gives a stratified infection process
with one row per stratum.

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

# The infection process's `ModelShape` implied by an observation model and a
# data value, at a given time-axis length. Shared by `IDProblem` (`time_steps`
# from `tspan`) and `forecast` (the fitted length plus the horizon), so the two
# build the shape of a data-driven model the same way. A plain vector gives the
# single-series `time_steps::Int` shape unchanged. A matrix or `NamedTuple` of
# streams gives `(n_strata, time_steps)`, with `n_strata` read from
# `infection_strata` applied to the observation stream count (the matrix's row
# count, or the number of NamedTuple entries).
# With `y_t === missing` there is no data to read a stream count from, so the
# shape falls back to the observation model alone: a `Split` with a weight
# `map` fixes the stratum count at `size(map, 2)`, a `Split` with named
# streams (no map) has as many strata as streams, and anything else is a
# single series.
_obs_data_shape(obs, y_t, time_steps) = time_steps
_obs_data_shape(obs, y_t::Missing, time_steps) = _obs_data_shape_missing(
    obs, time_steps
)
function _obs_data_shape(obs, y_t::AbstractMatrix, time_steps)
    return (infection_strata(obs, size(y_t, 1)), time_steps)
end
function _obs_data_shape(obs, y_t::NamedTuple, time_steps)
    return (infection_strata(obs, length(y_t)), time_steps)
end

_obs_data_shape_missing(obs, time_steps) = time_steps
function _obs_data_shape_missing(s::Split, time_steps)
    s.map === nothing || return (size(s.map, 2), time_steps)
    s.names === nothing || return (length(s.names), time_steps)
    return time_steps
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
