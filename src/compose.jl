# The headline composition: a full infections → observations model assembled from
# two components, each itself an `as_turing_model`. The latent (parameter) process
# is owned by the infection model, not a separate top-level slot.

@doc raw"
A composed epidemiological model linking an infection process and an observation
model.

The infection process owns its own latent (parameter) process internally, so
a composed model is just two parts: infections, then observations.

Sampling [`as_turing_model(model, y_t, n)`](@ref as_turing_model) runs the two
stages as submodels:

```math
I_t \;\xrightarrow{\text{infections}}\;
y_t \;\xrightarrow{\text{observations}}\; \text{data}
```

The returned generated quantities are `(; generated_y_t, expected_y_t, I_t, Z_t)`.
`generated_y_t` is the observation model's sampled `y_t` (the observed-or-simulated
series, or a `NamedTuple` of streams for a [`Split`](@ref)); `expected_y_t` is its
pre-error `expected` series (the uniform observation return contract). `Z_t` is the
infection model's internal latent draw (e.g. the (log) ``R_t`` path), kept
accessible as a generated quantity, or `nothing` for infection models with no
exposable latent (e.g. [`ODEProcess`](@ref)). Pass `y_t = missing` to simulate
from the prior, or a data vector to condition.

## Fields

  - `infection_model`: the infection process model generating ``I_t`` (and its
    internal latent ``Z_t``).
  - `observation_model`: the observation model mapping ``I_t`` to ``y_t``.

# Examples
```@example IDModel
using ComposableTuringIDModels, Distributions
model = IDModel(
    DirectInfections(; Z = RandomWalk(), initialisation = Normal()),
    PoissonError())
mdl = as_turing_model(model, missing, 20)
rand(mdl)
```
"
struct IDModel{I <: AbstractInfectionModel, O <: AbstractObservationModel} <:
    AbstractComposableModel
    "Infection process model generating ``I_t`` (and its internal latent ``Z_t``)."
    infection_model::I
    "Observation model mapping ``I_t`` to ``y_t``."
    observation_model::O
end

# A composed model's lead-in is its observation model's: the infection process
# runs the full `n` steps and the chain drops the head of the convolution.
observation_lead_in(model::IDModel) = observation_lead_in(model.observation_model)

@doc raw"
Narrow a data argument to a concrete element type before conditioning a model
on it, when nothing in it is actually `missing`.

Data simulated from a model's prior keeps a `Union{Missing, T}` element type
even once every entry is concrete. Conditioning on that makes DynamicPPL
`deepcopy` the argument on every evaluation, which Enzyme cannot differentiate
in reverse mode. Narrowing avoids both the copy and the failure.

A genuinely partially-missing vector (some entries observed, some `missing`)
is split into a [`MissingObservations`](@ref) carrier instead: a concrete
value vector plus a presence mask, neither of which admits `Missing` in its
type. A fully missing vector, a partially- or fully-missing matrix (the
strata data a data-driven [`Split`](@ref) reads shape from), and a scalar
`missing`, are returned unchanged, and so is an array whose element type is
already concrete. A `NamedTuple` of streams is narrowed field-wise, so each
per-stream submodel receives a concrete vector (or a `MissingObservations`
carrier).

[`as_turing_model(::IDModel, y_t, n)`](@ref as_turing_model) applies this to
its `y_t` automatically. It is exposed for data built elsewhere by simulating
from a prior.

# Arguments

  - `y`: the data argument to narrow — an `AbstractArray`, a `NamedTuple` of
    them, or a scalar (e.g. `missing`, returned unchanged).
"
function concrete_observations(y::AbstractArray)
    eltype(y) >: Missing || return y
    return any(ismissing, y) ? y : identity.(y)
end
function concrete_observations(y::AbstractVector)
    eltype(y) >: Missing || return y
    any(ismissing, y) || return identity.(y)
    all(ismissing, y) && return y
    T = nonmissingtype(eltype(y))
    present = .!ismissing.(y)
    value = identity.(coalesce.(y, zero(T)))
    return MissingObservations(value, present)
end
concrete_observations(y::NamedTuple) = map(concrete_observations, y)
concrete_observations(y) = y

# Narrowing must happen before `y_t` is stored on the `Model`: DynamicPPL's
# `hasmissing`/`deepcopy` check runs on the stored argument, so doing it inside
# this body would be too late. Hence the private `@model` plus the wrapper.
@model function _as_turing_model_idmodel(model::IDModel, y_t, n)
    infections ~ as_turing_submodel(model.infection_model, n)
    I_t = infections.I_t
    Z_t = infections.Z_t
    obs ~ as_turing_submodel(model.observation_model, y_t, I_t)
    # Uniform observation contract: sampled series and pre-error expected.
    generated_y_t = obs.y_t
    expected_y_t = obs.expected
    return (; generated_y_t, expected_y_t, I_t, Z_t)
end

function as_turing_model(model::IDModel, y_t, n)
    return _as_turing_model_idmodel(model, concrete_observations(y_t), n)
end

@doc raw"
Convenience 2-argument form: read the infection process's shape from the data.

The observation model and the data together fix the shape of the infection
process, so nothing about it needs to be supplied explicitly or stored on the
model. `as_turing_model(model, Y)` is `as_turing_model(model, Y, shape)` with
`shape` resolved via [`infection_strata`](@ref): the number of infection
strata the observation model consumes given the data's row count, paired with
the data's time length.

Three age strata observed as one hospitalisation stream is
`Split(NegativeBinomialError(), [1.0 1.0 1.0])`; a `1 x T` data matrix then
builds a 3-stratum infection process.

# Examples
```@example IDModel_shape
using ComposableTuringIDModels, Distributions
model = IDModel(
    DirectInfections(;
        Z = Stratify(RandomWalk(), Hierarchy(; across = IID(Normal(0, 0.5)))),
        initialisation = Normal(log(50), 0.2)),
    Split(PoissonError(), [1.0 1.0 1.0]))
# One observation stream, 12 time steps: the three infection strata come
# from the weight matrix, not from the data.
Ymiss = Matrix{Union{Missing, Float64}}(missing, 1, 12)
sim = as_turing_model(model, Ymiss)()
size(sim.I_t)
```
"
function as_turing_model(model::IDModel, Y::AbstractMatrix)
    return as_turing_model(model, Y, _shape(model, Y))
end

function _shape(model::IDModel, Y::AbstractMatrix)
    return (infection_strata(model.observation_model, size(Y, 1)), size(Y, 2))
end

@doc raw"
The number of infection strata an observation model consumes.

The seam an observation model uses to say how many infection strata it
consumes, given the number of observation streams in the data. The default
passes the observation stream count straight through (a one-to-one mapping).
[`Split`](@ref) with a weight `map` overrides it with the map's column count,
so a many-to-one or many-to-many mapping can build the right-shaped infection
process from the data alone.

# Arguments

  - `obs`: the observation model.
  - `n_obs_strata`: the number of observation streams in the data.

# Examples
```@example infection_strata
using ComposableTuringIDModels
ComposableTuringIDModels.infection_strata(
    Split(PoissonError(), [1.0 1.0 1.0]), 1)
```
"
infection_strata(obs, n_obs_strata::Int) = n_obs_strata
function infection_strata(s::Split, n_obs_strata::Int)
    return s.map === nothing ? n_obs_strata : size(s.map, 2)
end
