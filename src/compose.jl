# The latent parameter process is owned by the infection model rather than by a
# separate top-level slot.

@doc raw"
A composed epidemiological model linking an infection process and an observation
model.

The infection process owns its own latent (parameter) process internally, so
a composed model is just two parts: infections, then observations.

Sampling [`as_turing_model(model, y_t, n)`](@ref as_turing_model) runs the two
stages as submodels, with `n` the number of OBSERVATIONS:

```math
I_t \;\xrightarrow{\text{infections}}\;
y_t \;\xrightarrow{\text{observations}}\; \text{data}
```

The returned generated quantities are `(; generated_y_t, expected_y_t, I_t, Z_t)`,
with an `I_seed` for an infection model that exposes a seeding window (a
[`Renewal`](@ref) does).
`generated_y_t` is the observation model's sampled `y_t` (the observed-or-simulated
series, or a `NamedTuple` of streams for a [`Split`](@ref)); `expected_y_t` is its
pre-error `expected` series (the uniform observation return contract). `Z_t` is the
infection model's internal latent draw (e.g. the (log) ``R_t`` path), kept
accessible as a generated quantity, or `nothing` for infection models with no
exposable latent (e.g. [`ODEProcess`](@ref)). Pass `y_t = missing` to simulate
from the prior, or a data vector to condition.

An observation chain with delays in it consumes the head of the series it is
given, so the infection process is run over `n + observation_lead_in(model)`
time points and the observation model scores all `n` observations. That is
derived from the chain, not asked of the caller: `n` is always the number of
observations, and supplying a different number of them is an error rather than
a silent reinterpretation. [`data_requirements`](@ref) reports what a model
needs before it is built.

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

# A composed model's lead-in is its observation model's, and so is the chain a
# requirements report is read off.
observation_lead_in(model::IDModel) = observation_lead_in(model.observation_model)
_observation_chain(model::IDModel) = model.observation_model

@doc raw"
Narrow a data argument to a concrete element type before conditioning a model
on it, when nothing in it is actually `missing`.

Data simulated from a model's prior keeps a `Union{Missing, T}` element type
even once every entry is concrete. Conditioning on that makes DynamicPPL
`deepcopy` the argument on every evaluation, which Enzyme cannot differentiate
in reverse mode. Narrowing avoids both the copy and the failure.

A vector with any `missing` in it is split into a [`MissingObservations`](@ref)
carrier instead: a concrete value vector plus a presence mask, neither of which
admits `Missing` in its type. This is also what detaches the data from the
caller: the carrier is scored by reading it, so a blank's draw can never be
written back into the array that was handed over. Data reached through a
`NamedTuple` field or a struct field is not a model argument, so this is the
only thing standing between it and the `y_t[i] ~ …` sugar's `setindex!!`.

A partially- or fully-missing matrix (the strata data a data-driven
[`Split`](@ref) reads shape from) and a scalar `missing` are returned
unchanged, and so is an array whose element type is already concrete. A
`Vector{Missing}` is returned unchanged too: it carries no value type to build
a carrier from, and nothing can be written into it either (the sugar widens it
into a fresh array). A `NamedTuple` of streams is narrowed field-wise, so each
per-stream submodel receives a concrete vector (or a `MissingObservations`
carrier).

[`as_turing_model(::IDModel, y_t, n)`](@ref as_turing_model) applies this to
its `y_t` automatically, and so does every observation model that reads its
series out of a `NamedTuple` or a struct. It is exposed for data built
elsewhere by simulating from a prior.

A gap in a fitted series is marginalised out, as described under
[`MissingObservations`](@ref).

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
    T = nonmissingtype(eltype(y))
    # A `Vector{Missing}` has no value type to put in the carrier, and nothing
    # can be written into it either, so leave it.
    # Any other element type without a `zero` is copied instead, so the caller's
    # array is still out of reach.
    T === Union{} && return y
    (isconcretetype(T) && T <: Number) || return copy(y)
    present = .!ismissing.(y)
    value = identity.(coalesce.(y, zero(T)))
    return MissingObservations(value, present)
end
concrete_observations(y::NamedTuple) = map(concrete_observations, y)
concrete_observations(y) = y

# Narrowing must happen before `y_t` is stored on the `Model`, because
# DynamicPPL's `hasmissing`/`deepcopy` check runs on the stored argument.
# Hence the private `@model` plus the wrapper.
@model function _as_turing_model_idmodel(model::IDModel, y_t, n)
    infections ~ as_turing_submodel(model.infection_model, n)
    I_t = infections.I_t
    Z_t = infections.Z_t
    obs ~ as_turing_submodel(model.observation_model, y_t, I_t)
    # Uniform observation contract: sampled series and pre-error expected.
    generated_y_t = obs.y_t
    expected_y_t = obs.expected
    return merge(
        (; generated_y_t, expected_y_t, I_t, Z_t), _seed_quantity(infections)
    )
end

# The seeding window, passed on when the infection model exposes one.
# The names are a type parameter, so the branch is resolved at compile time and
# an infection model without a seeding window returns the tuple it always did.
function _seed_quantity(infections::NamedTuple{names}) where {names}
    return :I_seed in names ? (; I_seed = infections.I_seed) : (;)
end

function as_turing_model(model::IDModel, y_t, n)
    _check_observation_count(model, y_t, n)
    return _as_turing_model_idmodel(
        model, concrete_observations(y_t), _series_shape(model, n)
    )
end

# The shape the infection process is built at.
# `n` is the number of observations, so the series carries the chain's lead-in on
# top of it.
# The delays consume the head of the convolution, leaving exactly the `n`
# expected values the error model scores.
_series_shape(model, n::Int) = n + _series_lead_in(model)
function _series_shape(model, n::Dims{2})
    return (n[1], n[2] + _series_lead_in(model))
end

@doc raw"
Convenience 2-argument form: read the infection process's shape from the data.

The observation model and the data together fix the shape of the infection
process, so nothing about it needs to be supplied explicitly or stored on the
model. `as_turing_model(model, y_t)` is `as_turing_model(model, y_t, shape)`
with `shape` read from `y_t`. That is the data's time length, and for
stratified data the number of infection strata the observation model consumes,
resolved via [`infection_strata`](@ref).

The data may be a plain vector of observations, a `strata x time` matrix, or a
`NamedTuple` of per-stream series. Three age strata observed as one
hospitalisation stream is `Split(NegativeBinomialError(), [1.0 1.0 1.0])`, and
a `1 x T` data matrix then builds a 3-stratum infection process.

A scalar `missing` has no length to read, so simulating from the prior at a
chosen length is either `as_turing_model(model, missing, n)` or the two-argument
form over a blank series of that length.

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
function as_turing_model(model::IDModel, y_t)
    return as_turing_model(model, y_t, _data_shape(model, y_t))
end

_data_shape(model::IDModel, y_t) = _obs_data_shape(
    model.observation_model, y_t, _series_time_length(y_t)
)

# The two-argument `as_turing_model` is public, so a value with no time axis is
# refused by name rather than by a `MethodError` from a private helper.
# A `NamedTuple` of streams shares one time length, read off its first stream.
function _series_time_length(y_t)
    throw(
        ArgumentError(
            "cannot read a number of observations from data of type " *
                "$(typeof(y_t)); pass it explicitly, as " *
                "`as_turing_model(model, y_t, n)`"
        )
    )
end

function _series_time_length(::Missing)
    throw(
        ArgumentError(
            "a scalar `missing` carries no length to read a model shape from; \
            pass the number of observations, as `as_turing_model(model, missing, n)`, \
            or simulate over a blank series such as `Vector{Missing}(missing, n)`"
        )
    )
end
_series_time_length(y_t::AbstractVector) = length(y_t)
_series_time_length(y_t::AbstractMatrix) = size(y_t, 2)
_series_time_length(y_t::NamedTuple) = _series_time_length(first(y_t))

# The infection process's `ModelShape` implied by an observation model and a
# data value, shared by the two-argument `as_turing_model`, `IDProblem` and
# `forecast` so all three read a data-driven model's shape the same way.
# A scalar `missing` and a blank series have no stream axis to read, so the
# shape falls back to the observation model alone.
_obs_data_shape(obs, y_t, time_steps) = time_steps
_obs_data_shape(obs, y_t::Missing, time_steps) = _obs_data_shape_missing(
    obs, time_steps
)
function _obs_data_shape(obs, y_t::AbstractVector{Missing}, time_steps)
    return _obs_data_shape_missing(obs, time_steps)
end
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
