# Shared observation-error-model machinery: the error supertype, the generic
# `as_turing_model` loop, and the `observation_error` /
# `generate_observation_error_priors` interface.

@doc raw"
Internal supertype shared by simple observation-error models (Poisson, negative
binomial).

It exists only so that the generic observation-error `as_turing_model` loop —
which is identical across error families — can be written once and dispatch the
family-specific pieces ([`observation_error`](@ref) and
[`generate_observation_error_priors`](@ref)) on the concrete type. It is the
error sub-role of [`AbstractObservationModel`](@ref); the prototype keeps no
deeper hierarchy than this.
"
abstract type AbstractObservationErrorModel <: AbstractObservationModel end

@doc raw"
Generate observations from an observation-error model.

Supports missing observations (`y_t === missing`, simulating predictively) and
expected-observation vectors `Y_t` shorter than `y_t` (the expected values are
aligned to the last `length(Y_t)` entries). Expected values are nudged by a tiny
constant to avoid degenerate error distributions.

The error family supplies [`generate_observation_error_priors`](@ref) (sampled
as a submodel) and [`observation_error`](@ref) (the per-time-point distribution).

Returns the uniform `(; y_t, expected)` tuple: `y_t` is the observed (or
simulated) counts and `expected` is the pre-error series. Exposing `expected`
lets a [`Split`](@ref) thread one stream's expectation into another.
"
@model function as_turing_model(obs_model::AbstractObservationErrorModel, y_t, Y_t)
    priors ~ to_submodel(
        generate_observation_error_priors(obs_model, y_t, Y_t), false)

    # Extract the count series scored by this model (plain vector, `missing`, or
    # a NamedTuple carrying extra data). Rebinding `y_t` keeps DynamicPPL treating
    # the entries as conditioned observations.
    y_t = define_y_t(obs_model, y_t, Y_t)

    diff_t = length(y_t) - length(Y_t)
    @assert diff_t>=0 "The observation vector must be at least as long as the expected observation vector"

    pad_Y_t = Y_t .+ 1e-6
    for i in eachindex(Y_t)
        y_t[i + diff_t] ~ observation_error(obs_model, pad_Y_t[i], priors...)
    end
    return (; y_t, expected = Y_t)
end

@doc raw"
Unpack the observed count series an observation-error model scores from the data
`y_t`, dispatching on the model type.

The default method covers every count family (Poisson, negative binomial) and the
Gaussian family: it accepts a plain observation vector, a `missing` (replaced by a
length-`Y_t` vector of `missing` for predictive simulation), or a `NamedTuple`
carrying the counts in a `y` field alongside any extra per-time-point data (a
model that needs more than the counts — e.g. [`BinomialError`](@ref), which also
needs the number of trials — reads those extra fields itself). This keeps the
simple case ergonomic (a plain vector just works) while letting a model opt into a
richer `NamedTuple` data contract.

# Arguments

  - `obs_model`: the observation-error model.
  - `y_t`: the observed data — a vector, `missing`, or a `NamedTuple`.
  - `Y_t`: the expected-observation series (used to size a `missing` series).

# Examples
```@example define_y_t
using ComposableTuringIDModels
# A plain vector passes through; a NamedTuple's `y` field is unpacked.
define_y_t(PoissonError(), [1, 2, 3], fill(10.0, 3)),
define_y_t(PoissonError(), (y = [1, 2, 3],), fill(10.0, 3))
```
"
function define_y_t(::AbstractObservationErrorModel, y_t, Y_t)
    # A NamedTuple carries the counts in its `y` field; a plain value is the
    # counts directly. Either way, a `missing` count series becomes a length-`Y_t`
    # vector of `missing` for predictive simulation.
    y = y_t isa NamedTuple ? y_t.y : y_t
    return ismissing(y) ? Vector{Missing}(missing, length(Y_t)) : y
end

@doc raw"
Generate the priors required by an observation-error model. Returns a named
tuple consumed by [`observation_error`](@ref). The default is an empty tuple.

# Arguments

  - `obs_model`: the observation-error model whose priors are generated.
  - `y_t`: the observed series (or `missing` when simulating predictively).
  - `Y_t`: the expected-observation series.

# Examples
```@example generate_observation_error_priors
using ComposableTuringIDModels
m = generate_observation_error_priors(NegativeBinomialError(), missing, fill(10.0, 5))
rand(m)
```
"
@model function generate_observation_error_priors(
        obs_model::AbstractObservationErrorModel, y_t, Y_t)
    return NamedTuple()
end

@doc raw"
The per-time-point observation-error distribution given an expected value and
the sampled priors. Each error family implements its own method.

# Arguments

  - `obs_model`: the observation-error model.
  - `Y_t`: the expected observation at a single time point.
  - additional positional arguments: any sampled priors produced by
    [`generate_observation_error_priors`](@ref) for the family (e.g. the squared
    cluster factor for [`NegativeBinomialError`](@ref)).

# Examples
```@example observation_error
using ComposableTuringIDModels
observation_error(PoissonError(), 10.0)
```
"
function observation_error end
