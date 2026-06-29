# Shared observation-error-model machinery: the error supertype, the generic
# `as_turing_model` loop, and the `observation_error` /
# `generate_observation_error_priors` interface.

@doc raw"
Internal supertype shared by simple observation-error models (Poisson, negative
binomial).

It exists only so that the generic observation-error `as_turing_model` loop —
which is identical across error families — can be written once and dispatch the
family-specific pieces ([`observation_error`](@ref) and
[`generate_observation_error_priors`](@ref)) on the concrete type. It is a
subtype of [`AbstractEpiAwareModel`](@ref); the prototype keeps no deeper
hierarchy than this.
"
abstract type AbstractObservationErrorModel <: AbstractEpiAwareModel end

@doc raw"
Generate observations from an observation-error model.

Supports missing observations (`y_t === missing`, simulating predictively) and
expected-observation vectors `Y_t` shorter than `y_t` (the expected values are
aligned to the last `length(Y_t)` entries). Expected values are nudged by a tiny
constant to avoid degenerate error distributions.

The error family supplies [`generate_observation_error_priors`](@ref) (sampled
as a submodel) and [`observation_error`](@ref) (the per-time-point distribution).
"
@model function as_turing_model(obs_model::AbstractObservationErrorModel, y_t, Y_t)
    priors ~ to_submodel(
        generate_observation_error_priors(obs_model, y_t, Y_t), false)

    if ismissing(y_t)
        y_t = Vector{Missing}(missing, length(Y_t))
    end

    diff_t = length(y_t) - length(Y_t)
    @assert diff_t>=0 "The observation vector must be at least as long as the expected observation vector"

    pad_Y_t = Y_t .+ 1e-6
    for i in eachindex(Y_t)
        y_t[i + diff_t] ~ observation_error(obs_model, pad_Y_t[i], priors...)
    end
    return y_t
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
using EpiAwarePrototype
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
using EpiAwarePrototype
observation_error(PoissonError(), 10.0)
```
"
function observation_error end
