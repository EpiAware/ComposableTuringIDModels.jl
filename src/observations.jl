# Observation models. An observation model maps a path of expected observations
# `Y_t` (e.g. infections) to observed counts `y_t` via a single `as_turing_model`
# method. `y_t === missing` triggers prior/predictive simulation; a concrete
# `y_t` conditions the model on data.

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

@doc raw"
A Poisson observation-error model.

# Examples
```jldoctest PoissonError; output = false
using EpiAwarePrototype
poi = PoissonError()
mdl = as_turing_model(poi, missing, fill(10, 10))
rand(mdl)
nothing
# output
```
"
struct PoissonError <: AbstractObservationErrorModel end

observation_error(::PoissonError, Y_t) = SafePoisson(Y_t)

@doc raw"
A negative-binomial observation-error model with an inferred cluster factor.

The field `cluster_factor_prior` sets the prior distribution for the cluster
factor that is sampled and used to parameterise the negative-binomial error.

# Examples
```@example NegativeBinomialError
using EpiAwarePrototype, Distributions
nb = NegativeBinomialError()
mdl = as_turing_model(nb, missing, fill(10, 10))
rand(mdl)
```
"
@kwdef struct NegativeBinomialError{S <: Sampleable} <: AbstractObservationErrorModel
    "Prior distribution for the cluster factor."
    cluster_factor_prior::S = HalfNormal(0.01)
end

@model function generate_observation_error_priors(
        obs_model::NegativeBinomialError, y_t, Y_t)
    cluster_factor ~ obs_model.cluster_factor_prior
    sq_cluster_factor = cluster_factor^2
    return (; sq_cluster_factor)
end

function observation_error(::NegativeBinomialError, Y_t, sq_cluster_factor)
    return NegativeBinomialMeanClust(Y_t, sq_cluster_factor)
end

@doc raw"
Apply a reporting delay to an underlying observation model.

The expected observations are convolved with the (reversed) delay PMF before
being passed to the wrapped `model`. `LatentDelay` shortens the expected
observation vector by the length of the delay PMF to avoid fitting to partially
observed data.

## Constructors
  - `LatentDelay(model, pmf)` — from a delay PMF (non-negative, sums to 1).
  - `LatentDelay(model, distribution; D, Δd)` — discretise a continuous delay
    distribution via [`censored_pmf`](@ref).

## Fields

  - `model`: the wrapped observation model the delayed expected observations are
    passed to.
  - `rev_pmf`: the reversed delay PMF convolved with the expected observations.

# Examples
```@example LatentDelay
using EpiAwarePrototype, Distributions
obs = LatentDelay(NegativeBinomialError(), truncated(Normal(5.0, 2.0), 0.0, Inf))
mdl = as_turing_model(obs, missing, fill(10, 30))
mdl()
```
"
struct LatentDelay{M <: AbstractEpiAwareModel, T <: AbstractVector{<:Real}} <:
       AbstractEpiAwareModel
    model::M
    rev_pmf::T

    function LatentDelay(model::M, pmf::T) where {
            M <: AbstractEpiAwareModel, T <: AbstractVector{<:Real}}
        @assert all(pmf .>= 0) "Delay PMF must be non-negative"
        @assert isapprox(sum(pmf), 1) "Delay PMF must sum to 1"
        rev_pmf = reverse(pmf)
        new{typeof(model), typeof(rev_pmf)}(model, rev_pmf)
    end
end

function LatentDelay(model::M, distribution::C; D = nothing,
        Δd = 1.0) where {
        M <: AbstractEpiAwareModel, C <: ContinuousDistribution}
    pmf = censored_pmf(distribution; Δd = Δd, D = D)
    return LatentDelay(model, pmf)
end

@model function as_turing_model(obs_model::LatentDelay, y_t, Y_t)
    if ismissing(y_t)
        y_t = Vector{Missing}(missing, length(Y_t))
    end

    pmf_length = length(obs_model.rev_pmf)
    @assert pmf_length<=length(Y_t) "The delay PMF must be no longer than the observation vector"

    expected_obs = accumulate_scan(
        LDStep(obs_model.rev_pmf),
        (; val = 0, current = Y_t[1:pmf_length]),
        vcat(Y_t[(pmf_length + 1):end], 0.0))

    y_t ~ to_submodel(as_turing_model(obs_model.model, y_t, expected_obs), false)
    return y_t
end

@doc raw"
LatentDelay step for use with [`accumulate_scan`](@ref).
"
struct LDStep{D <: AbstractVector{<:Real}} <: AbstractAccumulationStep
    rev_pmf::D
end

function (ld::LDStep)(state, ϵ)
    val = dot(ld.rev_pmf, state.current)
    current = vcat(state.current[2:end], ϵ)
    return (; val, current)
end

get_state(::LDStep, initial_state, state) = state .|> x -> x.val
