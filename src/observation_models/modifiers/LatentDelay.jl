# Reporting-delay observation modifier (and its accumulation step).

@doc raw"
Apply a reporting delay to an underlying observation model.

The expected observations are convolved with the (reversed) delay PMF before
being passed to the wrapped `model`. `LatentDelay` shortens the expected
observation vector by the length of the delay PMF to avoid fitting to partially
observed data.

## Constructors
  - `LatentDelay(model, pmf)` — from a delay PMF (non-negative, sums to 1).
  - `LatentDelay(model, distribution; D, Δd)` — discretise a continuous delay
    distribution via double-interval censoring (CensoredDistributions.jl).

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
struct LatentDelay{M <: AbstractObservationModel, T <: AbstractVector{<:Real}} <:
       AbstractObservationModel
    model::M
    rev_pmf::T

    function LatentDelay(model::M, pmf::T) where {
            M <: AbstractObservationModel, T <: AbstractVector{<:Real}}
        @assert all(pmf .>= 0) "Delay PMF must be non-negative"
        @assert isapprox(sum(pmf), 1) "Delay PMF must sum to 1"
        rev_pmf = reverse(pmf)
        new{typeof(model), typeof(rev_pmf)}(model, rev_pmf)
    end
end

function LatentDelay(model::M, distribution::C; D = nothing,
        Δd = 1.0) where {
        M <: AbstractObservationModel, C <: ContinuousDistribution}
    pmf = _discretised_pmf(distribution; Δd = Δd, D = D)
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
