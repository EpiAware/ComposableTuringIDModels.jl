# Aggregation observation modifier (sum expected observations over reporting
# windows).

# Scatter the predicted observations for the present time points back into a
# length-`n` vector of expected observations (zeros where not present).
function _return_aggregate(pred_obs, present, n)
    agg_obs = zeros(eltype(pred_obs), n)
    agg_obs[findall(present)] = pred_obs
    return agg_obs
end

# Find a `LatentDelay` nested inside `model` by following the linear `.model`
# wrapper chain (the shape every modifier in this package uses). Returns
# `nothing` when no `.model` field is present (a leaf, or a branching
# composition such as `Split` that this conservative check does not look
# inside).
_find_latent_delay(m::LatentDelay) = m
function _find_latent_delay(m)
    hasproperty(m, :model) || return nothing
    return _find_latent_delay(m.model)
end

# The length of the (reversed or forward) delay PMF a `LatentDelay.delay`
# field yields, when that length is knowable without sampling. A fixed PMF
# and a deterministic per-time PMF sequence carry their length directly; an
# `UncertainDelay`'s PMF length is fixed by its `D`/`Δd` truncation horizon
# regardless of the (uncertain) distribution parameters. Any other
# `AbstractPriorModel` delay has no statically knowable length, so `nothing`
# is returned and the aggregation-alignment check is skipped for it.
_delay_pmf_length(pmf::AbstractVector{<:Real}) = length(pmf)
_delay_pmf_length(pmfs::AbstractVector{<:AbstractVector{<:Real}}) = length(first(pmfs))
function _delay_pmf_length(delay::UncertainDelay)
    return length(0.0:(delay.Δd):(delay.D - delay.Δd))
end
_delay_pmf_length(::AbstractPriorModel) = nothing

# Reject a `LatentDelay` nested inside an `Aggregate` when its PMF is longer
# than one bin. `LatentDelay` always trims that many lead-in windows off the
# *aggregated* window series to avoid fitting to partially observed data, so
# the inner model returns fewer predictions than there are reported windows
# and the scatter in `_return_aggregate` fails with an opaque
# `DimensionMismatch`. Aligning a partially covered window is a modelling
# decision (what does a partial window mean?) that this guard deliberately
# does not make; it only rejects the combination up front. Move the delay
# outside the `Aggregate` (applied to the pre-aggregation series) instead.
function _check_aggregate_delay(model)
    delay = _find_latent_delay(model)
    delay === nothing && return nothing
    d = _delay_pmf_length(delay.delay)
    (d === nothing || d <= 1) && return nothing
    throw(
        ArgumentError(
            "Aggregate cannot build: it wraps a LatentDelay with a " *
                "$d-long delay PMF. LatentDelay drops $(d - 1) lead-in " *
                "window(s) from the aggregated series to avoid fitting " *
                "partially observed windows, so the inner model would " *
                "return fewer predictions than reported windows. Move the " *
                "LatentDelay outside the Aggregate (apply it to the " *
                "pre-aggregation series instead)."
        )
    )
end

@doc raw"
Aggregate the expected observations of an underlying model over reporting windows.

Each entry of `aggregation` gives the window length to sum over at the
corresponding (broadcast) time point, and `present` (derived as
`aggregation .!= 0`) marks the time points that are reported. The aggregation and
presence vectors are broadcast to the observation length with
[`RepeatEach`](@ref), the expected observations are summed over each window, the
inner `model` is applied to the present windows, and the predictions are scattered
back into a full-length vector (zeros elsewhere).

# Arguments

  - `ag`: the [`Aggregate`](@ref) model.
  - `y_t`: the observed series (or `missing` when simulating predictively).
  - `Y_t`: the expected-observation series.

# Examples
```@example Aggregate
using ComposableTuringIDModels
obs = Aggregate(PoissonError(), [0, 0, 0, 0, 0, 0, 7])
mdl = as_turing_model(obs, missing, fill(10.0, 14))
rand(mdl)
```

## Fields

  - `model`: the underlying observation model applied to the aggregated windows.
  - `aggregation`: the per-period window lengths (`0` marks an unreported point).
  - `present`: the boolean presence mask (`aggregation .!= 0`).
"
struct Aggregate{
        M <: AbstractObservationModel, A <: AbstractVector{<:Int},
        P <: AbstractVector{<:Bool},
    } <: AbstractObservationModel
    "The underlying observation model."
    model::M
    "The per-period aggregation window lengths."
    aggregation::A
    "The boolean presence mask."
    present::P

    function Aggregate(
            model::M,
            aggregation::A
        ) where {
            M <: AbstractObservationModel, A <: AbstractVector{<:Int},
        }
        _check_aggregate_delay(model)
        present = aggregation .!= 0
        return new{M, A, typeof(present)}(model, aggregation, present)
    end
end

function Aggregate(; model::M, aggregation::A) where {
        M <: AbstractObservationModel, A <: AbstractVector{<:Int},
    }
    return Aggregate(model, aggregation)
end

@model function as_turing_model(ag::Aggregate, y_t, Y_t)
    if ismissing(y_t)
        y_t = Vector{Missing}(missing, length(Y_t))
    end
    n = length(y_t)
    m = length(ag.aggregation)
    aggregation = broadcast_rule(RepeatEach(), ag.aggregation, n, m)
    present = broadcast_rule(RepeatEach(), ag.present, n, m)
    agg_Y_t = map(findall(present)) do i
        sum(Y_t[max(1, i - aggregation[i] + 1):i])
    end
    inner ~ as_turing_submodel(ag.model, y_t[present], agg_Y_t)
    # Scatter counts and means back into length-`n` vectors (zeros where absent).
    y_t = _return_aggregate(inner.y_t, present, n)
    expected = _return_aggregate(inner.expected, present, n)
    return (; y_t, expected)
end
