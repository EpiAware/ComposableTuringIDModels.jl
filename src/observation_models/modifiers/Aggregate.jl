# Aggregation observation modifier (sum expected observations over reporting
# windows).

# Scatter the predictions for the present time points back into a length-`n`
# vector, zero where not present.
# A `LatentDelay` shortens the window series whichever side of the aggregation it
# sits on, so there can be fewer predictions than present windows.
# They are the last windows, matching the right-alignment every other observation
# chain uses, so the leading windows are left at zero.
function _return_aggregate(pred_obs, idx, n)
    k = length(pred_obs)
    @assert k <= length(idx) "The aggregated model returned more predictions ($k) than reporting windows ($(length(idx)))"
    agg_obs = zeros(eltype(pred_obs), n)
    agg_obs[idx[(end - k + 1):end]] = pred_obs
    return agg_obs
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

Because the outermost modifier is applied first, the nesting of a
[`LatentDelay`](@ref) decides the units the delay is measured in.
`Aggregate(LatentDelay(model, pmf), aggregation)` sums into windows and then
convolves, so the delay is in *windows* and the leading windows it consumes go
unpredicted, while `LatentDelay(Aggregate(model, aggregation), pmf)` delays the
daily series before summing, so the delay is in *time points*.

Either way the delay leaves the head of the series uncovered. A window with no
expected values left to sum is dropped rather than scored against a zero it
never measured, so the counts reported for it stay out of the likelihood. A
window the delay only partially uncovers still has values to sum and is kept.

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
        present = aggregation .!= 0
        return new{M, A, typeof(present)}(model, aggregation, present)
    end
end

# The form a field-wise rebuild supplies, which the public constructor has no
# signature for because it takes fewer arguments than the struct has fields.
# `present` is derived here too rather than taken. See `rewrap`'s docstring.
Aggregate(model, aggregation, _present) = Aggregate(model, aggregation)

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
    # An outer modifier can hand over a `Y_t` shorter than `y_t`.
    # It is right-aligned like every other chain, so it covers the last
    # `length(Y_t)` time points and the window indices shift by `offset`.
    # A window reaching back before the start of `Y_t` is clipped, as one
    # reaching before time 1 already is.
    offset = n - length(Y_t)
    idx = findall(present)
    # A window ending at or before `offset` covers no expected value, so its sum
    # would be an exact zero meaning "nothing here" rather than "we expect zero",
    # which is a silent and arbitrarily strong likelihood contribution.
    # Those leading windows are dropped from the window series instead, so the
    # error model's right-alignment leaves their counts unscored.
    # A window only partially covered still has expected values to sum and is
    # kept.
    scored = filter(>(offset), idx)
    @assert !isempty(scored) "Every reporting window ends before the start of the expected observations, so there is nothing to score. Shorten the delay applied outside the aggregation, or lengthen the series."
    agg_Y_t = map(scored) do i
        sum(Y_t[max(1, i - aggregation[i] + 1 - offset):(i - offset)])
    end
    inner ~ as_turing_submodel(ag.model, y_t[present], agg_Y_t)
    # Scatter counts and means back into length-`n` vectors (zeros where absent).
    y_t = _return_aggregate(inner.y_t, idx, n)
    expected = _return_aggregate(inner.expected, idx, n)
    return (; y_t, expected)
end
