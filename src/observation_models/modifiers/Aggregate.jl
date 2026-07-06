# Aggregation observation modifier (sum expected observations over reporting
# windows).

# Scatter the predicted observations for the present time points back into a
# length-`n` vector of expected observations (zeros where not present).
function _return_aggregate(pred_obs, present, n)
    agg_obs = zeros(eltype(pred_obs), n)
    agg_obs[findall(present)] = pred_obs
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

# Arguments

  - `ag`: the [`Aggregate`](@ref) model.
  - `y_t`: the observed series (or `missing` when simulating predictively).
  - `Y_t`: the expected-observation series.

# Examples
```@example Aggregate
using EpiAwarePrototype
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
    P <: AbstractVector{<:Bool}} <: AbstractObservationModel
    "The underlying observation model."
    model::M
    "The per-period aggregation window lengths."
    aggregation::A
    "The boolean presence mask."
    present::P

    function Aggregate(model::M,
            aggregation::A) where {
            M <: AbstractObservationModel, A <: AbstractVector{<:Int}}
        present = aggregation .!= 0
        return new{M, A, typeof(present)}(model, aggregation, present)
    end
end

function Aggregate(; model::M, aggregation::A) where {
        M <: AbstractObservationModel, A <: AbstractVector{<:Int}}
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
    inner ~ to_submodel(
        as_turing_model(ag.model, y_t[present], agg_Y_t), false)
    # Scatter both the sampled observations and the expected means back into
    # length-`n` vectors (zeros where not present) so `Aggregate` conforms to the
    # uniform `(; y_t, expected)` contract and can thread through a `Split`.
    y_t = _return_aggregate(inner.y_t, present, n)
    expected = _return_aggregate(inner.expected, present, n)
    return (; y_t, expected)
end
