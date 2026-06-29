# Broadcast helpers: day-of-week and weekly latent processes.

@doc raw"
Build a [`BroadcastLatentModel`](@ref) for a day-of-week effect: a transformed
inner model repeated across a 7-day period.

# Arguments

  - `model`: the inner latent model.
  - `link`: link applied before broadcasting (default `x -> 7 * softmax(x)`,
    constraining the week effects to sum to 7).

# Examples
```@example broadcast_dayofweek
using EpiAwarePrototype
broadcast_dayofweek(RandomWalk())
```
"
function broadcast_dayofweek(model::AbstractLatentModel; link = x -> 7 * softmax(x))
    return BroadcastLatentModel(TransformLatentModel(model, link), 7, RepeatEach())
end

@doc raw"
Build a [`BroadcastLatentModel`](@ref) for a piecewise-constant weekly process.

# Arguments

  - `model`: the inner latent model.

# Examples
```@example broadcast_weekly
using EpiAwarePrototype
broadcast_weekly(RandomWalk())
```
"
function broadcast_weekly(model::AbstractLatentModel)
    return BroadcastLatentModel(model, 7, RepeatBlock())
end
