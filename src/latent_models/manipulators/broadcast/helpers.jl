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
using ComposableTuringIDModels
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
using ComposableTuringIDModels
broadcast_weekly(RandomWalk())
```
"
function broadcast_weekly(model::AbstractLatentModel)
    return BroadcastLatentModel(model, 7, RepeatBlock())
end

@doc raw"
Build a [`BroadcastLatentModel`](@ref) for an effect confined to an index window.

The result is a latent model that equals `model` on `window` and zero everywhere else, so summing it onto a base path with [`CombineLatentModels`](@ref) gives an effect that is active only over the window.

The `model` slot is a PATH slot, which is what lets one helper serve both a fixed and an estimated effect.
A bare `Distribution` is wrapped in an [`Intercept`](@ref), giving one estimated level held across the window.
A [`FixedIntercept`](@ref) gives a known constant.
A process gives an effect that varies within the window, generated over `length(window)` steps rather than over the whole series.

Windows touching either end of the series are ordinary cases here.

# Arguments

  - `model`: the latent model generating the effect inside the window.
  - `window`: the index window the effect applies on.

# Examples
```@example broadcast_window
using ComposableTuringIDModels, Distributions
base = RandomWalk(; init = Normal(0, 0.1), ϵ_t = IID(Normal(0, 0.05)))
windowed = CombineLatentModels([base, broadcast_window(Normal(1, 0.2), 20:30)])
rand(as_turing_model(windowed, 40))
```
"
function broadcast_window(model, window::AbstractUnitRange{<:Integer})
    # `period` is meaningless for `InWindow`, which reads its length from the
    # window, so it is fixed at the smallest value `BroadcastLatentModel` accepts
    # rather than exposed as an argument.
    return BroadcastLatentModel(model, 1, InWindow(window))
end
