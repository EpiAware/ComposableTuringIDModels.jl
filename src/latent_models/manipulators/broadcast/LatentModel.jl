# Broadcasting of a shorter latent process to length `n`: the broadcast-rule
# supertype, the generic `broadcast_n`/`broadcast_rule` interface, and the
# `BroadcastLatentModel` wrapper.

@doc raw"
Abstract supertype for broadcast rules used by [`BroadcastLatentModel`](@ref).
A rule defines [`broadcast_n`](@ref) (how long an inner series to generate) and
[`broadcast_rule`](@ref) (how to expand it to length `n`).
"
abstract type AbstractBroadcastRule end

@doc raw"
Length of the inner series an [`AbstractBroadcastRule`](@ref) needs to produce a
length-`n` broadcasted series. Each rule implements its own method.

# Arguments

  - `rule`: the [`AbstractBroadcastRule`](@ref).
  - `n`: the length of the broadcasted output series.
  - `period`: the broadcast period.

# Examples
```@example broadcast_n
using ComposableTuringIDModels
broadcast_n(RepeatEach(), 10, 7), broadcast_n(RepeatBlock(), 10, 7)
```
"
function broadcast_n end

@doc raw"
Expand an inner latent series to length `n` under an
[`AbstractBroadcastRule`](@ref). Each rule implements its own method.

# Arguments

  - `rule`: the [`AbstractBroadcastRule`](@ref).
  - `latent`: the inner latent series to expand.
  - `n`: the length of the broadcasted output series.
  - `period`: the broadcast period.

# Examples
```@example broadcast_rule
using ComposableTuringIDModels
broadcast_rule(RepeatEach(), [1, 2], 5, 2)
```
"
function broadcast_rule end

@doc raw"
Broadcast a shorter latent process to length `n` under a broadcast rule.

The inner model is generated over the length the rule requires
([`broadcast_n`](@ref)), then expanded to length `n` ([`broadcast_rule`](@ref)).

# Arguments

  - `model`: the [`BroadcastLatentModel`](@ref).
  - `n`: the length of the broadcasted series to generate.

# Examples
```@example BroadcastLatentModel
using ComposableTuringIDModels
each = BroadcastLatentModel(RandomWalk(), 7, RepeatEach())
rand(as_turing_model(each, 10))
```

## Fields

  - `model`: the underlying latent model.
  - `period`: the broadcast period.
  - `broadcast_rule`: the [`AbstractBroadcastRule`](@ref) applied.
"
struct BroadcastLatentModel{
    M <: AbstractLatentModel, P <: Integer, B <: AbstractBroadcastRule} <:
       AbstractLatentModel
    "The underlying latent model."
    model::M
    "The period of the broadcast."
    period::P
    "The broadcast rule applied."
    broadcast_rule::B

    function BroadcastLatentModel(model::M,
            period::Integer,
            broadcast_rule::B) where {
            M <: AbstractLatentModel, B <: AbstractBroadcastRule}
        @assert period>0 "period must be greater than 0"
        new{typeof(model), typeof(period), typeof(broadcast_rule)}(
            model, period, broadcast_rule)
    end
end

function BroadcastLatentModel(model::M; period::Integer,
        broadcast_rule::B) where {
        M <: AbstractLatentModel, B <: AbstractBroadcastRule}
    return BroadcastLatentModel(model, period, broadcast_rule)
end

@model function as_turing_model(model::BroadcastLatentModel, n)
    m = broadcast_n(model.broadcast_rule, n, model.period)
    latent_period ~ to_submodel(as_turing_model(model.model, m), false)
    return broadcast_rule(model.broadcast_rule, latent_period, n, model.period)
end
