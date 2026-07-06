# Concrete broadcast rules: repeat-each and repeat-block.

@doc raw"
Broadcast rule that repeats the latent process at each position within a period
(e.g. a fixed day-of-week effect).

# Examples
```@example RepeatEach
using ComposableTuringIDModels
broadcast_rule(RepeatEach(), [1, 2], 10, 2)
```
"
struct RepeatEach <: AbstractBroadcastRule end

broadcast_n(::RepeatEach, n, period) = period

function broadcast_rule(::RepeatEach, latent, n, period)
    @assert length(latent)==period "length(latent) must equal period"
    broadcast_latent = repeat(latent; outer = ceil(Int, n / period))
    return broadcast_latent[1:n]
end

@doc raw"
Broadcast rule that repeats the latent process in blocks of length `period`
(e.g. a piecewise-constant weekly process).

# Examples
```@example RepeatBlock
using ComposableTuringIDModels
broadcast_rule(RepeatBlock(), [1, 2, 3, 4, 5], 10, 2)
```
"
struct RepeatBlock <: AbstractBroadcastRule end

broadcast_n(::RepeatBlock, n, period) = ceil(Int, n / period)

function broadcast_rule(::RepeatBlock, latent, n, period)
    @assert n<=period * length(latent) "n must be ≤ period * length(latent)"
    broadcast_latent = [latent[j] for j in 1:length(latent) for _ in 1:period]
    return broadcast_latent[1:n]
end
