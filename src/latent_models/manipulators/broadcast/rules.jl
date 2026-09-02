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
    @assert length(latent) == period "length(latent) must equal period"
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
    @assert n <= period * length(latent) "n must be ≤ period * length(latent)"
    broadcast_latent = [latent[j] for j in 1:length(latent) for _ in 1:period]
    return broadcast_latent[1:n]
end

@doc raw"
Broadcast rule that places the latent process on a declared index window and
returns zero at every index outside it.

The inner process is generated over `length(window)` steps, so it costs nothing
outside the window and its own dynamics run only across it.

# Examples
```@example InWindow
using ComposableTuringIDModels
broadcast_rule(InWindow(2:3), [1.0, 2.0], 5, 1)
```

## Fields

  - `window`: the index window the process is placed on.
"
struct InWindow{R <: AbstractUnitRange{<:Integer}} <: AbstractBroadcastRule
    "The index window the process is placed on."
    window::R

    function InWindow(window::AbstractUnitRange{<:Integer})
        @assert !isempty(window) "window must contain at least one index"
        @assert first(window) >= 1 "window must start at index 1 or later"
        return new{typeof(window)}(window)
    end
end

# The window is what distinguishes one `InWindow` from another, so it belongs in
# the label.
_rule_label(rule::InWindow) = string(nameof(typeof(rule)), "(", rule.window, ")")

broadcast_n(rule::InWindow, n, period) = length(rule.window)

function broadcast_rule(rule::InWindow, latent, n, period)
    w = rule.window
    @assert last(w) <= n "window $w must lie within 1:$n"
    @assert length(latent) == length(w) "length(latent) must equal length(window)"
    offset = first(w) - 1
    z = zero(eltype(latent))
    return [i in w ? latent[i - offset] : z for i in 1:n]
end
