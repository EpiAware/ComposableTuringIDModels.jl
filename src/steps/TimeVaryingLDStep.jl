# Time-varying reporting-delay accumulation step (used by the `LatentDelay`
# observation modifier when the delay is time-varying).

@doc raw"
Time-varying LatentDelay step for use with [`accumulate_scan`](@ref).

The time-invariant [`LDStep`](@ref) holds one reversed delay kernel and applies it
at every step. `TimeVaryingLDStep` instead reads a fresh reversed kernel from the
scan input at each step, so the delay can vary with time:

```math
\mu_t = \langle \text{rev\_pmf}_t,\, \text{window}_t \rangle,
```

with `rev_pmf_t` supplied through the driving tuple `input = (ϵ, rev_pmf_t)`. The
window slides exactly as in [`LDStep`](@ref) (drop the oldest entry, append `ϵ`),
so all kernels must share the same length `d` (a constant convolution window). When
every kernel is identical this reduces exactly to [`LDStep`](@ref).
"
struct TimeVaryingLDStep <: AbstractAccumulationStep end

function (::TimeVaryingLDStep)(state, input)
    ϵ, rev_pmf_t = input
    val = dot(rev_pmf_t, state.current)
    current = vcat(state.current[2:end], ϵ)
    return (; val, current)
end

get_state(::TimeVaryingLDStep, initial_state, state) = state .|> x -> x.val
