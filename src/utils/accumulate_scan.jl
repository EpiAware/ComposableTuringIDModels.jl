# `accumulate_scan` and its default `get_state`.

@doc raw"
Apply an [`AbstractAccumulationStep`](@ref) across an input sequence in a single
pass.

This is an optimised `accumulate`-based replacement for an explicit `for` loop.
`acc_step` is a callable step `(state, ϵ) -> new_state`, `initial_state` seeds
the scan, and `ϵ_t` is the driving sequence. The returned value is assembled by
[`get_state`](@ref) from the accumulated states.

# Arguments

  - `acc_step`: an [`AbstractAccumulationStep`](@ref), a callable
    `(state, ϵ) -> new_state` applied at each element of the sequence.
  - `initial_state`: the seed state passed to `accumulate` as `init`.
  - `ϵ_t`: the driving sequence accumulated over.

# Examples
```@example accumulate_scan
using EpiAwarePrototype
accumulate_scan(EpiAwarePrototype.RWStep(), 0.0, [1.0, 2.0, 3.0])
```
"
function accumulate_scan(acc_step::AbstractAccumulationStep, initial_state, ϵ_t)
    result = accumulate(acc_step, ϵ_t; init = initial_state)
    return get_state(acc_step, initial_state, result)
end

@doc raw"
Assemble the final sequence from the raw output of [`accumulate_scan`](@ref).

The default method prepends `initial_state` to the last element of each
accumulated state. Step structs whose state is a named tuple (e.g. the `MA` and
`LatentDelay` steps) override this method to extract the relevant field.

# Arguments

  - `acc_step`: the [`AbstractAccumulationStep`](@ref) used in the scan; the
    method dispatches on its concrete type.
  - `initial_state`: the seed state used by [`accumulate_scan`](@ref).
  - `state`: the raw accumulated output produced by `accumulate`.

# Examples
```@example get_state
using EpiAwarePrototype
accumulate_scan(EpiAwarePrototype.RWStep(), 0.0, [1.0, 2.0, 3.0])
```
"
function get_state(acc_step::AbstractAccumulationStep, initial_state, state)
    return vcat(initial_state, last.(state))
end
