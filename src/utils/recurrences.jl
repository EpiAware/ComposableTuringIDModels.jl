# Helpers for the ComposableRecurrences calls the latent and infection models
# share.

# A kernel of one unit weight on lag 1, so each value is the last plus its
# input: a running total.
const _UNIT_LAG = [1.0]

# A recurrence history: one past value is a length-1 window, and one past value
# per stratum is a `strata × 1` window.
# A vector is collected because a drawn parameter is often a view, and Enzyme
# reverse cannot differentiate the `vcat` a `prepend = true` call makes of a
# view (EpiAware/ComposableRecurrences.jl#299).
_history(x::Real) = [x]
_history(x::AbstractVector) = collect(x)
_strata_history(x::Real) = [x]
_strata_history(x::AbstractVector) = reshape(x, :, 1)

# The running total of `add` from `init`, returned after `init`.
function _cumulate(init, add)
    return Recurrence(_UNIT_LAG)(
        1.0; history = _history(init), add = add, prepend = true
    )
end
