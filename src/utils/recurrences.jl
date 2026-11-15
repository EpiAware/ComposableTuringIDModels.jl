# Helpers for the ComposableRecurrences calls the latent models share.

const _UNIT_LAG = [1.0]

# Collected because Enzyme cannot differentiate the `vcat` of a view, which a
# drawn vector often is, with the recurrence's output.
_history(x::Real) = [x]
_history(x::AbstractVector) = collect(x)

# The running total of `add` from `history`, returned after it.
function _cumulate(history, add)
    return vcat(history, Recurrence(_UNIT_LAG)(1.0; history, add))
end
