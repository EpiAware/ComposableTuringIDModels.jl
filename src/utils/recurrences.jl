# Helpers for the ComposableRecurrences calls the latent and infection models
# share.

const _UNIT_LAG = [1.0]

# Collected because Enzyme cannot differentiate the `vcat` a `prepend = true`
# call makes of a view, which a drawn vector often is.
_history(x::Real) = [x]
_history(x::AbstractVector) = collect(x)
_strata_history(x::Real) = [x]
_strata_history(x::AbstractVector) = reshape(x, :, 1)

# The running total of `add` from `history`, returned after it when
# `prepend` is set.
function _cumulate(history, add; prepend = true)
    return Recurrence(_UNIT_LAG)(
        1.0; history = history, add = add, prepend = prepend
    )
end
