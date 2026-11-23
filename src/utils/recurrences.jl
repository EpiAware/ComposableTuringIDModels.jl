# Helpers for the ComposableRecurrences calls the latent and infection models
# share.

const _UNIT_LAG = [1.0]

# Collected because Enzyme cannot differentiate the `vcat` of a view, which a
# drawn vector often is, with the recurrence's output.
_history(x::Real) = [x]
_history(x::AbstractVector) = collect(x)
_strata_history(x::Real) = [x]
_strata_history(x::AbstractVector) = reshape(x, :, 1)

# The running total of `add` from `history`, a unit lag-1 `Recurrence`.
_cumulate(history, add) = Recurrence(_UNIT_LAG)(1.0; history = history, add = add)
