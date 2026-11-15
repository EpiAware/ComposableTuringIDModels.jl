# Helpers for the ComposableRecurrences calls the latent models make.

const _UNIT_LAG = [1.0]

# Collected because Enzyme cannot differentiate the `vcat` of a view, which a
# drawn vector often is, with the recurrence's output.
_history(x::Real) = [x]
_history(x::AbstractVector) = collect(x)
