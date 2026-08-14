# Shape helpers: the two shapes a component's `n` argument can take, and the
# guard that stops a path model being asked for a strata axis it has none of.

@doc raw"
The shape a shape-aware component's `n` argument can take.

`n::Int` asks for a length-`n` path; `n::Dims{2}` asks for an
`(n_strata, n_time)` matrix, which is `size(I_t)` for a stratified infection
process. A bare `Distribution` is the one exception to the contract: it ignores
whichever shape it is asked for and draws a single scalar.

# Examples
```@example ModelShape
using ComposableTuringIDModels
ModelShape
```
"
const ModelShape = Union{Int, Dims{2}}

@doc raw"
The time-axis length implied by a [`ModelShape`](@ref): `n` itself for a
length-`n` path, the second dimension for an `(n_strata, n_time)` matrix.
"
_n_time(n::Int) = n
_n_time(n::Dims{2}) = n[2]

@doc raw"
The number of strata implied by a [`ModelShape`](@ref): `1` for a length-`n`
path, the first dimension for an `(n_strata, n_time)` matrix.
"
_n_strata(::Int) = 1
_n_strata(n::Dims{2}) = n[1]

@doc raw"
A path model has no strata axis of its own; ask for one explicitly rather than
guessing whether a shared path was meant to be pooled across strata.

Every [`AbstractPriorModel`](@ref) speaks the length-`n` path contract
(`as_turing_model(m, n::Int)`).
This guard rejects a `Dims{2}` shape for any model that has not opted into a
strata axis, so a stratum-shaped call to a plain path model fails at the point
it is made rather than returning something silently wrong.
Wrap the model in a `Stratify` (or a `Replicate`) to give it one.

# Examples
```@example ModelShape_guard
using ComposableTuringIDModels, Distributions
try
    as_turing_model(RandomWalk(), (3, 10))
catch e
    e
end
```
"
function as_turing_model(m::AbstractPriorModel, n::Dims{2})
    error(
        "$(nameof(typeof(m))) produces one path. Wrap it in a `Stratify` " *
            "to give it a strata axis."
    )
end
