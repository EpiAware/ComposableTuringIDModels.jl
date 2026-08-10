# Backend-agnostic accumulation-step supertype used with `accumulate_scan`.

@doc raw"
Abstract supertype for accumulation step structs used with
[`accumulate_scan`](@ref).

A concrete `AbstractAccumulationStep` is a callable `(step)(state, ϵ)` returning
the next state. It is backend-agnostic: it contains no `Turing`/`DynamicPPL`
machinery and is reused unchanged across model components (`RandomWalk`, `AR`,
`MA`, `LatentDelay`).
"
abstract type AbstractAccumulationStep end
