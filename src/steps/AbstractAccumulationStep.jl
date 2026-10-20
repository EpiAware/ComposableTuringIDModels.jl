# Backend-agnostic accumulation-step supertype used with `accumulate_scan`.

@doc raw"
Abstract supertype for accumulation step structs used with
[`accumulate_scan`](@ref).

A concrete `AbstractAccumulationStep` is a callable `(step)(state, ϵ)` returning
the next state. It is backend-agnostic: it contains no `Turing`/`DynamicPPL`
machinery.

The package's own recursions (the renewal, the latent processes and the
reporting delays) run on ComposableRecurrences.jl, which brings its own
gradient rules. `accumulate_scan` remains for a one-off recursion a model
author writes as a step.
"
abstract type AbstractAccumulationStep end
