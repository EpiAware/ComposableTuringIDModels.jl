# Random walk accumulation step (used by `RandomWalk`).

@doc raw"
Random walk step for use with [`accumulate_scan`](@ref).
"
struct RWStep <: AbstractAccumulationStep end
(::RWStep)(state, ϵ) = state + ϵ
