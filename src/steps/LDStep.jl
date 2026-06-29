# Reporting-delay accumulation step (used by the `LatentDelay` observation
# modifier).

@doc raw"
LatentDelay step for use with [`accumulate_scan`](@ref).
"
struct LDStep{D <: AbstractVector{<:Real}} <: AbstractAccumulationStep
    rev_pmf::D
end

function (ld::LDStep)(state, ϵ)
    val = dot(ld.rev_pmf, state.current)
    current = vcat(state.current[2:end], ϵ)
    return (; val, current)
end

get_state(::LDStep, initial_state, state) = state .|> x -> x.val
