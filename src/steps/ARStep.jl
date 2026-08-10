# Autoregressive accumulation step (used by `AR`).

@doc raw"
Autoregressive step for use with [`accumulate_scan`](@ref).
"
struct ARStep{D <: AbstractVector{<:Real}} <: AbstractAccumulationStep
    damp_AR::D
end

function (ar::ARStep)(state, ϵ)
    new_val = _wsum(ar.damp_AR, state) + ϵ
    # `collect` the window: Enzyme's type analysis rejects the lazy view.
    return vcat(collect(state[2:end]), new_val)
end
