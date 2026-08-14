# Autoregressive accumulation step (used by `AR`).

@doc raw"
Autoregressive step for use with [`accumulate_scan`](@ref).
"
struct ARStep{D <: AbstractVector{<:Real}} <: AbstractAccumulationStep
    damp_AR::D
end

function (ar::ARStep)(state, ϵ)
    new_val = dot(ar.damp_AR, state.window) + ϵ
    new_window = vcat(state.window[2:end], new_val)
    return (; val = new_val, window = new_window)
end

function get_state(::ARStep, initial_state, state)
    return vcat(collect(initial_state.window), state .|> x -> x.val)
end
