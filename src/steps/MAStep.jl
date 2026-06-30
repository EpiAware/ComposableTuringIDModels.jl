# Moving-average accumulation step (used by `MA`).

@doc raw"
Moving-average step for use with [`accumulate_scan`](@ref).
"
struct MAStep{C <: AbstractVector{<:Real}} <: AbstractAccumulationStep
    θ::C
end

function (ma::MAStep)(state, ϵ)
    new_val = ϵ + dot(ma.θ, state.state)
    new_state = vcat(ϵ, state.state[1:(end - 1)])
    return (; val = new_val, state = new_state)
end

function get_state(acc_step::MAStep, initial_state, state)
    init_vals = initial_state.state
    new_vals = state .|> x -> x.val
    return vcat(init_vals, new_vals)
end
