# Autoregressive accumulation step (used by `AR`).

@doc raw"
Autoregressive step for use with [`accumulate_scan`](@ref).
"
struct ARStep{D <: AbstractVector{<:Real}} <: AbstractAccumulationStep
    damp_AR::D
end

function (ar::ARStep)(state, ϵ)
    # A `mapreduce` inner product, not `LinearAlgebra.dot`: `damp_AR`/`state`
    # are lag-order-length (order `p`, typically 2-3), far too short for BLAS
    # `ddot` to buy anything, and Enzyme falls back to a slow BLAS
    # replacement for it (warns "Using fallback BLAS replacements for
    # cblas_ddot64_"). `mapreduce` is a plain generic reduction every
    # backend, including Enzyme, differentiates directly with no BLAS call.
    new_val = mapreduce(*, +, ar.damp_AR, state) + ϵ
    return vcat(state[2:end], new_val)
end
