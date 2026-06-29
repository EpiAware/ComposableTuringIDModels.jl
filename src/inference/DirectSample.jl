# Direct prior-sampling method.

@doc raw"
Direct sampling from a model's prior (no MCMC).

`apply_method(model, ::DirectSample)` samples the prior: with an integer
`n_samples` it draws that many times with `Turing.Prior()` (returning a chain),
and with `nothing` it draws once with `rand` (returning a `NamedTuple`).

## Fields

  - `n_samples`: number of prior draws, or `nothing` for a single `rand` draw.
"
@kwdef struct DirectSample <: AbstractEpiSamplingMethod
    "Number of prior draws, or `nothing` for a single `rand` draw."
    n_samples::Union{Int, Nothing} = nothing
end

function _apply_method(model::DynamicPPL.Model, method::DirectSample,
        prev_result = nothing; kwargs...)
    return _apply_direct_sample(model, method, method.n_samples; kwargs...)
end

function _apply_direct_sample(model, method, n_samples::Int; kwargs...)
    sample(
        model, Turing.Prior(), n_samples; kwargs...)
end
_apply_direct_sample(model, method, ::Nothing; kwargs...) = rand(model)
