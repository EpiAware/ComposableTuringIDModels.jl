# Epidemiological data container (generation interval + transformation).

@doc raw"
Epidemiological data shared by infection models: a discrete generation interval
and a transformation linking the unconstrained and constrained domains.

## Constructors

  - `IDData(gen_int, transformation)` — from a discrete generation interval
    vector (must be non-negative and sum to 1) and a transformation function.
  - `IDData(; gen_distribution, D_gen, Δd = 1.0, transformation = exp)` —
    discretise a continuous generation-interval distribution via
    double-interval censoring (CensoredDistributions.jl).

## Fields

  - `gen_int`: the discrete generation interval vector.
  - `len_gen_int`: the length of the discrete generation interval.
  - `transformation`: the transformation between unconstrained and constrained
    domains.

# Examples
```@example IDData
using ComposableTuringIDModels
data = IDData([0.2, 0.3, 0.5], exp)
```
"
struct IDData{T <: Real, F <: Function}
    "Discrete generation interval."
    gen_int::Vector{T}
    "Length of the discrete generation interval."
    len_gen_int::Integer
    "Transformation between unconstrained and constrained domains."
    transformation::F

    function IDData(gen_int, transformation::Function)
        @assert all(gen_int .>= 0) "Generation interval must be non-negative"
        @assert sum(gen_int)≈1 "Generation interval must sum to 1"
        new{eltype(gen_int), typeof(transformation)}(
            gen_int, length(gen_int), transformation)
    end
end

function IDData(; gen_distribution::ContinuousDistribution, D_gen = nothing,
        Δd = 1.0, transformation::Function = exp)
    # Drop the delay-0 bin (a generation interval has no mass at lag 0) and
    # renormalise, as the original EpiAware did.
    gen_int = _discretised_pmf(gen_distribution; Δd = Δd, D = D_gen) |>
              p -> p[2:end] ./ sum(p[2:end])
    return IDData(gen_int, transformation)
end
