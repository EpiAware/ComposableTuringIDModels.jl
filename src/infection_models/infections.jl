# Infection process models. An infection model maps a latent path `Z_t` to a
# path of unobserved infections `I_t` via a single `as_turing_model` method.

@doc raw"
Epidemiological data shared by infection models: a discrete generation interval
and a transformation linking the unconstrained and constrained domains.

## Constructors

  - `EpiData(gen_int, transformation)` — from a discrete generation interval
    vector (must be non-negative and sum to 1) and a transformation function.
  - `EpiData(; gen_distribution, D_gen, Δd = 1.0, transformation = exp)` —
    discretise a continuous generation-interval distribution via
    double-interval censoring (CensoredDistributions.jl).

## Fields

  - `gen_int`: the discrete generation interval vector.
  - `len_gen_int`: the length of the discrete generation interval.
  - `transformation`: the transformation between unconstrained and constrained
    domains.

# Examples
```@example EpiData
using EpiAwarePrototype
data = EpiData([0.2, 0.3, 0.5], exp)
```
"
struct EpiData{T <: Real, F <: Function}
    "Discrete generation interval."
    gen_int::Vector{T}
    "Length of the discrete generation interval."
    len_gen_int::Integer
    "Transformation between unconstrained and constrained domains."
    transformation::F

    function EpiData(gen_int, transformation::Function)
        @assert all(gen_int .>= 0) "Generation interval must be non-negative"
        @assert sum(gen_int)≈1 "Generation interval must sum to 1"
        new{eltype(gen_int), typeof(transformation)}(
            gen_int, length(gen_int), transformation)
    end
end

function EpiData(; gen_distribution::ContinuousDistribution, D_gen = nothing,
        Δd = 1.0, transformation::Function = exp)
    # Drop the delay-0 bin (a generation interval has no mass at lag 0) and
    # renormalise, as the original EpiAware did.
    gen_int = _discretised_pmf(gen_distribution; Δd = Δd, D = D_gen) |>
              p -> p[2:end] ./ sum(p[2:end])
    return EpiData(gen_int, transformation)
end

@doc raw"
Model unobserved infections as a direct transformation of a latent path.

```math
I_t = g\!\left(\hat I_0 + Z_t\right)
```

where ``g`` is `data.transformation` and the unconstrained initial infections
``\hat I_0`` are drawn from `initialisation_prior`.

## Fields

  - `data`: the [`EpiData`](@ref) object holding the generation interval and
    transformation.
  - `initialisation_prior`: the prior for the unconstrained initial infections.

# Examples
```@example DirectInfections
using EpiAwarePrototype, Distributions
data = EpiData([0.2, 0.3, 0.5], exp)
inf = DirectInfections(; data = data, initialisation_prior = Normal())
mdl = as_turing_model(inf, randn(10))
rand(mdl)
```
"
@kwdef struct DirectInfections{S <: Sampleable} <: AbstractEpiAwareModel
    "`EpiData` object."
    data::EpiData
    "Prior for the unconstrained initial infections."
    initialisation_prior::S = Normal()
end

@model function as_turing_model(model::DirectInfections, Z_t)
    init_incidence ~ model.initialisation_prior
    return model.data.transformation.(init_incidence .+ Z_t)
end
