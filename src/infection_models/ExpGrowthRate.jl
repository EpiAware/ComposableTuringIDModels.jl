# Exponential-growth-rate infection process model.

# `exp(y)` written through `LogExpFunctions.xexpy` to match the upstream
# numerics used by `ExpGrowthRate`.
_oneexpy(y::T) where {T} = xexpy(one(T), y)

@doc raw"
Model unobserved infections via a time-varying exponential growth rate.

```math
I_t = g(\hat I_0) \exp\!\left(\sum_{s \le t} r_s\right)
```

where the latent path supplies the log growth rates ``r_s``, ``g`` is
`data.transformation`, and the unconstrained initial infections ``\hat I_0``
come from `initialisation_prior`.

# Arguments

  - `model`: the [`ExpGrowthRate`](@ref) model.
  - `rt`: the latent path of (log) growth rates.

# Examples
```@example ExpGrowthRate
using EpiAwarePrototype, Distributions
data = EpiData([0.2, 0.3, 0.5], exp)
egr = ExpGrowthRate(; data = data, initialisation_prior = Normal())
rand(as_turing_model(egr, randn(10) * 0.05))
```

## Fields

  - `data`: the [`EpiData`](@ref) object.
  - `initialisation_prior`: prior for the unconstrained initial infections.
"
@kwdef struct ExpGrowthRate{S <: Sampleable} <: AbstractEpiAwareModel
    "`EpiData` object."
    data::EpiData
    "Prior for the unconstrained initial infections."
    initialisation_prior::S = Normal()
end

@model function as_turing_model(model::ExpGrowthRate, rt)
    init_incidence ~ model.initialisation_prior
    return _oneexpy.(init_incidence .+ cumsum(rt))
end
