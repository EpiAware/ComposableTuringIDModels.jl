# Direct-infections process model.

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
