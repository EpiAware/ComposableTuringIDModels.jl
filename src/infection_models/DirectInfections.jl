# Direct-infections process model.

@doc raw"
Model unobserved infections as a direct transformation of an internally
generated latent path.

```math
Z_t \sim \text{latent}, \qquad I_t = g\!\left(\hat I_0 + Z_t\right)
```

where the latent model `Z` supplies ``Z_t``, ``g`` is `transformation`, and the
unconstrained initial infections ``\hat I_0`` are drawn from
`initialisation_prior`. The latent process is generated *inside* the model rather
than threaded in from outside, so `as_turing_model` takes only the series length
`n` and returns the named tuple `(; I_t, Z_t)`.

This model carries no generation interval — it never uses one — so it takes a
`transformation` directly ([`Renewal`](@ref) is the only infection model that
carries a generation interval).

## Fields

  - `Z`: the latent process model (an [`AbstractLatentModel`](@ref)) generating
    ``Z_t``.
  - `transformation`: the link mapping the unconstrained sum to non-negative
    infections (default `exp`).
  - `initialisation_prior`: the prior for the unconstrained initial infections.

# Examples
```@example DirectInfections
using ComposableTuringIDModels, Distributions
inf = DirectInfections(; Z = RandomWalk(), initialisation_prior = Normal())
mdl = as_turing_model(inf, 10)
rand(mdl)
```
"
@kwdef struct DirectInfections{L <: AbstractLatentModel, F <: Function, S <: Sampleable} <:
              AbstractInfectionModel
    "Latent process model generating ``Z_t``."
    Z::L = RandomWalk()
    "Link mapping the unconstrained sum to non-negative infections."
    transformation::F = exp
    "Prior for the unconstrained initial infections."
    initialisation_prior::S = Normal()
end

@model function as_turing_model(model::DirectInfections, n)
    Z_t ~ to_submodel(as_turing_model(model.Z, n), false)
    init_incidence ~ model.initialisation_prior
    I_t = model.transformation.(init_incidence .+ Z_t)
    return (; I_t, Z_t)
end
