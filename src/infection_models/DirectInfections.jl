# Direct-infections process model.

@doc raw"
Model unobserved infections as a direct transformation of an internally
generated latent path.

```math
Z_t \sim \text{latent}, \qquad I_t = g\!\left(\hat I_0 + Z_t\right)
```

where the latent model `Z` supplies ``Z_t``, ``g`` is `transformation`, and the
unconstrained initial infections ``\hat I_0`` are drawn from the prior in
`initialisation`. The latent process is generated *inside* the model rather than
threaded in from outside, so `as_turing_model` takes only the series length `n`
and returns the named tuple `(; I_t, Z_t)`.

This model carries no generation interval — it never uses one — so it takes a
`transformation` directly ([`Renewal`](@ref) is the only infection model that
carries a generation interval).

## Fields

  - `Z`: the latent process model (an [`AbstractLatentModel`](@ref)) generating
    ``Z_t``. A length-`n` PATH slot: a bare `Distribution` here is auto-wrapped
    in an [`Intercept`](@ref), giving a constant path (one shared draw broadcast
    to length `n`); use [`IID`](@ref) for `n` independent draws.
  - `transformation`: the link mapping the unconstrained sum to non-negative
    infections (default `exp`).
  - `initialisation`: the prior for the unconstrained initial infections (a
    `Distribution` or prior model, sampled through [`as_turing_submodel`](@ref)).

# Examples
```@example DirectInfections
using ComposableTuringIDModels, Distributions
inf = DirectInfections(; Z = RandomWalk(), initialisation = Normal())
mdl = as_turing_model(inf, 10)
rand(mdl)
```
"
struct DirectInfections{L <: PriorLike, F <: Function, S <: PriorLike} <:
       AbstractInfectionModel
    "Latent process model generating ``Z_t``."
    Z::L
    "Link mapping the unconstrained sum to non-negative infections."
    transformation::F
    "Prior for the unconstrained initial infections."
    initialisation::S
end

function DirectInfections(; Z = RandomWalk(),
        transformation::Function = exp, initialisation = Normal())
    return DirectInfections(_path_prior(Z), transformation, initialisation)
end

@model function as_turing_model(model::DirectInfections, n)
    Z_t ~ as_turing_submodel(model.Z, n)
    init_incidence ~ as_turing_submodel(model.initialisation, 1; prefix = true)
    I_t = model.transformation.(only(init_incidence) .+ Z_t)
    return (; I_t, Z_t)
end
