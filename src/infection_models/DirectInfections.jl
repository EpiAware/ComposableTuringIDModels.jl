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
threaded in from outside, so `as_turing_model` takes a [`ModelShape`](@ref) `n`
and returns the named tuple `(; I_t, Z_t)`.

This model carries no generation interval — it never uses one — so it takes a
`transformation` directly ([`Renewal`](@ref) is the only infection model that
carries a generation interval).

Handing `Z` a [`Stratify`](@ref) (or [`Replicate`](@ref)) and calling
`as_turing_model(model, (n_strata, n_time))` gives one direct-infections series
per stratum, each with its own seed (from a vector-valued `initialisation`, or
the same scalar seed broadcast to every stratum). This model has no coupling
slot. A plain transformation of the latent path has no incidence window for a
mixing operator to act on. See [`Renewal`](@ref) for coupled strata.

## Fields

  - `Z`: the latent process model (an [`AbstractLatentModel`](@ref)) generating
    ``Z_t``. A length-`n` PATH slot: a bare `Distribution` here is auto-wrapped
    in an [`Intercept`](@ref), giving a constant path (one shared draw broadcast
    to length `n`); use [`IID`](@ref) for `n` independent draws, or
    [`Stratify`](@ref)/[`Replicate`](@ref) for a strata axis.
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

Stratified, sharing a random walk with partially pooled per-stratum deviations:

```@example DirectInfections
strat = DirectInfections(; Z = Stratify(RandomWalk(), Hierarchy()),
    initialisation = Normal())
size(as_turing_model(strat, (3, 10))().I_t)
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

    function DirectInfections(Z, transformation::Function, initialisation)
        # `Z` is a length-`n` PATH slot: a bare `Distribution` is wrapped in an
        # `Intercept` (a constant path), never left as a scalar.
        wrapped = path_prior(Z)
        return new{typeof(wrapped), typeof(transformation), typeof(initialisation)}(
            wrapped, transformation, initialisation
        )
    end
end

function DirectInfections(;
        Z = RandomWalk(),
        transformation::Function = exp, initialisation = Normal()
    )
    return DirectInfections(Z, transformation, initialisation)
end

@model function as_turing_model(model::DirectInfections, n::ModelShape)
    Z_t ~ as_turing_submodel(model.Z, n)
    init_incidence ~ as_turing_submodel(
        model.initialisation, _n_strata(n); prefix = true
    )
    I_t = model.transformation.(_seed(init_incidence, n) .+ Z_t)
    return (; I_t, Z_t)
end
