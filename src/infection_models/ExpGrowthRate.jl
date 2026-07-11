# Exponential-growth-rate infection process model.

# `exp(y)` written through `LogExpFunctions.xexpy` to match the upstream
# numerics used by `ExpGrowthRate`. This is the default `transformation`.
_oneexpy(y::T) where {T} = xexpy(one(T), y)

@doc raw"
Model unobserved infections via an internally generated time-varying exponential
growth rate.

```math
r_t \sim \text{latent}, \qquad I_t = g\!\left(\hat I_0 + \sum_{s \le t} r_s\right)
```

where the latent model `rt` supplies the (log) growth rates ``r_s``, ``g`` is
`transformation`, and the unconstrained initial infections ``\hat I_0`` come from
the prior in `initialisation`. The growth-rate process is generated *inside* the
model, so `as_turing_model` takes only the series length `n` and returns the
named tuple `(; I_t, Z_t)` with `Z_t` the growth-rate path.

This model carries no generation interval — it never uses one — so it takes a
`transformation` directly instead of an [`IDData`](@ref) object.

## Fields

  - `rt`: the latent process model (an [`AbstractLatentModel`](@ref)) generating
    the growth-rate path.
  - `transformation`: the link mapping the unconstrained cumulative sum to
    non-negative infections (default: numerically equivalent to `exp`,
    implemented via `LogExpFunctions.xexpy` for numerical stability).
  - `initialisation`: prior for the unconstrained initial infections (an
    [`AbstractPriorModel`](@ref); a bare `Distribution` is coerced via
    [`as_prior`](@ref)).

# Examples
```@example ExpGrowthRate
using ComposableTuringIDModels, Distributions
egr = ExpGrowthRate(; rt = RandomWalk(), initialisation = Normal())
rand(as_turing_model(egr, 10))
```
"
struct ExpGrowthRate{L <: AbstractLatentModel, F <: Function,
    S <: AbstractPriorModel} <: AbstractInfectionModel
    "Latent process model generating the growth-rate path."
    rt::L
    "Link mapping the unconstrained cumulative sum to non-negative infections."
    transformation::F
    "Prior for the unconstrained initial infections."
    initialisation::S
end

function ExpGrowthRate(; rt = RandomWalk(),
        transformation::Function = _oneexpy, initialisation = Normal())
    return ExpGrowthRate(
        as_prior(rt), transformation, as_prior(initialisation))
end

@model function as_turing_model(model::ExpGrowthRate, n)
    Z_t ~ to_submodel(as_turing_model(model.rt, n))
    init_incidence ~ to_submodel(
        as_turing_model(model.initialisation, 1))
    I_t = model.transformation.(only(init_incidence) .+ cumsum(Z_t))
    return (; I_t, Z_t)
end
