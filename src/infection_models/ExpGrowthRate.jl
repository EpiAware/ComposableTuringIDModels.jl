# Exponential-growth-rate infection process model.

# `exp(y)` written through `LogExpFunctions.xexpy` for numerical stability.
# This is the default `transformation`.
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
model, so `as_turing_model` takes a [`ModelShape`](@ref) `n` and returns the
named tuple `(; I_t, Z_t)` with `Z_t` the growth-rate path.

This model carries no generation interval — it never uses one — so it takes a
`transformation` directly ([`Renewal`](@ref) is the only infection model that
carries a generation interval).

Handing `rt` a [`Stratify`](@ref) (or [`Replicate`](@ref)) and calling
`as_turing_model(model, (n_strata, n_time))` gives one growth-rate series per
stratum, each cumulated (and seeded) independently. This model has no coupling
slot. The cumulative sum has no incidence window for a mixing operator to act
on. See [`Renewal`](@ref) for coupled strata.

## Fields

  - `rt`: the latent process model (an [`AbstractLatentModel`](@ref)) generating
    the growth-rate path. A length-`n` PATH slot: a bare `Distribution` here is
    auto-wrapped in an [`Intercept`](@ref), giving a constant path (one shared
    draw broadcast to length `n`); use [`IID`](@ref) for `n` independent draws,
    or [`Stratify`](@ref)/[`Replicate`](@ref) for a strata axis.
  - `transformation`: the link mapping the unconstrained cumulative sum to
    non-negative infections (default: numerically equivalent to `exp`,
    implemented via `LogExpFunctions.xexpy` for numerical stability).
  - `initialisation`: prior for the unconstrained initial infections (a
    `Distribution` or prior model, sampled through [`as_turing_submodel`](@ref)).

# Examples
```@example ExpGrowthRate
using ComposableTuringIDModels, Distributions
egr = ExpGrowthRate(; rt = RandomWalk(), initialisation = Normal())
rand(as_turing_model(egr, 10))
```

Stratified, each stratum an independent growth-rate path:

```@example ExpGrowthRate
strat = ExpGrowthRate(; rt = Replicate(RandomWalk()), initialisation = Normal())
size(as_turing_model(strat, (3, 10))().I_t)
```
"
struct ExpGrowthRate{L <: PriorLike, F <: Function, S <: PriorLike} <:
       AbstractInfectionModel
    "Latent process model generating the growth-rate path."
    rt::L
    "Link mapping the unconstrained cumulative sum to non-negative infections."
    transformation::F
    "Prior for the unconstrained initial infections."
    initialisation::S
end

function ExpGrowthRate(; rt = RandomWalk(),
        transformation::Function = _oneexpy, initialisation = Normal())
    return ExpGrowthRate(_path_prior(rt), transformation, initialisation)
end

# Cumulative sum along the time axis: the whole path for a single series, each
# stratum's row independently for a strata × time matrix.
_cumsum(Z_t::AbstractVector) = cumsum(Z_t)
_cumsum(Z_t::AbstractMatrix) = cumsum(Z_t; dims = 2)

@model function as_turing_model(model::ExpGrowthRate, n::ModelShape)
    Z_t ~ as_turing_submodel(model.rt, n)
    init_incidence ~ as_turing_submodel(
        model.initialisation, _n_strata(n); prefix = true)
    I_t = model.transformation.(_seed(init_incidence, n) .+ _cumsum(Z_t))
    return (; I_t, Z_t)
end
