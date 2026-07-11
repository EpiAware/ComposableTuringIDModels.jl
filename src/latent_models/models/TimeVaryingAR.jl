# Time-varying AR(1) latent process. Its accumulation step (`TVARStep`) lives in
# `src/steps/`.

@doc raw"
A first-order autoregressive process whose damping coefficient varies over time.

```math
z_t = \rho_t\, z_{t-1} + \epsilon_t, \qquad t = 2, \ldots, n,
```

with the initial value ``z_1`` from the prior in `init`, innovations from the
error model `ϵ_t`, and — unlike [`AR`](@ref), whose coefficient is constant — a
whole coefficient *path* ``\rho_t``. That path is `transform` applied to a
length-`(n-1)` draw from the prior process `damp`.

This is the opt-in, genuinely time-varying counterpart of [`AR`](@ref): `AR`
applies its damping as a constant length-`p` coefficient (a process supplied to
its `damp` slot only enriches the *prior* over that constant), whereas
`TimeVaryingAR` threads a per-step ``\rho_t`` through the recursion. Because the
coefficient is a component you choose, any latent process — a [`RandomWalk`](@ref)
for a smooth path, an [`AR`](@ref) for a mean-reverting one — drops into `damp`.

`damp` is drawn on the unconstrained scale and mapped by `transform` (default
`tanh`, giving the stationary band ``(-1, 1)``; pass `identity` to opt out) so an
unbounded process such as a [`RandomWalk`](@ref) does not produce an explosive
recursion.

Identifiability note: a single series informs each ``\rho_t`` through one
transition, so recovering the whole path leans on the smoothness of the `damp`
prior; a panel of series sharing one ``\rho_t`` draw sharpens it.

## Fields

  - `damp`: prior process for the raw (pre-`transform`) coefficient path (an
    [`AbstractPriorModel`](@ref); default [`RandomWalk`](@ref)).
  - `init`: prior for the initial value ``z_1`` (an [`AbstractPriorModel`](@ref)).
  - `ϵ_t`: error model for the innovations (default [`HierarchicalNormal`](@ref)).
  - `transform`: map from the raw path to the coefficient (default `tanh`).

Built with `as_turing_model(m, n)` it returns the numeric length-`n` path (like
every other latent model, so it drops straight into any latent slot — e.g. a
`Renewal`'s `rt` or a `DirectInfections`'s `Z`). The coefficient path ``\rho_t`` is
tracked as a generated quantity `ρ` (via `:=`), so it is recovered from the sampled
chain (`group(chain, :ρ)`) without changing the return contract.

# Examples
```@example TimeVaryingAR
using ComposableTuringIDModels, Distributions
tv = TimeVaryingAR()
length(as_turing_model(tv, 10)())
```
"
struct TimeVaryingAR{D <: AbstractPriorModel, I <: AbstractPriorModel,
    E <: AbstractLatentModel, F <: Function} <: AbstractLatentModel
    "Prior process for the raw (pre-`transform`) coefficient path."
    damp::D
    "Prior for the initial value ``z_1``."
    init::I
    "Error model for the innovations."
    ϵ_t::E
    "Map from the raw path to the coefficient (default `tanh`)."
    transform::F
end

function TimeVaryingAR(; damp = RandomWalk(), init = Normal(),
        ϵ_t::AbstractLatentModel = HierarchicalNormal(), transform = tanh)
    return TimeVaryingAR(
        as_prior(damp, :damp_tv), as_prior(init, :tvar_init), ϵ_t, transform)
end

@model function as_turing_model(model::TimeVaryingAR, n)
    @assert n>1 "n must be greater than 1"
    tvar_init ~ to_submodel(as_turing_model(model.init, 1), false)
    damp_tv ~ to_submodel(as_turing_model(model.damp, n - 1), false)
    # Track the coefficient path as a generated quantity so it is recoverable from
    # the chain, while the model still returns the numeric path (stays a drop-in
    # latent).
    ρ := model.transform.(damp_tv)
    ϵ_t ~ to_submodel(as_turing_model(model.ϵ_t, n - 1), false)
    z = accumulate_scan(TVARStep(), only(tvar_init), collect(zip(ρ, ϵ_t)))
    return z
end
