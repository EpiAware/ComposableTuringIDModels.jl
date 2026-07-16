# Time-varying AR(1): a thin constructor over `AR` with a `TimeVarying`-marked
# damping slot. The recursion lives in `AR`'s path-mode `as_turing_model` method
# (see `AR.jl`) and reuses the shared `TVARStep` (`src/steps/`); this file adds no
# new struct or recursion of its own.

@doc raw"
A first-order autoregressive process whose damping coefficient varies over time.

```math
z_t = \rho_t\, z_{t-1} + \epsilon_t, \qquad t = 2, \ldots, n,
```

with the initial value ``z_1`` from the prior in `init`, innovations from the
error model `ϵ_t`, and — unlike a constant [`AR`](@ref) — a whole coefficient
*path* ``\rho_t``.

`TimeVaryingAR` is not a distinct type: it is a convenience constructor that
returns an [`AR`](@ref) of order 1 whose `damp` slot is wrapped in a
[`TimeVarying`](@ref) marker. That marker is what tells `AR` to thread a per-step
coefficient path — drawn from the process in `damp` and mapped by `transform` —
through its recursion, instead of applying a constant coefficient. So a
time-varying AR is a genuine composition: `AR` plus a submodel on the coefficient.
A bare process in `AR`'s `damp` slot (without the marker) only enriches the *prior*
over a constant coefficient; the [`TimeVarying`](@ref) wrapper is the opt-in that
makes it a path.

Because the coefficient is a component you choose, any latent process — a
[`RandomWalk`](@ref) for a smooth path, an [`AR`](@ref) for a mean-reverting one —
drops into `damp`. `damp` is drawn on the unconstrained scale and mapped by
`transform` (default `tanh`, giving the stationary band ``(-1, 1)``; pass
`identity` to opt out) so an unbounded process such as a [`RandomWalk`](@ref) does
not produce an explosive recursion.

Identifiability note: a single series informs each ``\rho_t`` through one
transition, so recovering the whole path leans on the smoothness of the `damp`
prior; a panel of series sharing one ``\rho_t`` draw sharpens it.

# Keyword arguments

  - `damp`: prior process for the raw (pre-`transform`) coefficient path (a prior
    model or `Distribution`; default [`RandomWalk`](@ref)).
  - `init`: prior for the initial value ``z_1`` (a `Distribution` or prior model).
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
function TimeVaryingAR(; damp = RandomWalk(), init = Normal(),
        ϵ_t = HierarchicalNormal(), transform = tanh)
    return AR(; damp = TimeVarying(damp, transform), init = init, ϵ_t = ϵ_t)
end
