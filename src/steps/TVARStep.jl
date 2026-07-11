# Time-varying autoregressive accumulation step (used by `TimeVaryingAR`).

@doc raw"
Time-varying AR(1) step for use with [`accumulate_scan`](@ref).

Unlike [`ARStep`](@ref), whose damping coefficient is fixed across the series,
this step reads a *per-step* coefficient from the driving sequence. Each element
of that sequence is a `(ρ_t, ϵ_t)` pair, and the step applies

```math
z_t = \rho_t\, z_{t-1} + \epsilon_t,
```

so the coefficient ``\rho_t`` varies over time. The state is the scalar previous
value ``z_{t-1}``; the default [`get_state`](@ref) prepends the seed ``z_1``.
"
struct TVARStep <: AbstractAccumulationStep end

(::TVARStep)(state, ρϵ) = ρϵ[1] * state + ρϵ[2]
