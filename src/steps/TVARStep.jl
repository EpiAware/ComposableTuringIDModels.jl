# Order-1 autoregressive accumulation step (used by `AR` when `p == 1`).

@doc raw"
Order-1 AR step for use with [`accumulate_scan`](@ref).

The damping coefficient is read per step via [`_at`](@ref), so the *same* step
serves both a constant coefficient (a scalar `ρ`, drawn from a `Distribution`
prior) and a time-varying coefficient path (a vector `ρ`, drawn from a process
prior):

```math
z_t = \rho_t\, z_{t-1} + \epsilon_t, \qquad \rho_t = \texttt{\_at}(\rho, t).
```

Each element of the driving sequence is a `(t, ϵ_t)` pair; the state is the scalar
previous value ``z_{t-1}`` and the default [`get_state`](@ref) prepends the seed
``z_1``. A scalar `ρ` stays scalar (no per-step allocation), matching the
efficiency of a constant AR(1).
"
struct TVARStep{C} <: AbstractAccumulationStep
    ρ::C
end

(s::TVARStep)(state, tϵ) = _at(s.ρ, tϵ[1]) * state + tϵ[2]
