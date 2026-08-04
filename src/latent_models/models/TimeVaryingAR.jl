# Time-varying AR(1): a thin alias over [`AR`](@ref). Time-varying damping is no
# longer a bespoke struct — it is what an order-1 [`AR`](@ref) *is* when its `damp`
# slot holds a process rather than a `Distribution` (see `AR.jl` and the
# [time-varying-parameter mechanism](@ref as_turing_submodel)). This
# constructor is kept for discoverability and as the named entry point.

@doc raw"
A first-order autoregressive process whose damping coefficient varies over time.

```math
z_t = \rho_t\, z_{t-1} + \epsilon_t, \qquad \rho_t = \texttt{transform}(u_t),
```

`TimeVaryingAR(; damp, init, ϵ_t, transform)` is exactly
`AR(; damp, init, ϵ_t, transform)` with a process `damp` prior (default
[`RandomWalk`](@ref)): the coefficient is a length-`(n-1)` path drawn from `damp`
and mapped by `transform` (default `tanh`, the stationary band ``(-1, 1)``). There
is no separate type — `TimeVaryingAR()` and `AR(; damp = RandomWalk())` build the
same [`AR`](@ref), which returns the numeric length-`n` path and tracks the
coefficient path as the generated quantity `ρ`.

Because the damping is an ordinary prior slot, any latent process drops into it —
a [`RandomWalk`](@ref) for a smooth path, an [`AR`](@ref) for a mean-reverting one.

# Examples
```@example TimeVaryingAR
using ComposableTuringIDModels, Distributions
tv = TimeVaryingAR()
length(as_turing_model(tv, 10)())
```
"
function TimeVaryingAR(; damp = RandomWalk(), init = Normal(),
        ϵ_t = HierarchicalNormal(), transform = tanh)
    return AR(; damp = damp, init = init, ϵ_t = ϵ_t, transform = transform)
end
