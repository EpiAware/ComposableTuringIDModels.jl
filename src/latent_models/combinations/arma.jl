# ARMA(p, q) latent process combination.

@doc raw"
Build an ARMA(p, q) latent process: an [`AR`](@ref) whose innovation model is an
[`MA`](@ref).

# Arguments

  - `init`: prior(s) for the AR initial conditions.
  - `damp`: prior(s) for the AR damping coefficients.
  - `θ`: prior(s) for the MA coefficients.
  - `ϵ_t`: the innovation model (default [`HierarchicalNormal`](@ref)).

# Examples
```@example arma
using EpiAwarePrototype, Distributions
model = arma(; θ = [truncated(Normal(0.0, 0.02), -1, 1)],
    damp = [truncated(Normal(0.0, 0.02), 0, 1)])
rand(as_turing_model(model, 10))
```
"
function arma(; init = [Normal()], damp = [truncated(Normal(0.0, 0.05), 0, 1)],
        θ = [truncated(Normal(0.0, 0.05), -1, 1)], ϵ_t = HierarchicalNormal())
    ma = MA(; θ_priors = θ, ϵ_t = ϵ_t)
    return AR(; damp_priors = damp, init_priors = init, ϵ_t = ma)
end
