# ARIMA(p, d, q) latent process combination.

@doc raw"
Build an ARIMA(p, d, q) latent process: an [`arma`](@ref) wrapped in a
`d`-fold [`DiffLatentModel`](@ref).

# Arguments

  - `ar_init`: prior(s) for the AR initial conditions.
  - `diff_init`: prior(s) for the differencing initial conditions (sets `d`).
  - `damp`: prior(s) for the AR damping coefficients.
  - `θ`: prior(s) for the MA coefficients.
  - `ϵ_t`: the innovation model (default [`HierarchicalNormal`](@ref)).

# Examples
```@example arima
using ComposableTuringIDModels, Distributions
model = arima()
rand(as_turing_model(model, 10))
```
"
function arima(; ar_init = [Normal()], diff_init = [Normal()],
        damp = [truncated(Normal(0.0, 0.05), 0, 1)],
        θ = [truncated(Normal(0.0, 0.05), -1, 1)], ϵ_t = HierarchicalNormal())
    arma_model = arma(; init = ar_init, damp = damp, θ = θ, ϵ_t = ϵ_t)
    return DiffLatentModel(; model = arma_model, init_priors = diff_init)
end
