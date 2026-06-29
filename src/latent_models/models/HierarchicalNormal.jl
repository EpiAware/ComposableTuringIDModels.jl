# Non-centred hierarchical normal latent process model.

@doc raw"
A non-centred hierarchical normal latent process.

```math
\eta_t = \text{mean} + \sigma\, \epsilon_t, \quad \epsilon_t \sim
\mathrm{Normal}(0, 1), \quad \sigma \sim \text{std\_prior}
```

## Fields

  - `mean`: the mean of the normal process.
  - `std_prior`: the prior distribution for the standard deviation ``\sigma``.
  - `add_mean`: flag controlling whether `mean` is added (false when
    `mean == 0`).

# Examples
```@example HierarchicalNormal
using EpiAwarePrototype, Distributions
hn = HierarchicalNormal()
mdl = as_turing_model(hn, 10)
rand(mdl)
```
"
@kwdef struct HierarchicalNormal{R <: Real, D <: Sampleable, M <: Bool} <:
              AbstractLatentModel
    "Mean of the normal distribution."
    mean::R = 0.0
    "Prior distribution for the standard deviation."
    std_prior::D = truncated(Normal(0, 0.1), 0, Inf)
    "Flag controlling whether `mean` is added (false when `mean == 0`)."
    add_mean::M = mean != 0
end

HierarchicalNormal(std_prior::Distribution) = HierarchicalNormal(; std_prior = std_prior)
function HierarchicalNormal(mean::Real, std_prior::Distribution)
    return HierarchicalNormal(mean, std_prior, mean != 0)
end

@model function as_turing_model(model::HierarchicalNormal, n)
    std ~ model.std_prior
    ϵ_t ~ to_submodel(as_turing_model(IID(Normal()), n), false)
    η_t = model.add_mean ? model.mean .+ std * ϵ_t : std * ϵ_t
    return η_t
end
