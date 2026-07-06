# Non-centred hierarchical normal latent process model.

@doc raw"
A non-centred hierarchical normal latent process.

```math
\eta_t = \text{mean} + \sigma\, \epsilon_t, \quad \epsilon_t \sim
\mathrm{Normal}(0, 1), \quad \sigma \sim \text{std\_prior}
```

## Fields

  - `mean`: the mean of the normal process.
  - `std`: the prior for the standard deviation ``\sigma`` (an
    [`AbstractPriorModel`](@ref); a bare `Distribution` is coerced via
    [`as_prior`](@ref)).
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
struct HierarchicalNormal{R <: Real, D <: AbstractPriorModel, M <: Bool} <:
       AbstractLatentModel
    "Mean of the normal distribution."
    mean::R
    "Prior for the standard deviation."
    std::D
    "Flag controlling whether `mean` is added (false when `mean == 0`)."
    add_mean::M
end

function HierarchicalNormal(; mean::Real = 0.0,
        std = truncated(Normal(0, 0.1), 0, Inf), add_mean::Bool = mean != 0)
    return HierarchicalNormal(mean, as_prior(std, :std), add_mean)
end
HierarchicalNormal(std::Distribution) = HierarchicalNormal(; std = std)
function HierarchicalNormal(mean::Real, std::Distribution)
    return HierarchicalNormal(; mean = mean, std = std)
end

@model function as_turing_model(model::HierarchicalNormal, n)
    std ~ to_submodel(as_turing_model(model.std, 1), false)
    ϵ_t ~ to_submodel(as_turing_model(IID(Normal()), n), false)
    η_t = model.add_mean ? model.mean .+ only(std) * ϵ_t : only(std) * ϵ_t
    return η_t
end
