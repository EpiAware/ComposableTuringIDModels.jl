# Non-centred hierarchical normal latent process model.

@doc raw"
A non-centred hierarchical normal latent process.

```math
\eta_t = \text{mean} + \sigma\, \epsilon_t, \quad \epsilon_t \sim
\mathrm{Normal}(0, 1), \quad \sigma \sim \text{std\_prior}
```

## Fields

  - `mean`: the mean of the normal process.
  - `std`: the prior for the standard deviation ``\sigma`` — a `Distribution`
    (a constant ``\sigma``, one scalar RV) or a process (a length-`n`, e.g.
    time-varying, scale). Drawn through the single [`as_turing_submodel`](@ref)
    seam and broadcast over the innovations, so a process makes the scale
    time-varying (stochastic volatility) with no other change.

# Examples
```@example HierarchicalNormal
using ComposableTuringIDModels, Distributions
hn = HierarchicalNormal()
mdl = as_turing_model(hn, 10)
rand(mdl)
```
"
# `mean` is added only when it is non-zero, which was a stored flag read off
# `mean` at construction.
# A stored flag goes stale the moment `mean` is set, so the test is made where
# it is used instead: a field that is not stored cannot disagree with the field
# it came from.
struct HierarchicalNormal{R <: Real, S <: PriorLike} <: AbstractLatentModel
    "Mean of the normal distribution."
    mean::R
    "Prior for the standard deviation."
    std::S
end

function HierarchicalNormal(;
        mean::Real = 0.0, std = truncated(Normal(0, 0.1), 0, Inf)
    )
    return HierarchicalNormal(mean, std)
end
HierarchicalNormal(std::PriorLike) = HierarchicalNormal(; std = std)

@model function as_turing_model(model::HierarchicalNormal, n::Int)
    # The scale is drawn through the single seam, so it is a scalar or a
    # length-`n` path.
    # Broadcasting `std .* ϵ_t` consumes both uniformly.
    std ~ as_turing_submodel(model.std, n; prefix = true)
    ϵ_t ~ as_turing_submodel(IID(Normal()), n)
    η_t = iszero(model.mean) ? std .* ϵ_t : model.mean .+ std .* ϵ_t
    return η_t
end
