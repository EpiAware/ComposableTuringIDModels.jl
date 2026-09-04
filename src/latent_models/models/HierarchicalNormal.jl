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
  - `add_mean`: flag controlling whether `mean` is added, read off `mean` (false
    when `mean == 0`).

# Examples
```@example HierarchicalNormal
using ComposableTuringIDModels, Distributions
hn = HierarchicalNormal()
mdl = as_turing_model(hn, 10)
rand(mdl)
```
"
struct HierarchicalNormal{R <: Real, S <: PriorLike, M <: Bool} <:
    AbstractLatentModel
    "Mean of the normal distribution."
    mean::R
    "Prior for the standard deviation."
    std::S
    "Flag controlling whether `mean` is added, read off `mean`."
    add_mean::M
end

function HierarchicalNormal(;
        mean::Real = 0.0, std = truncated(Normal(0, 0.1), 0, Inf)
    )
    return HierarchicalNormal(mean, std, mean != 0)
end
HierarchicalNormal(std::PriorLike) = HierarchicalNormal(; std = std)
function HierarchicalNormal(mean::Real, std::PriorLike)
    return HierarchicalNormal(; mean = mean, std = std)
end

# `add_mean` says only whether `mean` is worth adding, so a rebuild reads it off
# `mean` rather than taking it and quietly dropping a new mean.
ConstructionBase.constructorof(::Type{<:HierarchicalNormal}) =
    _rebuild_hierarchical_normal

function _rebuild_hierarchical_normal(mean, std, _add_mean)
    return HierarchicalNormal(mean, std, mean != 0)
end

@model function as_turing_model(model::HierarchicalNormal, n::Int)
    # The scale is drawn through the single seam, so it is a scalar or a
    # length-`n` path.
    # Broadcasting `std .* ϵ_t` consumes both uniformly.
    std ~ as_turing_submodel(model.std, n; prefix = true)
    ϵ_t ~ as_turing_submodel(IID(Normal()), n)
    η_t = model.add_mean ? model.mean .+ std .* ϵ_t : std .* ϵ_t
    return η_t
end
