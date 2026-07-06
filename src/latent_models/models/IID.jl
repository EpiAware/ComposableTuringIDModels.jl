# IID latent process model.

@doc raw"
Model a latent process ``\epsilon_t`` as independent, identically distributed
draws from `ϵ_t`.

```math
\epsilon_t \sim \text{Prior}, \quad t = 1, \ldots, n
```

# Examples
```@example IID
using ComposableTuringIDModels, Distributions
model = IID(Normal(0, 1))
mdl = as_turing_model(model, 10)
rand(mdl)
```
"
@kwdef struct IID{D <: Sampleable} <: AbstractLatentModel
    ϵ_t::D = Normal(0, 1)
end

@model function as_turing_model(model::IID, n)
    ϵ_t ~ filldist(model.ϵ_t, n)
    return ϵ_t
end
