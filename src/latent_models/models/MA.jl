# Moving-average (MA) latent process model. Its accumulation step (`MAStep`)
# lives in `src/steps/`.

@doc raw"
A moving-average MA(`q`) latent process.

```math
Z_t = \epsilon_t + \sum_{i=1}^{q} \theta_i \epsilon_{t-i}
```

with coefficients ``\theta`` from the prior in `θ` and innovations from the
error model `ϵ_t`. The order `q` is the length of the coefficient prior.

# Examples
```@example MA
using EpiAwarePrototype, Distributions
ma = MA()
mdl = as_turing_model(ma, 10)
rand(mdl)
```
"
struct MA{C <: Sampleable, Q <: Int, E <: AbstractLatentModel} <:
       AbstractLatentModel
    "Prior distribution for the MA coefficients."
    θ::C
    "Order of the MA model."
    q::Q
    "Error model for the innovations."
    ϵ_t::E

    function MA(θ::Sampleable, q::Int, ϵ_t::AbstractLatentModel)
        @assert q>0 "q must be greater than 0"
        @assert q==length(θ) "q must equal the length of θ"
        new{typeof(θ), typeof(q), typeof(ϵ_t)}(θ, q, ϵ_t)
    end
end

function MA(θ::Distribution; q::Int = 1, ϵ_t::AbstractLatentModel = HierarchicalNormal())
    return MA(; θ_priors = fill(θ, q), ϵ_t = ϵ_t)
end

function MA(; θ_priors::Vector{C} = [truncated(Normal(0.0, 0.05), -1, 1)],
        ϵ_t::AbstractLatentModel = HierarchicalNormal()) where {C <: Distribution}
    q = length(θ_priors)
    return MA(_expand_dist(θ_priors), q, ϵ_t)
end

@model function as_turing_model(model::MA, n)
    q = model.q
    @assert n>q "n must be longer than the order of the moving average process"
    θ ~ model.θ
    ϵ_t ~ to_submodel(as_turing_model(model.ϵ_t, n), false)
    ma = accumulate_scan(MAStep(θ), (; val = 0, state = ϵ_t[1:q]), ϵ_t[(q + 1):end])
    return ma
end
