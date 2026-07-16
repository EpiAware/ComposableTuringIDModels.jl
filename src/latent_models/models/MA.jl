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
using ComposableTuringIDModels, Distributions
ma = MA()
mdl = as_turing_model(ma, 10)
rand(mdl)
```
"
struct MA{C <: PriorLike, Q <: Int, E <: PriorLike} <: AbstractLatentModel
    "Prior for the MA coefficients."
    θ::C
    "Order of the MA model."
    q::Q
    "Error model for the innovations."
    ϵ_t::E

    function MA(θ, q::Int, ϵ_t)
        @assert q>0 "q must be greater than 0"
        _assert_prior_length(θ, q, "θ")
        new{typeof(θ), typeof(q), typeof(ϵ_t)}(θ, q, ϵ_t)
    end
end

function MA(θ::Distribution; q::Int = 1, ϵ_t = HierarchicalNormal())
    return MA(; θ = fill(θ, q), ϵ_t = ϵ_t)
end

function MA(; θ = [truncated(Normal(0.0, 0.05), -1, 1)],
        ϵ_t = HierarchicalNormal())
    q = _prior_order(θ)
    return MA(θ, q, ϵ_t)
end

@model function as_turing_model(model::MA, n)
    q = model.q
    @assert n>q "n must be longer than the order of the moving average process"
    θ ~ as_turing_submodel(model.θ, q; prefix = true)
    ϵ_t ~ as_turing_submodel(model.ϵ_t, n)
    # `MAStep` keeps its innovation buffer newest-first (see `MAStep.jl`), so the
    # warm-up seed must be reversed: `reverse(ϵ_t[1:q]) = [ϵ_q, …, ϵ_1]` puts the
    # most recent innovation first, so `dot(θ, state)` pairs `θ[1]` with the most
    # recent innovation. Seeding oldest-first silently mis-ordered the first `q`
    # outputs for `q ≥ 2`. `get_state` reverses the seed back to natural order.
    ma = accumulate_scan(
        MAStep(θ), (; val = 0.0, state = reverse(ϵ_t[1:q])), ϵ_t[(q + 1):end])
    return ma
end
