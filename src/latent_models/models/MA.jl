# Moving-average (MA) latent process model. Its accumulation step (`MAStep`)
# lives in `src/steps/`.

@doc raw"
A moving-average MA(`q`) latent process.

```math
Z_t = \epsilon_t + \sum_{i=1}^{q} \theta_i \epsilon_{t-i}
```

with coefficients ``\theta`` from the prior in `θ` and innovations from the
error model `ϵ_t`. The order `q` is the length of the coefficient prior.

`ϵ_t` is a length-`n` path slot: a process gives time-varying innovations, while
a bare `Distribution` is auto-wrapped in an [`Intercept`](@ref), giving a
constant innovation path (one shared draw broadcast to every step). Use
[`IID`](@ref) for `n` independent innovations.

At order 1 the `θ` slot decides whether the coefficient is **constant or
time-varying**, through the same single-seam mechanism as [`AR`](@ref)'s damping:
`MA(θ = Normal(...))` is a constant coefficient (one scalar RV) while
`MA(θ = RandomWalk())` threads a per-step coefficient path. Higher-order (`q > 1`)
coefficients are constant.

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

    function MA(θ, _q, ϵ_t)
        # The order is read off `θ` rather than taken, which keeps every
        # construction path agreeing on it. See `rewrap`'s docstring.
        q = _slot_order(θ)
        @assert q > 0 "q must be greater than 0"
        # `ϵ_t` is a length-`n` PATH slot: a bare `Distribution` is wrapped in
        # an `Intercept` (a constant innovation path), never left as a scalar.
        wrapped = path_prior(ϵ_t)
        return new{typeof(θ), typeof(q), typeof(wrapped)}(θ, q, wrapped)
    end
end

function MA(θ::Distribution; q::Int = 1, ϵ_t = HierarchicalNormal())
    return MA(; θ = fill(θ, q), ϵ_t = ϵ_t)
end

function MA(;
        θ = [truncated(Normal(0.0, 0.05), -1, 1)],
        ϵ_t = HierarchicalNormal()
    )
    return MA(θ, nothing, ϵ_t)
end

@model function as_turing_model(model::MA, n::Int)
    q = model.q
    @assert n > q "n must be longer than the order of the moving average process"
    if q == 1
        # Order 1 draws the coefficient through the single seam.
        # MA is not recursive, so one broadcast over the lagged innovations
        # serves both a scalar and a path, giving `Z_t = ϵ_t + θ_{t-1} ϵ_{t-1}`.
        θ ~ as_turing_submodel(_order1_prior(model.θ), n - 1; prefix = true)
        ϵ_t ~ as_turing_submodel(model.ϵ_t, n)
        return vcat(ϵ_t[1], ϵ_t[2:n] .+ θ .* ϵ_t[1:(n - 1)])
    end
    θ ~ as_turing_submodel(model.θ, q; prefix = true)
    ϵ_t ~ as_turing_submodel(model.ϵ_t, n)
    # A causal convolution with kernel `[1; θ]` from lag 0, from time `q + 1`
    # on, reading the earlier innovations from the input itself.
    # The first `q` innovations pass through unchanged.
    ma = Convolution(vcat(one(eltype(θ)), θ))(ϵ_t; start = q + 1)
    return vcat(ϵ_t[1:q], ma)
end
