# Latent process models. Each struct subtypes `AbstractEpiAwareModel` and
# implements a single `as_turing_model` method returning a `DynamicPPL.Model`
# whose generated quantity is a length-`n` latent path `Z_t`.

@doc raw"
Model a latent process ``\epsilon_t`` as independent, identically distributed
draws from `ϵ_t`.

```math
\epsilon_t \sim \text{Prior}, \quad t = 1, \ldots, n
```

# Examples
```@example IID
using EpiAwarePrototype, Distributions
model = IID(Normal(0, 1))
mdl = as_turing_model(model, 10)
rand(mdl)
```
"
@kwdef struct IID{D <: Sampleable} <: AbstractEpiAwareModel
    ϵ_t::D = Normal(0, 1)
end

@model function as_turing_model(model::IID, n)
    ϵ_t ~ filldist(model.ϵ_t, n)
    return ϵ_t
end

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
              AbstractEpiAwareModel
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

@doc raw"
Model the latent process ``Z_t`` as a random walk.

```math
Z_t = Z_0 + \sum_{i=1}^{t} \epsilon_i
```

where ``Z_0`` is drawn from `init_prior` and the increments ``\epsilon_i`` come
from the error model `ϵ_t` (a `HierarchicalNormal` by default, giving an
inferred step standard deviation).

# Examples
```@example RandomWalk
using EpiAwarePrototype, Distributions
rw = RandomWalk()
mdl = as_turing_model(rw, 10)
rand(mdl)
```
"
@kwdef struct RandomWalk{D <: Sampleable, E <: AbstractEpiAwareModel} <:
              AbstractEpiAwareModel
    init_prior::D = Normal()
    ϵ_t::E = HierarchicalNormal()
end

@model function as_turing_model(model::RandomWalk, n)
    @assert n>0 "n must be greater than 0"
    rw_init ~ model.init_prior
    ϵ_t ~ to_submodel(as_turing_model(model.ϵ_t, n - 1), false)
    rw = accumulate_scan(RWStep(), rw_init, ϵ_t)
    return rw
end

@doc raw"
Random walk step for use with [`accumulate_scan`](@ref).
"
struct RWStep <: AbstractAccumulationStep end
(::RWStep)(state, ϵ) = state + ϵ

@doc raw"
An autoregressive AR(`p`) latent process.

```math
Z_t = \sum_{i=1}^{p} \rho_i Z_{t-i} + \epsilon_t
```

with damping coefficients ``\rho`` from `damp_prior`, initial conditions from
`init_prior`, and innovations from the error model `ϵ_t`. The order `p` is the
length of the damping/initial priors.

# Examples
```@example AR
using EpiAwarePrototype, Distributions
ar = AR()
mdl = as_turing_model(ar, 10)
rand(mdl)
```
"
struct AR{D <: Sampleable, I <: Sampleable, P <: Int, E <: AbstractEpiAwareModel} <:
       AbstractEpiAwareModel
    "Prior distribution for the damping coefficients."
    damp_prior::D
    "Prior distribution for the initial conditions."
    init_prior::I
    "Order of the AR model."
    p::P
    "Error model for the innovations."
    ϵ_t::E

    function AR(damp_prior::Sampleable, init_prior::Sampleable, p::Int,
            ϵ_t::AbstractEpiAwareModel)
        @assert p>0 "p must be greater than 0"
        @assert p==length(damp_prior)==length(init_prior) "p must equal the length of damp_prior and init_prior"
        new{typeof(damp_prior), typeof(init_prior), typeof(p), typeof(ϵ_t)}(
            damp_prior, init_prior, p, ϵ_t)
    end
end

function AR(damp_prior::Sampleable, init_prior::Sampleable; p::Int = 1,
        ϵ_t::AbstractEpiAwareModel = HierarchicalNormal())
    return AR(; damp_priors = fill(damp_prior, p), init_priors = fill(init_prior, p),
        ϵ_t = ϵ_t)
end

function AR(; damp_priors::Vector{D} = [truncated(Normal(0.0, 0.05), 0, 1)],
        init_priors::Vector{I} = [Normal()],
        ϵ_t::AbstractEpiAwareModel = HierarchicalNormal()) where {
        D <: Sampleable, I <: Sampleable}
    p = length(damp_priors)
    return AR(_expand_dist(damp_priors), _expand_dist(init_priors), p, ϵ_t)
end

@model function as_turing_model(model::AR, n)
    p = model.p
    @assert n>p "n must be longer than the order of the autoregressive process"
    ar_init ~ model.init_prior
    damp_AR ~ model.damp_prior
    ϵ_t ~ to_submodel(as_turing_model(model.ϵ_t, n - p), false)
    ar = accumulate_scan(ARStep(damp_AR), ar_init, ϵ_t)
    return ar
end

@doc raw"
Autoregressive step for use with [`accumulate_scan`](@ref).
"
struct ARStep{D <: AbstractVector{<:Real}} <: AbstractAccumulationStep
    damp_AR::D
end

function (ar::ARStep)(state, ϵ)
    new_val = dot(ar.damp_AR, state) + ϵ
    return vcat(state[2:end], new_val)
end

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
struct MA{C <: Sampleable, Q <: Int, E <: AbstractEpiAwareModel} <:
       AbstractEpiAwareModel
    "Prior distribution for the MA coefficients."
    θ::C
    "Order of the MA model."
    q::Q
    "Error model for the innovations."
    ϵ_t::E

    function MA(θ::Sampleable, q::Int, ϵ_t::AbstractEpiAwareModel)
        @assert q>0 "q must be greater than 0"
        @assert q==length(θ) "q must equal the length of θ"
        new{typeof(θ), typeof(q), typeof(ϵ_t)}(θ, q, ϵ_t)
    end
end

function MA(θ::Distribution; q::Int = 1, ϵ_t::AbstractEpiAwareModel = HierarchicalNormal())
    return MA(; θ_priors = fill(θ, q), ϵ_t = ϵ_t)
end

function MA(; θ_priors::Vector{C} = [truncated(Normal(0.0, 0.05), -1, 1)],
        ϵ_t::AbstractEpiAwareModel = HierarchicalNormal()) where {C <: Distribution}
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

@doc raw"
Moving-average step for use with [`accumulate_scan`](@ref).
"
struct MAStep{C <: AbstractVector{<:Real}} <: AbstractAccumulationStep
    θ::C
end

function (ma::MAStep)(state, ϵ)
    new_val = ϵ + dot(ma.θ, state.state)
    new_state = vcat(ϵ, state.state[1:(end - 1)])
    return (; val = new_val, state = new_state)
end

function get_state(acc_step::MAStep, initial_state, state)
    init_vals = initial_state.state
    new_vals = state .|> x -> x.val
    return vcat(init_vals, new_vals)
end

@doc raw"
Broadcast a single sampled intercept value to a length-`n` latent process.

The field `intercept_prior` sets the prior distribution the intercept is drawn
from.

# Examples
```@example Intercept
using EpiAwarePrototype, Distributions
int = Intercept(Normal(0, 1))
mdl = as_turing_model(int, 10)
rand(mdl)
```
"
@kwdef struct Intercept{D <: Sampleable} <: AbstractEpiAwareModel
    "Prior distribution for the intercept."
    intercept_prior::D
end

@model function as_turing_model(model::Intercept, n)
    intercept ~ model.intercept_prior
    return fill(intercept, n)
end

@doc raw"
A fixed (non-sampled) intercept broadcast to a length-`n` latent process.
"
@kwdef struct FixedIntercept{F <: Real} <: AbstractEpiAwareModel
    intercept::F
end

@model function as_turing_model(model::FixedIntercept, n)
    return fill(model.intercept, n)
end

@doc raw"
A null latent model that generates `nothing` (no latent variables).

# Examples
```jldoctest Null
using EpiAwarePrototype
null = Null()
mdl = as_turing_model(null, 10)
isnothing(mdl())

# output

true
```
"
struct Null <: AbstractEpiAwareModel end

@model function as_turing_model(model::Null, n)
    return nothing
end

@doc raw"
Model a latent process as a `d`-fold differenced version of an inner process.

If ``\tilde Z_t`` is the inner (undifferenced) latent path supplied via `model`,
then

```math
\Delta^{(d)} Z_t = \tilde Z_t,
```

and ``Z_t`` is recovered by applying `cumsum` `d` times. The `d` initial terms
are inferred from `init_prior`; `d` equals the length of `init_priors`.

Composing `DiffLatentModel` over an `AR` gives an ARIMA-style latent process.

# Examples
```@example DiffLatentModel
using EpiAwarePrototype, Distributions
diff = DiffLatentModel(; model = RandomWalk(), init_priors = [Normal(), Normal()])
mdl = as_turing_model(diff, 10)
rand(mdl)
```
"
struct DiffLatentModel{M <: AbstractEpiAwareModel, P <: Distribution} <:
       AbstractEpiAwareModel
    "Underlying (undifferenced) latent model."
    model::M
    "Prior distribution for the initial latent variables."
    init_prior::P
    "Number of times differenced."
    d::Int

    function DiffLatentModel(model::AbstractEpiAwareModel, init_prior::Distribution, d::Int)
        @assert d>0 "d must be greater than 0"
        @assert d==length(init_prior) "d must equal the length of init_prior"
        new{typeof(model), typeof(init_prior)}(model, init_prior, d)
    end
end

function DiffLatentModel(model::AbstractEpiAwareModel, init_prior::Distribution; d::Int)
    return DiffLatentModel(; model = model, init_priors = fill(init_prior, d))
end

function DiffLatentModel(; model::AbstractEpiAwareModel,
        init_priors::Vector{D} where {D <: Distribution} = [Normal()])
    d = length(init_priors)
    return DiffLatentModel(model, _expand_dist(init_priors), d)
end

@model function as_turing_model(model::DiffLatentModel, n)
    d = model.d
    @assert n>d "n must be longer than d"
    latent_init ~ model.init_prior
    diff_latent ~ to_submodel(as_turing_model(model.model, n - d), false)
    return _combine_diff(latent_init, diff_latent, d)
end

function _combine_diff(init, diff, d)
    combined = vcat(init, diff)
    for _ in 1:d
        combined = cumsum(combined)
    end
    return combined
end
