# Autoregressive (AR) latent process model. Its accumulation step (`ARStep`)
# lives in `src/steps/`.

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
using ComposableTuringIDModels, Distributions
ar = AR()
mdl = as_turing_model(ar, 10)
rand(mdl)
```
"
struct AR{D <: Sampleable, I <: Sampleable, P <: Int, E <: AbstractLatentModel} <:
       AbstractLatentModel
    "Prior distribution for the damping coefficients."
    damp_prior::D
    "Prior distribution for the initial conditions."
    init_prior::I
    "Order of the AR model."
    p::P
    "Error model for the innovations."
    ϵ_t::E

    function AR(damp_prior::Sampleable, init_prior::Sampleable, p::Int,
            ϵ_t::AbstractLatentModel)
        @assert p>0 "p must be greater than 0"
        @assert p==length(damp_prior)==length(init_prior) "p must equal the length of damp_prior and init_prior"
        new{typeof(damp_prior), typeof(init_prior), typeof(p), typeof(ϵ_t)}(
            damp_prior, init_prior, p, ϵ_t)
    end
end

function AR(damp_prior::Sampleable, init_prior::Sampleable; p::Int = 1,
        ϵ_t::AbstractLatentModel = HierarchicalNormal())
    return AR(; damp_priors = fill(damp_prior, p), init_priors = fill(init_prior, p),
        ϵ_t = ϵ_t)
end

function AR(; damp_priors::Vector{D} = [truncated(Normal(0.0, 0.05), 0, 1)],
        init_priors::Vector{I} = [Normal()],
        ϵ_t::AbstractLatentModel = HierarchicalNormal()) where {
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
    # `ARStep`'s state runs oldest→newest (`[Z_{t-p}, …, Z_{t-1}]`), so reverse
    # the damping coefficients so `damp_AR[i]` multiplies the lag-`i` term
    # `Z_{t-i}`, matching the documented recursion `Z_t = Σ ρ_i Z_{t-i}`. Without
    # the reversal `damp_AR[1]` was applied to the *longest* lag; identical for the
    # default i.i.d. priors, but wrong for heterogeneous per-lag priors.
    ar = accumulate_scan(ARStep(reverse(damp_AR)), ar_init, ϵ_t)
    return ar
end
