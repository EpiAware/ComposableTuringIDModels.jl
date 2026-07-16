# Autoregressive (AR) latent process model. Its accumulation step (`ARStep`)
# lives in `src/steps/`.

@doc raw"
An autoregressive AR(`p`) latent process.

```math
Z_t = \sum_{i=1}^{p} \rho_i Z_{t-i} + \epsilon_t
```

with damping coefficients ``\rho`` from the prior in `damp`, initial conditions
from the prior in `init`, and innovations from the error model `ϵ_t`. The order
`p` is the length of the damping/initial priors.

Each prior slot takes a raw prior: pass a bare `Distribution` (order 1), a vector
of them (order = its length), or a richer prior model (e.g. a latent process for a
time-varying coefficient). Each slot is sampled through
[`as_turing_submodel`](@ref).

# Examples
```@example AR
using ComposableTuringIDModels, Distributions
ar = AR()
mdl = as_turing_model(ar, 10)
rand(mdl)
```
"
struct AR{D <: PriorLike, I <: PriorLike, P <: Int, E <: PriorLike} <:
       AbstractLatentModel
    "Prior for the damping coefficients."
    damp::D
    "Prior for the initial conditions."
    init::I
    "Order of the AR model."
    p::P
    "Error model for the innovations."
    ϵ_t::E

    function AR(damp, init, p::Int, ϵ_t)
        @assert p>0 "p must be greater than 0"
        _assert_prior_length(damp, p, "damp")
        _assert_prior_length(init, p, "init")
        new{typeof(damp), typeof(init), typeof(p), typeof(ϵ_t)}(
            damp, init, p, ϵ_t)
    end
end

function AR(damp::Sampleable, init::Sampleable; p::Int = 1,
        ϵ_t = HierarchicalNormal())
    return AR(; damp = fill(damp, p), init = fill(init, p), ϵ_t = ϵ_t)
end

function AR(; damp = [truncated(Normal(0.0, 0.05), 0, 1)], init = [Normal()],
        ϵ_t = HierarchicalNormal())
    p = _prior_order(damp)
    return AR(damp, init, p, ϵ_t)
end

@model function as_turing_model(model::AR, n)
    p = model.p
    @assert n>p "n must be longer than the order of the autoregressive process"
    ar_init ~ as_turing_submodel(model.init, p; prefix = true)
    damp_AR ~ as_turing_submodel(model.damp, p; prefix = true)
    ϵ_t ~ as_turing_submodel(model.ϵ_t, n - p)
    # `ARStep`'s state runs oldest→newest (`[Z_{t-p}, …, Z_{t-1}]`), so reverse
    # the damping coefficients so `damp_AR[i]` multiplies the lag-`i` term
    # `Z_{t-i}`, matching the documented recursion `Z_t = Σ ρ_i Z_{t-i}`. Without
    # the reversal `damp_AR[1]` was applied to the *longest* lag; identical for the
    # default i.i.d. priors, but wrong for heterogeneous per-lag priors.
    ar = accumulate_scan(ARStep(reverse(damp_AR)), ar_init, ϵ_t)
    return ar
end
