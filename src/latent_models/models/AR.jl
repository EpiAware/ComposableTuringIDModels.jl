# Autoregressive (AR) latent process model. Its constant-coefficient accumulation
# step (`ARStep`) and the time-varying step it reuses when its `damp` slot is
# wrapped in a `TimeVarying` marker (`TVARStep`) both live in `src/steps/`.

@doc raw"
An autoregressive AR(`p`) latent process.

```math
Z_t = \sum_{i=1}^{p} \rho_i Z_{t-i} + \epsilon_t
```

with damping coefficients ``\rho`` from the prior in `damp`, initial conditions
from the prior in `init`, and innovations from the error model `ϵ_t`. The order
`p` is fixed by the `damp` prior (a length-`k` vector ⇒ order `k`, a single
distribution ⇒ order 1); the `init` prior is sized to match.

Each prior slot takes a raw prior: pass a bare `Distribution` (order 1), a vector
of them (order = its length), or a richer prior model (a latent process supplying a
structured *prior* over the constant coefficient). Each slot is sampled through
[`as_turing_submodel`](@ref).

Wrapping the `damp` prior in a [`TimeVarying`](@ref) marker switches an AR(1) to a
genuinely time-varying coefficient: `AR(damp = TimeVarying(RandomWalk()))` threads a
per-step coefficient *path* ``\rho_t`` (drawn from the marked process and mapped by
its `transform`) through the AR(1) recursion, rather than a constant. The named
constructor [`TimeVaryingAR`](@ref) builds exactly this. Only order 1 supports a
time-varying coefficient; higher-order time-varying AR is tracked in #113.

# Examples
```@example AR
using ComposableTuringIDModels, Distributions
ar = AR()
mdl = as_turing_model(ar, 10)
rand(mdl)
```
"
struct AR{D <: Union{PriorLike, TimeVarying}, I <: PriorLike, P <: Int,
    E <: PriorLike} <: AbstractLatentModel
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

function AR(; damp = truncated(Normal(0.0, 0.05), 0, 1), init = Normal(),
        ϵ_t = HierarchicalNormal())
    # Order `p` is fixed by the damping prior (a length-`k` vector ⇒ order `k`, a
    # single distribution / process ⇒ order 1). The initial-conditions prior is
    # sized to match: a single distribution is sampled at length `p` (`filldist`),
    # while an explicitly-passed vector must already have length `p`.
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

@doc raw"
Time-varying AR(1): the same first-order recursion as [`AR`](@ref), but with the
constant damping coefficient replaced by a per-step path ``\rho_t`` drawn from the
[`TimeVarying`](@ref)-marked `damp` slot.

```math
z_t = \rho_t\, z_{t-1} + \epsilon_t, \qquad \rho_t = \texttt{transform}(u_t),
```

where the raw path ``u_t`` is a length-`(n-1)` draw from the marked process and
`transform` (default `tanh`) maps it into the stationary band. The recursion reuses
the shared [`TVARStep`](@ref) with [`accumulate_scan`](@ref), so time-varying AR is
expressed through AR's own machinery plus a submodel on the coefficient rather than
a separate process. The coefficient path is tracked as a generated quantity `ρ`
(via `:=`) so it is recovered from the sampled chain, while the model still returns
the numeric length-`n` path (staying a drop-in latent).
"
@model function as_turing_model(model::AR{<:TimeVarying}, n)
    @assert n>1 "n must be greater than 1 for a time-varying AR"
    ar_init ~ as_turing_submodel(model.init, 1; prefix = true)
    damp_AR ~ as_turing_submodel(model.damp.prior, n - 1; prefix = true)
    # Track the coefficient path as a generated quantity so it is recoverable from
    # the chain, while the model still returns the numeric path.
    ρ := model.damp.transform.(damp_AR)
    ϵ_t ~ as_turing_submodel(model.ϵ_t, n - 1)
    ar = accumulate_scan(TVARStep(), only(ar_init), collect(zip(ρ, ϵ_t)))
    return ar
end
