# Autoregressive (AR) latent process model. Its constant higher-order step
# (`ARStep`) and the order-1 step it uses for a scalar or time-varying coefficient
# (`TVARStep`) both live in `src/steps/`.

# Default coefficient transform: a process `damp` prior is unbounded, so `tanh`
# maps its draws into the stationary band ``(-1, 1)``; a `Distribution` prior is
# already bounded by the user, so it is used as-is. Overridable via the `transform`
# keyword.
_default_transform(::Distribution) = identity
_default_transform(::AbstractVector{<:Distribution}) = identity
_default_transform(::AbstractPriorModel) = tanh

# Order-1 prior: unwrap a length-1 vector slot to its single element so an order-1
# `damp` is drawn as a scalar coefficient regardless of how it was written.
_order1_prior(prior) = prior
_order1_prior(prior::AbstractVector) = only(prior)

@doc raw"
An autoregressive AR(`p`) latent process.

```math
Z_t = \sum_{i=1}^{p} \rho_i Z_{t-i} + \epsilon_t
```

with damping coefficients ``\rho`` from the prior in `damp`, initial conditions
from the prior in `init`, and innovations from the error model `ϵ_t`. The order
`p` is fixed by the `damp` prior (a length-`k` vector ⇒ order `k`, a single
distribution or a process ⇒ order 1); the `init` prior is sized to match.

Each prior slot takes a raw prior: pass a bare `Distribution` (order 1), a vector
of them (order = its length), or a process (a latent model). The `ϵ_t`
innovations are a length-`n` PATH slot, so a bare `Distribution` there is
auto-wrapped in an [`Intercept`](@ref) — a **constant** innovation path (one
shared draw broadcast to every step); use [`IID`](@ref) for `n` independent
innovations. At order 1 the `damp` slot decides whether the coefficient is
**constant or time-varying**, through the same
[time-varying-parameter mechanism](@ref as_turing_submodel) any component can
use:

  - `AR(damp = Normal(...))` — a `Distribution` gives a **constant** coefficient,
    drawn as a single scalar RV (efficient, no length-`n` allocation);
  - `AR(damp = RandomWalk())` — a process gives a **time-varying** coefficient
    *path* ``\rho_t``, drawn at length `n-1` and threaded per step.

The coefficient is mapped through `transform` (default `tanh` for a process, so an
unbounded path stays in the stationary band; `identity` for a bounded
`Distribution`) and tracked as the generated quantity `ρ`, recoverable from the
chain (`group(chain, :ρ)`). Higher-order (`p > 1`) coefficients are constant;
time-varying higher-order AR is tracked in
[#113](https://github.com/EpiAware/ComposableTuringIDModels.jl/issues/113).

# Examples
```@example AR
using ComposableTuringIDModels, Distributions
ar = AR()
mdl = as_turing_model(ar, 10)
rand(mdl)
```
"
struct AR{D <: PriorLike, I <: PriorLike, P <: Int, E <: PriorLike,
    F <: Function} <: AbstractLatentModel
    "Prior for the damping coefficients."
    damp::D
    "Prior for the initial conditions."
    init::I
    "Order of the AR model."
    p::P
    "Error model for the innovations."
    ϵ_t::E
    "Map from the raw coefficient to the damping (default `tanh` for a process
    prior, `identity` for a bounded `Distribution`)."
    transform::F

    function AR(damp, init, p::Int, ϵ_t, transform)
        @assert p>0 "p must be greater than 0"
        _assert_prior_length(damp, p, "damp")
        _assert_prior_length(init, p, "init")
        new{typeof(damp), typeof(init), typeof(p), typeof(ϵ_t),
            typeof(transform)}(damp, init, p, ϵ_t, transform)
    end
end

function AR(damp::Sampleable, init::Sampleable; p::Int = 1,
        ϵ_t = HierarchicalNormal())
    return AR(; damp = fill(damp, p), init = fill(init, p), ϵ_t = ϵ_t)
end

function AR(; damp = truncated(Normal(0.0, 0.05), 0, 1), init = Normal(),
        ϵ_t = HierarchicalNormal(), transform = _default_transform(damp))
    # Order `p` is fixed by the damping prior (a length-`k` vector ⇒ order `k`, a
    # single distribution / process ⇒ order 1). The order-`p` initial-conditions
    # slot needs `p` values, so a bare `Distribution` is sized to a length-`p`
    # vector (one i.i.d. draw per lag); an explicit vector must already be length
    # `p`, and a process supplies its own length.
    p = _prior_order(damp)
    init = (p > 1 && init isa Distribution) ? fill(init, p) : init
    return AR(damp, init, p, _path_prior(ϵ_t), transform)
end

@model function as_turing_model(model::AR, n)
    p = model.p
    @assert n>p "n must be longer than the order of the autoregressive process"
    ar_init ~ as_turing_submodel(model.init, p; prefix = true)
    if p == 1
        # Order 1: draw the coefficient through the single seam. A `Distribution`
        # gives a scalar (constant, no length-`n` allocation); a process gives a
        # length-`(n-1)` path. `TVARStep` reads it per step with `_at`, so one
        # recursion serves both.
        damp_AR ~ as_turing_submodel(_order1_prior(model.damp), n - 1;
            prefix = true)
        # Track the (possibly time-varying) coefficient as a generated quantity so
        # it is recoverable from the chain; `transform` broadcasts over a scalar or
        # a path.
        ρ := model.transform.(damp_AR)
        ϵ_t ~ as_turing_submodel(model.ϵ_t, n - 1)
        z = accumulate_scan(TVARStep(ρ), only(ar_init),
            collect(zip(1:(n - 1), ϵ_t)))
        return z
    end
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
