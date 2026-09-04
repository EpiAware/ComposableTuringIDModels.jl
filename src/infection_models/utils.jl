# Growth-rate / reproduction-number conversion utilities shared by the
# infection models.

@doc raw"
Negative moment generating function of a discrete generation interval `w` at
rate `r`: ``\sum_i w_i e^{-r i}``.

# Arguments

  - `r`: the exponential growth rate.
  - `w`: the discrete generation interval weights.

# Examples
```@example neg_MGF
using ComposableTuringIDModels
ComposableTuringIDModels.neg_MGF(0.1, [0.2, 0.3, 0.5])
```
"
function neg_MGF(r, w::AbstractVector)
    return sum(w[i] * exp(-r * i) for i in 1:length(w))
end

# Derivative of `neg_MGF` with respect to `r`, used by the Newton step in
# `R_to_r`.
function _dneg_MGF_dr(r, w::AbstractVector)
    return -sum(w[i] * i * exp(-r * i) for i in 1:length(w))
end

@doc raw"
Approximate the exponential growth rate `r` implied by a reproduction number
`R₀` and discrete generation interval `w`.

Solves ``R_0 \sum_i w_i e^{-r i} = 1`` by a small-`r` initial guess refined with
`newton_steps` Newton iterations.

# Arguments

  - `R₀`: the reproduction number.
  - `w`: the discrete generation interval weights (or a [`Renewal`](@ref) model,
    whose generation interval is used).

# Keyword Arguments

  - `newton_steps`: number of Newton refinement steps (default `2`).
  - `Δd`: generation-interval discretisation width (default `1.0`).

# Examples
```@example R_to_r
using ComposableTuringIDModels
R_to_r(1.5, [0.2, 0.3, 0.5])
```
"
function R_to_r(R₀, w::AbstractVector{T}; newton_steps = 2, Δd = 1.0) where {T <: Real}
    mean_gen_time = dot(w, 1:length(w)) * Δd
    r_approx = (R₀ - 1) / (R₀ * mean_gen_time)
    for _ in 1:newton_steps
        r_approx -= (R₀ * neg_MGF(r_approx, w) - 1) /
            (R₀ * _dneg_MGF_dr(r_approx, w))
    end
    return r_approx
end

# The two renewal infection models, which are the ones carrying a generation
# interval.
# The deterministic scan and the centred stochastic loop take the same
# arguments, so the deterministic summaries below serve both.
const _RenewalModel = Union{Renewal, StochasticRenewal}

# The fixed generation interval of a renewal model, or a clear error when it is
# inferred.
# An uncertain interval varies per draw, so there is no single interval for these
# deterministic summaries to use.
function _fixed_gen_int(infection::_RenewalModel)
    infection.gen_int isa AbstractVector && return infection.gen_int
    throw(
        ArgumentError(
            "`R_to_r`/`expected_Rt` need a fixed generation interval, but this " *
                "model has an inferred (uncertain) generation interval that varies " *
                "per draw. Summarise the sampled interval per posterior draw instead."
        )
    )
end

# Only the renewal models carry a generation interval, so the model-typed
# method dispatches on them specifically (the other infection models have no
# `gen_int`).
function R_to_r(R₀, infection::_RenewalModel; newton_steps = 2, Δd = 1.0)
    return R_to_r(
        R₀, _fixed_gen_int(infection); newton_steps = newton_steps,
        Δd = Δd
    )
end

@doc raw"
Expected reproduction number ``R_t`` from a discrete generation interval and an
infection series.

```math
R_t = \frac{I_t}{\sum_{i=1}^{n} I_{t-i} g_i}
```

# Arguments

  - `gen_int`: the discrete generation interval weights (or a [`Renewal`](@ref)
    model, whose generation interval is used).
  - `infections`: the infection series (longer than the generation interval).

# Examples
```@example expected_Rt
using ComposableTuringIDModels
expected_Rt([0.2, 0.3, 0.5], [100.0, 200, 300, 400, 500])
```
"
function expected_Rt(gen_int::AbstractVector, infections::Vector{<:Real})
    n = length(gen_int)
    @assert n < length(infections) "Infections vector must be longer than the generation interval"
    denom_Rt = [
        dot(reverse(gen_int), infections[(t - n):(t - 1)])
            for t in (n + 1):length(infections)
    ]
    return infections[(n + 1):end] ./ denom_Rt
end

function expected_Rt(infection::_RenewalModel, infections::Vector{<:Real})
    return expected_Rt(_fixed_gen_int(infection), infections)
end

@doc raw"
Reproduction number implied by an exponential growth rate `r` and discrete
generation interval `w`: ``1 / \sum_i w_i e^{-r i}``.

# Arguments

  - `r`: the exponential growth rate.
  - `w`: the discrete generation interval weights.

# Examples
```@example r_to_R
using ComposableTuringIDModels
r_to_R(0.1, [0.2, 0.3, 0.5])
```
"
function r_to_R(r, w::AbstractVector)
    return 1 / neg_MGF(r, w)
end

# Shape helpers shared by `Renewal`, `DirectInfections` and `ExpGrowthRate`: how
# the initial-infections seed and the latent path are read at a `ModelShape`.

# Drive a scan one time step at a time: a path by its elements, a
# strata × time matrix by its columns.
_steps(Rt::AbstractVector) = Rt
_steps(Rt::AbstractMatrix) = collect.(eachcol(Rt))

# The initial-infections seed is one value for a single series and one per
# stratum otherwise.
# Broadcast against `n::Dims{2}` gives the same seed to every stratum.
# The slot is a level, so a draw of any other length is an error rather than a
# silently collapsed path, and `SeedingPath` is the way to estimate a seeding
# window.
function _seed(x, ::Int)
    length(x) == 1 || throw(
        ArgumentError(
            "`initialisation` is the level of infections at t₀, so it must " *
                "draw one value for a single series, but it drew " *
                "$(length(x)). A `Renewal` can estimate a whole seeding " *
                "window instead: wrap the process in a `SeedingPath`."
        )
    )
    return only(x)
end
_seed(x::Real, n::Dims{2}) = fill(x, n[1])
function _seed(x::AbstractVector, n::Dims{2})
    @assert length(x) == n[1] "`initialisation` drew $(length(x)) values but the model has $(n[1]) strata"
    return x
end

# The implied initial growth rate from a seed `R_t` and generation interval.
# A shared interval broadcasts across strata with `Ref(w)`, and a per-stratum
# interval is read row by row.
_init_rate(Rt₀::Real, w::AbstractVector) = R_to_r(Rt₀, w)
_init_rate(Rt₀::AbstractVector, w::AbstractVector) = R_to_r.(Rt₀, Ref(w))
function _init_rate(Rt₀::AbstractVector, w::AbstractMatrix)
    return [R_to_r(Rt₀[g], view(w, g, :)) for g in eachindex(Rt₀)]
end

# Everything both renewal models do before their recursions diverge.
# That is the `R_t` path, the generation interval and the step it bakes, the
# seeding window and the initial state.
# `Renewal` scans it and `StochasticRenewal` loops over it drawing as it goes.
#
# Called through `to_submodel(..., false)`, which keeps `Z_t`, `gen`,
# `init_incidence` and `scan_step` at the names they carry when drawn inline.
@model function _renewal_setup(infection::_RenewalModel, n::ModelShape)
    Z_t ~ as_turing_submodel(infection.rt, n)
    Rt = infection.transformation.(Z_t)

    # An inferred generation interval is sampled through the single seam, with
    # the lag-0 bin dropped and the remainder renormalised per draw.
    # The step is rebuilt per draw so the gradient flows through the
    # discretisation.
    if infection.gen_int isa AbstractPriorModel
        gen ~ as_turing_submodel(
            infection.gen_int, _gen_int_shape(n)...; prefix = true
        )
        gen_int = _drop_lag_zero(gen)
        step = _bake_step(
            _renewal_step(gen_int, infection.mixing), infection.modifiers
        )
    else
        gen_int = infection.gen_int
        step = infection.recurrent_step
    end

    # The `initialisation` slot is either a level at ``t_0`` or a whole seeding
    # path drawn at the incidence window's own shape.
    # Both branches are decided by the slot's type, so nothing is chosen at run
    # time.
    if infection.initialisation isa SeedingPath
        init_incidence ~ as_turing_submodel(
            infection.initialisation.model,
            _seeding_shape(n, _n_lags(gen_int)); prefix = true
        )
    else
        init_incidence ~ as_turing_submodel(
            infection.initialisation, _n_strata(n); prefix = true
        )
    end

    # Resolve the step before the recursion, so a modifier or a `mixing` model
    # carrying priors draws them here and hands back what the recursion uses.
    # A purely deterministic step returns itself, so nothing here branches on
    # what the step holds.
    scan_step ~ as_turing_submodel(step, n)

    Rts = _steps(Rt)
    # A seeding path is the incidence window itself; a level is decayed into one
    # at the growth rate implied by ``\mathcal R_1``.
    if infection.initialisation isa SeedingPath
        seed_window = infection.transformation.(init_incidence)
        n_seed = _n_lags(seed_window)
        @assert n_seed == _n_lags(gen_int) "a `SeedingPath` must draw one " *
            "value per generation-interval lag, but it drew $(n_seed) for " *
            "$(_n_lags(gen_int)) lags"
        init = renewal_init_state(scan_step, seed_window)
    else
        I₀ = infection.transformation.(_seed(init_incidence, n))
        init = _make_renewal_init(scan_step, gen_int, I₀, first(Rts))
    end
    return (; Z_t, Rts, scan_step, init)
end
