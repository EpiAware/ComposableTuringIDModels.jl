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
function R_to_r(R₀, w::Vector{T}; newton_steps = 2, Δd = 1.0) where {T <: AbstractFloat}
    mean_gen_time = dot(w, 1:length(w)) * Δd
    r_approx = (R₀ - 1) / (R₀ * mean_gen_time)
    for _ in 1:newton_steps
        r_approx -= (R₀ * neg_MGF(r_approx, w) - 1) /
                    (R₀ * _dneg_MGF_dr(r_approx, w))
    end
    return r_approx
end

# Only `Renewal` carries a generation interval, so the model-typed method
# dispatches on it specifically (the other infection models have no `gen_int`).
function R_to_r(R₀, infection::Renewal; newton_steps = 2, Δd = 1.0)
    return R_to_r(R₀, infection.gen_int; newton_steps = newton_steps, Δd = Δd)
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
    @assert n<length(infections) "Infections vector must be longer than the generation interval"
    denom_Rt = [dot(reverse(gen_int), infections[(t - n):(t - 1)])
                for t in (n + 1):length(infections)]
    return infections[(n + 1):end] ./ denom_Rt
end

function expected_Rt(infection::Renewal, infections::Vector{<:Real})
    return expected_Rt(infection.gen_int, infections)
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
