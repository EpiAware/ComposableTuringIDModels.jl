# Renewal accumulation steps: the force-of-infection core, the coupling seam it
# dispatches through, and the window arithmetic shared by every renewal step.

@doc raw"
Abstract supertype for renewal accumulation steps (constant generation interval,
with or without susceptible depletion).
"
abstract type AbstractConstantRenewalStep <: AbstractAccumulationStep end

@doc raw"
Generation-time-weighted incidence pressure on each stratum, before ``R_t``.

```math
\Lambda_{g,t} = \sum_h K_{gh} \sum_i g_i I_{h,t-i}
```

This is the one dispatch point for coupling between strata. It sits *between*
the generation-time convolution and ``R_t``.
That is why coupling is a slot on the renewal core rather than a renewal
modifier: a modifier transforms the finished force of infection, and
``K(R \odot v) \neq R \odot Kv``.

The `mixing` argument chooses the method.

  - `I` (`LinearAlgebra.I`), the default, leaves the strata uncoupled. For a
    single series the window is a vector and this is a dot product, so a
    one-stratum model runs exactly the arithmetic it would with no strata axis
    at all.
  - A `strata × strata` matrix mixes the convolved histories, so off-diagonal
    mass is infection arriving from another stratum. Row `g` says where stratum
    `g`'s force of infection comes from.
  - A `strata × lags` generation interval gives each stratum its own interval,
    with or without a mixing matrix on top.
  - A `strata × strata × lags` array gives a per-pair generation interval, so a
    between-stratum transmission can carry a longer effective interval. Build
    one with [`pairwise_gen_int`](@ref). Its lags are indexed forwards, so
    `K[:, :, 1]` is the lag-1 weight.

Nothing about a mixing matrix is checked or normalised. A zero column
`K[:, h]` means stratum `h` infects nobody, though its own incidence still
evolves through `K[h, h]`. A zero row `K[g, :]` means stratum `g` receives no
force at all, so its incidence falls to zero and stays there.

Extend it by adding a method: any operator on the convolved window is a valid
coupling.

# Arguments

  - `mixing`: the coupling operator.
  - `g`: the reversed generation interval, a vector or a `strata × lags` matrix.
  - `window`: the recent incidence, oldest to newest.

# Examples
```@example renewal_pressure
using ComposableTuringIDModels, LinearAlgebra
window = [10.0 20.0; 5.0 8.0]        # 2 strata, 2 lags
g = [0.4, 0.6]
renewal_pressure(I, g, window)       # uncoupled
```
"
function renewal_pressure end

renewal_pressure(::UniformScaling, g, window::AbstractVector) = sum(window .* g)

function renewal_pressure(::UniformScaling, g, window::AbstractMatrix)
    vec(sum(window .* reshape(g, 1, :); dims = 2))
end

# A `strata × lags` generation interval: each stratum convolves its own row of
# the window with its own interval.
function renewal_pressure(
        ::UniformScaling, g::AbstractMatrix, window::AbstractMatrix)
    return vec(sum(window .* g; dims = 2))
end

# A mixing matrix redistributes the convolved histories, whichever generation
# interval produced them.
function renewal_pressure(K::AbstractMatrix, g, window::AbstractMatrix)
    return vec(sum(K .* reshape(renewal_pressure(I, g, window), 1, :); dims = 2))
end

# A per-pair generation interval, `K[g, h, i]` being the weight stratum `h`
# puts on stratum `g` at lag `i`. The window runs oldest to newest, so lag `i`
# is column `end - i + 1`; the reversed generation interval `g` is unused,
# because the array carries the intervals itself.
function renewal_pressure(
        K::AbstractArray{<:Any, 3}, g, window::AbstractMatrix)
    last_lag = size(window, 2)
    return sum(view(K, :, :, i) * view(window, :, last_lag - i + 1)
    for i in axes(K, 3))
end

@doc raw"
Force of infection for a constant-generation-interval renewal step: the
reproduction number ``R_t`` times the [`renewal_pressure`](@ref) on each
stratum,

```math
R_t \sum_{i=1}^{n-1} I_{t-i} g_i.
```

This is the raw new-incidence term before any modifier (e.g. susceptible
depletion) is applied. `renewal_foi` is shared by the internal
`ConstantRenewalStep` core and the composable [`RenewalStep`](@ref) so the two
cannot drift. It reads the coupling operator off the step and defers to
[`renewal_pressure`](@ref), so a new coupling needs no change here.

# Arguments

  - `step`: the renewal step supplying the generation interval and the mixing.
  - `window`: the recent incidence, oldest to newest.
  - `Rt`: the reproduction number, a scalar or one value per stratum.

# Examples
```@example renewal_foi
using ComposableTuringIDModels
step = ComposableTuringIDModels.ConstantRenewalStep(reverse([0.2, 0.3, 0.5]))
ComposableTuringIDModels.renewal_foi(step, [10.0, 20.0, 30.0], 1.5)
```
"
function renewal_foi(step::AbstractConstantRenewalStep, window, Rt)
    return Rt .* renewal_pressure(step.mixing, step.rev_gen_int, window)
end

@doc raw"
Renewal step with a constant generation interval (stored reversed) and a
coupling operator.

```math
I_t = R_t \sum_h K_{gh} \sum_i I_{h,t-i} g_i
```

`rev_gen_int` is a vector for one shared generation interval and a
`strata × lags` matrix for one interval per stratum. `mixing` defaults to `I`,
which leaves the strata uncoupled; see [`renewal_pressure`](@ref) for what else
it accepts.

## Fields

  - `rev_gen_int`: the reversed generation interval.
  - `mixing`: the coupling operator applied to the convolved window.
"
struct ConstantRenewalStep{T, K} <: AbstractConstantRenewalStep
    "The reversed generation interval."
    rev_gen_int::T
    "The coupling operator applied to the convolved window."
    mixing::K
end

ConstantRenewalStep(rev_gen_int) = ConstantRenewalStep(rev_gen_int, I)

# --- window arithmetic ------------------------------------------------------
#
# One incidence window serves both shapes. A single series is a vector of the
# last `lags` incidences. Several strata are a `strata × lags` matrix, one row
# per stratum. Each helper has one method per shape, so the recursion, the
# state assembly and the seeding are written once.

# The newest entry of the window is the incidence just committed.
_newest(window::AbstractVector) = last(window)
_newest(window::AbstractMatrix) = window[:, end]

# Assemble the scanned states into the returned series: a path stays a vector,
# per-stratum columns become a `strata × time` matrix, which is `size(I_t)`.
_series(v::AbstractVector{<:Real}) = v
_series(v::AbstractVector{<:AbstractVector}) = reduce(hcat, v)

# Drop the oldest entry and commit the new incidence.
_advance(window::AbstractVector, new) = vcat(window[2:end], new)
function _advance(window::AbstractMatrix, new)
    return hcat(collect(view(window, :, 2:size(window, 2))), new)
end

# The number of lags a generation interval covers.
_n_lags(g::AbstractVector) = length(g)
_n_lags(g::AbstractMatrix) = size(g, 2)

# Reverse a generation interval along its lag axis, whatever its shape.
_reverse_lags(g::AbstractVector) = reverse(g)
_reverse_lags(g::AbstractMatrix) = reverse(g; dims = 2)

function (recurrent_step::ConstantRenewalStep)(window, Rt)
    new_incidence = renewal_foi(recurrent_step, window, Rt)
    return _advance(window, new_incidence)
end

# One series: a window decaying at the implied rate. Several strata: one row
# each, at that stratum's own seed and implied rate. `r` broadcasts, so it is a
# scalar for a shared generation interval and a vector for per-stratum ones.
function _renewal_init_state(::ConstantRenewalStep, I₀::Real, r, len_gen_int)
    return I₀ * [exp(-r * t) for t in (len_gen_int - 1):-1:0]
end

function _renewal_init_state(
        ::ConstantRenewalStep, I₀::AbstractVector, r, len_gen_int)
    return I₀ .* exp.(.-r .* permutedims((len_gen_int - 1):-1:0))
end

function get_state(::ConstantRenewalStep, initial_state, state)
    return _series(_newest.(state))
end

# `ConstantRenewalStep` is the force-of-infection primitive. The renewal step
# users build through the [`Renewal`](@ref) helper is `RenewalStep` (see
# `RenewalStep.jl`), which wraps this core and composes modifiers on top.
