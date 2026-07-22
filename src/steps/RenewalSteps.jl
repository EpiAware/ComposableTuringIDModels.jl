# Renewal accumulation steps (constant generation interval, with/without
# susceptible depletion).

@doc raw"
Abstract supertype for renewal accumulation steps (constant generation interval,
with or without susceptible depletion).
"
abstract type AbstractConstantRenewalStep <: AbstractAccumulationStep end

@doc raw"
Force of infection for a constant-generation-interval renewal step: the
reproduction number ``R_t`` times the convolution of the recent incidence window
with the (reversed) generation interval,

```math
R_t \sum_{i=1}^{n-1} I_{t-i} g_i.
```

This is the raw new-incidence term before any modifier (e.g. susceptible
depletion) is applied. It is shared by the internal `ConstantRenewalStep` core
and the composable [`RenewalStep`](@ref) so the two cannot drift.
"
function renewal_foi(step::AbstractConstantRenewalStep, recent_incidence, Rt)
    return Rt * dot(recent_incidence, step.rev_gen_int)
end

@doc raw"
Renewal step with a constant generation interval (stored reversed).

```math
I_t = R_t \sum_{i=1}^{n-1} I_{t-i} g_i
```
"
struct ConstantRenewalStep{T} <: AbstractConstantRenewalStep
    rev_gen_int::Vector{T}
end

function (recurrent_step::ConstantRenewalStep)(recent_incidence, Rt)
    new_incidence = renewal_foi(recurrent_step, recent_incidence, Rt)
    return vcat(recent_incidence[2:end], new_incidence)
end

function _renewal_init_state(::ConstantRenewalStep, I₀, r_approx, len_gen_int)
    return I₀ * [exp(-r_approx * t) for t in (len_gen_int - 1):-1:0]
end

get_state(::ConstantRenewalStep, initial_state, state) = last.(state)

# `ConstantRenewalStep` is the internal force-of-infection primitive. The public
# renewal step is `RenewalStep` (see `RenewalStep.jl`), which wraps this core and
# composes modifiers (e.g. `SusceptibleDepletion`) on top.
