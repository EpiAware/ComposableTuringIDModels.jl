# Composable renewal accumulation step (#48 Phase 2, proof-of-concept).
#
# The renewal family is the one place where composing accumulation steps has a
# naturally uniform state contract: every mechanism reads and writes the same
# shared incidence window. `ComposedRenewalStep` builds a renewal recurrence from
# a `ConstantRenewalStep` force-of-infection core plus an ordered tuple of
# modifiers that transform the proposed incidence and carry their own substate
# (susceptible depletion, and — later — waning immunity, seasonality, …).
#
# This is deliberately NOT the naive `state = step(state, ϵ)` sequential-threading
# composite sketched upstream: that double-advances time for steps whose call
# commits a new value. The contract here separates the *contribution* to the new
# incidence (each modifier transforms it) from the single shared-window *advance*
# the composite performs once per step. AR/MA step-fusion (a different, non-shared
# state contract) is out of scope and stays as model nesting; see #48.

@doc raw"
Abstract supertype for renewal modifiers used with [`ComposedRenewalStep`](@ref).

A modifier transforms the proposed new incidence and carries its own substate
across the scan. Concrete modifiers implement

  - `modifier_init_state(mod)` — the modifier's initial substate.
  - `apply_modifier(mod, incidence, substate)` — return
    `(new_incidence, new_substate)`.
"
abstract type AbstractRenewalModifier end

@doc raw"
Susceptible-depletion modifier for [`ComposedRenewalStep`](@ref).

Scales the proposed incidence by the available susceptible fraction and depletes
the susceptible pool,

```math
I_t = \frac{S_{t-1}}{N} \, R_t \sum_{i=1}^{n-1} I_{t-i} g_i, \qquad
S_t = S_{t-1} - I_t,
```

with population size ``N`` = `pop_size`. Its substate is the current susceptible
count ``S``. Composing it over a [`ConstantRenewalStep`](@ref) reproduces
[`ConstantRenewalWithPopulationStep`](@ref) exactly.
"
struct SusceptibleDepletion{T} <: AbstractRenewalModifier
    pop_size::T
end

modifier_init_state(mod::SusceptibleDepletion) = mod.pop_size

function apply_modifier(mod::SusceptibleDepletion, incidence, S)
    new_incidence = max(S / mod.pop_size, 1e-6) * incidence
    return new_incidence, S - new_incidence
end

@doc raw"
A renewal accumulation step assembled from a force-of-infection `core` (a
[`ConstantRenewalStep`](@ref)) and a tuple of `modifiers`
([`AbstractRenewalModifier`](@ref)s), sharing one incidence window.

Its state is `[window, substate₁, …, substateₙ]`: the incidence window followed
by one substate per modifier. Each step computes the core force of infection,
threads it through the modifiers (each transforming the incidence and updating
its own substate), then advances the shared window once with the final incidence.

`ComposedRenewalStep(ConstantRenewalStep(g), (SusceptibleDepletion(N),))`
reproduces `ConstantRenewalWithPopulationStep(g, N)`; this equivalence is the
acceptance test for the composable contract (#48).
"
struct ComposedRenewalStep{R <: AbstractConstantRenewalStep, M <: Tuple} <:
       AbstractConstantRenewalStep
    core::R
    modifiers::M
end

# Thread the proposed incidence through the modifier tuple, collecting each
# modifier's updated substate. Recursive over the tuple to stay type-stable and
# AD-friendly (no mutation of tracked state).
_thread_modifiers(::Tuple{}, incidence, ::Tuple{}) = (incidence, ())
function _thread_modifiers(mods::Tuple, incidence, substates::Tuple)
    inc, s = apply_modifier(first(mods), incidence, first(substates))
    rest_inc, rest_states = _thread_modifiers(
        Base.tail(mods), inc, Base.tail(substates))
    return rest_inc, (s, rest_states...)
end

function (step::ComposedRenewalStep)(state, Rt)
    window = state[1]
    substates = ntuple(i -> state[i + 1], length(step.modifiers))
    foi = renewal_foi(step.core, window, Rt)
    new_incidence, new_substates = _thread_modifiers(step.modifiers, foi, substates)
    new_window = vcat(window[2:end], new_incidence)
    return [new_window, new_substates...]
end

function _renewal_init_state(
        step::ComposedRenewalStep, I₀, r_approx, len_gen_int)
    window = _renewal_init_state(step.core, I₀, r_approx, len_gen_int)
    return [window, map(modifier_init_state, step.modifiers)...]
end

get_state(::ComposedRenewalStep, initial_state, state) = state .|> st -> last(st[1])
