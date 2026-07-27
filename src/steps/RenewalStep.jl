# The renewal accumulation step (#48).
#
# `RenewalStep` is *the* renewal step: a constant-generation-interval force of
# infection with an ordered tuple of modifiers composing on top, sharing one
# incidence window. With no modifiers it is a plain renewal process; a
# `SusceptibleDepletion` modifier makes it a renewal with a fixed population, and
# further renewal-family mechanisms (waning immunity, seasonality, …) compose the
# same way. The plain force-of-infection core (`ConstantRenewalStep`) is an
# internal primitive; users build `RenewalStep`s through the [`Renewal`](@ref)
# helper.
#
# This is deliberately NOT the naive `state = step(state, ϵ)` sequential-threading
# composite sketched upstream: that double-advances time for steps whose call
# commits a new value. The contract here separates the *contribution* to the new
# incidence (each modifier transforms it) from the single shared-window *advance*
# performed once per step. AR/MA step-fusion (a different, non-shared state
# contract) is out of scope and stays as model nesting; see #48.

@doc raw"
Abstract supertype for renewal modifiers composed onto a [`RenewalStep`](@ref).

A modifier transforms the proposed new incidence and carries its own substate
across the scan. Concrete modifiers implement

  - `modifier_init_state(mod)` — the modifier's initial substate.
  - `apply_modifier(mod, incidence, substate)` — return
    `(new_incidence, new_substate)`.
"
abstract type AbstractRenewalModifier end

@doc raw"
Susceptible-depletion modifier for [`RenewalStep`](@ref).

Scales the proposed incidence by the available susceptible fraction and depletes
the susceptible pool,

```math
I_t = \frac{S_{t-1}}{N} \, R_t \sum_{i=1}^{n-1} I_{t-i} g_i, \qquad
S_t = S_{t-1} - I_t,
```

with population size ``N`` = `pop_size`. Its substate is the current susceptible
count ``S``. Adding it to a renewal step gives a renewal process with a fixed
population and susceptible depletion, e.g. `Renewal(data, SusceptibleDepletion(N))`.
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
The renewal accumulation step: a force-of-infection `core` (a constant
generation interval by default) with a tuple of `modifiers`
([`AbstractRenewalModifier`](@ref)s) composing on top, sharing one incidence
window.

With no modifiers it is a plain renewal recurrence. With modifiers its state is
`[window, substate₁, …, substateₙ]` — the incidence window followed by one
substate per modifier; each step computes the core force of infection, threads it
through the modifiers (each transforming the incidence and updating its own
substate), then advances the shared window once with the final incidence.

`RenewalStep(core, (SusceptibleDepletion(N),))` is a renewal process with a fixed
population `N` and susceptible depletion; a [`Renewal`](@ref) built with a
`SusceptibleDepletion(N)` modifier uses exactly this step.
"
struct RenewalStep{R <: AbstractConstantRenewalStep, M <: Tuple} <:
       AbstractConstantRenewalStep
    core::R
    modifiers::M
end

RenewalStep(core::AbstractConstantRenewalStep) = RenewalStep(core, ())

# The renewal force of infection is defined by the core primitive.
renewal_foi(step::RenewalStep, window, Rt) = renewal_foi(step.core, window, Rt)

# No modifiers: behave exactly as the plain force-of-infection core (bare-window
# state, no per-step overhead), so a modifier-free renewal is unchanged.
const _PlainRenewalStep = RenewalStep{<:AbstractConstantRenewalStep, Tuple{}}

(step::_PlainRenewalStep)(state, Rt) = step.core(state, Rt)

function _renewal_init_state(step::_PlainRenewalStep, I₀, r_approx, len_gen_int)
    return _renewal_init_state(step.core, I₀, r_approx, len_gen_int)
end

function get_state(step::_PlainRenewalStep, initial_state, state)
    return get_state(step.core, initial_state, state)
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

function (step::RenewalStep)(state, Rt)
    window = state[1]
    substates = ntuple(i -> state[i + 1], length(step.modifiers))
    foi = renewal_foi(step.core, window, Rt)
    new_incidence, new_substates = _thread_modifiers(step.modifiers, foi, substates)
    new_window = vcat(window[2:end], new_incidence)
    return [new_window, new_substates...]
end

function _renewal_init_state(step::RenewalStep, I₀, r_approx, len_gen_int)
    window = _renewal_init_state(step.core, I₀, r_approx, len_gen_int)
    return [window, map(modifier_init_state, step.modifiers)...]
end

get_state(::RenewalStep, initial_state, state) = state .|> st -> last(st[1])
