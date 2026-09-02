# The renewal accumulation step.
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
# This deliberately avoids the naive `state = step(state, ϵ)` sequential-threading
# composite: that double-advances time for steps whose call commits a new
# value. The contract here separates the *contribution* to the new
# incidence (each modifier transforms it) from the single shared-window *advance*
# performed once per step. AR/MA step-fusion (a different, non-shared state
# contract) is out of scope and stays as model nesting.

@doc raw"
Abstract supertype for renewal modifiers composed onto a [`RenewalStep`](@ref).

The type covers two shapes. A **scan** modifier is called by the scan itself: it
transforms the proposed new incidence and carries its own substate. It
implements

  - `modifier_init_state(mod, window)` — the modifier's initial substate, given
    the step's initial incidence window. The window is passed because a
    substate that tracks the incidence has to match its shape: a scalar for one
    series, one value per stratum for a stratified renewal.
  - `apply_modifier(mod, incidence, substate)` — return
    `(new_incidence, new_substate)`.

A **prior-carrying** modifier instead needs parameters sampled *before* the
scan — a per-time importation rate, say — which the scan cannot draw, because a
scan step is a plain deterministic function. It implements the other part of the
interface:

  - `as_turing_model(mod, n)` — a `DynamicPPL.Model` returning the modifier
    used in the scan. The default method samples nothing and returns `mod`
    unchanged, so a purely deterministic modifier (e.g.
    `SusceptibleDepletion(1000.0)`) is a scan modifier and implements nothing
    extra. A prior-carrying modifier implements this method, draws its slots
    through [`as_turing_submodel`](@ref), and returns a *resolved* scan
    modifier holding the drawn values — e.g. [`ImportedCases`](@ref) resolves
    to an [`ImportedRate`](@ref) and implements no scan interface of its own.

[`RenewalStep`](@ref) resolves its whole modifier tuple through this one seam
(see its `as_turing_model` method), so a sampling modifier needs no special
handling anywhere in the renewal model. A prior-carrying modifier reaching the
scan means the step was built by hand and never resolved, so the scan interface
errors with that message rather than a bare `MethodError`.

Each modifier's sampled variables are prefixed by its **position** in the
modifier tuple, `modifier_<i>`, so two modifiers of the same kind cannot
collide on a variable name — in
`Renewal(gen_int, SusceptibleDepletion(N), ImportedCases(Normal()))` the
importation rate is `modifier_2.import_rates`. Inserting or reordering a
modifier therefore renames the variables of every modifier after it.
"
abstract type AbstractRenewalModifier end

# A prior-carrying modifier resolves to a separate scan modifier and has no
# scan interface of its own, so reaching the scan with one means the step was
# never resolved. Name the modifier and point at the seam instead of failing
# with a bare `MethodError` from inside the recursion.
function _unresolved_modifier(mod)
    return error(
        "$(typeof(mod)) has no scan interface. A modifier " *
            "carrying priors must be resolved before the scan through " *
            "the step's `as_turing_model` seam, i.e. " *
            "`as_turing_submodel(step, n)` (what `Renewal` does); a " *
            "deterministic modifier must implement " *
            "`modifier_init_state` and `apply_modifier`."
    )
end

function modifier_init_state(mod::AbstractRenewalModifier, window)
    return _unresolved_modifier(mod)
end

function apply_modifier(mod::AbstractRenewalModifier, incidence, substate)
    return _unresolved_modifier(mod)
end

# The default pre-scan seam: a modifier with no parameters of its own samples
# nothing and scans as itself. Kept as a `@model` so every modifier resolves
# through the same submodel call, with no branch on whether it samples.
@model function as_turing_model(mod::AbstractRenewalModifier, n)
    return mod
end

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
population and susceptible depletion, e.g.
`Renewal(gen_int, SusceptibleDepletion(N))`.

`pop_size` is a scalar for one series, and a per-stratum vector for a
stratified renewal.
Each stratum then depletes its own pool.
A scalar given to a stratified renewal is shared by every stratum, so each
depletes a separate pool of the same size.

It samples nothing, so it is a plain scan modifier: the pre-scan seam (see
[`AbstractRenewalModifier`](@ref)) returns it unchanged.

## Fields

  - `pop_size`: the population size, one value or one per stratum.
"
struct SusceptibleDepletion{T} <: AbstractRenewalModifier
    "The population size, one value or one per stratum."
    pop_size::T
end

# The susceptible pool has to match the incidence it depletes, so a scalar pool
# given to a stratified renewal is spread over the strata rather than left as a
# scalar that the first step would silently widen.
function modifier_init_state(mod::SusceptibleDepletion, window)
    return _match_strata(mod.pop_size, window)
end

_match_strata(pop_size, ::AbstractVector) = pop_size
_match_strata(pop_size::Real, window::AbstractMatrix) = fill(pop_size, size(window, 1))
function _match_strata(pop_size::AbstractVector, window::AbstractMatrix)
    @assert length(pop_size) == size(window, 1) "`pop_size` has " *
        "$(length(pop_size)) entries " *
        "but the renewal has " *
        "$(size(window, 1)) strata"
    return pop_size
end

# The susceptible fraction is floored because a large force of infection can
# take more than the pool holds, leaving `S` negative for the rest of the run.
# The floor keeps the recursion going there. Every operation is broadcast, so
# one line serves a scalar and a per-stratum pool.
function apply_modifier(mod::SusceptibleDepletion, incidence, S)
    new_incidence = max.(S ./ mod.pop_size, 1.0e-6) .* incidence
    return new_incidence, S .- new_incidence
end

@doc raw"
The renewal accumulation step: a force-of-infection `core` (a constant
generation interval by default) with a tuple of `modifiers`
([`AbstractRenewalModifier`](@ref)s) composing on top, sharing one incidence
window.

With no modifiers it is a plain renewal recurrence. With modifiers its state is
`(; val, window, substates)` — the newest incidence, the shared window, and one
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

function renewal_init_state(step::_PlainRenewalStep, I₀, r_approx, len_gen_int)
    return renewal_init_state(step.core, I₀, r_approx, len_gen_int)
end

function renewal_init_state(step::_PlainRenewalStep, window::AbstractArray)
    return renewal_init_state(step.core, window)
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
        Base.tail(mods), inc, Base.tail(substates)
    )
    return rest_inc, (s, rest_states...)
end

# `exp_val` is the force of infection reaching the modifier tuple; `val` is that
# same incidence after every modifier has transformed it.
function (step::RenewalStep)(state, Rt)
    foi = renewal_foi(step.core, state.window, Rt)
    new_incidence, new_substates = _thread_modifiers(
        step.modifiers, foi, state.substates
    )
    new_window = _advance(state.window, new_incidence)
    return (;
        val = new_incidence, exp_val = foi, window = new_window,
        substates = new_substates,
    )
end

function renewal_init_state(step::RenewalStep, I₀, r_approx, len_gen_int)
    return renewal_init_state(
        step, renewal_init_window(step.core, I₀, r_approx, len_gen_int)
    )
end

function renewal_init_state(step::RenewalStep, window::AbstractArray)
    substates = map(mod -> modifier_init_state(mod, window), step.modifiers)
    return (;
        val = _newest(window), exp_val = _newest(window), window = window,
        substates = substates,
    )
end

function get_state(::RenewalStep, initial_state, state)
    return _series(state .|> x -> x.val)
end

# Assemble the noise-free expectation series. A `_PlainRenewalStep` delegates
# its scan to the bare core, whose states carry no `exp_val`; there is nothing
# to transform, so the expectation is the draw itself.
function get_expected_state(step::_PlainRenewalStep, initial_state, state)
    return get_state(step, initial_state, state)
end

@doc raw"
Assemble the noise-free renewal expectation series from a scan's accumulated
states, mirroring [`get_state`](@ref).

Each step's expectation is the force of infection the core computed, before any
modifier transformed it. Without modifiers it equals the committed draw.

# Arguments

  - `step`: the renewal step scanned (the same step `accumulate_scan` was
    called with).
  - `initial_state`: the seed state the scan started from.
  - `state`: the raw accumulated output of `accumulate`.

# Examples
```@example get_expected_state
using ComposableTuringIDModels
core = ComposableTuringIDModels.ConstantRenewalStep(reverse([0.2, 0.3, 0.5]))
step = RenewalStep(core, (SusceptibleDepletion(100.0),))
init = ComposableTuringIDModels.renewal_init_state(
    step, [1.0], 0.0, 3)
result = accumulate(step, [1.0, 2.0, 3.0]; init = init)
ComposableTuringIDModels.get_expected_state(step, init, result)
```
"
function get_expected_state(::RenewalStep, initial_state, state)
    return _series(state .|> x -> x.exp_val)
end

# --- the pre-scan seam ------------------------------------------------------
#
# A scan step is a deterministic function, so a modifier that needs sampled
# parameters (an importation rate, say) cannot draw them inside the scan. The
# step is resolved ONCE before the scan: every modifier is sampled through its
# own `as_turing_model` and the step is rebuilt from the resolved modifiers.
# Modifiers that sample nothing return themselves, so one call covers both and
# nothing here ever tests what a modifier is.

@doc raw"
Resolve an accumulation step ahead of the scan, sampling any parameters its
parts carry.

The default method samples nothing and returns the step unchanged.
A [`RenewalStep`](@ref) resolves both its core and its
[`AbstractRenewalModifier`](@ref)s through their own `as_turing_model` methods
and rebuilds itself from the resolved parts. [`Renewal`](@ref) draws its step
through this one seam, so neither a modifier with priors nor a drawn coupling
operator (a [`MixingStep`](@ref)) needs special handling in the infection
model.

Each modifier is prefixed by its position in the tuple, so the ``i``th
modifier's variables are namespaced `modifier_<i>` (see
[`AbstractRenewalModifier`](@ref)).

# Arguments

  - `step`: the accumulation step to resolve.
  - `n`: the series length the scan will run over.
"
@model function as_turing_model(step::AbstractAccumulationStep, n)
    return step
end

@model function as_turing_model(step::RenewalStep, n)
    # The core composes flat, as it does when a modifier-free renewal scans the
    # core directly, so a drawn coupling operator is named the same way whether
    # or not the step carries modifiers. The modifiers below carry their own
    # positional prefixes, so nothing here can collide.
    core ~ as_turing_submodel(step.core, n)
    modifiers ~ to_submodel(_resolve_modifiers(step.modifiers, n, 1), false)
    return RenewalStep(core, modifiers)
end

# Resolve a modifier tuple, one submodel per modifier. The recursion is
# structural (over the tuple, not an index into it) so each level has a
# concrete tuple type. Each modifier is prefixed by its position, so two
# modifiers of the same kind cannot collide on a variable name.
@model function _resolve_modifiers(mods::Tuple{}, n, index)
    return ()
end

@model function _resolve_modifiers(mods::Tuple, n, index)
    modifier ~ to_submodel(
        prefix(as_turing_model(first(mods), n), Symbol(:modifier_, index)),
        false
    )
    rest ~ to_submodel(
        _resolve_modifiers(Base.tail(mods), n, index + 1), false
    )
    return (modifier, rest...)
end
