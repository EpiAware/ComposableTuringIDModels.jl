# Imported (externally seeded) cases as a renewal modifier (#189).
#
# Unlike `SusceptibleDepletion`, which is a pure function of the proposed
# incidence and its own substate, importation adds a per-time rate that is
# *sampled* — it is a prior slot, not a constant carried on the step. So the
# modifier itself is a no-op in the `apply_modifier` chain and only marks the
# step as importing; `Renewal` samples the rate as a length-`n` submodel and
# zips it into the scan input, and `RenewalStep`'s tuple-input method adds it
# after the modifier chain. That keeps importation additive on the *committed*
# incidence rather than on the force of infection, so a later multiplicative
# modifier cannot scale imports away.
#
# The slot is on the UNCONSTRAINED scale, like `Renewal`'s `rt` and
# `initialisation`: `Renewal` pushes the sampled value through its
# `transformation` before adding it. An importation rate is a count per unit
# time and cannot be negative, and a slot that accepts any latent process has to
# make that hold structurally — an unconstrained process (a `RandomWalk`, an
# `AR`, a GP) spends much of its time below zero, and subtracting that from the
# incidence is not importation.

@doc raw"
Imported-cases modifier for [`RenewalStep`](@ref).

Adds an externally seeded importation rate to the renewal incidence,

```math
\iota_t = g(\tilde\iota_t), \qquad
I_t = \iota_t + \mathcal R_t \sum_{i=1}^{n-1} I_{t-i} g_i,
```

so infections can arrive from outside the modelled population — the mechanism
behind a renewal process that would otherwise die out from a zero initial
incidence, and behind reintroduction after local elimination.

`importation_rate` is a length-`n` prior slot holding the **unconstrained**
rate ``\tilde\iota_t``, exactly as [`Renewal`](@ref)'s `rt` slot holds the
unconstrained ``Z_t``. [`Renewal`](@ref) maps it to the positive rate
``\iota_t`` through its `transformation` ``g`` (default `exp`) before adding it,
so importation is positive by construction and *any* latent process can drive
it — a bare `Distribution` gives a constant rate shared across time, and a
process (e.g. a [`RandomWalk`](@ref)) gives a time-varying ``\iota_t`` that
stays positive even as the underlying path crosses zero.

The rate is added to the committed incidence after the modifier chain, so it
composes with modifiers such as [`SusceptibleDepletion`](@ref) without being
scaled by them. Because it is added *after* that chain, imported infections do
not themselves deplete the susceptible pool.

The resulting incidence is floored at a small positive value, keeping the
renewal recursion well defined for observation models that require a positive
mean.

At most one `ImportedCases` may be composed onto a step — a second one would
never be sampled, so [`Renewal`](@ref) rejects it at construction.

# Arguments

  - `importation_rate`: the unconstrained importation-rate prior — a
    `Distribution`, a `Vector{<:Distribution}`, or any prior/latent process.

# Examples

A constant rate with median ``\exp(-1) \approx 0.37`` imports per unit time:

```@example ImportedCases
using ComposableTuringIDModels, Distributions
Renewal([0.2, 0.3, 0.5], ImportedCases(Normal(-1.0, 0.5));
    rt = RandomWalk(), initialisation = Normal())
```

A time-varying rate, positive at every time despite the walk being
unconstrained:

```@example ImportedCases
Renewal([0.2, 0.3, 0.5], ImportedCases(RandomWalk());
    rt = RandomWalk(), initialisation = Normal())
```
"
struct ImportedCases{I <: PriorLike} <: AbstractRenewalModifier
    importation_rate::I
    function ImportedCases(importation_rate::PriorLike)
        return new{typeof(importation_rate)}(importation_rate)
    end
end

# The modifier carries no state of its own and leaves the incidence untouched:
# the sampled rate is added by `RenewalStep`'s tuple-input method (see above).
modifier_init_state(::ImportedCases) = 0.0
apply_modifier(::ImportedCases, inc, _) = (inc, 0.0)
