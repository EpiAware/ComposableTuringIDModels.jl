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

@doc raw"
Imported-cases modifier for [`RenewalStep`](@ref).

Adds an externally seeded importation rate to the renewal incidence,

```math
I_t = \iota_t + R_t \sum_{i=1}^{n-1} I_{t-i} g_i,
```

so infections can arrive from outside the modelled population — the mechanism
behind a renewal process that would otherwise die out from a zero initial
incidence, and behind reintroduction after local elimination.

`importation_rate` is a length-`n` prior slot: a bare `Distribution` gives a
constant rate ``\iota`` shared across time, and a process (e.g. a
[`RandomWalk`](@ref)) gives a time-varying ``\iota_t``. The rate is sampled by
[`Renewal`](@ref) and added to the committed incidence after the modifier
chain, so it composes with modifiers such as [`SusceptibleDepletion`](@ref)
without being scaled by them.

The resulting incidence is floored at a small positive value, keeping the
renewal recursion well defined for observation models that require a positive
mean.

# Arguments

  - `importation_rate`: the importation-rate prior — a `Distribution`, a
    `Vector{<:Distribution}`, or any prior/latent process.

# Examples
```@example
using ComposableTuringIDModels, Distributions
Renewal([0.2, 0.3, 0.5], ImportedCases(Exponential(1.0));
    rt = RandomWalk(), initialisation = Normal())
```
"
struct ImportedCases{I} <: AbstractRenewalModifier
    importation_rate::I
    function ImportedCases(importation_rate)
        return new{typeof(importation_rate)}(importation_rate)
    end
end

# The modifier carries no state of its own and leaves the incidence untouched:
# the sampled rate is added by `RenewalStep`'s tuple-input method (see above).
modifier_init_state(::ImportedCases) = 0.0
apply_modifier(::ImportedCases, inc, _) = (inc, 0.0)
