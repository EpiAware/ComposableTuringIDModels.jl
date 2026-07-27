# Imported (externally seeded) cases as a renewal modifier (#189).
#
# `ImportedCases` is the worked example of a modifier that carries priors. The
# rate it adds is *sampled*, and a scan step is a deterministic function, so the
# rate is drawn once before the scan through the modifier seam
# (`as_turing_model(mod, n)`, see `AbstractRenewalModifier`). That returns an
# `ImportedRate` — the same modifier with the drawn path in place of the prior —
# which then behaves like any other modifier in the scan: it transforms the
# proposed incidence and carries its own substate.
#
# The substate is the step counter, which is how a per-time quantity is read
# inside a scan whose step sees no clock. Nothing about this is specific to
# importation: any modifier needing a per-time input resolves it the same way.
#
# The prior is on the UNCONSTRAINED scale, like `Renewal`'s `rt` and
# `initialisation` slots, and the modifier maps it through its own
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

`importation_rate` is a per-step parameter slot holding the **unconstrained**
rate ``\tilde\iota_t``, read at step ``t`` with [`_at`](@ref): a bare
`Distribution` is one unknown constant shared across time, a
`Vector{<:Distribution}` or a latent process (e.g. a [`RandomWalk`](@ref)) is a
length-`n` path. The modifier maps whatever it gets onto the positive rate
``\iota_t`` with its own `transformation` ``g`` (default `exp`), so importation
is positive by construction and *any* latent process can drive it — a
time-varying ``\iota_t`` stays positive even as the underlying path crosses
zero.

The rate is drawn before the scan through the modifier seam (see
[`AbstractRenewalModifier`](@ref)), giving an [`ImportedRate`](@ref) that adds
``\iota_t`` at step ``t``. Where it sits in the modifier tuple therefore decides
how it composes: placed *after* a [`SusceptibleDepletion`](@ref) the imports are
added to the depleted incidence, so they are not scaled by the susceptible
fraction and do not themselves deplete the pool; placed *before* it they are
treated as part of the incidence the pool depletes.

The drawn rate is named `import_rates`, prefixed by the modifier's **position**
in the renewal step's modifier tuple, so its posterior is read from a chain as

```julia
model = Renewal(gen_int, SusceptibleDepletion(N), ImportedCases(Normal());
    rt = RandomWalk(), initialisation = Normal())
# `modifier_2` because importation is the second modifier; `exp` puts the
# draws back on the scale of imports per unit time.
exp.(vec(chain[@varname(modifier_2.import_rates)]))
```

Inserting or reordering modifiers renames it.

The resulting incidence is floored at a small positive value, keeping the
renewal recursion well defined for observation models that require a positive
mean.

# Arguments

  - `importation_rate`: the unconstrained importation-rate prior — a
    `Distribution`, a `Vector{<:Distribution}`, or any prior/latent process.

# Keyword Arguments

  - `transformation`: the map from the unconstrained rate onto the positive
    importation rate (default `exp`).

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

## Fields

  - `importation_rate`: the unconstrained importation-rate prior.
  - `transformation`: the map onto the positive importation rate.
"
struct ImportedCases{I <: PriorLike, F <: Function} <: AbstractRenewalModifier
    "The unconstrained importation-rate prior."
    importation_rate::I
    "The map from the unconstrained rate onto the positive rate."
    transformation::F

    function ImportedCases(
            importation_rate::PriorLike; transformation::Function = exp)
        return new{typeof(importation_rate), typeof(transformation)}(
            importation_rate, transformation)
    end
end

@doc raw"
A resolved [`ImportedCases`](@ref): the drawn importation-rate path, ready to
scan.

This is what [`ImportedCases`](@ref) returns from its pre-scan
`as_turing_model` seam — the prior has been sampled and transformed, so the scan
sees a plain deterministic modifier. Its substate is the step counter, so step
``t`` adds ``\iota_t``, and the incidence is floored at a small positive value.

## Fields

  - `rate`: the positive importation rate at each time.
"
struct ImportedRate{V <: AbstractVector} <: AbstractRenewalModifier
    "The positive importation rate at each time."
    rate::V
end

# The substate is the step counter: a scan step has no clock of its own, so a
# per-time modifier carries the index it has reached.
modifier_init_state(::ImportedRate) = 0

function apply_modifier(mod::ImportedRate, incidence, t)
    return max(incidence + mod.rate[t + 1], 1e-6), t + 1
end

@doc raw"
Sample the importation rate ahead of the scan.

Draws the unconstrained rate slot through [`as_turing_submodel`](@ref) — a bare
`Distribution` giving one constant rate, a process giving a length-`n` path —
maps it onto the positive scale with the modifier's `transformation`, and
returns the [`ImportedRate`](@ref) the scan uses.

# Arguments

  - `mod`: the [`ImportedCases`](@ref) modifier.
  - `n`: the length of the renewal series.
"
@model function as_turing_model(mod::ImportedCases, n)
    import_rates ~ as_turing_submodel(mod.importation_rate, n; prefix = true)
    rate = [mod.transformation(_at(import_rates, t)) for t in 1:n]
    return ImportedRate(rate)
end
