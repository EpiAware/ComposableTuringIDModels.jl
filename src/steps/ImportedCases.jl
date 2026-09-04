# A scan step is a deterministic function, so the rate is drawn once before the
# scan through the modifier seam (`as_turing_model(mod, n)`), which returns an
# `ImportedRate` holding the drawn path.
#
# The prior is on the UNCONSTRAINED scale and is mapped through the modifier's
# own `transformation` before being added.
# An importation rate cannot be negative, and a slot that accepts any latent
# process has to make that hold structurally, because an unconstrained process
# spends much of its time below zero.

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

For infection arriving from another *modelled* stratum, rather than from
outside the system altogether, use the renewal `mixing` slot instead (see
[`renewal_pressure`](@ref)).

On a stratified renewal the rate is read per step with [`at`](@ref) whatever
shape it has, and the two shapes mean different things.
A bare `Distribution` is one constant shared across strata, which lifts every
stratum by the same amount at each step.
A [`Stratify`](@ref) or [`Replicate`](@ref) rate is a `strata × time` rate, one
exogenous stream per stratum.

`importation_rate` is a per-step parameter slot holding the **unconstrained**
rate ``\tilde\iota_t``, read at step ``t`` with [`at`](@ref): a bare
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

Nothing clamps the incidence afterwards: the renewal recursion stays positive
because ``\iota_t`` is, which is `transformation`'s job. Passing
`transformation = identity` hands that job to you, and a rate that goes negative
then subtracts from the incidence.

# Arguments

  - `importation_rate`: the unconstrained importation-rate prior — a
    `Distribution`, a `Vector{<:Distribution}`, or any prior/latent process.

# Keyword Arguments

  - `transformation`: the map from the unconstrained rate onto the positive
    importation rate (default `exp`).

# Examples

Drawing the modifier through its pre-scan seam shows the rate the scan will
add. A bare `Distribution` gives one constant, here with median
``\exp(-1) \approx 0.37`` imports per unit time:

```@example ImportedCases
using ComposableTuringIDModels, Distributions, Random
Random.seed!(189)
as_turing_model(ImportedCases(Normal(-1.0, 0.5)), 5)().rate
```

A latent process gives a path, positive at every time despite the walk itself
being unconstrained:

```@example ImportedCases
Random.seed!(189)
as_turing_model(ImportedCases(RandomWalk()), 5)().rate
```

Either goes onto a renewal process as a positional modifier, where the rate
joins the model's parameters under the modifier's position:

```@example ImportedCases
r = Renewal([0.2, 0.3, 0.5], ImportedCases(Normal(-1.0, 0.5));
    rt = FixedIntercept(0.0), initialisation = Normal())
keys(rand(as_turing_model(r, 5)))
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
            importation_rate::PriorLike; transformation::Function = exp
        )
        return new{typeof(importation_rate), typeof(transformation)}(
            importation_rate, transformation
        )
    end
end

@doc raw"
A resolved [`ImportedCases`](@ref): the drawn importation rate, ready to scan.

This is what [`ImportedCases`](@ref) returns from its pre-scan
`as_turing_model` seam — the prior has been sampled and transformed, so the scan
sees a plain deterministic modifier. Its substate is the step counter, so step
``t`` adds ``\iota_t`` read with [`at`](@ref) (a constant rate stays a scalar).

## Fields

  - `rate`: the positive importation rate, constant or one value per time.
"
struct ImportedRate{V} <: AbstractRenewalModifier
    "The positive importation rate: one constant, or one value per time."
    rate::V
end

# A scan step has no clock of its own, so a per-time modifier carries the index
# it has reached as its substate.
# The window is unused, because a step counter has no shape to match.
modifier_init_state(::ImportedRate, window) = 0

# Broadcast, so a shared constant rate lifts every stratum of a stratified
# renewal rather than failing on a scalar added to a vector.
# A single series adds two scalars and allocates nothing.
function apply_modifier(mod::ImportedRate, incidence, t)
    return incidence .+ at(mod.rate, t + 1), t + 1
end

@doc raw"
Sample the importation rate ahead of the scan.

Draws the unconstrained rate slot through [`as_turing_submodel`](@ref) — a bare
`Distribution` giving one constant rate, a process giving a length-`n` path —
maps it onto the positive scale with the modifier's `transformation`, and
returns the [`ImportedRate`](@ref) the scan uses. The map is broadcast, so a
constant stays a scalar (no length-`n` allocation) and the scan reads either
shape with [`at`](@ref).

# Arguments

  - `mod`: the [`ImportedCases`](@ref) modifier.
  - `n`: the length of the renewal series.
"
@model function as_turing_model(mod::ImportedCases, n)
    import_rates ~ as_turing_submodel(mod.importation_rate, n; prefix = true)
    return ImportedRate(mod.transformation.(import_rates))
end
