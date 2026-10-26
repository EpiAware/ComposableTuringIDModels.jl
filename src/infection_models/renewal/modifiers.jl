# Renewal modifiers: the model-level description of what a renewal does to
# each step's new infections, resolved before the scan into the
# ComposableRecurrences modifiers the recurrence runs.

@doc raw"
Abstract supertype for renewal modifiers composed onto a [`Renewal`](@ref).

A renewal modifier transforms each step's new infections after the
generation-interval convolution and ``R_t``, and may carry state of its own
(a susceptible pool, say).
The recursion itself is run by ComposableRecurrences.jl, so a modifier is
described here and **resolved** before the scan into a ComposableRecurrences
modifier, which is what the recurrence applies at each step.

A modifier implements one method,

```julia
as_turing_model(mod::MyModifier, n)  # ⇒ a DynamicPPL.Model returning the
                                     #   ComposableRecurrences modifier
```

where `n` is the renewal's [`ModelShape`](@ref).
A deterministic modifier samples nothing and returns its ComposableRecurrences
form, as [`SusceptibleDepletion`](@ref) returns a floored `Depletion`.
A modifier carrying priors draws them here, through
[`as_turing_submodel`](@ref), and returns a modifier holding the drawn values,
as [`ImportedCases`](@ref) returns an `Add` of the drawn importation rate.
A new kind of step-level behaviour is a ComposableRecurrences modifier (a type
with a `ComposableRecurrences.forward!` method for the `Step` role), wrapped
by a renewal modifier that returns it.

Modifiers apply in the order given.
Each modifier's sampled variables are prefixed by its **position** in the
modifier tuple, `modifier_<i>`, so two modifiers of the same kind cannot
collide on a variable name — in
`Renewal(gen_int, SusceptibleDepletion(N), ImportedCases(Normal()))` the
importation rate is `modifier_2.import_rates`. Inserting or reordering a
modifier therefore renames the variables of every modifier after it.
"
abstract type AbstractRenewalModifier end

@doc raw"
Susceptible-depletion modifier for a [`Renewal`](@ref).

Scales the proposed incidence by the available susceptible fraction and depletes
the susceptible pool,

```math
I_t = \max\left(\frac{S_{t-1}}{N}, 10^{-6}\right) \, R_t \sum_{i=1}^{n-1} I_{t-i} g_i,
\qquad S_t = S_{t-1} - I_t,
```

with population size ``N`` = `pop_size` and ``S_0 = N``.
The fraction is floored because a large force of infection can take more than
the pool holds, leaving ``S`` negative for the rest of the run.
Adding it to a renewal gives a renewal process with a fixed population and
susceptible depletion, e.g. `Renewal(gen_int, SusceptibleDepletion(N))`.

`pop_size` is a scalar for one series, and a per-stratum vector for a
stratified renewal.
Each stratum then depletes its own pool.
A scalar given to a stratified renewal is shared by every stratum, so each
depletes a separate pool of the same size.

It samples nothing, and resolves to ComposableRecurrences'
`Depletion(N; form = Floor())`.

## Fields

  - `pop_size`: the population size, one value or one per stratum.
"
struct SusceptibleDepletion{T} <: AbstractRenewalModifier
    "The population size, one value or one per stratum."
    pop_size::T
end

# A modifier with no resolution would otherwise fail with a bare `MethodError`
# inside the pre-scan seam.
function as_turing_model(mod::AbstractRenewalModifier, n)
    throw(
        ArgumentError(
            "$(typeof(mod)) must implement `as_turing_model(mod, n)` " *
                "returning a ComposableRecurrences modifier"
        )
    )
end

@doc raw"
Resolve a [`SusceptibleDepletion`](@ref) into the floored
ComposableRecurrences `Depletion` the recurrence applies.

# Arguments

  - `mod`: the [`SusceptibleDepletion`](@ref) modifier.
  - `n`: the renewal's [`ModelShape`](@ref).
"
@model function as_turing_model(mod::SusceptibleDepletion, n)
    return CR.Depletion(_strata_param(mod.pop_size, n); form = CR.Floor())
end

# A population is one value shared by every stratum or one value per stratum.
_strata_param(x::Real, n) = x
function _strata_param(x::AbstractVector, n)
    _check_strata(x, n)
    return PerStratum(x)
end

_check_strata(x, ::Int) = throw(
    ArgumentError(
        "a per-stratum `pop_size` needs a stratified renewal, but the " *
            "renewal has one series"
    )
)
function _check_strata(x, n::Dims{2})
    length(x) == n[1] || throw(
        DimensionMismatch(
            "`pop_size` has $(length(x)) entries but the renewal has " *
                "$(n[1]) strata"
        )
    )
    return nothing
end

# --- the pre-scan seam ------------------------------------------------------
#
# A scan is a deterministic function, so a modifier that needs sampled
# parameters cannot draw them inside it.
# The modifiers are resolved once beforehand, in order.

# The recursion is structural, over the tuple rather than an index into it, so
# each level has a concrete tuple type.
# Each modifier is prefixed by its position, so two modifiers of the same kind
# cannot collide on a variable name.
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
