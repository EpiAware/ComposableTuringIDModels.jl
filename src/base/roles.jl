# Role supertypes: the shallow layer beneath `AbstractComposableModel` that encodes
# what role a component plays (latent / infection / observation). Each role fixes
# the `as_turing_model` signature its members must implement; the composer and the
# manipulators dispatch and constrain on these so a wrong-role component fails at
# construction. This is the only type structure beyond the single root supertype —
# there is no deeper hierarchy and no per-concept `generate_*` functions.

@doc raw"
Supertype for **prior** models — a parameter prior expressed as a length-`n`
submodel rather than a bare `Distribution`.

A prior model maps a length `n` to a length-`n` vector of parameter values via the
same `as_turing_model` protocol every other component speaks:

```julia
as_turing_model(prior::AbstractPriorModel, n)  # ⇒ a length-`n` vector
```

A raw `Distribution` (or vector of them) is *not* a prior model but flows through
the same [`as_turing_submodel`](@ref) seam: `as_turing_model` has `Distribution`
and `Vector{<:Distribution}` methods, so a bare distribution composes as a
length-`n` prior submodel exactly like a model does. This is the single role for
every parameter *process*: a latent process (a `RandomWalk` for a time-varying
parameter, an [`AR`](@ref) process, …) satisfies the same
`as_turing_model(m, n) ⇒ length-n` contract, so it drops into any prior slot
directly. The former `AbstractLatentModel` role has been folded into this one (it
survives as a deprecated alias). A genuinely scalar parameter is drawn with a
native tilde (`σ ~ model.std`), keeping the chain as small as a bare `~ dist`.

This delivers issue #37 (priors as length-`n` submodels).
"
abstract type AbstractPriorModel <: AbstractComposableModel end

@doc raw"
Deprecated alias for [`AbstractPriorModel`](@ref).

A latent process and a parameter prior share one role — both map a length `n` to a
length-`n` vector via `as_turing_model(m, n)` — so the separate
`AbstractLatentModel` type has been collapsed into [`AbstractPriorModel`](@ref).
`AbstractLatentModel` remains as a `const` alias for one release for backwards
compatibility (`AbstractLatentModel === AbstractPriorModel`); prefer
[`AbstractPriorModel`](@ref) in new code.
"
const AbstractLatentModel = AbstractPriorModel

@doc raw"
Supertype for **infection process** models.

An infection model maps a series length `n` to a path of unobserved infections
`I_t`. It owns its own latent (parameter) process internally — generating, e.g.,
a (log) reproduction number or growth-rate path — and maps that to infections, so
no external latent path is threaded in. Its role interface is

```julia
as_turing_model(model::AbstractInfectionModel, n)  # ⇒ (; I_t, Z_t)
```

where the returned named tuple carries the infection path `I_t` and the model's
internal latent draw `Z_t` (the (log) ``R_t`` / growth-rate path, or `nothing`
for models with no exposable latent such as [`ODEProcess`](@ref)). Exposing
`Z_t` keeps the latent recoverable as a generated quantity downstream.

Members include [`DirectInfections`](@ref), [`ExpGrowthRate`](@ref),
[`Renewal`](@ref) and [`ODEProcess`](@ref). Only [`Renewal`](@ref) carries a
generation interval ([`IDData`](@ref)); the others take a `transformation`
directly.
"
abstract type AbstractInfectionModel <: AbstractComposableModel end

@doc raw"
Supertype for **observation** models.

An observation model maps a path of expected observations `Y_t` to observed
counts `y_t`. Its role interface is

```julia
as_turing_model(model::AbstractObservationModel, y_t, Y_t)  # ⇒ observed counts y_t
```

`y_t === missing` triggers prior/predictive simulation; a concrete `y_t`
conditions the model on data. Observation *modifiers* (e.g. [`LatentDelay`](@ref),
[`Ascertainment`](@ref), [`Aggregate`](@ref)) are themselves
`AbstractObservationModel`s: wrapping an observation model yields another
observation model. Their inner-model slots are typed `AbstractObservationModel`,
so only observation components can be wrapped.

[`AbstractObservationErrorModel`](@ref) is the sub-role for the simple
error families (Poisson, negative binomial).
"
abstract type AbstractObservationModel <: AbstractComposableModel end
