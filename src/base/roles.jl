# Role supertypes: the shallow layer beneath `AbstractEpiAwareModel` that encodes
# what role a component plays (latent / infection / observation). Each role fixes
# the `as_turing_model` signature its members must implement; the composer and the
# manipulators dispatch and constrain on these so a wrong-role component fails at
# construction. This is the only type structure beyond the single root supertype —
# there is no deeper hierarchy and no per-concept `generate_*` functions.

@doc raw"
Supertype for **latent process** models.

A latent model maps a series length `n` to a length-`n` latent path. Its role
interface is

```julia
as_turing_model(model::AbstractLatentModel, n)  # ⇒ a length-`n` latent path
```

Latent *modifiers* and *manipulators* (e.g. [`DiffLatentModel`](@ref),
[`CombineLatentModels`](@ref), [`BroadcastLatentModel`](@ref)) are themselves
`AbstractLatentModel`s: wrapping a latent model yields another latent model, so
they compose freely. Their inner-model slots are typed `AbstractLatentModel`, so
only latent components can be wrapped.
"
abstract type AbstractLatentModel <: AbstractEpiAwareModel end

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
generation interval ([`EpiData`](@ref)); the others take a `transformation`
directly.
"
abstract type AbstractInfectionModel <: AbstractEpiAwareModel end

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
abstract type AbstractObservationModel <: AbstractEpiAwareModel end
