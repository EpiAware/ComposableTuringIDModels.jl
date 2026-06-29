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

An infection model maps a latent path `Z_t` to a path of unobserved infections
`I_t`. Its role interface is

```julia
as_turing_model(model::AbstractInfectionModel, Z_t)  # ⇒ an infection path I_t
```

Members include [`DirectInfections`](@ref), [`ExpGrowthRate`](@ref),
[`Renewal`](@ref) and [`ODEProcess`](@ref).
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
