# Per-role interface contracts. Each role's members must implement the
# `as_turing_model` signature fixed by its supertype:
#
#   latent       as_turing_model(m, n)        ⇒ a DynamicPPL.Model
#   infection    as_turing_model(m, n)        ⇒ a DynamicPPL.Model ⇒ (; I_t, Z_t)
#   observation  as_turing_model(m, y_t, Y_t) ⇒ a DynamicPPL.Model
#
# The helpers below encode those contracts in a form usable from tests: each
# returns `true` when `model` is in the given role *and* its role-specific
# `as_turing_model` call returns a `DynamicPPL.Model`. They construct the model
# but do not sample it, so they are cheap conformance checks.

@doc raw"
Deprecated alias for [`implements_prior_interface`](@ref).

The latent and prior roles are one ([`AbstractLatentModel`](@ref) `===`
[`AbstractPriorModel`](@ref)), so this simply forwards to
[`implements_prior_interface`](@ref).

# Examples
```@example
using ComposableTuringIDModels
implements_latent_interface(RandomWalk())
```
"
implements_latent_interface(model; kwargs...) = implements_prior_interface(model; kwargs...)

@doc raw"
Check that `model` satisfies the [`AbstractPriorModel`](@ref) interface: it is a
prior model and `as_turing_model(model, n)` returns a `DynamicPPL.Model`.

Every [`AbstractLatentModel`](@ref) is also an `AbstractPriorModel`, so this holds
for latent models (a latent process used directly as a prior) as well as for the
[`BroadcastPrior`](@ref) wrapper and any bespoke prior submodel.

# Arguments

  - `model`: the component to check.

# Keyword Arguments

  - `n`: the prior length used for the construction check (default `10`).

# Examples
```@example
using ComposableTuringIDModels, Distributions
implements_prior_interface(BroadcastPrior(Normal()))
```
"
function implements_prior_interface(model; n::Int = 10)
    model isa AbstractPriorModel || return false
    return as_turing_model(model, n) isa DynamicPPL.Model
end

@doc raw"
Check that `model` satisfies the [`AbstractInfectionModel`](@ref) interface: it is
an infection model and `as_turing_model(model, n)` returns a `DynamicPPL.Model`.

The infection model owns its latent process internally, so the construction check
passes only a series length `n` (no external latent path).

# Arguments

  - `model`: the component to check.

# Keyword Arguments

  - `n`: the infection series length used for the construction check (default
    `10`).

# Examples
```@example
using ComposableTuringIDModels, Distributions
implements_infection_interface(
    DirectInfections(; Z = RandomWalk(), initialisation = Normal()))
```
"
function implements_infection_interface(model; n::Int = 10)
    model isa AbstractInfectionModel || return false
    return as_turing_model(model, n) isa DynamicPPL.Model
end

@doc raw"
Check that `model` satisfies the [`AbstractObservationModel`](@ref) interface: it
is an observation model and `as_turing_model(model, y_t, Y_t)` returns a
`DynamicPPL.Model`.

# Arguments

  - `model`: the component to check.

# Keyword Arguments

  - `y_t`: the observed series for the construction check (default `missing`,
    i.e. predictive simulation).
  - `Y_t`: the expected-observation series (default `fill(10.0, 10)`).

# Examples
```@example
using ComposableTuringIDModels
implements_observation_interface(PoissonError())
```
"
function implements_observation_interface(model; y_t = missing, Y_t = fill(10.0, 10))
    model isa AbstractObservationModel || return false
    return as_turing_model(model, y_t, Y_t) isa DynamicPPL.Model
end
