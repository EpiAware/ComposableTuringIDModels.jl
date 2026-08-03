# The headline composition: a full infections → observations model assembled from
# two components, each itself an `as_turing_model`. The latent (parameter) process
# is owned by the infection model, not a separate top-level slot.

@doc raw"
A composed epidemiological model linking an infection process and an observation
model.

The infection process owns its own latent (parameter) process internally — it is
no longer a separate top-level component — so a composed model is just two parts:
infections, then observations.

Sampling [`as_turing_model(model, y_t, n)`](@ref as_turing_model) runs the two
stages as submodels:

```math
I_t \;\xrightarrow{\text{infections}}\;
y_t \;\xrightarrow{\text{observations}}\; \text{data}
```

The returned generated quantities are `(; generated_y_t, expected_y_t, I_t, Z_t)`.
`generated_y_t` is the observation model's sampled `y_t` (the observed-or-simulated
series, or a `NamedTuple` of streams for a [`Split`](@ref)); `expected_y_t` is its
pre-error `expected` series (the uniform observation return contract). `Z_t` is the
infection model's internal latent draw (e.g. the (log) ``R_t`` path), kept
accessible as a generated quantity even though it is no longer a top-level
component — or `nothing` for infection models with no exposable latent (e.g.
[`ODEProcess`](@ref)). Pass `y_t = missing` to simulate from the prior, or a data
vector to condition.

## Fields

  - `infection_model`: the infection process model generating ``I_t`` (and its
    internal latent ``Z_t``).
  - `observation_model`: the observation model mapping ``I_t`` to ``y_t``.

# Examples
```@example IDModel
using ComposableTuringIDModels, Distributions
model = IDModel(
    DirectInfections(; Z = RandomWalk(), initialisation = Normal()),
    PoissonError())
mdl = as_turing_model(model, missing, 20)
rand(mdl)
```
"
struct IDModel{I <: AbstractInfectionModel, O <: AbstractObservationModel} <:
       AbstractComposableModel
    "Infection process model generating ``I_t`` (and its internal latent ``Z_t``)."
    infection_model::I
    "Observation model mapping ``I_t`` to ``y_t``."
    observation_model::O
end

@doc raw"
Build a grouped/panel `IDModel`: one shared infection process observed by
several groups, each reporting it at its own partially pooled level.

This is the shared constructor issue #180 asked for — the same `IDModel(...)`
call that builds a single-group model also builds a panel, distinguished only
by passing a `group_effect` prior between the infection and observation
arguments. It composes [`GroupedInfections`](@ref) (the shared curve replicated
per group) with `Split(observation_model)` (every group observed through the
same observation model, namespaced by group):

```julia
IDModel(infection_model, group_effect, observation_model) ==
    IDModel(GroupedInfections(infection_model, group_effect),
        Split(observation_model))
```

Sample it with `as_turing_model(model, Y, (n_time, n_groups))`, or the
convenience `as_turing_model(model, Y)` that reads `n_time` and `n_groups`
from the shape of the data matrix `Y` (`n_groups x n_time`, one row per group).

# Examples
```@example IDModel_grouped
using ComposableTuringIDModels, Distributions
model = IDModel(
    DirectInfections(; Z = RandomWalk(), initialisation = Normal(log(50.0), 0.2)),
    Hierarchy(; mean = Normal(0.0, 0.5), across = IID(Normal(0.0, 0.5))),
    PoissonError())
Ymiss = Matrix{Union{Missing, Float64}}(missing, 3, 12)   # 3 groups, 12 time steps
sim = as_turing_model(model, Ymiss)()
size(sim.I_t)
```

## Arguments

  - `infection_model`: the shared infection process, drawn once for all groups.
  - `group_effect`: the prior over the grouping axis (a [`Hierarchy`](@ref), a
    latent process, or a bare `Distribution`).
  - `observation_model`: the observation model each group is observed through.

## Keyword Arguments

  - `combiner`: the function `(I_t, level_g)` mapping the shared curve and a
    group's effect to that group's row (default: multiplicative on the
    exponential scale).
"
function IDModel(infection_model::AbstractInfectionModel, group_effect,
        observation_model::AbstractObservationModel;
        combiner = _grouped_combiner)
    return IDModel(GroupedInfections(infection_model, group_effect; combiner),
        Split(observation_model))
end

@doc raw"
Lift an existing `IDModel` (a shared infection + observation) to a grouped/panel
model by adding a per-group effect over the grouping axis.

Equivalent to `IDModel(idmodel.infection_model, group_effect,
idmodel.observation_model; combiner)`; see the other [`IDModel`](@ref) grouped
constructor for the full grouped-model behaviour.
"
function IDModel(idmodel::IDModel, group_effect; combiner = _grouped_combiner)
    return IDModel(idmodel.infection_model, group_effect,
        idmodel.observation_model; combiner)
end

@model function as_turing_model(model::IDModel, y_t, n)
    infections ~ as_turing_submodel(model.infection_model, n)
    I_t = infections.I_t
    Z_t = infections.Z_t
    obs ~ as_turing_submodel(model.observation_model, y_t, I_t)
    # Uniform observation contract: sampled series and pre-error expected.
    generated_y_t = obs.y_t
    expected_y_t = obs.expected
    return (; generated_y_t, expected_y_t, I_t, Z_t)
end

@doc raw"
Convenience 2-argument form for a grouped `IDModel` (one built with
[`GroupedInfections`](@ref)): reads `n_time` and `n_groups` from the shape of
the data matrix `Y` (`n_groups x n_time`, one row per group) instead of
requiring them to be supplied explicitly.
"
function as_turing_model(model::IDModel{<:GroupedInfections}, Y::AbstractMatrix)
    n_groups, n_time = size(Y)
    return as_turing_model(model, Y, (n_time = n_time, n_groups = n_groups))
end
