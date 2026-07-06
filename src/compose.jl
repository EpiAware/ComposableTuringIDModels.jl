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
```@example EpiAwareModel
using EpiAwarePrototype, Distributions
model = EpiAwareModel(
    DirectInfections(; Z = RandomWalk(), initialisation_prior = Normal()),
    PoissonError())
mdl = as_turing_model(model, missing, 20)
rand(mdl)
```
"
struct EpiAwareModel{I <: AbstractInfectionModel, O <: AbstractObservationModel} <:
       AbstractEpiAwareModel
    "Infection process model generating ``I_t`` (and its internal latent ``Z_t``)."
    infection_model::I
    "Observation model mapping ``I_t`` to ``y_t``."
    observation_model::O
end

@model function as_turing_model(model::EpiAwareModel, y_t, n)
    infections ~ to_submodel(as_turing_model(model.infection_model, n), false)
    I_t = infections.I_t
    Z_t = infections.Z_t
    obs ~ to_submodel(
        as_turing_model(model.observation_model, y_t, I_t), false)
    # Uniform observation contract: sampled series and pre-error expected.
    generated_y_t = obs.y_t
    expected_y_t = obs.expected
    return (; generated_y_t, expected_y_t, I_t, Z_t)
end
