# The headline composition: a full latent → infections → observations model
# assembled from three components, each itself an `as_turing_model`.

@doc raw"
A composed epidemiological model linking a latent process, an infection process,
and an observation model.

Sampling [`as_turing_model(model, y_t, n)`](@ref as_turing_model) runs the three
stages as submodels:

```math
Z_t \;\xrightarrow{\text{latent}}\; I_t \;\xrightarrow{\text{infections}}\;
y_t \;\xrightarrow{\text{observations}}\; \text{data}
```

The returned generated quantities are `(; generated_y_t, I_t, Z_t)`. Pass
`y_t = missing` to simulate from the prior, or a data vector to condition.

## Fields

  - `latent_model`: the latent process model generating ``Z_t``.
  - `epi_model`: the infection process model mapping ``Z_t`` to ``I_t``.
  - `observation_model`: the observation model mapping ``I_t`` to ``y_t``.

# Examples
```@example EpiAwareModel
using EpiAwarePrototype, Distributions
data = EpiData([0.2, 0.3, 0.5], exp)
model = EpiAwareModel(
    RandomWalk(),
    DirectInfections(; data = data, initialisation_prior = Normal()),
    PoissonError())
mdl = as_turing_model(model, missing, 20)
rand(mdl)
```
"
struct EpiAwareModel{
    L <: AbstractLatentModel, I <: AbstractInfectionModel,
    O <: AbstractObservationModel} <: AbstractEpiAwareModel
    "Latent process model generating ``Z_t``."
    latent_model::L
    "Infection process model mapping ``Z_t`` to ``I_t``."
    epi_model::I
    "Observation model mapping ``I_t`` to ``y_t``."
    observation_model::O
end

@model function as_turing_model(model::EpiAwareModel, y_t, n)
    Z_t ~ to_submodel(as_turing_model(model.latent_model, n), false)
    I_t ~ to_submodel(as_turing_model(model.epi_model, Z_t), false)
    generated_y_t ~ to_submodel(
        as_turing_model(model.observation_model, y_t, I_t), false)
    return (; generated_y_t, I_t, Z_t)
end
