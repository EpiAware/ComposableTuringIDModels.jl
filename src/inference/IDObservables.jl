# `IDObservables` container and `generated_observables` wrapper.

@doc raw"
Container for the outputs of an inference run: the model, the data, the posterior
samples, and any generated quantities.

## Fields

  - `model`: the model that was sampled.
  - `data`: the data the model was conditioned on.
  - `samples`: the posterior samples (or optimiser result).
  - `generated`: generated quantities, or `missing` if not computed.
"
struct IDObservables{M, D, S, G}
    "The model that was sampled."
    model::M
    "The data the model was conditioned on."
    data::D
    "The posterior samples (or optimiser result)."
    samples::S
    "Generated quantities, or `missing`."
    generated::G
end

@doc raw"
Wrap a model, data, and inference solution into an [`IDObservables`](@ref).

# Arguments

  - `model`: the model that was sampled.
  - `data`: the data the model was conditioned on.
  - `solution`: the inference solution (samples or optimiser result).

# Examples
```@example generated_observables
using ComposableTuringIDModels, Distributions
m = as_turing_model(
    IDModel(
        DirectInfections(; Z = RandomWalk(), initialisation_prior = Normal()),
        PoissonError()), missing, 10)
generated_observables(m, (; y_t = missing), rand(m))
```
"
function generated_observables(model, data, solution)
    return IDObservables(model, data, solution, missing)
end
