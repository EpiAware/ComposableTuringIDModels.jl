# `EpiAwareObservables` container and `generated_observables` wrapper.

@doc raw"
Container for the outputs of an inference run: the model, the data, the posterior
samples, and any generated quantities.

## Fields

  - `model`: the model that was sampled.
  - `data`: the data the model was conditioned on.
  - `samples`: the posterior samples (or optimiser result).
  - `generated`: generated quantities, or `missing` if not computed.
"
struct EpiAwareObservables{M, D, S, G}
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
Wrap a model, data, and inference solution into an [`EpiAwareObservables`](@ref).

When `solution` is an MCMC `Chains`, the model is re-run over the draws with
`DynamicPPL.returned` to recover the model's returned generated quantities (e.g.
`(; generated_y_t, I_t, Z_t)`) per sample, stored in the `generated` field. For
any other solution (an optimiser result, a prior draw, …) there are no per-draw
generated quantities, so `generated` is `missing`.

# Arguments

  - `model`: the model that was sampled (a conditioned `DynamicPPL.Model`).
  - `data`: the data the model was conditioned on.
  - `solution`: the inference solution (samples or optimiser result).

# Examples
```@example generated_observables
using EpiAwarePrototype, Distributions
m = as_turing_model(
    EpiAwareModel(
        DirectInfections(; Z = RandomWalk(), initialisation_prior = Normal()),
        PoissonError()), missing, 10)
generated_observables(m, (; y_t = missing), rand(m))
```
"
function generated_observables(model, data, solution)
    return EpiAwareObservables(model, data, solution, _generated_quantities(model, solution))
end

# Re-run the (conditioned) model over the posterior draws to recover its returned
# generated quantities per draw. `solution` may be an `MCMCChains.Chains` or a
# FlexiChains chain (Turing's sampler output), so we don't dispatch on a single
# concrete chain type — we ask `returned` to consume it and treat anything it
# cannot (an optimiser result, a prior draw, …) as "no generated quantities".
_generated_quantities(model, solution) = missing
# A single prior/predictive draw (a `NamedTuple` of parameter values, e.g. from
# `rand(model)`) is not a chain of posterior draws, so there are no per-draw
# generated quantities to collect.
_generated_quantities(model::DynamicPPL.Model, solution::NamedTuple) = missing
function _generated_quantities(model::DynamicPPL.Model, solution)
    return try
        returned(model, solution)
    catch
        missing
    end
end
