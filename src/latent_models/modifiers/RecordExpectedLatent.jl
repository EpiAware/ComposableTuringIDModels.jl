# Record-the-expected-latent modifier.

@doc raw"
Record the inner latent vector as a tracked generated quantity (`exp_latent`).

# Arguments

  - `model`: the inner latent model whose output is recorded.
  - `n`: the length of the latent series to generate.

# Examples
```@example RecordExpectedLatent
using ComposableTuringIDModels
rm = RecordExpectedLatent(FixedIntercept(0.1))
rand(as_turing_model(rm, 1))
```

The `model` slot is an [`AbstractPriorModel`](@ref): a bare `Distribution` (or a
vector of them) is coerced via [`as_prior`](@ref), as at the top-level slots.

## Fields

  - `model`: the latent model whose expected latent vector is recorded.
"
struct RecordExpectedLatent{M <: AbstractLatentModel} <: AbstractLatentModel
    "The latent model whose expected latent vector is recorded."
    model::M
end

# Coerce a bare `Distribution` (or vector) member to the prior interface so it is
# accepted alongside a process, matching the top-level slots and Combine/Concat.
RecordExpectedLatent(model) = RecordExpectedLatent(as_prior(model))

@model function as_turing_model(model::RecordExpectedLatent, n)
    latent ~ to_submodel(as_turing_model(model.model, n))
    exp_latent := latent
    return latent
end
