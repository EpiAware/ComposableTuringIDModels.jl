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

The `model` slot takes a raw component: a latent model, or a `Distribution` (or a
vector of them), composed through [`as_turing_submodel`](@ref).

## Fields

  - `model`: the latent model whose expected latent vector is recorded.
"
struct RecordExpectedLatent{M <: PriorLike} <: AbstractLatentModel
    "The latent model whose expected latent vector is recorded."
    model::M
end

@model function as_turing_model(model::RecordExpectedLatent, n)
    latent ~ as_turing_submodel(model.model, n)
    exp_latent := latent
    return latent
end
