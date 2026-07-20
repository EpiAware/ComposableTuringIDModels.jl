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

The `model` slot is a length-`n` PATH slot: a bare `Distribution` there is
auto-wrapped in an [`Intercept`](@ref), giving a constant inner path; a process,
an [`IID`](@ref), or a vector passes through. Use [`IID`](@ref) for `n`
independent draws. It is composed through [`as_turing_submodel`](@ref).

## Fields

  - `model`: the latent model whose expected latent vector is recorded.
"
struct RecordExpectedLatent{M <: PriorLike} <: AbstractLatentModel
    "The latent model whose expected latent vector is recorded."
    model::M

    function RecordExpectedLatent(model)
        # `model` is a length-`n` PATH slot: a bare `Distribution` is wrapped in
        # an `Intercept` (a constant inner path), never left as a scalar.
        wrapped = _path_prior(model)
        new{typeof(wrapped)}(wrapped)
    end
end

@model function as_turing_model(model::RecordExpectedLatent, n)
    latent ~ as_turing_submodel(model.model, n)
    exp_latent := latent
    return latent
end
