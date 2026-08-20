# Record-the-expected-latent modifier.

@doc raw"
Record the inner latent vector as a tracked generated quantity (`exp_latent`).

# Arguments

  - `model`: the inner latent model whose output is recorded.
  - `n`: the shape to generate — a length or an `(n_strata, n_time)` shape,
    whatever the inner model accepts.

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
        wrapped = path_prior(model)
        return new{typeof(wrapped)}(wrapped)
    end
end

# `RecordExpectedLatent` wraps whatever its inner model returns, so it serves
# both a length-`n` path and an `(n_strata, n_time)` shape. The two methods
# below delegate to one shared `@model`, which avoids both the duplication and
# a dispatch ambiguity against the `AbstractPriorModel`/`Dims{2}` guard.
@model function _record_expected_latent(model::RecordExpectedLatent, n)
    latent ~ as_turing_submodel(model.model, n)
    exp_latent := latent
    return latent
end
function as_turing_model(model::RecordExpectedLatent, n::Int)
    return _record_expected_latent(model, n)
end
function as_turing_model(model::RecordExpectedLatent, n::Dims{2})
    return _record_expected_latent(model, n)
end
