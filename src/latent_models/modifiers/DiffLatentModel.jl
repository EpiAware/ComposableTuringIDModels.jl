# Differenced latent process modifier.

@doc raw"
Model a latent process as a `d`-fold differenced version of an inner process.

If ``\tilde Z_t`` is the inner (undifferenced) latent path supplied via `model`,
then

```math
\Delta^{(d)} Z_t = \tilde Z_t,
```

and ``Z_t`` is recovered by applying `cumsum` `d` times. The `d` initial terms
are inferred from the prior in `init`; `d` equals the length of that prior.

The `init` slot is an [`AbstractPriorModel`](@ref): pass a vector of
`Distribution`s (coerced via [`as_prior`](@ref)) as before, or a richer prior
model.

Composing `DiffLatentModel` over an `AR` gives an ARIMA-style latent process.

# Examples
```@example DiffLatentModel
using ComposableTuringIDModels, Distributions
diff = DiffLatentModel(; model = RandomWalk(), init = [Normal(), Normal()])
mdl = as_turing_model(diff, 10)
rand(mdl)
```
"
struct DiffLatentModel{M <: AbstractLatentModel, P <: AbstractPriorModel} <:
       AbstractLatentModel
    "Underlying (undifferenced) latent model."
    model::M
    "Prior for the initial latent variables."
    init::P
    "Number of times differenced."
    d::Int

    function DiffLatentModel(
            model::AbstractLatentModel, init::AbstractPriorModel, d::Int)
        @assert d>0 "d must be greater than 0"
        _assert_prior_length(init, d, "init")
        new{typeof(model), typeof(init)}(model, init, d)
    end
end

function DiffLatentModel(model, init::Distribution; d::Int)
    return DiffLatentModel(; model = model, init = fill(init, d))
end

function DiffLatentModel(; model, init = [Normal()])
    init_prior = as_prior(init)
    d = _prior_order(init_prior)
    return DiffLatentModel(as_prior(model), init_prior, d)
end

@model function as_turing_model(model::DiffLatentModel, n)
    d = model.d
    @assert n>d "n must be longer than d"
    latent_init ~ to_submodel(as_turing_model(model.init, d))
    diff_latent ~ to_submodel(as_turing_model(model.model, n - d))
    return _combine_diff(latent_init, diff_latent, d)
end

function _combine_diff(init, diff, d)
    combined = vcat(init, diff)
    for _ in 1:d
        combined = cumsum(combined)
    end
    return combined
end
