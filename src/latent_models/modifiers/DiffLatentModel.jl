# Differenced latent process modifier.

@doc raw"
Model a latent process as a `d`-fold differenced version of an inner process.

If ``\tilde Z_t`` is the inner (undifferenced) latent path supplied via `model`,
then

```math
\Delta^{(d)} Z_t = \tilde Z_t,
```

and ``Z_t`` is recovered by applying `cumsum` `d` times. The `d` initial terms
are inferred from `init_prior`; `d` equals the length of `init_priors`.

Composing `DiffLatentModel` over an `AR` gives an ARIMA-style latent process.

# Examples
```@example DiffLatentModel
using EpiAwarePrototype, Distributions
diff = DiffLatentModel(; model = RandomWalk(), init_priors = [Normal(), Normal()])
mdl = as_turing_model(diff, 10)
rand(mdl)
```
"
struct DiffLatentModel{M <: AbstractEpiAwareModel, P <: Distribution} <:
       AbstractEpiAwareModel
    "Underlying (undifferenced) latent model."
    model::M
    "Prior distribution for the initial latent variables."
    init_prior::P
    "Number of times differenced."
    d::Int

    function DiffLatentModel(model::AbstractEpiAwareModel, init_prior::Distribution, d::Int)
        @assert d>0 "d must be greater than 0"
        @assert d==length(init_prior) "d must equal the length of init_prior"
        new{typeof(model), typeof(init_prior)}(model, init_prior, d)
    end
end

function DiffLatentModel(model::AbstractEpiAwareModel, init_prior::Distribution; d::Int)
    return DiffLatentModel(; model = model, init_priors = fill(init_prior, d))
end

function DiffLatentModel(; model::AbstractEpiAwareModel,
        init_priors::Vector{D} where {D <: Distribution} = [Normal()])
    d = length(init_priors)
    return DiffLatentModel(model, _expand_dist(init_priors), d)
end

@model function as_turing_model(model::DiffLatentModel, n)
    d = model.d
    @assert n>d "n must be longer than d"
    latent_init ~ model.init_prior
    diff_latent ~ to_submodel(as_turing_model(model.model, n - d), false)
    return _combine_diff(latent_init, diff_latent, d)
end

function _combine_diff(init, diff, d)
    combined = vcat(init, diff)
    for _ in 1:d
        combined = cumsum(combined)
    end
    return combined
end
