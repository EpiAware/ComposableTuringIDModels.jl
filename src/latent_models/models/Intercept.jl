# Intercept latent process models (sampled and fixed).

@doc raw"
Broadcast a single sampled intercept value to a length-`n` latent process.

The field `intercept` sets the prior the intercept is drawn from (an
[`AbstractPriorModel`](@ref); a bare `Distribution` is coerced via
[`as_prior`](@ref)).

# Examples
```@example Intercept
using ComposableTuringIDModels, Distributions
int = Intercept(Normal(0, 1))
mdl = as_turing_model(int, 10)
rand(mdl)
```
"
struct Intercept{D <: AbstractPriorModel} <: AbstractLatentModel
    "Prior for the intercept."
    intercept::D
    Intercept(intercept::AbstractPriorModel) = new{typeof(intercept)}(intercept)
end

Intercept(intercept) = Intercept(as_prior(intercept, :intercept))
Intercept(; intercept) = Intercept(as_prior(intercept, :intercept))

@model function as_turing_model(model::Intercept, n)
    intercept ~ to_submodel(as_turing_model(model.intercept, 1), false)
    return fill(only(intercept), n)
end

@doc raw"
A fixed (non-sampled) intercept broadcast to a length-`n` latent process.
"
@kwdef struct FixedIntercept{F <: Real} <: AbstractLatentModel
    intercept::F
end

@model function as_turing_model(model::FixedIntercept, n)
    return fill(model.intercept, n)
end
