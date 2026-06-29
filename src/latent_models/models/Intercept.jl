# Intercept latent process models (sampled and fixed).

@doc raw"
Broadcast a single sampled intercept value to a length-`n` latent process.

The field `intercept_prior` sets the prior distribution the intercept is drawn
from.

# Examples
```@example Intercept
using EpiAwarePrototype, Distributions
int = Intercept(Normal(0, 1))
mdl = as_turing_model(int, 10)
rand(mdl)
```
"
@kwdef struct Intercept{D <: Sampleable} <: AbstractEpiAwareModel
    "Prior distribution for the intercept."
    intercept_prior::D
end

@model function as_turing_model(model::Intercept, n)
    intercept ~ model.intercept_prior
    return fill(intercept, n)
end

@doc raw"
A fixed (non-sampled) intercept broadcast to a length-`n` latent process.
"
@kwdef struct FixedIntercept{F <: Real} <: AbstractEpiAwareModel
    intercept::F
end

@model function as_turing_model(model::FixedIntercept, n)
    return fill(model.intercept, n)
end
