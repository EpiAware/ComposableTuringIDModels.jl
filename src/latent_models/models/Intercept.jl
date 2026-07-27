# Intercept latent process models (sampled and fixed).

@doc raw"
Broadcast a single sampled intercept value to a length-`n` latent process.

The field `intercept` sets the prior the intercept is drawn from — a
`Distribution`, drawn with a native tilde (a single scalar draw broadcast to
length `n`).

# Examples
```@example Intercept
using ComposableTuringIDModels, Distributions
int = Intercept(Normal(0, 1))
mdl = as_turing_model(int, 10)
rand(mdl)
```
"
struct Intercept{D <: Distribution} <: AbstractLatentModel
    "Prior for the intercept."
    intercept::D
end

Intercept(; intercept) = Intercept(intercept)

@model function as_turing_model(model::Intercept, n)
    intercept ~ model.intercept
    return fill(intercept, n)
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
