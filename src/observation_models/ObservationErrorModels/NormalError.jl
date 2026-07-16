# Normal (Gaussian) observation-error model.

@doc raw"
A normal (Gaussian) observation-error model with an inferred standard deviation.

Unlike [`PoissonError`](@ref) and [`NegativeBinomialError`](@ref), which model
**count** observations, `NormalError` models **continuous** observations: each
observed value is normally distributed about its expected value,

```math
y_t \sim \mathrm{Normal}(Y_t, \sigma),
```

with the standard deviation ``\sigma`` drawn from the prior in `std`. It is the
minimal non-count observation error, useful for already-aggregated or transformed
quantities (e.g. log-incidence, prevalence proportions, wastewater
concentrations) where a Gaussian likelihood is appropriate.

The field `std` sets the prior for ``\sigma`` — a `Distribution`, drawn with a
native tilde (a plain scalar draw, no submodel).

# Examples
```@example NormalError
using ComposableTuringIDModels, Distributions
ne = NormalError()
mdl = as_turing_model(ne, missing, fill(10.0, 10))
rand(mdl)
```
"
struct NormalError{S <: Distribution} <: AbstractObservationErrorModel
    "Prior for the observation standard deviation."
    std::S
end

NormalError(; std = HalfNormal(0.1)) = NormalError(std)

@model function generate_observation_error_priors(obs_model::NormalError, y_t, Y_t)
    σ ~ obs_model.std
    return (; σ = σ)
end

observation_error(::NormalError, Y_t, σ) = Normal(Y_t, σ)
