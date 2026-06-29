# Normal (Gaussian) observation-error model.

@doc raw"
A normal (Gaussian) observation-error model with an inferred standard deviation.

Unlike [`PoissonError`](@ref) and [`NegativeBinomialError`](@ref), which model
**count** observations, `NormalError` models **continuous** observations: each
observed value is normally distributed about its expected value,

```math
y_t \sim \mathrm{Normal}(Y_t, \sigma),
```

with the standard deviation ``\sigma`` drawn from `std_prior`. It is the minimal
non-count observation error, useful for already-aggregated or transformed
quantities (e.g. log-incidence, prevalence proportions, wastewater
concentrations) where a Gaussian likelihood is appropriate.

The field `std_prior` sets the prior distribution for ``\sigma``.

# Examples
```@example NormalError
using EpiAwarePrototype, Distributions
ne = NormalError()
mdl = as_turing_model(ne, missing, fill(10.0, 10))
rand(mdl)
```
"
@kwdef struct NormalError{S <: Sampleable} <: AbstractObservationErrorModel
    "Prior distribution for the observation standard deviation."
    std_prior::S = HalfNormal(0.1)
end

@model function generate_observation_error_priors(obs_model::NormalError, y_t, Y_t)
    σ ~ obs_model.std_prior
    return (; σ)
end

observation_error(::NormalError, Y_t, σ) = Normal(Y_t, σ)
