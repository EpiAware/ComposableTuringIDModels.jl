# Negative-binomial observation-error model.

@doc raw"
A negative-binomial observation-error model with an inferred cluster factor.

The field `cluster_factor_prior` sets the prior distribution for the cluster
factor that is sampled and used to parameterise the negative-binomial error.

# Examples
```@example NegativeBinomialError
using ComposableTuringIDModels, Distributions
nb = NegativeBinomialError()
mdl = as_turing_model(nb, missing, fill(10, 10))
rand(mdl)
```
"
@kwdef struct NegativeBinomialError{S <: Sampleable} <: AbstractObservationErrorModel
    "Prior distribution for the cluster factor."
    cluster_factor_prior::S = HalfNormal(0.01)
end

@model function generate_observation_error_priors(
        obs_model::NegativeBinomialError, y_t, Y_t)
    cluster_factor ~ obs_model.cluster_factor_prior
    sq_cluster_factor = cluster_factor^2
    return (; sq_cluster_factor)
end

function observation_error(::NegativeBinomialError, Y_t, sq_cluster_factor)
    return NegativeBinomialMeanClust(Y_t, sq_cluster_factor)
end
