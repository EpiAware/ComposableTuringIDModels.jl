# Negative-binomial observation-error model.

@doc raw"
A negative-binomial observation-error model with an inferred cluster factor.

The field `cluster_factor` sets the prior for the cluster factor that is sampled
and used to parameterise the negative-binomial error (an
[`AbstractPriorModel`](@ref); a bare `Distribution` is coerced via
[`as_prior`](@ref)).

# Examples
```@example NegativeBinomialError
using ComposableTuringIDModels, Distributions
nb = NegativeBinomialError()
mdl = as_turing_model(nb, missing, fill(10, 10))
rand(mdl)
```
"
struct NegativeBinomialError{S <: AbstractPriorModel} <:
       AbstractObservationErrorModel
    "Prior for the cluster factor."
    cluster_factor::S
    function NegativeBinomialError(cluster_factor::AbstractPriorModel)
        return new{typeof(cluster_factor)}(cluster_factor)
    end
end

function NegativeBinomialError(cluster_factor)
    return NegativeBinomialError(as_prior(cluster_factor))
end
function NegativeBinomialError(; cluster_factor = HalfNormal(0.01))
    return NegativeBinomialError(as_prior(cluster_factor))
end

@model function generate_observation_error_priors(
        obs_model::NegativeBinomialError, y_t, Y_t)
    cluster_factor ~ to_submodel(
        as_turing_model(obs_model.cluster_factor, 1), false)
    sq_cluster_factor = only(cluster_factor)^2
    return (; sq_cluster_factor)
end

function observation_error(::NegativeBinomialError, Y_t, sq_cluster_factor)
    return NegativeBinomialMeanClust(Y_t, sq_cluster_factor)
end
