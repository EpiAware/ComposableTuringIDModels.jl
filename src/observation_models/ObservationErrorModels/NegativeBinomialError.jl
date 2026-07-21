# Negative-binomial observation-error model.

@doc raw"
A negative-binomial observation-error model with an inferred cluster factor.

The field `cluster_factor` sets the prior for the cluster factor — a
`Distribution` (a constant, one scalar RV) or a process (a length-`n`, e.g.
time-varying, overdispersion). It is drawn through the single
[`as_turing_submodel`](@ref) seam and read per time point via `_at`, so a process
makes the overdispersion time-varying with no other change.

# Examples
```@example NegativeBinomialError
using ComposableTuringIDModels, Distributions
nb = NegativeBinomialError()
mdl = as_turing_model(nb, missing, fill(10, 10))
rand(mdl)
```
"
struct NegativeBinomialError{S <: PriorLike} <: AbstractObservationErrorModel
    "Prior for the cluster factor."
    cluster_factor::S
end

function NegativeBinomialError(; cluster_factor = HalfNormal(0.01))
    return NegativeBinomialError(cluster_factor)
end

@model function generate_observation_error_priors(
        obs_model::NegativeBinomialError, y_t, Y_t)
    cluster_factor ~ as_turing_submodel(
        obs_model.cluster_factor, length(Y_t); prefix = true)
    sq_cluster_factor = cluster_factor .^ 2
    return (; sq_cluster_factor)
end

function observation_error(::NegativeBinomialError, Y_t, sq_cluster_factor)
    return NegativeBinomialMeanClust(Y_t, sq_cluster_factor)
end
