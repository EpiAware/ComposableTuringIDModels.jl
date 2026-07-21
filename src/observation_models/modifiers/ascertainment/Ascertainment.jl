# Ascertainment observation modifier (scale expected observations by a latent
# process).

@doc raw"
Scale the expected observations of an underlying observation model by an
ascertainment prior process.

The `latent_model` slot takes a latent model for a **time-varying** ascertainment
effect, or a bare `Distribution` for a single **constant** ascertainment factor
shared across the series (wrapped in an [`Intercept`](@ref), so one value is drawn
and broadcast). Whatever is passed generates a length-`length(Y_t)`
series which is combined with the expected observations `Y_t` through `transform`
before being passed to the inner observation `model`. The default `transform`
applies a multiplicative effect on the exponential scale
(`(Y_t, x) -> xexpy.(Y_t, x)`), so a value `x` multiplies the expected
observation by `exp(x)`. The prior is prefixed with `latent_prefix` (a latent
model via [`PrefixLatentModel`](@ref), a distribution via its sampled-variable
name) unless the prefix is the empty string.

# Arguments

  - `obs_model`: the [`Ascertainment`](@ref) model.
  - `y_t`: the observed series (or `missing` when simulating predictively).
  - `Y_t`: the expected-observation series.

# Examples
```@example Ascertainment
using ComposableTuringIDModels, Distributions
# A latent model gives a time-varying ascertainment effect. The default
# transform reads the effect on the log scale, so `log(0.1)` is a 10%
# ascertainment rate.
obs = Ascertainment(PoissonError(), FixedIntercept(log(0.1)))
rand(as_turing_model(obs, missing, fill(10.0, 5)))
# A bare Distribution / prior gives a single constant ascertainment factor.
obs_const = Ascertainment(PoissonError(), Normal(0.0, 0.1))
rand(as_turing_model(obs_const, missing, fill(10.0, 5)))
```

## Fields

  - `model`: the underlying observation model the ascertained expected
    observations are passed to.
  - `latent_model`: the prior model generating the ascertainment effect (a latent
    model for a time-varying effect, a distribution/prior for a constant factor),
    prefixed unless `latent_prefix` is empty.
  - `transform`: the function `(Y_t, x)` combining expected observations with the
    ascertainment effect.
  - `latent_prefix`: the prefix applied to the ascertainment prior's variables.
"
struct Ascertainment{
    M <: AbstractObservationModel, L <: AbstractPriorModel, F <: Function,
    P <: String} <: AbstractObservationModel
    "The underlying observation model."
    model::M
    "The prior model generating the ascertainment effect."
    latent_model::L
    "The function combining expected observations with the ascertainment effect."
    transform::F
    "The prefix applied to the ascertainment prior's variables."
    latent_prefix::P

    function Ascertainment(model::M, latent_model, transform::F,
            latent_prefix::P) where {M <: AbstractObservationModel,
            F <: Function, P <: String}
        @assert hasmethod(transform, Tuple{Vector, Vector}) "transform must have a method for (Vector, Vector)"
        # A latent model is used as-is (a time-varying effect); a bare
        # `Distribution` is wrapped in an `Intercept` (a single constant factor
        # drawn once and broadcast). Namespace it under `latent_prefix` via
        # `PrefixLatentModel` so the ascertainment variables stay distinct; an
        # empty prefix opts out, leaving the prior unprefixed.
        coerced = _ascertainment_prior(latent_model)
        prior = latent_prefix == "" ? coerced :
                PrefixLatentModel(coerced, latent_prefix)
        return new{M, typeof(prior), F, P}(
            model, prior, transform, latent_prefix)
    end
end

_ascertainment_prior(latent_model::AbstractPriorModel) = latent_model
_ascertainment_prior(dist::Distribution) = Intercept(dist)

function Ascertainment(model::AbstractObservationModel, latent_model;
        transform = (Y_t, x) -> xexpy.(Y_t, x),
        latent_prefix::String = "Ascertainment")
    return Ascertainment(model, latent_model, transform, latent_prefix)
end

function Ascertainment(; model::AbstractObservationModel, latent_model,
        transform = (Y_t, x) -> xexpy.(Y_t, x),
        latent_prefix::String = "Ascertainment")
    return Ascertainment(model, latent_model, transform, latent_prefix)
end

@model function as_turing_model(obs_model::Ascertainment, y_t, Y_t)
    expected_obs_mod ~ as_turing_submodel(obs_model.latent_model, length(Y_t))
    expected_obs = obs_model.transform(Y_t, expected_obs_mod)
    inner ~ as_turing_submodel(obs_model.model, y_t, expected_obs)
    return (; y_t = inner.y_t, expected = inner.expected)
end
