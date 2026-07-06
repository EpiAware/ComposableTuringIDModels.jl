# Ascertainment observation modifier (scale expected observations by a latent
# process).

@doc raw"
Scale the expected observations of an underlying observation model by a latent
ascertainment process.

A latent model generates a length-`length(Y_t)` series which is combined with the
expected observations `Y_t` through `transform` before being passed to the inner
observation `model`. The default `transform` applies a multiplicative effect on
the exponential scale (`(Y_t, x) -> xexpy.(Y_t, x)`), so a latent value `x`
multiplies the expected observation by `exp(x)`. The latent model is wrapped in a
[`PrefixLatentModel`](@ref) (with prefix `latent_prefix`) unless the prefix is the
empty string.

# Arguments

  - `obs_model`: the [`Ascertainment`](@ref) model.
  - `y_t`: the observed series (or `missing` when simulating predictively).
  - `Y_t`: the expected-observation series.

# Examples
```@example Ascertainment
using ComposableTuringIDModels, Distributions
obs = Ascertainment(PoissonError(), FixedIntercept(0.1))
mdl = as_turing_model(obs, missing, fill(10.0, 5))
rand(mdl)
```

## Fields

  - `model`: the underlying observation model the ascertained expected
    observations are passed to.
  - `latent_model`: the latent model generating the ascertainment effect
    (prefix-wrapped unless `latent_prefix` is empty).
  - `transform`: the function `(Y_t, x)` combining expected observations with the
    latent effect.
  - `latent_prefix`: the prefix applied to the latent model's variables.
"
struct Ascertainment{
    M <: AbstractObservationModel, L <: AbstractLatentModel, F <: Function,
    P <: String} <: AbstractObservationModel
    "The underlying observation model."
    model::M
    "The latent model generating the ascertainment effect."
    latent_model::L
    "The function combining expected observations with the latent effect."
    transform::F
    "The prefix applied to the latent model's variables."
    latent_prefix::P

    function Ascertainment(model::M, latent_model::L, transform::F,
            latent_prefix::P) where {M <: AbstractObservationModel,
            L <: AbstractLatentModel, F <: Function, P <: String}
        @assert hasmethod(transform, Tuple{Vector, Vector}) "transform must have a method for (Vector, Vector)"
        wrapped_latent_model = latent_prefix == "" ? latent_model :
                               PrefixLatentModel(latent_model, latent_prefix)
        return new{M, typeof(wrapped_latent_model), F, P}(
            model, wrapped_latent_model, transform, latent_prefix)
    end
end

function Ascertainment(model::M, latent_model::L;
        transform = (Y_t, x) -> xexpy.(Y_t, x),
        latent_prefix::String = "Ascertainment") where {
        M <: AbstractObservationModel, L <: AbstractLatentModel}
    return Ascertainment(model, latent_model, transform, latent_prefix)
end

function Ascertainment(; model::M, latent_model::L,
        transform = (Y_t, x) -> xexpy.(Y_t, x),
        latent_prefix::String = "Ascertainment") where {
        M <: AbstractObservationModel, L <: AbstractLatentModel}
    return Ascertainment(model, latent_model, transform, latent_prefix)
end

@model function as_turing_model(obs_model::Ascertainment, y_t, Y_t)
    expected_obs_mod ~ to_submodel(
        as_turing_model(obs_model.latent_model, length(Y_t)), false)
    expected_obs = obs_model.transform(Y_t, expected_obs_mod)
    y_t ~ to_submodel(as_turing_model(obs_model.model, y_t, expected_obs), false)
    return y_t
end
