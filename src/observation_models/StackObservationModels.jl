# Stack several observation models, each applied to a named data component.

@doc raw"
Stack several observation models, each applied to a named component of the data.

Each inner model is wrapped in a [`PrefixObservationModel`](@ref) keyed by its
name, so the stacked variables stay distinct. The model is constructed either from
parallel vectors of models and names, or from a `NamedTuple` of models (the keys
supply the names). When sampled, each component model is applied to the matching
entry of the `y_t` / `Y_t` named tuples; a single expected-observation vector is
broadcast across all components.

# Arguments

  - `obs_model`: the [`StackObservationModels`](@ref) model.
  - `y_t`: a `NamedTuple` of observed series, one per stacked model.
  - `Y_t`: a `NamedTuple` of expected-observation series (or a single vector
    broadcast across the components).

# Examples
```@example StackObservationModels
using EpiAwarePrototype
obs = StackObservationModels((cases = PoissonError(), deaths = PoissonError()))
mdl = as_turing_model(obs, (cases = missing, deaths = missing), fill(10.0, 5))
rand(mdl)
```

## Fields

  - `models`: the vector of observation models (each prefix-wrapped by its name).
  - `model_names`: the names identifying each stacked model.
"
struct StackObservationModels{
    M <: AbstractVector, N <: AbstractVector{<:AbstractString}} <:
       AbstractObservationModel
    "The vector of observation models (each prefix-wrapped by its name)."
    models::M
    "The names identifying each stacked model."
    model_names::N

    function StackObservationModels(
            models::M, model_names::N) where {
            M <: AbstractVector, N <: AbstractVector{<:AbstractString}}
        @assert length(models)==length(model_names) "The number of models and model names must be equal"
        prefix_models = [PrefixObservationModel(models[i], model_names[i])
                         for i in eachindex(models)]
        return new{typeof(prefix_models), N}(prefix_models, model_names)
    end
end

function StackObservationModels(models::NamedTuple)
    model_names = keys(models) .|> string |> collect
    return StackObservationModels(collect(values(models)), model_names)
end

@model function as_turing_model(
        obs_model::StackObservationModels, y_t::NamedTuple, Y_t::NamedTuple)
    @assert length(obs_model.models)==length(y_t) "The number of models must match the number of observed series"
    @assert obs_model.model_names==(keys(y_t) .|> string |> collect) "The model names must match the keys of the observed series"
    @assert keys(y_t)==keys(Y_t) "The keys of the observed and expected series must match"
    obs = Vector{Any}(undef, length(obs_model.models))
    for i in eachindex(obs_model.models)
        name = obs_model.model_names[i]
        obs_i ~ to_submodel(
            as_turing_model(
                obs_model.models[i], y_t[Symbol(name)], Y_t[Symbol(name)]),
            false)
        obs[i] = obs_i
    end
    return obs
end

@model function as_turing_model(
        obs_model::StackObservationModels, y_t::NamedTuple, Y_t::AbstractVector)
    tuple_Y_t = NamedTuple{keys(y_t)}(fill(Y_t, length(y_t)))
    obs ~ to_submodel(as_turing_model(obs_model, y_t, tuple_Y_t), false)
    return obs
end
