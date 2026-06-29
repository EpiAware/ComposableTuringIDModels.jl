# Combine latent models by summing their outputs.

@doc raw"
Combine several latent models of the same length by summing their outputs.

Each component is generated over the full length `n` and the results are added.
When a non-empty prefix is supplied for a component it is wrapped in a
[`PrefixLatentModel`](@ref) so its variables stay distinct.

# Arguments

  - `latent_models`: the [`CombineLatentModels`](@ref) collection.
  - `n`: the length of the latent series to generate.

# Examples
```@example CombineLatentModels
using EpiAwarePrototype, Distributions
combined = CombineLatentModels([Intercept(Normal(2, 0.2)), AR()])
rand(as_turing_model(combined, 10))
```

## Fields

  - `models`: the vector of latent models (prefix-wrapped where a prefix is set).
  - `prefixes`: the vector of prefixes, one per model.
"
struct CombineLatentModels{
    M <: AbstractVector{<:AbstractEpiAwareModel}, P <: AbstractVector{<:String}} <:
       AbstractEpiAwareModel
    "A vector of latent models."
    models::M
    "A vector of prefixes for the latent models."
    prefixes::P

    function CombineLatentModels(models::M,
            prefixes::P) where {
            M <: AbstractVector{<:AbstractEpiAwareModel},
            P <: AbstractVector{<:String}}
        @assert length(models)>1 "At least two models are required"
        @assert length(models)==length(prefixes) "The number of models and prefixes must be equal"
        prefix_models = [prefixes[i] == "" ? models[i] :
                         PrefixLatentModel(models[i], prefixes[i])
                         for i in eachindex(models)]
        return new{AbstractVector{<:AbstractEpiAwareModel},
            AbstractVector{<:String}}(prefix_models, prefixes)
    end
end

function CombineLatentModels(models::M) where {
        M <: AbstractVector{<:AbstractEpiAwareModel}}
    prefixes = "Combine." .* string.(1:length(models))
    return CombineLatentModels(models, prefixes)
end

@model function as_turing_model(latent_models::CombineLatentModels, n)
    final_latent ~ to_submodel(
        _accumulate_latents(latent_models.models, 1, fill(0.0, n), n,
            length(latent_models.models)), false)
    return final_latent
end

@model function _accumulate_latents(models, index, acc_latent, n, n_models)
    if index > n_models
        return acc_latent
    else
        latent ~ to_submodel(as_turing_model(models[index], n), false)
        updated_latent ~ to_submodel(
            _accumulate_latents(models, index + 1, acc_latent .+ latent, n,
                n_models), false)
        return updated_latent
    end
end
