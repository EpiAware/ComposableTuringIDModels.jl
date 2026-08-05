# Combine latent models by summing their outputs.

@doc raw"
Combine several latent models of the same length by summing their outputs.

Each component is generated over the full length `n` and the results are added.
When a non-empty prefix is supplied for a component it is wrapped in a
[`PrefixLatentModel`](@ref) so its variables stay distinct.

# Arguments

  - `latent_models`: the [`CombineLatentModels`](@ref) collection.
  - `n`: the shape to generate — a length or an `(n_strata, n_time)` shape,
    whatever the component models accept.

# Examples
```@example CombineLatentModels
using ComposableTuringIDModels, Distributions
combined = CombineLatentModels([Intercept(Normal(2, 0.2)), AR()])
rand(as_turing_model(combined, 10))
```

## Fields

  - `models`: the vector of latent models (prefix-wrapped where a prefix is set).
  - `prefixes`: the vector of prefixes, one per model.
"
struct CombineLatentModels{M <: AbstractVector, P <: AbstractVector{<:String}} <:
       AbstractLatentModel
    "A vector of latent models."
    models::M
    "A vector of prefixes for the latent models."
    prefixes::P

    function CombineLatentModels(models::AbstractVector,
            prefixes::P) where {P <: AbstractVector{<:String}}
        @assert length(models)>1 "At least two models are required"
        @assert length(models)==length(prefixes) "The number of models and prefixes must be equal"
        # Each member is a length-`n` PATH slot, so a bare `Distribution` is
        # wrapped in an `Intercept` (a constant path) before it is namespaced;
        # then non-empty prefixes get a `PrefixLatentModel` so variables stay
        # distinct. A process / `IID` / vector member passes through unchanged.
        prefix_models = [prefixes[i] == "" ? _path_prior(models[i]) :
                         PrefixLatentModel(_path_prior(models[i]), prefixes[i])
                         for i in eachindex(models)]
        return new{AbstractVector, AbstractVector{<:String}}(
            prefix_models, prefixes)
    end
end

function CombineLatentModels(models::AbstractVector)
    prefixes = "Combine." .* string.(1:length(models))
    return CombineLatentModels(models, prefixes)
end

# `CombineLatentModels` wraps whatever its member models return, so it serves
# both a length-`n` path and an `(n_strata, n_time)` shape (`fill(0.0, n)`
# builds a zero of either shape). The two methods below delegate to one shared
# `@model`, which avoids both the duplication and a dispatch ambiguity against
# the `AbstractPriorModel`/`Dims{2}` guard.
@model function _combine_latents(latent_models::CombineLatentModels, n)
    final_latent ~ to_submodel(
        _accumulate_latents(latent_models.models, 1, fill(0.0, n), n,
            length(latent_models.models)), false)
    return final_latent
end
function as_turing_model(latent_models::CombineLatentModels, n::Int)
    return _combine_latents(latent_models, n)
end
function as_turing_model(latent_models::CombineLatentModels, n::Dims{2})
    return _combine_latents(latent_models, n)
end

@model function _accumulate_latents(models, index, acc_latent, n, n_models)
    if index > n_models
        return acc_latent
    else
        latent ~ as_turing_submodel(models[index], n)
        updated_latent ~ to_submodel(
            _accumulate_latents(models, index + 1, acc_latent .+ latent, n,
                n_models), false)
        return updated_latent
    end
end
