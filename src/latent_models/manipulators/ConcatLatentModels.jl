# Concatenate latent models along time.

@doc raw"
Concatenate several latent models along time into one length-`n` series.

The length `n` is partitioned across the component models by `dimension_adaptor`
(default [`equal_dimensions`](@ref)); each component generates its own segment
and the segments are concatenated.

# Arguments

  - `latent_models`: the [`ConcatLatentModels`](@ref) collection.
  - `n`: the total length of the latent series to generate.

# Examples
```@example ConcatLatentModels
using ComposableTuringIDModels, Distributions
combined = ConcatLatentModels([Intercept(Normal(2, 0.2)), AR()])
rand(as_turing_model(combined, 10))
```

## Fields

  - `models`: the vector of latent models (prefix-wrapped where a prefix is set).
  - `no_models`: the number of models in the collection.
  - `dimension_adaptor`: maps `(n, no_models)` to a vector of segment lengths.
  - `prefixes`: the vector of prefixes, one per model.
"
struct ConcatLatentModels{
    M <: AbstractVector{<:AbstractLatentModel}, N <: Int, F <: Function,
    P <: AbstractVector{<:String}} <: AbstractLatentModel
    "A vector of latent models."
    models::M
    "The number of models in the collection."
    no_models::N
    "Maps `(n, no_models)` to a vector of per-model segment lengths."
    dimension_adaptor::F
    "A vector of prefixes for the latent models."
    prefixes::P

    function ConcatLatentModels(models::M, no_models::I, dimension_adaptor::F,
            prefixes::P) where {M <: AbstractVector{<:AbstractLatentModel},
            I <: Int, F <: Function, P <: AbstractVector{<:String}}
        @assert length(models)>1 "At least two models are required"
        @assert length(models)==no_models "no_models must be equal to the number of models"
        check_dim = dimension_adaptor(no_models, no_models)
        @assert typeof(check_dim)<:AbstractVector{Int} "Output of dimension_adaptor must be a vector of integers"
        @assert length(check_dim)==no_models "The vector of dimensions must have the same length as the number of models"
        @assert length(prefixes)==no_models "The number of models and prefixes must be equal"
        prefix_models = [prefixes[i] == "" ? models[i] :
                         PrefixLatentModel(models[i], prefixes[i])
                         for i in eachindex(models)]
        return new{AbstractVector{<:AbstractLatentModel}, Int, Function,
            AbstractVector{<:String}}(
            prefix_models, no_models, dimension_adaptor, prefixes)
    end
end

function ConcatLatentModels(models::M, dimension_adaptor::Function;
        prefixes = nothing) where {M <: AbstractVector{<:AbstractLatentModel}}
    no_models = length(models)
    if isnothing(prefixes)
        prefixes = "Concat." .* string.(1:no_models)
    end
    return ConcatLatentModels(models, no_models, dimension_adaptor, prefixes)
end

function ConcatLatentModels(models::M;
        dimension_adaptor::Function = equal_dimensions, prefixes = nothing) where {
        M <: AbstractVector{<:AbstractLatentModel}}
    return ConcatLatentModels(models, dimension_adaptor; prefixes = prefixes)
end

function ConcatLatentModels(; models::M,
        dimension_adaptor::Function = equal_dimensions, prefixes = nothing) where {
        M <: AbstractVector{<:AbstractLatentModel}}
    return ConcatLatentModels(models, dimension_adaptor; prefixes = prefixes)
end

@doc raw"
Partition `n` elements into `m` segments of as-equal-as-possible length.

Each segment gets `floor(n / m)`, and the `n mod m` leftover elements are handed
out one apiece to the leading segments. The segment lengths therefore always sum
to exactly `n` (differing by at most one), which is required by
[`ConcatLatentModels`](@ref)'s `@assert sum(dims) == n` check. This is the
default `dimension_adaptor` for [`ConcatLatentModels`](@ref).

# Arguments

  - `n`: the total number of elements.
  - `m`: the number of segments.

# Examples
```@example equal_dimensions
using ComposableTuringIDModels
ComposableTuringIDModels.equal_dimensions(10, 3)
```
"
function equal_dimensions(n::Int, m::Int)::Vector{Int}
    base, r = divrem(n, m)
    return [base + (i <= r ? 1 : 0) for i in 1:m]
end

@model function as_turing_model(latent_models::ConcatLatentModels, n)
    @assert latent_models.no_models<n "The number of latent variables must be greater than the number of models"
    dims = latent_models.dimension_adaptor(n, latent_models.no_models)
    @assert all(x -> x > 0, dims) "Non-positive dimensions are not allowed"
    @assert sum(dims)==n "Sum of dimensions must equal the latent dimension"
    final_latent ~ to_submodel(
        _concat_latents(latent_models.models, 1, nothing, dims,
            latent_models.no_models), false)
    return final_latent
end

@model function _concat_latents(
        models, index::Int, acc_latent, dims::AbstractVector{<:Int}, n_models::Int)
    if index > n_models
        return acc_latent
    else
        latent ~ to_submodel(as_turing_model(models[index], dims[index]), false)
        acc_latent = isnothing(acc_latent) ? latent : vcat(acc_latent, latent)
        updated_latent ~ to_submodel(
            _concat_latents(models, index + 1, acc_latent, dims, n_models), false)
        return updated_latent
    end
end
