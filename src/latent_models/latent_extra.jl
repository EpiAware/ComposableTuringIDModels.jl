# Latent-model modifiers, manipulators, combinations, and broadcasting. Each
# struct subtypes `AbstractEpiAwareModel` and implements one `as_turing_model`
# method. Prefixing — previously the `prefix_submodel` helper — is done here with
# `DynamicPPL.prefix` applied to the inner model before `to_submodel(..., false)`.

# --- modifiers --------------------------------------------------------------

@doc raw"
Apply a transformation function to the output of an inner latent model.

# Arguments

  - `model`: the inner latent model whose output is transformed.
  - `n`: the length of the latent series to generate.

# Examples
```@example TransformLatentModel
using EpiAwarePrototype, Distributions
trans = TransformLatentModel(Intercept(Normal(2, 0.2)), x -> exp.(x))
rand(as_turing_model(trans, 5))
```

## Fields

  - `model`: the latent model to transform.
  - `transform`: the transformation function applied to the latent vector.
"
@kwdef struct TransformLatentModel{M <: AbstractEpiAwareModel, F <: Function} <:
              AbstractEpiAwareModel
    "The latent model to transform."
    model::M
    "The transformation function."
    transform::F
end

@model function as_turing_model(model::TransformLatentModel, n)
    untransformed ~ to_submodel(as_turing_model(model.model, n), false)
    return model.transform(untransformed)
end

@doc raw"
Wrap an inner latent model so its sampled variables are prefixed with `prefix`.

This replaces the original `prefix_submodel` helper: the inner model is prefixed
with `DynamicPPL.prefix` before being sampled as a submodel, so its variables
appear as `prefix.varname`.

# Arguments

  - `model`: the inner latent model.
  - `n`: the length of the latent series to generate.

# Examples
```@example PrefixLatentModel
using EpiAwarePrototype
pm = PrefixLatentModel(; model = HierarchicalNormal(), prefix = \"Test\")
rand(as_turing_model(pm, 10))
```

## Fields

  - `model`: the latent model to prefix.
  - `prefix`: the string prefix applied to the inner model's variables.
"
@kwdef struct PrefixLatentModel{M <: AbstractEpiAwareModel, P <: String} <:
              AbstractEpiAwareModel
    "The latent model."
    model::M
    "The prefix for the latent model."
    prefix::P
end

@model function as_turing_model(model::PrefixLatentModel, n)
    submodel ~ to_submodel(
        prefix(as_turing_model(model.model, n), Symbol(model.prefix)), false)
    return submodel
end

@doc raw"
Record the inner latent vector as a tracked generated quantity (`exp_latent`).

# Arguments

  - `model`: the inner latent model whose output is recorded.
  - `n`: the length of the latent series to generate.

# Examples
```@example RecordExpectedLatent
using EpiAwarePrototype
rm = RecordExpectedLatent(FixedIntercept(0.1))
rand(as_turing_model(rm, 1))
```

## Fields

  - `model`: the latent model whose expected latent vector is recorded.
"
struct RecordExpectedLatent{M <: AbstractEpiAwareModel} <: AbstractEpiAwareModel
    "The latent model whose expected latent vector is recorded."
    model::M
end

@model function as_turing_model(model::RecordExpectedLatent, n)
    latent ~ to_submodel(as_turing_model(model.model, n), false)
    exp_latent := latent
    return latent
end

# --- manipulators -----------------------------------------------------------

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
using EpiAwarePrototype, Distributions
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
    M <: AbstractVector{<:AbstractEpiAwareModel}, N <: Int, F <: Function,
    P <: AbstractVector{<:String}} <: AbstractEpiAwareModel
    "A vector of latent models."
    models::M
    "The number of models in the collection."
    no_models::N
    "Maps `(n, no_models)` to a vector of per-model segment lengths."
    dimension_adaptor::F
    "A vector of prefixes for the latent models."
    prefixes::P

    function ConcatLatentModels(models::M, no_models::I, dimension_adaptor::F,
            prefixes::P) where {M <: AbstractVector{<:AbstractEpiAwareModel},
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
        return new{AbstractVector{<:AbstractEpiAwareModel}, Int, Function,
            AbstractVector{<:String}}(
            prefix_models, no_models, dimension_adaptor, prefixes)
    end
end

function ConcatLatentModels(models::M, dimension_adaptor::Function;
        prefixes = nothing) where {M <: AbstractVector{<:AbstractEpiAwareModel}}
    no_models = length(models)
    if isnothing(prefixes)
        prefixes = "Concat." .* string.(1:no_models)
    end
    return ConcatLatentModels(models, no_models, dimension_adaptor, prefixes)
end

function ConcatLatentModels(models::M;
        dimension_adaptor::Function = equal_dimensions, prefixes = nothing) where {
        M <: AbstractVector{<:AbstractEpiAwareModel}}
    return ConcatLatentModels(models, dimension_adaptor; prefixes = prefixes)
end

function ConcatLatentModels(; models::M,
        dimension_adaptor::Function = equal_dimensions, prefixes = nothing) where {
        M <: AbstractVector{<:AbstractEpiAwareModel}}
    return ConcatLatentModels(models, dimension_adaptor; prefixes = prefixes)
end

@doc raw"
Partition `n` elements into `m` segments of as-equal-as-possible length.

The first segment gets `ceil(n / m)` and the rest `floor(n / m)`. This is the
default `dimension_adaptor` for [`ConcatLatentModels`](@ref).

# Arguments

  - `n`: the total number of elements.
  - `m`: the number of segments.

# Examples
```@example equal_dimensions
using EpiAwarePrototype
EpiAwarePrototype.equal_dimensions(10, 3)
```
"
function equal_dimensions(n::Int, m::Int)::AbstractVector{Int}
    return vcat(ceil(n / m), fill(floor(n / m), m - 1))
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

# --- broadcasting -----------------------------------------------------------

@doc raw"
Abstract supertype for broadcast rules used by [`BroadcastLatentModel`](@ref).
A rule defines [`broadcast_n`](@ref) (how long an inner series to generate) and
[`broadcast_rule`](@ref) (how to expand it to length `n`).
"
abstract type AbstractBroadcastRule end

@doc raw"
Length of the inner series an [`AbstractBroadcastRule`](@ref) needs to produce a
length-`n` broadcasted series. Each rule implements its own method.

# Arguments

  - `rule`: the [`AbstractBroadcastRule`](@ref).
  - `n`: the length of the broadcasted output series.
  - `period`: the broadcast period.

# Examples
```@example broadcast_n
using EpiAwarePrototype
broadcast_n(RepeatEach(), 10, 7), broadcast_n(RepeatBlock(), 10, 7)
```
"
function broadcast_n end

@doc raw"
Expand an inner latent series to length `n` under an
[`AbstractBroadcastRule`](@ref). Each rule implements its own method.

# Arguments

  - `rule`: the [`AbstractBroadcastRule`](@ref).
  - `latent`: the inner latent series to expand.
  - `n`: the length of the broadcasted output series.
  - `period`: the broadcast period.

# Examples
```@example broadcast_rule
using EpiAwarePrototype
broadcast_rule(RepeatEach(), [1, 2], 5, 2)
```
"
function broadcast_rule end

@doc raw"
Broadcast a shorter latent process to length `n` under a broadcast rule.

The inner model is generated over the length the rule requires
([`broadcast_n`](@ref)), then expanded to length `n` ([`broadcast_rule`](@ref)).

# Arguments

  - `model`: the [`BroadcastLatentModel`](@ref).
  - `n`: the length of the broadcasted series to generate.

# Examples
```@example BroadcastLatentModel
using EpiAwarePrototype
each = BroadcastLatentModel(RandomWalk(), 7, RepeatEach())
rand(as_turing_model(each, 10))
```

## Fields

  - `model`: the underlying latent model.
  - `period`: the broadcast period.
  - `broadcast_rule`: the [`AbstractBroadcastRule`](@ref) applied.
"
struct BroadcastLatentModel{
    M <: AbstractEpiAwareModel, P <: Integer, B <: AbstractBroadcastRule} <:
       AbstractEpiAwareModel
    "The underlying latent model."
    model::M
    "The period of the broadcast."
    period::P
    "The broadcast rule applied."
    broadcast_rule::B

    function BroadcastLatentModel(model::M,
            period::Integer,
            broadcast_rule::B) where {
            M <: AbstractEpiAwareModel, B <: AbstractBroadcastRule}
        @assert period>0 "period must be greater than 0"
        new{typeof(model), typeof(period), typeof(broadcast_rule)}(
            model, period, broadcast_rule)
    end
end

function BroadcastLatentModel(model::M; period::Integer,
        broadcast_rule::B) where {
        M <: AbstractEpiAwareModel, B <: AbstractBroadcastRule}
    return BroadcastLatentModel(model, period, broadcast_rule)
end

@model function as_turing_model(model::BroadcastLatentModel, n)
    m = broadcast_n(model.broadcast_rule, n, model.period)
    latent_period ~ to_submodel(as_turing_model(model.model, m), false)
    return broadcast_rule(model.broadcast_rule, latent_period, n, model.period)
end

@doc raw"
Broadcast rule that repeats the latent process at each position within a period
(e.g. a fixed day-of-week effect).

# Examples
```@example RepeatEach
using EpiAwarePrototype
broadcast_rule(RepeatEach(), [1, 2], 10, 2)
```
"
struct RepeatEach <: AbstractBroadcastRule end

broadcast_n(::RepeatEach, n, period) = period

function broadcast_rule(::RepeatEach, latent, n, period)
    @assert length(latent)==period "length(latent) must equal period"
    broadcast_latent = repeat(latent; outer = ceil(Int, n / period))
    return broadcast_latent[1:n]
end

@doc raw"
Broadcast rule that repeats the latent process in blocks of length `period`
(e.g. a piecewise-constant weekly process).

# Examples
```@example RepeatBlock
using EpiAwarePrototype
broadcast_rule(RepeatBlock(), [1, 2, 3, 4, 5], 10, 2)
```
"
struct RepeatBlock <: AbstractBroadcastRule end

broadcast_n(::RepeatBlock, n, period) = ceil(Int, n / period)

function broadcast_rule(::RepeatBlock, latent, n, period)
    @assert n<=period * length(latent) "n must be ≤ period * length(latent)"
    broadcast_latent = [latent[j] for j in 1:length(latent) for _ in 1:period]
    return broadcast_latent[1:n]
end

@doc raw"
Build a [`BroadcastLatentModel`](@ref) for a day-of-week effect: a transformed
inner model repeated across a 7-day period.

# Arguments

  - `model`: the inner latent model.
  - `link`: link applied before broadcasting (default `x -> 7 * softmax(x)`,
    constraining the week effects to sum to 7).

# Examples
```@example broadcast_dayofweek
using EpiAwarePrototype
broadcast_dayofweek(RandomWalk())
```
"
function broadcast_dayofweek(model::AbstractEpiAwareModel; link = x -> 7 * softmax(x))
    return BroadcastLatentModel(TransformLatentModel(model, link), 7, RepeatEach())
end

@doc raw"
Build a [`BroadcastLatentModel`](@ref) for a piecewise-constant weekly process.

# Arguments

  - `model`: the inner latent model.

# Examples
```@example broadcast_weekly
using EpiAwarePrototype
broadcast_weekly(RandomWalk())
```
"
function broadcast_weekly(model::AbstractEpiAwareModel)
    return BroadcastLatentModel(model, 7, RepeatBlock())
end

# --- combinations -----------------------------------------------------------

@doc raw"
Build an ARMA(p, q) latent process: an [`AR`](@ref) whose innovation model is an
[`MA`](@ref).

# Arguments

  - `init`: prior(s) for the AR initial conditions.
  - `damp`: prior(s) for the AR damping coefficients.
  - `θ`: prior(s) for the MA coefficients.
  - `ϵ_t`: the innovation model (default [`HierarchicalNormal`](@ref)).

# Examples
```@example arma
using EpiAwarePrototype, Distributions
model = arma(; θ = [truncated(Normal(0.0, 0.02), -1, 1)],
    damp = [truncated(Normal(0.0, 0.02), 0, 1)])
rand(as_turing_model(model, 10))
```
"
function arma(; init = [Normal()], damp = [truncated(Normal(0.0, 0.05), 0, 1)],
        θ = [truncated(Normal(0.0, 0.05), -1, 1)], ϵ_t = HierarchicalNormal())
    ma = MA(; θ_priors = θ, ϵ_t = ϵ_t)
    return AR(; damp_priors = damp, init_priors = init, ϵ_t = ma)
end

@doc raw"
Build an ARIMA(p, d, q) latent process: an [`arma`](@ref) wrapped in a
`d`-fold [`DiffLatentModel`](@ref).

# Arguments

  - `ar_init`: prior(s) for the AR initial conditions.
  - `diff_init`: prior(s) for the differencing initial conditions (sets `d`).
  - `damp`: prior(s) for the AR damping coefficients.
  - `θ`: prior(s) for the MA coefficients.
  - `ϵ_t`: the innovation model (default [`HierarchicalNormal`](@ref)).

# Examples
```@example arima
using EpiAwarePrototype, Distributions
model = arima()
rand(as_turing_model(model, 10))
```
"
function arima(; ar_init = [Normal()], diff_init = [Normal()],
        damp = [truncated(Normal(0.0, 0.05), 0, 1)],
        θ = [truncated(Normal(0.0, 0.05), -1, 1)], ϵ_t = HierarchicalNormal())
    arma_model = arma(; init = ar_init, damp = damp, θ = θ, ϵ_t = ϵ_t)
    return DiffLatentModel(; model = arma_model, init_priors = diff_init)
end
