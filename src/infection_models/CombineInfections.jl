# Stack several (possibly different) infection processes into one
# infection-strata x time `I_t` matrix — the "different many" many-to-many
# infection-side mapping (issue #180). Complements `Split`/`StrataMap`, which
# already carry the observation-side mapping from an `inf_strata x time`
# matrix onto observation streams.

@doc raw"
Combine several infection processes into one `n_strata x n_time` `I_t` matrix.

Each model in `models` is drawn **independently** over the same series length
`n` — a genuinely different infection process per stratum (a different region,
a different pathogen variant), each with its own internal latent. Stratum `k`
is row `k` of the returned `I_t`, the same `inf_strata x time` orientation
[`StrataMap`](@ref) uses, so `CombineInfections` composes directly with
[`Split`](@ref) on the observation side: `IDModel(CombineInfections(models),
Split(template, W))` maps several distinct infection processes onto
observation streams through a weight matrix `W`, covering many-to-one (a
single aggregation row) and many-to-many (a general `W`).

Each sub-model is **deliberately prefixed** by its `names` entry before being
sampled as a submodel (like [`Split`](@ref)): two independent infection
processes would otherwise collide on variable names (e.g. two [`RandomWalk`](@ref)s
both naming their innovation `ϵ_t`), so this is one of the components that
prefixes on purpose rather than following the package's flat, prefix-off
default.

## Fields

  - `models`: the vector of infection models, one per stratum.
  - `names`: the stratum names, used both as the submodel prefix and as the
    `Z_t` NamedTuple keys.

# Examples
```@example CombineInfections
using ComposableTuringIDModels, Distributions
north = DirectInfections(; Z = RandomWalk(), initialisation = Normal(log(50.0), 0.2))
south = DirectInfections(; Z = RandomWalk(), initialisation = Normal(log(20.0), 0.2))
model = CombineInfections([north, south], [\"north\", \"south\"])
sim = as_turing_model(model, 12)()
size(sim.I_t)   # 2 strata x 12 time points
```
"
struct CombineInfections{M <: AbstractVector, P <: AbstractVector{<:String}} <:
       AbstractInfectionModel
    "The vector of infection models, one per stratum."
    models::M
    "The stratum names (submodel prefix and `Z_t` keys)."
    names::P

    function CombineInfections(
            models::M, names::P) where {M <: AbstractVector, P <: AbstractVector{<:String}}
        @assert !isempty(models) "CombineInfections needs at least one model"
        @assert length(models)==length(names) "The number of models ($(length(models))) and names ($(length(names))) must be equal"
        return new{M, P}(models, names)
    end
end

function CombineInfections(models::AbstractVector)
    names = "inf" .* string.(1:length(models))
    return CombineInfections(models, names)
end

@model function as_turing_model(model::CombineInfections, n)
    n_strata = length(model.models)
    Is = Vector{Any}(undef, n_strata)
    Zs = Vector{Any}(undef, n_strata)
    for i in 1:n_strata
        res ~ to_submodel(
            prefix(as_turing_model(model.models[i], n), Symbol(model.names[i])),
            false)
        Is[i] = res.I_t
        Zs[i] = res.Z_t
    end
    I_t = permutedims(reduce(hcat, Is))
    keysyms = Tuple(Symbol.(model.names))
    Z_t = NamedTuple{keysyms}(Tuple(Zs))
    return (; I_t, Z_t)
end
