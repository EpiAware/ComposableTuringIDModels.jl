# A gravity coupling operator, fixed or inferred.

@doc raw"
Build a gravity coupling operator from populations and distances.

```math
K_{gh} \propto \frac{n_g^{\alpha} n_h^{\beta}}{d_{gh}^{\gamma}},
\qquad g \neq h, \qquad K_{gg} \propto \text{within},
\qquad \textstyle\sum_h K_{gh} = 1.
```

Off-diagonal entry ``K_{gh}`` is the share of stratum `g`'s force of infection
that comes from stratum `h`; the diagonal is the share from itself. `within`
sets the self-weight relative to a typical pairwise term, so the default of
`1.0` makes a stratum's own history count about as much as a neighbour's.

Populations are taken in units of their mean and each row is normalised to sum
to one. Both are required rather than optional: [`renewal_pressure`](@ref)
computes ``R_g \sum_h K_{gh} \Lambda_h``, so any overall scale on `K` is
indistinguishable from scaling ``R_t``. Normalising leaves `K` carrying only
where a stratum's force comes from and ``R_t`` carrying its size. `pop` can
therefore be raw counts.

``\alpha`` is identifiable through the split between the diagonal and the
off-diagonals, so it needs `within` to be non-zero. At `within = 0` a common
row factor divides out and ``\alpha`` cancels.

# Arguments

  - `pop`: the population of each stratum, in any units.
  - `dist`: the `strata × strata` distance matrix. Only its off-diagonal
    entries are read.

# Keyword Arguments

  - `α`: the exponent on the destination population (default `1.0`).
  - `β`: the exponent on the origin population (default `1.0`).
  - `γ`: the exponent on distance (default `2.0`).
  - `within`: the own-stratum weight, relative to a typical pairwise term
    (default `1.0`).

# Examples
```@example gravity
using ComposableTuringIDModels
pop = [1.0e6, 2.0e5]
dist = [0.0 50.0; 50.0 0.0]
gravity(pop, dist; α = 0.0, β = 1.0, γ = 2.0)
```
"
function gravity(pop, dist; α = 1.0, β = 1.0, γ = 2.0, within = 1.0)
    n = length(pop)
    @assert size(dist)==(n, n) "`dist` must be $n x $n for $n strata"
    p = pop ./ (sum(pop) / n)
    K = [g == h ? within : p[g]^α * p[h]^β / dist[g, h]^γ
         for g in 1:n, h in 1:n]
    return K ./ sum(K; dims = 2)
end

@doc raw"
A gravity coupling operator with inferred exponents.

The exponents of [`gravity`](@ref) are prior slots, so the operator is
rebuilt from sampled values on every draw.
The coupling between strata is therefore estimated rather than assumed.
Handing one to [`Renewal`](@ref)'s `mixing` keyword is the only difference
from a fixed matrix.
The renewal draws it through the step seam before the scan.

Each exponent slot takes a bare `Distribution` (one scalar draw), a
[`FixedIntercept`](@ref) to hold it fixed, or any other prior model, so a
partly-fixed operator needs no separate type.

Both the fixed and the inferred path call the same [`gravity`](@ref) function,
so they cannot drift.

## Fields

  - `pop`: the population of each stratum.
  - `dist`: the `strata × strata` distance matrix.
  - `α`: prior for the exponent on the destination population.
  - `β`: prior for the exponent on the origin population.
  - `γ`: prior for the exponent on distance.
  - `within`: the own-stratum weight, relative to a typical pairwise term.

# Examples
```@example Gravity
using ComposableTuringIDModels, Distributions
pop = [1.0e6, 2.0e5]
dist = [0.0 50.0; 50.0 0.0]
g = Gravity(pop, dist; α = Normal(0, 0.5), β = Normal(0, 0.5),
    γ = truncated(Normal(2, 0.5), 0, Inf))
size(as_turing_model(g, (2, 20))())
```
"
struct Gravity{P, D, A <: PriorLike, B <: PriorLike, C <: PriorLike, W} <:
       AbstractMixingModel
    "The population of each stratum."
    pop::P
    "The `strata × strata` distance matrix."
    dist::D
    "Prior for the exponent on the destination population."
    α::A
    "Prior for the exponent on the origin population."
    β::B
    "Prior for the exponent on distance."
    γ::C
    "The own-stratum weight, relative to a typical pairwise term."
    within::W
end

function Gravity(pop, dist; α = Normal(0, 0.5), β = Normal(0, 0.5),
        γ = truncated(Normal(2, 0.5), 0, Inf), within = 1.0)
    return Gravity(pop, dist, α, β, γ, within)
end

@model function as_turing_model(m::Gravity, n)
    α ~ as_turing_submodel(m.α, 1; prefix = true)
    β ~ as_turing_submodel(m.β, 1; prefix = true)
    γ ~ as_turing_submodel(m.γ, 1; prefix = true)
    return gravity(m.pop, m.dist; α = only(α), β = only(β), γ = only(γ),
        within = m.within)
end

@doc raw"
Fold per-pair generation intervals into a coupling operator.

Returns the `strata × strata × lags` array [`renewal_pressure`](@ref) takes,
with entry `[g, h, i]` equal to `K[g, h] * G[g, h, i]`: the weight stratum `h`
puts on stratum `g` at lag `i`. A between-stratum transmission can therefore
carry a longer effective interval than a within-stratum one, which a single
shared interval cannot express.

Lags are indexed forwards, so `[:, :, 1]` is lag 1.

# Arguments

  - `K`: the `strata × strata` coupling weights.
  - `G`: the per-pair generation intervals, either a `strata × strata × lags`
    array whose `[g, h, :]` is the interval from `h` to `g`, or a single vector
    used for every pair.

# Examples
```@example pairwise_gen_int
using ComposableTuringIDModels
K = [0.9 0.1; 0.05 0.95]
size(pairwise_gen_int(K, [0.2, 0.3, 0.5]))
```
"
function pairwise_gen_int(K::AbstractMatrix, G::AbstractArray{<:Any, 3})
    @assert size(G)[1:2]==size(K) "`G` must be $(size(K)) x lags to match `K`"
    for idx in CartesianIndices(K)
        _assert_pmf(view(G, idx[1], idx[2], :))
    end
    return K .* G
end

function pairwise_gen_int(K::AbstractMatrix, g::AbstractVector)
    _assert_pmf(g)
    return pairwise_gen_int(
        K, repeat(reshape(g, 1, 1, :), size(K, 1), size(K, 2)))
end

function _assert_pmf(g)
    @assert all(>=(0), g) "A generation interval must be non-negative"
    @assert sum(g)≈1 "A generation interval must sum to 1"
    return nothing
end
