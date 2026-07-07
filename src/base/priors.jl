# The default prior wrapper and the `as_prior` coercion seam. This is the
# foundation for issue #37 (priors as length-`n` submodels): it lets every prior
# join the same `as_turing_model` protocol as the other components, fronted by a
# wrapper that turns a plain `Distribution` into a length-`n` prior submodel.

@doc raw"
Default [`AbstractPriorModel`](@ref) wrapper: turn a plain `Distribution` (or a
vector of `Distribution`s) into a length-`n` prior submodel.

This is the drop-in replacement for a bare-`Distribution` prior field. It has two
modes, dispatching on what it wraps:

  - **A single `Distribution` — *repeat-one*.** One value is drawn and repeated to
    length `n`: `as_turing_model(BroadcastPrior(d), n)` samples a single `θ ~ d`
    and returns `fill(θ, n)`. A single global coefficient stays a single random
    variable — it is *not* silently turned into `n` i.i.d. draws. For a genuinely
    scalar parameter use `n == 1` and read the element with `only(...)`, so the
    chain stays as small as a bare `~ dist`.
  - **A vector of `Distribution`s — one i.i.d. draw per element.** The length is
    fixed by the vector (`n` must match). A homogeneous vector uses `filldist` and
    a heterogeneous one `arraydist` — reproducing the eager `_expand_dist` helper,
    but as a submodel. This is the explicit way to ask for `n` independent draws.

Because `BroadcastPrior` is itself an [`AbstractPriorModel`](@ref), anything richer
a user writes (a partially-pooled prior, a time-varying prior wrapping a
[`RandomWalk`](@ref)) drops into the same prior slot with no struct changes.

An optional `name` fixes the sampled variable's name in the chain. When a
component samples a prior as a submodel via `to_submodel(..., false)`, the
component sets `name` to its own parameter name (e.g. `:damp_AR`, `:std`,
`:init_incidence`) so the flat chain names are unchanged when the default wrapper
replaces a bare `~ dist`. It also keeps two broadcast priors in one model from
colliding on a shared internal name. When `name === nothing` (the default) the
variable is sampled under the generic name `θ`.

## Fields

  - `dist`: a `Distribution`, or a `Vector{<:Distribution}`.
  - `name`: the sampled variable name (a `Symbol`), or `nothing` for the generic
    `θ`.

# Examples
```@example BroadcastPrior
using ComposableTuringIDModels, Distributions
# scalar parameter: length-1, read back with `only`
only(as_turing_model(BroadcastPrior(Normal()), 1)())
```
"
struct BroadcastPrior{D, N <: Union{Symbol, Nothing}} <: AbstractPriorModel
    "A `Distribution`, or a `Vector{<:Distribution}`."
    dist::D
    "The sampled variable name, or `nothing` for the generic `θ`."
    name::N
end

BroadcastPrior(dist) = BroadcastPrior(dist, nothing)

# Sample from `dist` under the wrapper's fixed name, or the generic `θ` when no
# name is set. `NamedDist` fixes the chain variable name regardless of the
# tilde's left-hand side.
_named(dist, ::Nothing) = dist
_named(dist, name::Symbol) = NamedDist(dist, name)

@model function as_turing_model(prior::BroadcastPrior{<:Distribution}, n)
    @assert n>0 "n must be greater than 0"
    # Repeat-one: a single random variable repeated to length `n`, so a global
    # parameter is not expanded into `n` i.i.d. draws.
    θ ~ _named(prior.dist, prior.name)
    return fill(θ, n)
end

@model function as_turing_model(prior::BroadcastPrior{<:AbstractVector}, n)
    @assert length(prior.dist)==n "BroadcastPrior of a length-$(length(prior.dist)) vector cannot produce a length-$n prior"
    # One i.i.d. draw per element; `filldist` for a homogeneous vector (as
    # `_expand_dist` did), `arraydist` otherwise.
    product_dist = all(first(prior.dist) .== prior.dist) ?
                   filldist(first(prior.dist), n) : arraydist(prior.dist)
    θ ~ _named(product_dist, prior.name)
    return θ
end

@doc raw"
Coerce a user-supplied prior into an [`AbstractPriorModel`](@ref).

This is the seam that keeps constructors ergonomic: a prior field is typed
`::AbstractPriorModel` and the constructor calls `as_prior` on whatever the user
passed, so a bare `Distribution` (or a vector of them) is wrapped in a
[`BroadcastPrior`](@ref) while a prior submodel — including any
[`AbstractLatentModel`](@ref), e.g. `RandomWalk()` for a time-varying parameter —
is accepted unchanged.

# Arguments

  - `p`: an `AbstractPriorModel`, a `Distribution`, or a `Vector{<:Distribution}`.

# Examples
```@example as_prior
using ComposableTuringIDModels, Distributions
as_prior(Normal())        # a BroadcastPrior
as_prior(RandomWalk())    # a latent model is already a prior model
```
"
as_prior(p::AbstractPriorModel) = p
as_prior(d::Distribution) = BroadcastPrior(d)
as_prior(v::AbstractVector{<:Distribution}) = BroadcastPrior(v)

# Name-carrying coercion used by component constructors: a bare `Distribution`
# (or vector) is wrapped in a `BroadcastPrior` that samples under the component's
# own parameter `name`, so the flat chain names are unchanged. A non-latent prior
# model (e.g. a `BroadcastPrior`) already carries its own naming, so it is passed
# through unchanged.
as_prior(p::AbstractPriorModel, ::Symbol) = p
# A latent model used directly as another component's prior carries its own inner
# variable names (`std`, `ϵ_t`, `rw_init`, …). Under the prefix-off submodel
# convention these collide with the host component's own latent, so a linked
# log-density mis-threads the flattened parameter vector (issue #80): e.g. a bare
# `AR(damp = RandomWalk())` samples via `rand` but errors as a linked log-density.
# Auto-prefix the latent-model prior with the component's parameter `name` inside
# the coercion seam, so the bare form threads without a manual `PrefixLatentModel`.
as_prior(p::AbstractLatentModel, name::Symbol) = PrefixLatentModel(p, String(name))
as_prior(d::Distribution, name::Symbol) = BroadcastPrior(d, name)
function as_prior(v::AbstractVector{<:Distribution}, name::Symbol)
    return BroadcastPrior(v, name)
end

# Order (p / q / d) implied by a prior: a vector wrapper fixes it to the vector
# length; a single-distribution (repeat-one) or richer prior defaults to order 1.
_prior_order(p::BroadcastPrior{<:AbstractVector}) = length(p.dist)
_prior_order(::AbstractPriorModel) = 1

# Assert a vector wrapper's length matches the required order `k`. A
# single-distribution or richer prior broadcasts to `k` and imposes no
# constraint.
function _assert_prior_length(p::BroadcastPrior{<:AbstractVector}, k, what)
    @assert length(p.dist)==k "$what prior length $(length(p.dist)) must equal $k"
    return nothing
end
_assert_prior_length(::AbstractPriorModel, k, what) = nothing
