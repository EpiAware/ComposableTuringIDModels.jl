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
    a heterogeneous one `arraydist`. This is the explicit way to ask for `n`
    independent draws.

Because `BroadcastPrior` is itself an [`AbstractPriorModel`](@ref), anything richer
a user writes (a partially-pooled prior, a time-varying prior wrapping a
[`RandomWalk`](@ref)) drops into the same prior slot with no struct changes.

The sampled variable is always the generic `θ`; the component that samples the
prior namespaces it by turning submodel prefixing on at the slot (the left-hand
name becomes the chain prefix), so two priors in one model never collide.

## Fields

  - `dist`: a `Distribution`, or a `Vector{<:Distribution}`.

# Examples
```@example BroadcastPrior
using ComposableTuringIDModels, Distributions
# scalar parameter: length-1, read back with `only`
only(as_turing_model(BroadcastPrior(Normal()), 1)())
```
"
struct BroadcastPrior{D} <: AbstractPriorModel
    "A `Distribution`, or a `Vector{<:Distribution}`."
    dist::D
end

@model function as_turing_model(prior::BroadcastPrior{<:Distribution}, n)
    @assert n>0 "n must be greater than 0"
    # Repeat-one: a single random variable repeated to length `n`, so a global
    # parameter is not expanded into `n` i.i.d. draws.
    θ ~ prior.dist
    return fill(θ, n)
end

@model function as_turing_model(prior::BroadcastPrior{<:AbstractVector}, n)
    @assert length(prior.dist)==n "BroadcastPrior of a length-$(length(prior.dist)) vector cannot produce a length-$n prior"
    # One i.i.d. draw per element; `filldist` for a homogeneous vector, `arraydist`
    # otherwise.
    product_dist = all(first(prior.dist) .== prior.dist) ?
                   filldist(first(prior.dist), n) : arraydist(prior.dist)
    θ ~ product_dist
    return θ
end

@doc raw"
Coerce a user-supplied prior into an [`AbstractPriorModel`](@ref).

This is the seam that keeps constructors ergonomic: a prior field is typed
`::AbstractPriorModel` and the constructor calls `as_prior` on whatever the user
passed, so a bare `Distribution` (or a vector of them) is wrapped in a
[`BroadcastPrior`](@ref) while a prior submodel — including any latent process,
e.g. `RandomWalk()` for a time-varying parameter — is accepted unchanged.

Collision-free naming is handled at the call site rather than here: a component
turns submodel prefixing on when it samples the prior, so the slot's left-hand
name namespaces the whole prior submodel (a process prior's inner variables can
never collide with the host's own — issue #80).

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
