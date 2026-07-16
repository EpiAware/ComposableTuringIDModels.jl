# The composition seam. A component threads its sub-components and its prior
# slots through Turing submodels; `as_turing_submodel` names that one pattern, and
# `as_turing_model` gains `Distribution` / vector-of-`Distribution` methods so a
# raw prior flows through the seam identically to a full model (issue #37: priors
# as length-`n` submodels).

@doc raw"
Compose a component as a Turing submodel: `to_submodel(as_turing_model(m,
args...), prefix)`.

This is the single public composition seam of the package. Every composition point
— a manipulator wrapping an inner model, an infection model owning its latent
process, a component sampling a vector/process-valued prior slot — threads its
sub-component through here, and third-party component authors use it as *the* way
to compose an `as_turing_model` inside their own `@model` body:

```julia
latent ~ as_turing_submodel(inner_model, n)
```

Because `as_turing_model` also has `Distribution` and
`Vector{<:Distribution}` methods, the same call composes a raw prior:

```julia
damp ~ as_turing_submodel(model.damp, p; prefix = true)   # a Distribution or a process
```

`prefix` defaults to `false` — the package standard, keeping the submodel's
variable names flat. Two kinds of call site pass `prefix = true`:

  - a **prior slot** (a component's `damp` / `init` / `θ` etc.), so the slot's
    left-hand name namespaces the whole prior submodel and a process-valued prior
    can never collide with the host's own variables (issue #80);
  - the **deliberately-prefixing** components ([`PrefixLatentModel`](@ref),
    [`Split`](@ref)), which stream their children under an explicit name.

# Arguments

  - `m`: the component (or raw prior) to compose.
  - `args...`: positional arguments forwarded to `as_turing_model` (e.g. the
    series length `n`, or the observed/expected series for an observation model).

# Keyword Arguments

  - `prefix`: whether to prefix the submodel's variables with the tilde
    left-hand name (default `false`).

# Examples

Inside a component's `@model` body it is used on the right of a `~`:

```julia
latent ~ as_turing_submodel(inner_model, n)
damp ~ as_turing_submodel(model.damp, p; prefix = true)
```

It returns a Turing submodel; the underlying prior submodel returns a length-`n`
value:

```@example as_turing_submodel
using ComposableTuringIDModels, Distributions
length(as_turing_model(Normal(), 4)())
```
"
function as_turing_submodel(m, args...; prefix::Bool = false)
    return to_submodel(as_turing_model(m, args...), prefix)
end

@doc raw"
The types accepted in a prior / process slot: a raw `Distribution`, a vector of
`Distribution`s, or an [`AbstractPriorModel`](@ref) (a latent process used as a
prior). Bounding a widened slot to `PriorLike` keeps the fail-fast role guard — a
wrong-role component (an observation or infection model) is rejected at
construction — while accepting a bare distribution alongside a process.
"
const PriorLike = Union{Distribution, AbstractVector{<:Distribution},
    AbstractPriorModel}

@doc raw"
Sample a raw prior `Distribution` as a length-`n` prior submodel.

Giving `as_turing_model` a `Distribution` method lets a **bare distribution** flow
through [`as_turing_submodel`](@ref) exactly like a full model: a component's prior
slot samples `θ ~ as_turing_submodel(model.slot, n)` whether the slot holds a bare
distribution or a latent process. The draw is `n` independent values from `prior`
(`filldist`); for a genuinely scalar parameter sample the distribution with a
native tilde instead (`σ ~ model.std`), which is a plain scalar draw with no
submodel. To broadcast a **single shared** value to length `n` use
[`Intercept`](@ref).

# Arguments

  - `prior`: the prior distribution.
  - `n`: the number of independent draws.

# Examples
```@example as_turing_model_distribution
using ComposableTuringIDModels, Distributions
length(as_turing_model(Normal(), 3)())
```
"
@model function as_turing_model(prior::Distribution, n::Int)
    θ ~ filldist(prior, n)
    return θ
end

@doc raw"
Sample a vector of prior `Distribution`s as a length-`n` prior submodel, one
independent draw per element.

The length is fixed by the vector, so `n` must match it. A homogeneous vector uses
`filldist` and a heterogeneous one `arraydist`. This is the explicit way a prior
slot asks for `n` independent draws with per-element priors (e.g. an `AR`'s
per-lag damping coefficients).

# Arguments

  - `prior`: the vector of prior distributions.
  - `n`: the required length (must equal `length(prior)`).

# Examples
```@example as_turing_model_vector
using ComposableTuringIDModels, Distributions
as_turing_model([Normal(0, 1), Normal(5, 0.1)], 2)()
```
"
@model function as_turing_model(prior::AbstractVector{<:Distribution}, n::Int)
    @assert length(prior)==n "a length-$(length(prior)) prior vector cannot produce a length-$n prior"
    # One i.i.d. draw per element; `filldist` for a homogeneous vector, `arraydist`
    # otherwise.
    product_dist = all(first(prior) .== prior) ? filldist(first(prior), n) :
                   arraydist(prior)
    θ ~ product_dist
    return θ
end

@doc raw"
Mark a prior slot as supplying a *time-varying* coefficient path rather than a
constant.

Wrapping a prior in `TimeVarying` is the opt-in that turns a slot's draw from a
single constant into a whole length-`(n-1)` path: the inner `prior` is sampled at
that length and mapped by `transform` before it is threaded through the recursion.
It is currently understood by [`AR`](@ref)'s `damp` slot — `AR(damp =
TimeVarying(RandomWalk()))` is a first-order AR whose damping coefficient changes
over the series (see [`TimeVaryingAR`](@ref) for the named constructor) — whereas a
bare process in `damp` only enriches the *prior* over a constant coefficient.

# Fields

  - `prior`: the prior process for the raw (pre-`transform`) path (a `Distribution`
    or a prior model such as a [`RandomWalk`](@ref)).
  - `transform`: map from the raw path to the coefficient (default `tanh`, giving
    the stationary band ``(-1, 1)``; pass `identity` to opt out).
"
struct TimeVarying{P <: PriorLike, F <: Function}
    "Prior process for the raw (pre-`transform`) path."
    prior::P
    "Map from the raw path to the coefficient (default `tanh`)."
    transform::F
end

TimeVarying(prior; transform = tanh) = TimeVarying(prior, transform)

# Order (p / q / d) implied by a prior: a vector fixes it to the vector length; a
# single distribution or a richer prior model defaults to order 1. A `TimeVarying`
# slot is always order 1 (only AR(1) supports a time-varying coefficient; higher
# orders are tracked in #113).
_prior_order(p::AbstractVector{<:Distribution}) = length(p)
_prior_order(::Distribution) = 1
_prior_order(::AbstractPriorModel) = 1
_prior_order(::TimeVarying) = 1

# Assert a vector prior's length matches the required order `k`. A single
# distribution, a richer prior model, or a `TimeVarying` marker broadcasts to `k`
# and imposes no constraint.
function _assert_prior_length(p::AbstractVector{<:Distribution}, k, what)
    @assert length(p)==k "$what prior length $(length(p)) must equal $k"
    return nothing
end
_assert_prior_length(_, k, what) = nothing
