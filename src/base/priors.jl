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

# --- Single-seam specialisations (clean names for the constant case) ---------
#
# A bare `Distribution` returns the distribution itself, so a component's
# `θ ~ as_turing_submodel(model.slot, n)` is a plain native scalar draw named `θ`
# (a single constant RV, zero submodel overhead, no `.θ` namespace). A vector of
# distributions returns the product distribution (a native, per-element draw). A
# process (or any other model) returns a namespaced submodel via the generic
# method above — the length-`n`, e.g. time-varying / hierarchical, path. A
# component consumes whichever it gets with [`_at`](@ref), so supplying a process
# makes the parameter vary with no rewiring while a `Distribution` keeps its clean
# constant name. `n` is ignored for the scalar case.

as_turing_submodel(d::Distribution, ::Int; prefix::Bool = false) = d

function as_turing_submodel(
        v::AbstractVector{<:Distribution}, n::Int; prefix::Bool = false)
    @assert length(v)==n "a length-$(length(v)) prior vector cannot produce a length-$n prior"
    # `filldist` for a homogeneous vector, `arraydist` otherwise.
    return all(first(v) .== v) ? filldist(first(v), n) : arraydist(v)
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
Sample a raw prior `Distribution` as a **single scalar** RV.

Giving `as_turing_model` a `Distribution` method lets a **bare distribution** flow
through [`as_turing_submodel`](@ref) exactly like a full model, so a component's
parameter slot samples `θ ~ as_turing_submodel(model.slot, n)` uniformly whether
the slot holds a bare distribution or a process. A bare distribution draws ONE
scalar value (a constant, no length-`n` allocation); `n` is ignored. A component
then reads a possibly-time-varying parameter per step with [`_at`](@ref), so the
scalar stays constant while a process-valued slot varies — this is the single seam
behind [`AR`](@ref)'s optionally-time-varying damping and the other per-step
parameters.

For `n` **independent** draws (a white-noise process) use the explicit
[`IID`](@ref) component; for a **single shared** value broadcast to length `n` use
[`Intercept`](@ref); for **per-element** priors use a `Vector{<:Distribution}`.

# Arguments

  - `prior`: the prior distribution.
  - `n`: accepted for a uniform seam signature; ignored (the draw is scalar).

# Examples
```@example as_turing_model_distribution
using ComposableTuringIDModels, Distributions
as_turing_model(Normal(), 3)()   # a single scalar draw
```
"
@model function as_turing_model(prior::Distribution, n::Int)
    θ ~ prior
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

# --- Time-varying-capable parameters ---------------------------------------
#
# The single seam above draws a parameter slot as EITHER a scalar (a bare
# `Distribution` ⇒ one RV, a constant, no length-`n` allocation) OR a length-`n`
# path (an `AbstractPriorModel` process ⇒ the process's own draw). A component
# reads the result per step with `_at`, so ONE recursion serves both the constant
# and the time-varying (or hierarchical) case with no per-component special-casing
# and no efficiency loss when the parameter is constant. This is general: any
# per-step parameter is widened to optionally-time-varying just by drawing its slot
# through [`as_turing_submodel`](@ref) and consuming it with [`_at`](@ref) — see
# [`AR`](@ref)'s damping for the worked example.

@doc raw"
Read a possibly-time-varying parameter at step `t`.

A scalar (a constant parameter, drawn from a `Distribution` prior through the
single [`as_turing_submodel`](@ref) seam) is returned unchanged at every step; a
vector (a per-step path, drawn from a process prior) is indexed at `t`. A
component's recursion writes `_at(ρ, t) * …` so the *same* code serves a constant
and a time-varying parameter — the scalar branch is zero-cost (no per-step
allocation).
"
_at(p::Number, t) = p
_at(p::AbstractVector, t) = p[t]

# Order (p / q / d) implied by a prior: a vector fixes it to the vector length; a
# single distribution or a richer prior model defaults to order 1.
_prior_order(p::AbstractVector{<:Distribution}) = length(p)
_prior_order(::Distribution) = 1
_prior_order(::AbstractPriorModel) = 1

# Assert a vector prior's length matches the required order `k`. A single
# distribution or a richer prior model broadcasts to `k` and imposes no
# constraint.
function _assert_prior_length(p::AbstractVector{<:Distribution}, k, what)
    @assert length(p)==k "$what prior length $(length(p)) must equal $k"
    return nothing
end
_assert_prior_length(_, k, what) = nothing
