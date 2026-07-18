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

# --- Time-varying-capable parameters ---------------------------------------
#
# A parameter that may vary over time is drawn from its prior slot as EITHER a
# scalar (a bare `Distribution` ⇒ one native-tilde draw, a constant, no length-n
# allocation) OR a length-`n` path (an `AbstractPriorModel` process ⇒ drawn
# through `as_turing_submodel`). A component then reads it per step with `_at`, so
# a single recursion serves both the constant and the time-varying case with no
# per-component special-casing and no efficiency loss when the parameter is
# constant. This mechanism is general: any component with a scalar parameter can
# widen it to a time-varying one just by drawing the slot through
# [`as_timevarying_submodel`](@ref) and consuming it with [`_at`](@ref) (see
# [`AR`](@ref) for the worked example).

@doc raw"
Read a possibly-time-varying parameter at step `t`.

A scalar (a constant parameter, drawn from a `Distribution` prior) is returned
unchanged at every step; a vector (a per-step path, drawn from a process prior) is
indexed at `t`. A component's recursion writes `_at(ρ, t) * …` so the *same* code
serves a constant and a time-varying parameter — the scalar branch is zero-cost
(no per-step allocation), matching [`as_timevarying_submodel`](@ref).
"
_at(p::Number, t) = p
_at(p::AbstractVector, t) = p[t]

# The drawing half of the time-varying mechanism. A bare `Distribution` is a
# single scalar RV (constant, efficient); a process prior is a length-`n` path.
# Both are exposed through `as_timevarying_submodel` so a component draws the slot
# uniformly regardless of which was supplied.
@model function _timevarying_prior(prior::Distribution, n::Int)
    θ ~ prior
    return θ
end

@model function _timevarying_prior(prior::AbstractPriorModel, n::Int)
    θ ~ as_turing_submodel(prior, n)
    return θ
end

@doc raw"
Draw a possibly-time-varying parameter from its prior slot as a submodel.

The companion of [`_at`](@ref): a bare `Distribution` slot yields ONE scalar RV (a
constant parameter, no length-`n` allocation) while an [`AbstractPriorModel`](@ref)
process slot yields a length-`n` path. A component draws its slot through this seam
and consumes the result per step with [`_at`](@ref), so one recursion serves both
the constant and the time-varying case. This is the general way any component
widens a scalar parameter to an optionally-time-varying one; [`AR`](@ref) is the
worked example (its `damp` coefficient).

# Arguments

  - `prior`: the parameter's prior — a `Distribution` (constant) or a process.
  - `n`: the path length used when the prior is a process.

# Keyword Arguments

  - `prefix`: whether to namespace the submodel's variables under the tilde's
    left-hand name (default `false`; a prior slot passes `true`).
"
function as_timevarying_submodel(prior, n; prefix::Bool = false)
    return to_submodel(_timevarying_prior(prior, n), prefix)
end

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
