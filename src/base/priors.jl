# The composition seam.
# `as_turing_model` gains `Distribution` and vector-of-`Distribution` methods so
# a raw prior flows through it identically to a full model.

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
    can never collide with the host's own variables;
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
# `θ ~ as_turing_submodel(model.slot, n)` is a plain native scalar draw named
# `θ` with no submodel overhead and no `.θ` namespace.
# A component consumes whichever it gets with [`at`](@ref), so a `Distribution`
# keeps its clean constant name while a process makes the parameter vary.

as_turing_submodel(d::Distribution, ::ModelShape; prefix::Bool = false) = d

function as_turing_submodel(
        v::AbstractVector{<:Distribution}, n::Int; prefix::Bool = false
    )
    @assert length(v) == n "a length-$(length(v)) prior vector cannot produce a length-$n prior"
    # `product_distribution` unconditionally, because `arraydist`'s deprecated
    # `Distributions.Product` path breaks Enzyme reverse and choosing `filldist`
    # at runtime returns a `Union` its type analysis rejects.
    return product_distribution(v)
end

@doc raw"
Widen a raw prior into a well-defined length-`n` PATH prior.

A bare `Distribution` given to a length-`n` PATH slot (an innovation, or a
latent process such as [`Hierarchy`](@ref)'s `across`) is wrapped in an
[`Intercept`](@ref) so it is a **constant path** (one shared draw broadcast to
length `n`), not a scalar. A process, an explicit [`IID`](@ref)/`Intercept`,
or a vector passes through unchanged.

A component author writing a new PATH-slot constructor calls this on the raw
argument, exactly as [`AR`](@ref), [`MA`](@ref) and [`Hierarchy`](@ref) do for
their innovation/across slots. **Per-step PARAMETER slots** (`damp`, `θ`,
`std`, …) keep the bare `Distribution` — a scalar constant — and must **not**
use this; use [`IID`](@ref) for `n` independent draws instead.

# Examples
```@example path_prior
using ComposableTuringIDModels, Distributions
ComposableTuringIDModels.path_prior(Normal()) isa Intercept
```
"
path_prior(p::Distribution) = Intercept(p)
path_prior(p) = p

# Deprecated alias for `path_prior`; the same generic function.
const _path_prior = path_prior

@doc raw"
The types accepted in a prior / process slot: a raw `Distribution`, a vector of
`Distribution`s, or an [`AbstractPriorModel`](@ref) (a latent process used as a
prior). Bounding a widened slot to `PriorLike` keeps the fail-fast role guard — a
wrong-role component (an observation or infection model) is rejected at
construction — while accepting a bare distribution alongside a process.
"
const PriorLike = Union{
    Distribution, AbstractVector{<:Distribution},
    AbstractPriorModel,
}

@doc raw"
Sample a raw prior `Distribution` as a **single scalar** RV.

Giving `as_turing_model` a `Distribution` method lets a **bare distribution** flow
through [`as_turing_submodel`](@ref) exactly like a full model, so a component's
parameter slot samples `θ ~ as_turing_submodel(model.slot, n)` uniformly whether
the slot holds a bare distribution or a process. A bare distribution draws ONE
scalar value (a constant, no length-`n` allocation) whatever shape it is asked
for — `n` is ignored, whether it is a length or an `(n_strata, n_time)` shape.
A component then reads a possibly-time-varying parameter per step with
[`at`](@ref), so the scalar stays constant while a process-valued slot
varies — this is the single seam behind [`AR`](@ref)'s optionally-time-varying
damping and the other per-step parameters.

For `n` **independent** draws (a white-noise process) use the explicit
[`IID`](@ref) component; for a **single shared** value broadcast to length `n` use
[`Intercept`](@ref); for **per-element** priors use a `Vector{<:Distribution}`.

# Arguments

  - `prior`: the prior distribution.
  - `n`: accepted for a uniform seam signature; ignored (the draw is scalar
    whatever shape is asked for).

# Examples
```@example as_turing_model_distribution
using ComposableTuringIDModels, Distributions
as_turing_model(Normal(), 3)()   # a single scalar draw
```
"
@model function as_turing_model(prior::Distribution, n::ModelShape)
    θ ~ prior
    return θ
end

@doc raw"
Sample a vector of prior `Distribution`s as a length-`n` prior submodel, one
independent draw per element.

The length is fixed by the vector, so `n` must match it. This is the explicit
way a prior slot asks for `n` independent draws with per-element priors (e.g. an
`AR`'s per-lag damping coefficients).

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
    @assert length(prior) == n "a length-$(length(prior)) prior vector cannot produce a length-$n prior"
    # One i.i.d. draw per element. See `as_turing_submodel`'s vector method
    # for why this is `product_distribution` and not `arraydist`.
    product_dist = product_distribution(prior)
    θ ~ product_dist
    return θ
end

# --- Time-varying-capable parameters ---------------------------------------
#
# The seam above draws a parameter slot as either a scalar or a length-`n` path.
# A component reads the result per step with `at`, so one recursion serves both
# the constant and the time-varying case with no per-component special-casing.

@doc raw"
Read a possibly-time-varying parameter at step `t`.

A scalar (a constant parameter, drawn from a `Distribution` prior through the
single [`as_turing_submodel`](@ref) seam) is returned unchanged at every step; a
vector (a per-step path, drawn from a process prior) is indexed at `t`; a
matrix (a strata × time parameter) is indexed at column `t`. A component's
recursion writes `at(ρ, t) * …` so the *same* code serves a constant and a
time-varying (or strata-varying) parameter — the scalar branch is zero-cost
(no per-step allocation).

This is the read side of the widening seam a custom component uses to make any
per-step parameter optionally time-varying: draw the slot through
[`as_turing_submodel`](@ref), then read it per step with `at`. See [Time-varying
damping in an AR process](@ref tutorial-tvdamp) for the worked example.

# Examples
```@example at
using ComposableTuringIDModels
(ComposableTuringIDModels.at(0.5, 3), ComposableTuringIDModels.at([0.1, 0.2, 0.3], 2))
```
"
at(p::Number, t) = p
at(p::AbstractVector, t) = p[t]
# A strata × time parameter read at step `t` is that step's column.
at(p::AbstractMatrix, t) = view(p, :, t)

# Deprecated alias for `at`; the same generic function.
const _at = at

@doc raw"
The order (`p`/`q`/`d`) implied by a prior slot.

A vector of `Distribution`s fixes the order to the vector length (one
independent per-lag/per-element prior); a single `Distribution` or a richer
prior model (an [`AbstractPriorModel`](@ref)) defaults to order 1. This is
how [`AR`](@ref) and [`MA`](@ref) infer their order from `damp`/`θ` without a
separate order argument; a component author defining a similar per-lag slot
uses it the same way.

# Examples
```@example prior_order
using ComposableTuringIDModels, Distributions
(ComposableTuringIDModels.prior_order(Normal()),
    ComposableTuringIDModels.prior_order([Normal(), Normal()]))
```
"
prior_order(p::AbstractVector{<:Distribution}) = length(p)
prior_order(::Distribution) = 1
prior_order(::AbstractPriorModel) = 1

const _prior_order = prior_order

# The order of a slot that may not hold a prior at all. A non-prior is passed
# through to the slot's own `PriorLike` bound, which rejects it by role, rather
# than raising a `MethodError` about reading an order it never had.
_slot_order(p::PriorLike) = prior_order(p)
_slot_order(_) = 1

@doc raw"
Assert that a vector-of-`Distribution`s prior has exactly `k` elements.

Pairs with [`prior_order`](@ref): once a slot has fixed the order `k` (e.g.
from `damp`), a second per-lag/per-element slot given as a vector (e.g.
`init`) must match it. A single `Distribution` or a richer prior model
broadcasts to `k` and imposes no constraint, so it always passes.

# Arguments

  - `p`: the prior for the slot being checked. A vector of `Distribution`s is
    length-checked; anything else broadcasts and always passes.
  - `k`: the required number of elements, fixed earlier by another slot.
  - `what`: a short description of the slot, used in the assertion message.

# Examples
```@example assert_prior_length
using ComposableTuringIDModels, Distributions
ComposableTuringIDModels.assert_prior_length([Normal(), Normal()], 2, :damp)
```
"
function assert_prior_length(p::AbstractVector{<:Distribution}, k, what)
    @assert length(p) == k "$what prior length $(length(p)) must equal $k"
    return nothing
end
assert_prior_length(_, k, what) = nothing

const _assert_prior_length = assert_prior_length
