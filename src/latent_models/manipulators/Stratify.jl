# Give a shared latent path a strata axis: a per-stratum offset (or
# time-varying deviation) broadcast against a path drawn once over time.

@doc raw"
Add a strata axis to a shared latent path.

`Stratify(shared, across)` draws `shared` once over the time axis and `across`
over the strata axis, then broadcasts the two together with `combine` into a
`n_strata × n_time` matrix,

```math
Z_{g,t} = \text{combine}(a_{g}, s_t), \qquad
g = 1, \ldots, n_{\text{strata}}, \quad t = 1, \ldots, n_{\text{time}},
```

with ``s`` the `shared` draw and ``a`` the `across` draw. Built with
`as_turing_model(m, (n_strata, n_time))`, so the stratum count is **not** a
field of the struct: it arrives through the shape argument, exactly the way a
series length `n` is passed to a plain path model. One `Stratify` therefore
serves any panel width.

`across` is drawn at [`across_shape`](@ref)`(m.across, n)`, which is
`n_strata` by default: `across` is a length-`n_strata` vector, and `combine`
broadcasts it against a `1 × n_time` row, `(g,) .+ (1, t) ⇒ (g, t)` — a
constant offset per stratum. An `across` model that itself spans both axes
(e.g. [`Replicate`](@ref)) overrides `across_shape` to draw the full
`(n_strata, n_time)` matrix instead, and the same broadcast,
`(g, t) .+ (1, t) ⇒ (g, t)`, now gives a time-varying deviation per stratum.
One `combine` line therefore serves both a constant per-stratum offset and a
time-varying one, with no branch on which `across` was supplied.

The `across` slot is where the pooling choice lives:

  - [`Hierarchy`](@ref) partially pools the per-stratum deviations towards a
    shared level;
  - [`IID`](@ref) draws each deviation independently, giving no pooling;
  - [`FixedIntercept`](@ref)`(0.0)` fixes every deviation to zero, giving full
    pooling — every stratum then shares `shared` exactly;
  - a [`RandomWalk`](@ref) correlates neighbouring strata, useful when the
    strata are ordered (e.g. adjacent age bands).

Both slots take a raw prior (a bare `Distribution`, or a latent/prior model),
sampled through [`as_turing_submodel`](@ref); a bare `Distribution` is
auto-wrapped in an [`Intercept`](@ref) (a constant path), the same PATH
convention used elsewhere (e.g. [`RandomWalk`](@ref)'s `ϵ_t`).

## Fields

  - `shared`: the path drawn once over the time axis.
  - `across`: the cross-stratum relationship generating the per-stratum
    offsets (or deviations, for an `across` that spans both axes).
  - `combine`: how `across` and `shared` are broadcast together (default `+`,
    a multiplicative effect on a log scale).

# Examples
```@example Stratify
using ComposableTuringIDModels, Distributions
# A shared random walk in R_t-space with partially pooled per-stratum
# deviations; the strata count (3) arrives through the shape, not a field.
strat = Stratify(RandomWalk(), Hierarchy())
size(as_turing_model(strat, (3, 20))())
```
"
struct Stratify{S <: PriorLike, A <: PriorLike, F <: Function} <:
    AbstractPriorModel
    "The path drawn once over the time axis."
    shared::S
    "The cross-stratum relationship generating the per-stratum offsets."
    across::A
    "How `across` and `shared` are broadcast together."
    combine::F
end

function Stratify(shared, across; combine = +)
    return Stratify(path_prior(shared), path_prior(across), combine)
end

@model function as_turing_model(m::Stratify, n::Dims{2})
    n_strata, n_time = n
    shared ~ as_turing_submodel(m.shared, n_time)
    across ~ as_turing_submodel(
        m.across, across_shape(m.across, n); prefix = true
    )
    return m.combine.(across, permutedims(shared))
end

@doc raw"
The shape `Stratify` draws its `across` slot at.

The default method returns `n[1]`, the stratum count alone, so `across` draws
one value per stratum: a per-stratum offset, constant over time. An `across`
model that spans both axes (e.g. [`Replicate`](@ref)) overrides this method to
return the full shape `n`, so it draws a `(n_strata, n_time)` matrix instead —
a per-stratum deviation that varies over time. This is the extension point: a
new `across` model that needs a shape other than a bare stratum count adds a
method here.

# Arguments

  - `across`: the model in the `across` slot.
  - `n`: the `(n_strata, n_time)` shape the `Stratify` is being built at.

# Examples
```@example across_shape
using ComposableTuringIDModels
ComposableTuringIDModels.across_shape(Hierarchy(), (3, 20))
```
"
across_shape(across, n::Dims{2}) = n[1]
