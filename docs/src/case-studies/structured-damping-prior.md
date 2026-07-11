# [A structured prior on the AR damping](@id case-study-ar-structured-damping)

Every parameter slot of a component is an [`AbstractPriorModel`](@ref): as well as
a bare `Distribution`, a slot accepts a *prior model* — including a latent process
— so a coefficient can carry a structured prior instead of a single global draw.
This page uses that on the damping of an autoregressive process, the workhorse
inside [`AR`](@ref), [`arma`](@ref) and [`arima`](@ref), entirely through the
struct's own constructor with no custom model code.

This is early-stage, actively developed software; the API may change.

!!! note "Scope: a structured prior, not a time-varying coefficient"
    This structures the *prior* over the damping coefficient; the coefficient is
    still **constant in time**. An [`AR`](@ref) applies its damping as a fixed
    length-`p` coefficient vector (`accumulate_scan(ARStep(damp), …)`), so a
    process supplied to `damp` enriches the prior over that vector rather than
    making the coefficient vary over the series. Genuinely time-varying AR
    coefficients (a length-`n` ``\rho_t`` path threaded through the recursion) are
    a separate feature tracked in [#113](https://github.com/EpiAware/ComposableTuringIDModels.jl/issues/113).

## The damping slot takes any prior model

Before the priors weave, a component's damping was a bare `Distribution` sampled
once. Now the [`AR`](@ref) `damp` slot is an [`AbstractPriorModel`](@ref): a bare
`Distribution` is coerced to the default wrapper, while a prior model — such as a
[`HierarchicalNormal`](@ref) whose scale is itself inferred — is accepted
unchanged and gives the coefficient a hierarchical (adaptive-scale) prior. Both
forms go through the same constructor:

```@example damping
using ComposableTuringIDModels, Distributions, Turing, Random, Statistics
using Turing: to_submodel, returned
using DynamicPPL: LogDensityFunction, VarInfo, link, getlogjoint
import LogDensityProblems as LDP
Random.seed!(80)

flat = AR(; damp = truncated(Normal(0.4, 0.1), 0, 1))   # a plain-Distribution prior
structured = AR(; damp = HierarchicalNormal())          # a hierarchical prior
(flat_order = flat.p, structured_order = structured.p)
```

The same works through [`arima`](@ref) — the damping prior threads down the
[`arima`](@ref) → [`arma`](@ref) → [`AR`](@ref) nesting stack unchanged:

```@example damping
mdl = arima(; damp = HierarchicalNormal())
sampled = rand(as_turing_model(mdl, 12))
# the damping carries the hierarchical structure (its own inferred scale)
filter(v -> occursin("damp", string(v)), collect(keys(sampled)))
```

## It threads under a gradient sampler

A latent-model prior is auto-prefixed inside the [`as_prior`](@ref) coercion seam,
so its internal variables (`std`, `ϵ_t`) cannot collide with the AR innovation's
own latent under the prefix-off submodel convention. Without that prefixing the
bare form sampled via `rand` but errored when evaluated as a *linked* log-density
(the target a gradient sampler differentiates):

```@example damping
m = as_turing_model(AR(; damp = HierarchicalNormal()), 40)
vi = link(VarInfo(m), m)
ldf = LogDensityFunction(m, getlogjoint, vi)
(rand_ok = rand(m) !== nothing,
    linked_logdensity_finite = isfinite(
        LDP.logdensity(ldf, zeros(LDP.dimension(ldf)))))
```

So the structured-damping model draws a value *and* evaluates as a linked
log-density, and samples under NUTS end-to-end:

```@example damping
chain = sample(m, NUTS(0.8; adtype = Turing.AutoForwardDiff()), 300;
    progress = false)
size(chain, 1)
```

The damping coefficient is recovered alongside its inferred prior scale, all from
the [`AR`](@ref) struct's own [`as_turing_model`](@ref) with no bespoke recursion.
The same slot takes any prior process: swapping the prior over the damping (a
different scale prior, an [`IID`](@ref) draw, a [`RandomWalk`](@ref)) is a one-line
change to the `damp` argument, and the AR recursion is untouched.
