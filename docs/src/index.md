# ComposableTuringIDModels

!!! warning "Prototype"
    This package is an exploratory **prototype** for composable probabilistic
    infectious disease modelling. Treat everything here as exploratory and
    subject to change.

`ComposableTuringIDModels` builds epidemiological models from small, reusable
components — infection processes (each owning its own latent process) and
observation models — and turns each one into a
[Turing](https://turinglang.org) / `DynamicPPL` model through the single generic
constructor [`as_turing_model`](@ref). Components compose by sampling one another
as submodels, so a full model is *assembled* from parts rather than written by
hand.

## Getting started

```@example index
using ComposableTuringIDModels, Distributions

# Compose infections -> observations. The infection model owns its latent
# process (here a random walk in its `Z` slot); the latent is folded in rather
# than threaded as a separate top-level component.
model = IDModel(
    DirectInfections(; Z = RandomWalk(), initialisation_prior = Normal()),
    PoissonError())

# Build a Turing model; `missing` observations simulate from the prior.
turing_model = as_turing_model(model, missing, 20)
draw = rand(turing_model)
nothing # hide
```

The composed model exposes its generated quantities directly:

```@example index
(; generated_y_t, I_t, Z_t) = turing_model()
length(generated_y_t), length(I_t), length(Z_t)
```

See [Composable design](@ref) for how the pieces fit together, the
[case studies](@ref case-studies-overview) for worked end-to-end examples, and the
[Public API](@ref public-api) for the full surface.

## [Automatic differentiation backend](@id ad-backend)

Hamiltonian samplers such as NUTS need gradients of the log density, computed by
automatic differentiation (AD). The recommended default for this package is
[Mooncake](https://chalk-lab.github.io/Mooncake.jl/), a reverse-mode backend that
handles the `DynamicPPL` `@model` + `to_submodel` composition used throughout
well. Select it by passing an `adtype` to the sampler:

```julia
using Turing, Mooncake
using ADTypes: AutoMooncake
chain = sample(model, NUTS(; adtype = AutoMooncake(; config = nothing)), 1_000)
```

The examples in [Composable design](@ref) and the renewal/delay
[case studies](@ref case-studies-overview) all use Mooncake. Other backends —
`ForwardDiff`, `ReverseDiff`, and `Enzyme` — are supported and exercised by the
package's AD test suite; swap them in with `AutoForwardDiff()`, `AutoReverseDiff()`,
or `AutoEnzyme()`. ForwardDiff is the right choice for the ODE models in the
[SIR case study](@ref case-study-sir): Mooncake-driven NUTS through an ODE solver
is not available yet (see that page).

## Adapted from

The modelling code in this package is **ported and adapted** from the
open-source, Apache-2.0 licensed `EpiAware` package
([CDCgov/Rt-without-renewal](https://github.com/CDCgov/Rt-without-renewal),
ported from the fork
[seabbs/Rt-without-renewal](https://github.com/seabbs/Rt-without-renewal)).
`ComposableTuringIDModels` is a modified, derived work: renamed, re-architected around
`as_turing_model`, and upgraded to the latest Turing. See the `NOTICE` file for
attribution and the list of changes.
