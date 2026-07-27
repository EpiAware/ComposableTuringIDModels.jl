# [Getting started](@id getting-started)

`ComposableTuringIDModels` assembles an infectious disease model from
interchangeable latent, infection, and observation parts, each turned into a
[Turing](https://turinglang.org) model by the single [`as_turing_model`](@ref)
constructor.
This package is early-stage and under active development; expect breaking
changes.

## Installation

Install the package from the Julia General registry:

```julia
using Pkg
Pkg.add("ComposableTuringIDModels")
```

Then load it alongside `Distributions` and `Turing`, which supply the priors and
the samplers:

```julia
using ComposableTuringIDModels, Distributions, Turing
```

## Where to go next

  - The [home page](../index.md) composes a model, simulates from its prior, and
    fits it, end to end.
  - The [Overview](@ref overview) introduces the three roles and the one
    interface they share.
  - [Composable design](@ref) explains how parts nest as submodels.
  - The [case studies](@ref case-studies-overview) fit complete models to real
    surveillance data.
  - The [Public API](@ref public-api) lists every component you can compose.
