# EpiAwarePrototype.jl

> **Prototype.** This package is an exploratory prototype for composable
> probabilistic infectious disease modelling in Julia. Expect rough edges and
> breaking changes.

`EpiAwarePrototype` builds epidemiological models from small, reusable
components — latent processes, infection processes, and observation models —
and turns each one into a [Turing](https://turinglang.org) /
[DynamicPPL](https://github.com/TuringLang/DynamicPPL.jl) model through a single
generic constructor, `as_turing_model`. Components compose by sampling one
another as submodels, so a full model is *assembled* from parts rather than
hand-written.

## Installation

This package is not registered. Install it directly from the repository:

```julia
using Pkg
Pkg.add(url = "https://github.com/EpiAware/EpiAwarePrototype.jl")
```

## A composable example

```julia
using EpiAwarePrototype, Distributions, Turing

# Generation interval and an EpiData container.
data = EpiData([0.2, 0.3, 0.5], exp)

# Compose a model: an ARIMA-style latent process (differenced AR) feeds a
# direct-infections process, observed with Poisson error.
model = EpiAwareModel(
    DiffLatentModel(; model = AR(), init_priors = [Normal(), Normal()]),
    DirectInfections(; data = data, initialisation_prior = Normal()),
    PoissonError())

# Build a Turing model. `missing` observations simulate from the prior.
n = 30
prior_model = as_turing_model(model, missing, n)

# Sample from the prior, then inspect the generated quantities.
draw = rand(prior_model)
(; generated_y_t, I_t, Z_t) = prior_model()

# Condition on data and run inference.
y = generated_y_t
posterior_model = as_turing_model(model, y, n)
chain = sample(posterior_model, NUTS(), 1_000)
```

Every component is itself an `as_turing_model`, so the same swap-in/swap-out
pattern applies throughout: replace `PoissonError()` with
`NegativeBinomialError()`, wrap it in a `LatentDelay`, or change the latent
process, without touching the rest of the model.

## Adapted from

The modelling code in this package is **ported and adapted** from the
open-source, Apache-2.0 licensed `EpiAware` package
([CDCgov/Rt-without-renewal](https://github.com/CDCgov/Rt-without-renewal),
ported from the fork [seabbs/Rt-without-renewal](https://github.com/seabbs/Rt-without-renewal)).
EpiAwarePrototype is a modified, derived work: it has been renamed, re-architected
around the generic `as_turing_model` constructor, and upgraded to build against
the latest Turing.jl. See the [`NOTICE`](NOTICE) file for full attribution and a
summary of the changes, and [`LICENSE`](LICENSE) for the Apache-2.0 terms.

## License

Apache License 2.0. See [`LICENSE`](LICENSE) and [`NOTICE`](NOTICE).
