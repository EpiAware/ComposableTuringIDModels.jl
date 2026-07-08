# ComposableTuringIDModels.jl <img src="docs/src/assets/logo.svg" width="150" alt="ComposableTuringIDModels logo" align="right">

<!-- badges:start -->
| **Documentation** | **Build Status** | **Code Quality** | **License & DOI** | **Downloads** |
|:-----------------:|:----------------:|:----------------:|:-----------------:|:-------------:|
| [![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://composableturingidmodels.epiaware.org/stable/) [![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://composableturingidmodels.epiaware.org/dev/) | [![Test](https://github.com/EpiAware/ComposableTuringIDModels.jl/actions/workflows/test.yaml/badge.svg?branch=main)](https://github.com/EpiAware/ComposableTuringIDModels.jl/actions/workflows/test.yaml) [![codecov](https://codecov.io/gh/EpiAware/ComposableTuringIDModels.jl/graph/badge.svg)](https://codecov.io/gh/EpiAware/ComposableTuringIDModels.jl) [![AD](https://github.com/EpiAware/ComposableTuringIDModels.jl/actions/workflows/ad.yaml/badge.svg?branch=main)](https://github.com/EpiAware/ComposableTuringIDModels.jl/actions/workflows/ad.yaml) | [![SciML Code Style](https://img.shields.io/static/v1?label=code%20style&message=SciML&color=9558b2&labelColor=389826)](https://github.com/SciML/SciMLStyle) [![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl) [![JET](https://img.shields.io/badge/%E2%9C%88%EF%B8%8F%20tested%20with%20-%20JET.jl%20-%20red)](https://github.com/aviatesk/JET.jl) | [![License: Apache-2.0](https://img.shields.io/badge/License-Apache--2.0-blue.svg)](https://opensource.org/licenses/Apache-2.0) | [![Downloads](https://img.shields.io/badge/dynamic/json?url=http%3A%2F%2Fjuliapkgstats.com%2Fapi%2Fv1%2Ftotal_downloads%2FComposableTuringIDModels&query=total_requests&label=Downloads)](https://juliapkgstats.com/pkg/ComposableTuringIDModels) [![Downloads](https://img.shields.io/badge/dynamic/json?url=http%3A%2F%2Fjuliapkgstats.com%2Fapi%2Fv1%2Fmonthly_downloads%2FComposableTuringIDModels&query=total_requests&suffix=%2Fmonth&label=Downloads)](https://juliapkgstats.com/pkg/ComposableTuringIDModels) |

| ForwardDiff | ReverseDiff (tape) | Enzyme forward | Enzyme reverse | Mooncake reverse | Mooncake forward |
|:---:|:---:|:---:|:---:|:---:|:---:|
| [![cov ForwardDiff](https://codecov.io/gh/EpiAware/ComposableTuringIDModels.jl/graph/badge.svg?flag=ad-forwarddiff)](https://app.codecov.io/gh/EpiAware/ComposableTuringIDModels.jl?flags%5B0%5D=ad-forwarddiff) | [![cov ReverseDiff](https://codecov.io/gh/EpiAware/ComposableTuringIDModels.jl/graph/badge.svg?flag=ad-reversediff)](https://app.codecov.io/gh/EpiAware/ComposableTuringIDModels.jl?flags%5B0%5D=ad-reversediff) | [![cov Enzyme forward](https://codecov.io/gh/EpiAware/ComposableTuringIDModels.jl/graph/badge.svg?flag=ad-enzyme-forward)](https://app.codecov.io/gh/EpiAware/ComposableTuringIDModels.jl?flags%5B0%5D=ad-enzyme-forward) | [![cov Enzyme reverse](https://codecov.io/gh/EpiAware/ComposableTuringIDModels.jl/graph/badge.svg?flag=ad-enzyme-reverse)](https://app.codecov.io/gh/EpiAware/ComposableTuringIDModels.jl?flags%5B0%5D=ad-enzyme-reverse) | [![cov Mooncake reverse](https://codecov.io/gh/EpiAware/ComposableTuringIDModels.jl/graph/badge.svg?flag=ad-mooncake-reverse)](https://app.codecov.io/gh/EpiAware/ComposableTuringIDModels.jl?flags%5B0%5D=ad-mooncake-reverse) | [![cov Mooncake forward](https://codecov.io/gh/EpiAware/ComposableTuringIDModels.jl/graph/badge.svg?flag=ad-mooncake-forward)](https://app.codecov.io/gh/EpiAware/ComposableTuringIDModels.jl?flags%5B0%5D=ad-mooncake-forward) |
<!-- badges:end -->

*A toolkit for composable probabilistic infectious disease modelling in Julia.*

> This package is in early development. Expect rough edges and breaking changes.

## Why ComposableTuringIDModels?

- **Composable models**: Assemble a model from interchangeable infection and
  observation parts — each infection model owning its own latent process —
  instead of writing one monolithic model.
- **Swap a part to change an assumption**: Change the latent process, infection
  process, or observation model on its own to test how each assumption shapes
  your conclusions.
- **One interface**: Every part becomes a [Turing](https://turinglang.org) /
  [DynamicPPL](https://github.com/TuringLang/DynamicPPL.jl) model through the
  single `as_turing_model` constructor, so parts nest freely as submodels. The
  full Turing inference toolbox (NUTS, Pathfinder, prior simulation) applies.
- **Simulate and infer**: Generate synthetic data from any model, then run
  Bayesian inference on the same model with real data.
- **A library of parts**: Random walks, AR/MA/ARIMA latent processes, renewal
  and exponential-growth infection models, ODE (SIR/SEIR) processes, Poisson and
  negative-binomial observations, reporting delays, ascertainment, and
  aggregation — all interchangeable.

## Getting started

You assemble a model from parts, and `ComposableTuringIDModels` turns the assembly into
a single Turing model you can simulate from and fit.
Each part is itself a model, joined through the generic `as_turing_model`
constructor.

```julia
using ComposableTuringIDModels, Distributions, Turing

# Compose a model: an ARIMA-style latent process (a differenced AR) folded
# into a direct-infections process, observed with Poisson error.
model = IDModel(
    DirectInfections(;
        Z = DiffLatentModel(; model = AR(), init_priors = [Normal(), Normal()]),
        initialisation_prior = Normal()),
    PoissonError())

# Build a Turing model; `missing` observations simulate from the prior.
n = 30
prior_model = as_turing_model(model, missing, n)

# Sample from the prior and inspect the generated quantities.
draw = rand(prior_model)
(; generated_y_t, I_t, Z_t) = prior_model()

# Condition on data and run inference.
posterior_model = as_turing_model(model, generated_y_t, n)
chain = sample(posterior_model, NUTS(), 1_000)
```

## Installation

This package is not registered. Install it directly from the repository:

```julia
using Pkg
Pkg.add(url = "https://github.com/EpiAware/ComposableTuringIDModels.jl")
```

## Swap a part

Because every part is interchangeable, the same swap-in/swap-out pattern applies
throughout.
Replace `PoissonError()` with `NegativeBinomialError()`, wrap the observation in
a `LatentDelay`, or change the latent process — without touching the rest of
the model.

```julia
# Same infection process, two observation assumptions.
latent = DiffLatentModel(; model = AR(), init_priors = [Normal(), Normal()])
infections = DirectInfections(; Z = latent, initialisation_prior = Normal())

poisson_model = IDModel(infections, PoissonError())
negbin_model = IDModel(infections, NegativeBinomialError())
```

That is the point: you compare modelling assumptions by swapping parts, not by
rewriting models.

## Where to learn more

- New here? Start with the
  [Overview](https://composableturingidmodels.epiaware.org/dev/overview) and the
  [Composable design](https://composableturingidmodels.epiaware.org/dev/design) page.
- Want worked examples? See the
  [case studies](https://composableturingidmodels.epiaware.org/dev/case-studies), which
  fit complete models to real surveillance data.
- Want the full interface? Browse the
  [Public API](https://composableturingidmodels.epiaware.org/dev/lib/public).
- Want to see the code or report a problem? Check the
  [GitHub repository](https://github.com/EpiAware/ComposableTuringIDModels.jl).

## Adapted from

The modelling code in this package is **ported and adapted** from the
open-source, Apache-2.0 licensed `EpiAware` package
([CDCgov/Rt-without-renewal](https://github.com/CDCgov/Rt-without-renewal),
ported from the fork [seabbs/Rt-without-renewal](https://github.com/seabbs/Rt-without-renewal)).
ComposableTuringIDModels is a modified, derived work: it has been renamed, re-architected
around the generic `as_turing_model` constructor, and upgraded to build against
the latest Turing.jl. See the [`NOTICE`](NOTICE) file for full attribution and a
summary of the changes, and [`LICENSE`](LICENSE) for the Apache-2.0 terms.

## License

Apache License 2.0. See [`LICENSE`](LICENSE) and [`NOTICE`](NOTICE).
