# [Tutorials](@id tutorials-overview)

These worked examples build complete models from the package's components and fit them with [Turing](https://turinglang.org).
Most use real epidemic surveillance data and recreate published analyses.
The Gaussian-process study instead fits data simulated from the model itself, where the point is to recover a known truth.
Each one is self-contained and runs when the documentation is built, so the numbers you see are produced by the code on the page.

!!! note "Sampling settings"
    To keep the documentation build to a sensible time, these tutorials draw moderate NUTS samples, two chains of 200-300 draws each, run in parallel with `MCMCThreads()`.
    The target acceptance rate (`adapt_delta`) has been raised across these tutorials to avoid divergent transitions.
    That is enough to demonstrate the models and produce stable figures.
    A real analysis would use more draws and check convergence diagnostics carefully.

To fit a renewal model to reported cases, start with [Renewal model with negative-binomial reporting](@ref tutorial-renewal).
Follow it with [Reporting delays and day-of-week effects](@ref tutorial-delays) for a realistic observation process.

The table groups the pages by theme.
The infection, observation, and latent themes are the package's three component roles.
Those names are also used on the [Composable design](@ref) page and in the [Public API](@ref public-api).
A page covering more than one theme appears under each.

| Theme | Tutorials |
|:--|:--|
| End-to-end worked examples | [Renewal with negative-binomial reporting](@ref tutorial-renewal), [Reporting delays and day-of-week effects](@ref tutorial-delays), [An SIR compartmental model](@ref tutorial-sir) |
| Infection processes | [Renewal with negative-binomial reporting](@ref tutorial-renewal), [Renewal modifiers](@ref renewal-modifiers), [An SIR compartmental model](@ref tutorial-sir), [Declarative compartmental models](@ref tutorial-catalyst), [Coupled patch models](@ref tutorial-patches) |
| Observation models | [Reporting delays and day-of-week effects](@ref tutorial-delays), [Real-time nowcasting](@ref tutorial-nowcast), [Multiple observation streams](@ref tutorial-split) |
| Latent processes | [A Gaussian-process latent process](@ref tutorial-gp), [Time-varying damping in an AR process](@ref tutorial-tvdamp), [Renewal with negative-binomial reporting](@ref tutorial-renewal) |
| Multiple strata and groups | [Multiple observation streams](@ref tutorial-split), [Partial pooling across groups](@ref tutorial-hierarchy), [Coupled patch models](@ref tutorial-patches) |
| Mechanistic and ODE models | [An SIR compartmental model](@ref tutorial-sir), [Declarative compartmental models](@ref tutorial-catalyst) |
| Nowcasting and forecasting | [Real-time nowcasting](@ref tutorial-nowcast), [Renewal with negative-binomial reporting](@ref tutorial-renewal) |

## The tutorials

  - [Renewal model with negative-binomial reporting](@ref tutorial-renewal).
    A time-varying reproduction number ``R_t`` driven by an autoregressive latent process, mapped to infections through the renewal equation and observed with overdispersed counts.
    This is the canonical renewal model of [cori2013new](@citet) and [mishra2020derivation](@citet).
  - [Reporting delays and day-of-week effects](@ref tutorial-delays).
    The same renewal core wrapped in an observation model that convolves infections through reporting delays and modulates them with a day-of-week reporting pattern, in the style of real-time estimation tools [abbott2020estimating](@citep).
  - [Real-time nowcasting: correcting right-truncation](@ref tutorial-nowcast).
    The same renewal core fit to a right-truncated real-time snapshot, contrasting a naive fit (which shows the artefactual recent-``R_t`` down-turn) with a [`RightTruncate`](@ref)-corrected fit that removes it, again following real-time estimation practice [abbott2020estimating](@citep).
  - [Multiple observation streams: cases, deaths, and strata](@ref tutorial-split).
    One renewal infection process observed through several named streams with a single [`Split`](@ref) construct, covering parallel streams (cases and deaths off shared infections), a cascade (deaths downstream of reported cases, achieved by placing the split lower in the pipeline), and data-driven strata (one stream per age band), motivated by the differing biases of surveillance streams [sherratt2021surveillance](@citep).
  - [An SIR compartmental model](@ref tutorial-sir).
    An alternative infection process where dynamics come from an ordinary differential equation solved by the SciML stack [rackauckas2017differentialequations](@citep), following the Bayesian compartmental-inference example of [chatzilena2019contemporary](@citet).
  - [Declarative compartmental models with Catalyst](@ref tutorial-catalyst).
    The same SIR dynamics declared as a [Catalyst.jl](https://github.com/SciML/Catalyst.jl) reaction network [loman2023catalyst](@citep) instead of a hand-written vector field, so the ODE system and its Jacobian are generated rather than maintained by hand.
  - [A Gaussian-process latent process](@ref tutorial-gp).
    Plugs a Hilbert-space approximate Gaussian process [riutortmayol2023practical](@citep) into the renewal model as the latent ``\log R_t`` process, fits it under NUTS with Mooncake, and checks it against an exact GP and against the simulated latent.
  - [Time-varying damping in an AR process](@ref tutorial-tvdamp).
    The same autoregressive latent process with its damping coefficient itself allowed to vary over time, so the amount of mean reversion in ``\log R_t`` is learned rather than fixed.
  - [Partial pooling across groups](@ref tutorial-hierarchy).
    A panel of groups sharing one infection process structure but with group-level parameters drawn from a common population distribution, so information is pooled across groups without forcing them to be identical.
  - [Coupled patch models](@ref tutorial-patches).
    Several patches over one time axis, built a layer at a time, from independent patches to a shared partially pooled process to a mixing matrix that moves infection pressure between them.
    The finished model is fit to data simulated from itself, and the coupling is then swapped for a [`Gravity`](@ref) model and for one written as a single method.
  - [Renewal modifiers: depletion and importation](@ref renewal-modifiers).
    One delayed renewal process compared against the same process extended with [`SusceptibleDepletion`](@ref) and [`ImportedCases`](@ref), so each modifier's contribution is visible against a shared baseline.

Every example uses the same recipe.
Assemble components into a model, call [`as_turing_model`](@ref) on it, simulate by passing `missing` data, and fit by passing observed data and sampling.
Because the components share one interface, you swap a modelling assumption by swapping a struct.
The [Composable design](@ref) page explains the mechanism.

## References

The methods these tutorials recreate and adapt are described in the following works.
Individual pages link back to the relevant entries here.

```@bibliography
```
