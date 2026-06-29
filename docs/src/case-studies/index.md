# [Case studies](@id case-studies-overview)

These worked examples build complete models from the package's components,
simulate data from them, and fit them back with [Turing](https://turinglang.org).
Each one is self-contained and runs when the documentation is built, so the
numbers you see are produced by the code on the page.

They progress from a single renewal model to a layered observation process and
then to a mechanistic compartmental model:

  - [Renewal model with negative-binomial reporting](@ref case-study-renewal) —
    a time-varying reproduction number ``R_t`` driven by an autoregressive
    latent process, mapped to infections through the renewal equation and
    observed with overdispersed counts. This is the canonical renewal model of
    [cori2013new](@citet) and [mishra2020derivation](@citet).
  - [Reporting delays and day-of-week effects](@ref case-study-delays) — the
    same renewal core wrapped in an observation model that convolves infections
    through reporting delays and modulates them with a day-of-week reporting
    pattern, in the style of real-time estimation tools [abbott2020estimating](@citep).
  - [Joint cases and deaths from shared infections](@ref case-study-cases-deaths) —
    one renewal infection process observed through two parallel streams at once,
    each with its own delay, ascertainment, and overdispersion, composed with
    [`StackObservationModels`](@ref) and fit jointly to use the signal in both
    [abbott2020estimating](@citep).
  - [An SIR compartmental model](@ref case-study-sir) — an alternative infection
    process where dynamics come from an ordinary differential equation solved by
    the SciML stack [rackauckas2017differentialequations](@citep), following the
    Bayesian compartmental-inference example of [chatzilena2019contemporary](@citet).

Every example uses the same recipe: assemble components into a model, call
[`as_turing_model`](@ref) (directly or through [`EpiAwareModel`](@ref) /
[`EpiProblem`](@ref)), simulate by passing `missing` data, and fit by passing
observed data and sampling. Because the components share one interface, you swap
a modelling assumption by swapping a struct — the
[Composable design](@ref) page explains the mechanism.

!!! note "Prototype"
    These examples are written to be didactic and fast to build. They use short
    sampler runs and simulated data; for real analyses you would use more
    iterations, check convergence diagnostics, and supply observed data.

## References

The methods these case studies recreate and adapt are described in the following
works. Individual pages link back to the relevant entries here.

```@bibliography
```
