# [Case studies](@id case-studies-overview)

These worked examples build complete models from the package's components and fit
them to real epidemic surveillance data with [Turing](https://turinglang.org),
recreating published analyses. Each one is self-contained and runs when the
documentation is built, so the numbers you see are produced by the code on the
page.

!!! note "Interim sampling settings"
    The case studies currently draw fewer NUTS samples than we would use for a
    real analysis.
    This is a documented interim measure, not a modelling choice.
    The prefix-on model composition this package relies on is slow to evaluate on
    DynamicPPL 0.41, and the fix is the nested-submodel type-inference improvement
    in DynamicPPL 0.42.1 (Turing 0.46).
    We cannot yet move to those versions because Pathfinder, a dependency, has no
    Turing 0.46 release.
    The draws will be restored once Pathfinder supports Turing 0.46 and the
    already-widened compat bounds resolve the update automatically.
    See [issue #124](https://github.com/EpiAware/ComposableTuringIDModels.jl/issues/124)
    and the upstream fix [mlcolab/Pathfinder.jl#400](https://github.com/mlcolab/Pathfinder.jl/pull/400).

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
  - [Real-time nowcasting: correcting right-truncation](@ref case-study-nowcast) —
    the same renewal core fit to a right-truncated real-time snapshot, contrasting
    a naive fit (which shows the artefactual recent-``R_t`` down-turn) with a
    [`RightTruncate`](@ref)-corrected fit that removes it, again following
    real-time estimation practice [abbott2020estimating](@citep).
  - [Multiple observation streams: cases, deaths, and strata](@ref case-study-split) —
    one renewal infection process observed through several named streams with a
    single [`Split`](@ref) construct, covering parallel streams (cases and deaths
    off shared infections), a cascade (deaths downstream of reported cases,
    achieved by placing the split lower in the pipeline), and data-driven strata
    (one stream per age band), motivated by the differing biases of surveillance
    streams [sherratt2021surveillance](@citep).
  - [An SIR compartmental model](@ref case-study-sir) — an alternative infection
    process where dynamics come from an ordinary differential equation solved by
    the SciML stack [rackauckas2017differentialequations](@citep), following the
    Bayesian compartmental-inference example of [chatzilena2019contemporary](@citet).

Every example uses the same recipe: assemble components into a model, call
[`as_turing_model`](@ref) (directly or through [`IDModel`](@ref) /
[`IDProblem`](@ref)), simulate by passing `missing` data, and fit by passing
observed data and sampling. Because the components share one interface, you swap
a modelling assumption by swapping a struct — the
[Composable design](@ref) page explains the mechanism.

## References

The methods these case studies recreate and adapt are described in the following
works. Individual pages link back to the relevant entries here.

```@bibliography
```
