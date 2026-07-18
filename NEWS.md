# News

## Unreleased

### Breaking

- **Priors are now length-`n` submodels and the latent/prior roles are one.**
  Every component parameter (damping, initial conditions, innovation model,
  reproduction number, ascertainment, …) accepts a bare `Distribution`, a
  `Vector{<:Distribution}`, or any process (e.g. a `RandomWalk`) uniformly.
  Components store the raw prior in a parametric field and sample it through the
  public `as_turing_submodel` seam; scalar slots draw natively (`σ ~ m.std`,
  zero-allocation), and `as_turing_model` gains `Distribution` /
  `Vector{<:Distribution}` methods so a bare prior threads through as a length-`n`
  submodel exactly like a model. This widens the `rt` / `Z` / `ϵ_t` slots
  and the manipulator members (`CombineLatentModels` / `ConcatLatentModels` /
  `DiffLatentModel`), so constant reproduction numbers (`Renewal(data; rt =
  Normal())`), constant innovations (`AR(; ϵ_t = Normal())`) and mixed
  `[Distribution, process]` collections are now valid.
- **`AbstractLatentModel` is collapsed into `AbstractPriorModel`.** It remains for
  one release as a deprecated alias (`const AbstractLatentModel =
  AbstractPriorModel`); prefer `AbstractPriorModel`. `implements_latent_interface`
  forwards to `implements_prior_interface`.
- **Chain variable names gain a namespace segment.** Prior variables are
  namespaced at the component's call site (`as_turing_submodel(…; prefix = true)`)
  instead of by a carried name, so e.g. `damp_AR` becomes `damp_AR.θ` and a process
  prior nests deeper. Update any code that reads exact flat variable names.
- **A bare `Distribution` in a length-`n` slot now draws `n` i.i.d. values**
  (`filldist`), not one shared value repeated. The one load-bearing constant case
  (`Ascertainment` with a bare `Distribution`) is preserved by wrapping it in an
  `Intercept` (a single shared draw broadcast across the series).
- **`as_prior`, `BroadcastPrior`, and `sample_prior` are removed** — components
  store raw priors in parametric fields and the `as_turing_submodel` seam replaces
  the coercion — along with the internal `NamedDist`/`_named` naming and the dead
  `_expand_dist` helper.
- **`Renewal` takes one `generation_time` keyword.** The positional
  `Renewal(gen_int; …)` and the keyword `Renewal(; gen_distribution = …)`
  constructors are replaced by a single `Renewal(; generation_time = …, …)` that
  dispatches on the value: a discrete probability vector is used directly, and a
  continuous `Distribution` is discretised internally (double-interval censoring,
  `D_gen`/`Δd`). Update calls to `Renewal(; generation_time = …)`.

### Removed

- **The Pathfinder-based `ManyPathfinder` pre-sampler and the Pathfinder NUTS
  warm-start are removed**, and `Pathfinder` is dropped as a dependency. It was
  the only dependency without a Turing 0.46 release, so it blocked the move to
  DynamicPPL 0.42.1 and its nested-submodel type-inference fix. The now-unused
  `IDMethod` optimisation-then-sampler combinator and the `AbstractIDOptMethod`
  supertype are removed with it; a variational warm-start can be re-added later
  in an extension or a separate package. Resolves
  [#124](https://github.com/EpiAware/ComposableTuringIDModels.jl/issues/124).

### Changed

- **Moved to Turing 0.46 / DynamicPPL 0.42.1.** This tracks the latest Turing
  and picks up DynamicPPL 0.42.1's nested-submodel type-inference fix. For this
  package's composed models the fix is a modest evaluation speedup (~1.2x on a
  representative renewal model, not the order-of-magnitude seen for pathological
  pure-nesting), so it does not by itself make the heavy case-study builds cheap.
- **Docs case studies sample moderate NUTS draws** (250 draws x 2 chains for the
  multi-chain fits; 200-300 for the single-chain ones) rather than research-grade
  counts. This keeps the documentation build to a sensible time while remaining
  statistically valid for the demonstrations; a production analysis would use
  more draws.

### Added

- **`forecast(model, y, chain, horizon)`** — out-of-sample forecasting. Fits at
  length `T`, then predicts the observations over a future horizon
  `t = T+1 … T+h`, carrying each posterior draw forward: the fitted parameters and
  in-sample latent path are held fixed while the non-centred latent process is
  continued over the horizon with fresh prior innovations. Returns a chain of the
  predicted `y_t[T+1] … y_t[T+h]`; also accepts an `IDProblem`. Demonstrated in a
  short section of the renewal case study.
- **General time-varying parameters.** Any component parameter can be constant or
  time-varying through *which prior* fills its slot, with no per-component
  special-casing: a bare `Distribution` draws one scalar RV (constant, no
  length-`n` allocation) while a process draws a length-`n` path, and the component
  reads it per step with the accessor `_at(p, t)` (`_at(::Number)` is the constant,
  `_at(::AbstractVector)` indexes the path). [`AR`](@ref) is the worked example:
  `AR(; damp = Normal())` is a constant coefficient (unchanged) and
  `AR(; damp = RandomWalk())` a per-step time-varying path (`ρ_t = tanh.(damp)`),
  tracked as the generated quantity `ρ`. `TimeVaryingAR` is now a thin alias for
  `AR(; damp = <process>)` rather than a separate type, and `TVARStep` is the
  shared order-1 step for both the constant and time-varying cases.
