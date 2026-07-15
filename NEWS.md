# News

## Unreleased

### Breaking

- **Priors are now length-`n` submodels and the latent/prior roles are one.**
  Every component parameter (damping, initial conditions, innovation model,
  reproduction number, ascertainment, …) is now an `AbstractPriorModel`, accepting
  a bare `Distribution`, a `Vector{<:Distribution}`, or any process (e.g. a
  `RandomWalk`) uniformly via `as_prior`. This widens the `rt` / `Z` / `ϵ_t` slots
  and the manipulator members (`CombineLatentModels` / `ConcatLatentModels` /
  `DiffLatentModel`), so constant reproduction numbers (`Renewal(data; rt =
  Normal())`), constant innovations (`AR(; ϵ_t = Normal())`) and mixed
  `[Distribution, process]` collections are now valid.
- **`AbstractLatentModel` is collapsed into `AbstractPriorModel`.** It remains for
  one release as a deprecated alias (`const AbstractLatentModel =
  AbstractPriorModel`); prefer `AbstractPriorModel`. `implements_latent_interface`
  forwards to `implements_prior_interface`.
- **Chain variable names gain a namespace segment.** Prior variables are
  namespaced at the component's call site (prefix-on `to_submodel`) instead of by
  a carried name, so e.g. `damp_AR` becomes `damp_AR.θ` and a process prior nests
  deeper. Update any code that reads exact flat variable names.
- **`as_prior(p, name)` and `BroadcastPrior`'s `name` field are removed**, along
  with the internal `NamedDist`/`_named` naming and the dead `_expand_dist` helper.

### Added

- **`forecast(model, y, chain, horizon)`** — out-of-sample forecasting. Fits at
  length `T`, then predicts the observations over a future horizon
  `t = T+1 … T+h`, carrying each posterior draw forward: the fitted parameters and
  in-sample latent path are held fixed while the non-centred latent process is
  continued over the horizon with fresh prior innovations. Returns a chain of the
  predicted `y_t[T+1] … y_t[T+h]`; also accepts an `IDProblem`. See the
  out-of-sample forecasting case study.
- **`TimeVaryingAR`** — a first-order AR whose damping coefficient is a per-step
  path drawn from a latent process (`ρ_t = tanh.(damp)`), threaded through the new
  `TVARStep`. It returns the numeric path (drops into any latent slot) and tracks
  the coefficient path as the generated quantity `ρ`.
