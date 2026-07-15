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

- **Moved to Turing 0.46 / DynamicPPL 0.42.1.** DynamicPPL 0.42.1's
  nested-submodel type-inference fix makes the package's prefix-on composition
  fast to evaluate and differentiate again, so the docs case studies run their
  full NUTS draw counts.

### Added

- **`TimeVaryingAR`** — a first-order AR whose damping coefficient is a per-step
  path drawn from a latent process (`ρ_t = tanh.(damp)`), threaded through the new
  `TVARStep`. It returns the numeric path (drops into any latent slot) and tracks
  the coefficient path as the generated quantity `ρ`.
