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

### Changed

- **Turing/DynamicPPL compat widened to allow 0.46/0.42**, so the update resolves
  automatically once the upstream block clears. The move is currently held back by
  Pathfinder, which has no Turing 0.46 release; DynamicPPL 0.42.1 carries the
  nested-submodel type-inference fix that makes prefix-on composition fast again.
  Until then the docs case studies draw fewer NUTS samples than a real analysis
  would as a documented interim, restored when Pathfinder supports Turing 0.46.
  See [#124](https://github.com/EpiAware/ComposableTuringIDModels.jl/issues/124).

### Added

- **`TimeVaryingAR`** — a first-order AR whose damping coefficient is a per-step
  path drawn from a latent process (`ρ_t = tanh.(damp)`), threaded through the new
  `TVARStep`. It returns the numeric path (drops into any latent slot) and tracks
  the coefficient path as the generated quantity `ρ`.
