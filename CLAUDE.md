# ComposableTuringIDModels.jl

Composable probabilistic infectious disease modelling in Julia.
Goal: reach a state external collaborators can install and test. This is
early-stage software under active development; expect breaking changes.

> **Read this file before starting any work.** It is the standing brief for every
> agent and contributor on this repo.

## Provenance, attribution & licensing (IMPORTANT)

- The modelling code in this package is **ported and adapted** from the
  open-source, **Apache-2.0** licensed `EpiAware` package developed at CDC:
  - Upstream: https://github.com/CDCgov/Rt-without-renewal (`EpiAware/`)
  - Fork ported from: https://github.com/seabbs/Rt-without-renewal
- Because we incorporate Apache-2.0 code, **this package is licensed Apache-2.0**
  for compatibility. We must:
  - Keep the Apache-2.0 `LICENSE`.
  - Ship a `NOTICE` file attributing the upstream authors and stating that this
    is a modified/derived work, with a link back to the upstream repo.
  - Add a short **"Adapted from"** disclaimer in the `README` and docs index.
  - State that significant changes have been made (renamed, re-architected
    around `as_turing_model`, upgraded to the latest Turing).
- Attribute **only** the open Apache-2.0 code above. Do **not** reference any
  non-public repositories or unpublished material in package-facing files, commit
  messages, issues, or docs. All package documentation must be original to this
  package.

## Tooling & template

- This package is scaffolded from **`EpiAwarePackageTools.jl`** (the
  `EpiAwareTestUtils` kit): https://github.com/EpiAware/EpiAwarePackageTools.jl
- Use its `scaffold(pkgdir)` / `update(pkgdir)` and its test helpers
  (`test_aqua`, `test_jet`, `test_doctest`, `test_formatting`, AD harness, etc.)
  **instead of re-inventing** test/CI/dev scaffolding.
- **Adopt the template; adapt our ported code to fit it** — not the other way
  round.
- If the template or `scaffold`/`update` does not work, is missing something, or
  forces an awkward workaround, **file an issue against
  `EpiAware/EpiAwarePackageTools.jl`** describing the gap, rather than silently
  patching around it locally.

## Architecture directive: `as_turing_model`

Replace the upstream abstract **type hierarchy** + per-concept generate
functions (`generate_latent`, `generate_observations`, `generate_latent_infs`,
`generate_epiaware`, dispatching on `AbstractLatentModel` /
`AbstractObservationModel` / `AbstractEpiModel`) with a **single generic
constructor**:

```julia
as_turing_model(model, args...; kwargs...)  # returns a DynamicPPL.Model
```

- Every model struct implements one `@model function as_turing_model(m::MyModel, ...)`.
- Compose via submodels of `as_turing_model(component, ...)`.
- Keep backend-agnostic pieces backend-agnostic (`accumulate_scan`,
  `AbstractAccumulationStep` step structs, distribution/utility helpers).
- Collapse the deep abstract hierarchy. A single light supertype (e.g.
  `AbstractComposableModel`) for shared behaviour/printing is fine; the deep
  `AbstractTuring*` tree is not needed.

## Turing / DynamicPPL: target the latest

- Build against the **latest** released `Turing.jl` / `DynamicPPL.jl`.
- The `@submodel` macro is **removed**. Use the tilde + `to_submodel` form:

  ```julia
  # old:  @submodel ϵ_t = generate_latent(m.ϵ_t, n - p)
  # new:  ϵ_t ~ to_submodel(as_turing_model(m.ϵ_t, n - p), false)
  ```

  **Prefix off is the standard here.** The current DynamicPPL default for
  `to_submodel` is **prefix = true** (it prefixes the submodel's variables with
  the left-hand name), which differs from upstream's flat `@submodel` behaviour.
  Pass `false` as the standard on *every* submodel conversion to preserve the
  existing variable names/behaviour. The only exceptions are the components that
  upstream deliberately prefixed (`PrefixLatentModel`, `PrefixObservationModel`,
  `StackObservationModels`, the old `prefix_submodel` call sites) — implement
  their prefixing explicitly via `to_submodel`'s prefix argument.

## Naming

- Package name and top-level module: **`ComposableTuringIDModels`** (with its own
  fresh UUID in `Project.toml`, distinct from the upstream EpiAware package).
- User-facing model wrappers use the `ID*` prefix: **`IDModel`** (the composed
  infection + observation model) and **`IDProblem`** (the inference wrapper that
  ties latent, infection, and observation models to a dataset).
- The root supertype for shared behaviour/printing is
  **`AbstractComposableModel`** (`IDModel <: AbstractComposableModel`).

## Docs

- Built with **DocumenterVitepress** (the EpiAware org standard).
- Make it clear throughout that this is **composable infectious disease
  modelling** that is early-stage and under active development.
- Document the composable-modelling design (the component DSL idea + the
  Turing.jl backend / `as_turing_model` API) with **original, package-specific
  documentation** written for this package.
- Keep a focused, honest set of pages (getting started, the composable design, a
  worked example / case studies, API reference). Don't add stub/placeholder pages
  for features we are not shipping yet.

## Workflow

- `main` is **branch-protected**: all changes land through **review PRs**. Do not
  commit directly to `main`; do not merge your own PRs (the human gates merges).
