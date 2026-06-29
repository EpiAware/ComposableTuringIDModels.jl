# EpiAwarePackageTools template gaps

Template/scaffold gaps found while adopting
[EpiAwarePackageTools.jl](https://github.com/EpiAware/EpiAwarePackageTools.jl) in
this package have all been filed upstream. This file is now just a pointer to
those issues and a note of what we carry locally until each is fixed.

| Upstream issue | What it is | Status here |
|---|---|---|
| [#14](https://github.com/EpiAware/EpiAwarePackageTools.jl/issues/14) | `scaffold` wrote `LICENSE` as a *managed* MIT file, so `update()` reverted a deliberate Apache-2.0 licence | **Fixed upstream & adopted.** `LICENSE` is now package-owned; `update()` no longer touches it. Our Apache-2.0 `LICENSE`/`NOTICE` survive a template sync. |
| [#16](https://github.com/EpiAware/EpiAwarePackageTools.jl/issues/16) | JET runner could not analyse a DynamicPPL `@model` package (false `UndefVarErrorReport`s for every `~`-assigned local) | **Fixed upstream & adopted.** The managed `test/jet/runtests.jl` now reads a package-owned `test/jet/jet_config.jl`; ours sets `JET_REPORT_FILTER = dynamicppl_model_filter` (shipped by the kit). The previous local runner override is gone. |
| [#33](https://github.com/EpiAware/EpiAwarePackageTools.jl/issues/33) | `scaffold` set up plain Documenter, not DocumenterVitepress (the org docs standard) | **Fixed upstream.** This package already migrated its docs to DocumenterVitepress, so no further change was needed here. |
| [#17](https://github.com/EpiAware/EpiAwarePackageTools.jl/issues/17) | `test_explicit_imports` forwards `ignore` only to `check_all_explicit_imports_are_public`, not `check_no_implicit_imports`, so `@reexport` cannot pass | **Open.** Moot for us: the package does not blanket-reexport (users `using EpiAwarePrototype, Distributions, Turing`), matching the upstream EpiAware. The helper limitation still stands for packages that *do* reexport. |
| [#18](https://github.com/EpiAware/EpiAwarePackageTools.jl/issues/18) | `test_doctest` + `@meta CurrentModule` fails under `TestItemRunner` isolation (`Main` has no package binding) | **Open.** Local workaround: docs pages omit `@meta CurrentModule`; `docs/make.jl` uses exported names plus `setdocmeta!` for the full build's cross-references. |

When #17 or #18 is resolved upstream, adopt the fix and update the row above.
