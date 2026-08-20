<!-- epiaware-standards:start MANAGED by EpiAwarePackageTools.scaffold -->

# ComposableTuringIDModels

Standards and how-to live in the docs, not here.

- [Ecosystem design](https://epiaware.org/approaches/) — why composition, what is wanted from it, and the approaches being explored.
- [Package standards](https://epiawarepackagetools.epiaware.org/stable/standards)
- [Test infrastructure](https://epiawarepackagetools.epiaware.org/stable/getting-started/test-infrastructure)
- [Infrastructure and template sync](https://epiawarepackagetools.epiaware.org/stable/getting-started/infrastructure)
- [ComposableTuringIDModels documentation](https://composableturingidmodels.epiaware.org)
- [EpiAware](https://epiaware.org)
- [ColPrac](https://github.com/SciML/ColPrac)
<!-- epiaware-standards:end -->

## Docs

Documentation pages do not use `@assert`.
A figure a page claims is computed and rendered, so a regression is visible to a reader rather than breaking the build.
Checks that must fail belong in the test suite, because gating the docs build on sampler output makes it a version-coupled tripwire.
