Changes are documented in Github releases.

## Unreleased

  - `Aggregate` no longer fails to build when a `LatentDelay` is nested inside
    it (#292). The delay shortens the window series, and the inner model's
    predictions are now scattered onto the reporting windows that survive it,
    dropping the leading ones as every other observation chain does. An
    expected series shorter than the observed series is likewise right-aligned,
    so a delay applied outside the aggregation composes too. `Aggregate`'s
    docstring now states that the nesting decides whether the delay is measured
    in windows or in time points.

  - Infection models can now generate several strata at once. `Stratify` puts
    a stratum axis on a shared process (a partially pooled panel of
    reproduction numbers, say), and `Renewal`'s `mixing` slot couples the
    resulting patches, from a fixed weight matrix to an inferred `Gravity`
    model. `Split` projects any number of infection strata onto any number of
    observation streams through a weight matrix, so `IDModel` builds
    many-to-one and many-to-many infection↔observation mappings from a plain
    two-argument constructor with no separate panel API.

    **Breaking**: `GroupedInfections` and the grouped `IDModel` constructors
    are removed; a `Stratify`-based infection model observed through `Split`
    replaces them.
