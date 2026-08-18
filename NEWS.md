Changes are documented in Github releases.

## Unreleased

  - Observation data supplied as a `NamedTuple` (`BinomialError`'s
    `(y = successes, N = trials)`, and any error model given a `y` field) no
    longer has its `missing` entries written back into the caller's array. A
    blank is now sampled as a tracked latent `y_t[i]` on every evaluation, as
    it already was for a plain observation vector, so it appears in the chain.
    The same applies to a `ReportTriangle` handed a ready-built
    `ReportingTriangle`. `forecast` still does not work with a
    `BinomialError`-shaped `(y = …, N = …)` series: it reads the two fields as
    two observation streams. That is a separate data-contract defect, tracked
    in #273.

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
