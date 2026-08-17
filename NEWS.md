Changes are documented in Github releases.

## Unreleased

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

  - Fixed: discretising a continuous delay distribution without an explicit
    `D` now takes the horizon from the distribution's own support when it is
    finite (e.g. a distribution the caller has already `truncated`), instead
    of silently re-deriving a shorter horizon from a quantile. `LatentDelay`,
    `UncertainDelay`, `ReportingCDF`, `ReportingPMF` and `Renewal`'s
    generation-interval slot all pick this up automatically; a caller can now
    express "use this distribution's own bound" by truncating the
    distribution, without a separate `D` keyword.
