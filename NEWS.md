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

  - `SafePoisson` and `SafeNegativeBinomial` now raise a clear `DomainError`
    naming the offending value when a non-finite (`Inf`/`NaN`) rate or shape
    parameter reaches an integer conversion, instead of an opaque
    `InexactError`. The guard covers every entry point that floors a
    parameter to an integer: `rand` (via `_safe_int_floor`, including the
    Gamma-mixing draw inside `SafeNegativeBinomial.rand`), `quantile` on both
    types, and `SafeNegativeBinomial`'s `mode`.
