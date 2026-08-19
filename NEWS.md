Changes are documented in Github releases.

## Unreleased

  - The prior-slot widening helpers a custom component's constructor and
    recursion call directly are now `public` (documented, not exported):
    `at`, `path_prior`, `prior_order`, `assert_prior_length`. The previous
    private names (`_at`, `_path_prior`, `_prior_order`,
    `_assert_prior_length`) still work as aliases of the same functions, so
    this is non-breaking.

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
