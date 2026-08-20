Changes are documented in Github releases.

## Unreleased

  - An observation model now right-aligns its data against the expected series
    in both directions. An expected series longer than the data is read as
    unobserved run-in and its head is left unscored, where it used to be
    rejected. This is what a `Split` produces whenever its streams consume
    different lead-ins, so a caller passes each stream the observations it has
    rather than padding the shorter-lead-in stream with leading `missing`s.

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
