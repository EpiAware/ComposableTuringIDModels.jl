Changes are documented in Github releases.

## Unreleased

  - **Breaking**: `GroupedIDModel` is removed. Grouping/panel modelling is now
    part of `IDModel` itself, via `IDModel(infection_model, group_effect,
    observation_model)` and `IDModel(idmodel, group_effect)`. Two new
    infection-side components, `GroupedInfections` and `CombineInfections`,
    together with an extended `Split`, cover the full range of
    infection↔observation mapping cardinalities (one-to-one, one-to-many,
    many-to-one, many-to-many). See issue #180 and the "Partial pooling across
    groups" case study for the migration.

    **The observed data matrix is transposed relative to `GroupedIDModel`.**
    It is now `n_groups × n_time`, where it was `n_time × n_groups`. The
    constructor argument lists are otherwise unchanged, so a mechanical rename
    compiles and runs: a `24 × 8` matrix builds a 24-group, 8-step panel
    without error. Transpose your data when migrating.
