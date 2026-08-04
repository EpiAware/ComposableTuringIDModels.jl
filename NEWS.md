Changes are documented in Github releases.

## Unreleased

  - **Breaking**: `TimeVaryingAR` is removed. It was a thin wrapper over
    `AR(; damp = <process>)`; replace `TimeVaryingAR()` with
    `AR(; damp = RandomWalk())`. See issue #182.

    **`TimeVaryingAR` always defaulted `transform` to `tanh`; `AR` picks
    its default from `damp`'s type instead** (`identity` for a bare
    `Distribution` or vector of them, `tanh` for a process such as
    `RandomWalk`). The migration above is exact for the default,
    process-valued `damp`. Anyone who passed a bare `Distribution` as
    `damp` (e.g. `TimeVaryingAR(; damp = Normal(0, 0.05))`) must add
    `transform = tanh` explicitly to `AR` to keep the old behaviour — a
    mechanical rename silently switches them to `identity`, with no error.
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
