Changes are documented in Github releases.

## Unreleased

  - **Breaking**: `GroupedIDModel` is removed. Grouping/panel modelling is now
    part of `IDModel` itself, via `IDModel(infection_model, group_effect,
    observation_model)` and `IDModel(idmodel, group_effect)` (same argument
    lists as the old `GroupedIDModel` constructors). Two new infection-side
    components, `GroupedInfections` and `CombineInfections`, together with an
    extended `Split`, cover the full range of infection↔observation mapping
    cardinalities (one-to-one, one-to-many, many-to-one, many-to-many). See
    issue #180 and the "Partial pooling across groups" case study for the
    migration.
