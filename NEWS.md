Changes are documented in Github releases.

## Unreleased

  - Rebuilding an observation component no longer corrupts it. `LatentDelay`
    stores its delay PMF reversed, `Ascertainment` stores its prior already
    wrapped in a `PrefixLatentModel`, and `Aggregate` derives its presence
    mask, so feeding the stored fields back to the public constructor reversed,
    prefixed or rejected them a second time. The first two did so without
    throwing, and an `Accessors.@set` of any other field returned a model that
    generated different data. Each now points `ConstructionBase.constructorof`
    at a constructor taking its fields as stored, so `Accessors` operations are
    safe on these types.

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
