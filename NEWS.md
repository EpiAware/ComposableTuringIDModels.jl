Changes are documented in Github releases.

## Unreleased

  - An observation chain can now be inspected and rewritten through a supported
    API rather than by field path. `wrapped_models` reports what a single
    component wraps, `observation_components` walks a whole chain outermost
    first with a `Split`'s branches followed in stream order, and `rewrap`
    rebuilds a component around new wrapped models. All three are public but
    not exported.

  - Rebuilding an observation component no longer corrupts it. `LatentDelay`
    stores its delay PMF reversed, `Ascertainment` stores its prior already
    wrapped in a `PrefixLatentModel`, and `Aggregate` derives its presence
    mask, so feeding the stored fields back to the public constructor reversed,
    prefixed or rejected them a second time. The first two did so without
    throwing, and an `Accessors.@set` of any other field returned a model that
    generated different data. Each now points `ConstructionBase.constructorof`
    at a constructor taking its fields as stored, so `Accessors` operations are
    safe on these types.

  - `observation_lead_in` reports how many leading time points an observation
    chain drops, so a caller can write
    `n = length(y) + observation_lead_in(model)` instead of summing
    `length(pmf) - 1` over the chain's delays by hand.
    `observation_coverage` reports how many of the supplied observations a
    given `n` actually scores.

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
