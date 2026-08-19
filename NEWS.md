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
    given `n` actually scores, counting them through the contract of the model
    that scores them: the `y` field of a `BinomialError`'s `(y, N)` data, the
    reference days of a `ReportTriangle`'s reporting triangle, and the
    reporting windows of an `Aggregate` rather than the time points it is
    given.

  - `forecast` takes the series length the fit used as `n`, so a model fitted
    with `n = length(y) + observation_lead_in(model)` forecasts from the same
    length rather than being rebuilt at `length(y) + horizon`. A chain fitted
    over a longer series than the horizon model covers is now an error rather
    than a silent misalignment.

    **Breaking**: `forecast(::IDProblem, ...)` now takes that length from the
    problem's `tspan` rather than from the length of `y`. The two differ
    exactly when the span was set to cover an observation chain's lead-in, and
    for such a problem the forecast is now made at a different length than
    before. Pass `n` explicitly to keep the old behaviour.

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
