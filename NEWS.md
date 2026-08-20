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

  - **Breaking**: the `n` passed to `as_turing_model(model, y_t, n)` is now the
    number of OBSERVATIONS, not the length of the infection series.

    An observation chain with delays in it consumes the head of the series it
    is handed: each `LatentDelay` returns a series shorter by
    `length(pmf) - 1`, and the observation-error loop right-aligns the data
    against what is left. Under the old meaning
    `as_turing_model(model, y, length(y))` therefore dropped the leading
    observations from the likelihood without saying so, and a caller had to
    write the arithmetic out to avoid it:

    ```julia
    # before
    n = length(y) + observation_lead_in(model)
    mdl = as_turing_model(model, y, n)

    # after
    mdl = as_turing_model(model, y, length(y))
    ```

    The model now derives the lead-in from the chain and runs the infection
    process over `n + observation_lead_in(model)` time points, so all `n`
    observations are scored. `observation_lead_in` remains, for reading back
    why a model runs longer than its data.

    The old call is **detected, not reinterpreted**: passing
    `length(y) + observation_lead_in(model)` now asks for more observations
    than were supplied, short by exactly the lead-in, and raises an
    `ArgumentError` naming the change. Any other mismatch between `n` and the
    data supplied is an error too.

    An `IDProblem`'s `tspan` is likewise the span of the observations, and
    `forecast` takes no `n`: it extends the observations by the horizon and
    derives the rest. Simulating from the prior with `y_t = missing` now
    returns `n` observations, rather than a longer series with the lead-in at
    its head coming back as `missing`.

  - `data_requirements(model, n)` reports what a caller must supply: one entry
    per observation stream, naming the stream, the length its series must have,
    how many of those entries are scored, the lead-in its chain consumes, and
    the shape the data takes (a plain series, a `(y, N)` pair for
    `BinomialError`, a reference-day × delay matrix for `ReportTriangle`, a
    full series of which an `Aggregate` scores one entry per reporting window).
    Pass the data too and each stream reports what was supplied alongside;
    `data_fits` is the same question as a yes or no. The report is printable,
    indexable by stream name, and iterable.

  - `observation_lead_in` reports how many leading time points an observation
    chain consumes, walking the chain through the `wrapped_models` traversal
    seam. A `Split`'s streams run in parallel, so unequal streams report one
    lead-in each and the series covers the deepest. A delay nested inside an
    `Aggregate` consumes reporting windows rather than time points, which is
    not a series length, so that raises rather than guessing.

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
