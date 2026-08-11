# Documentation navigation. `lib/public.md` and `lib/internals.md` are generated
# by make.jl from the module's documented bindings.
pages = [
    "Home" => "index.md",
    # Installation is an authored page rather than a README section, so every
    # ```julia block on the README-derived home page can run.
    "Installation" => "getting-started/installation.md",
    "Overview" => "overview.md",
    "Composable design" => "design.md",
    "Tutorials" => [
        "Overview" => "tutorials/index.md",
        "Renewal model with negative-binomial reporting" => "tutorials/renewal-negbin.md",
        "Reporting delays and day-of-week effects" => "tutorials/delays-dayofweek.md",
        "Real-time nowcasting: correcting right-truncation" => "tutorials/realtime-nowcast.md",
        "An SIR compartmental model" => "tutorials/sir-ode.md",
        "Declarative compartmental models with Catalyst" => "tutorials/catalyst-ode.md",
        "Multiple observation streams: cases, deaths, and strata" => "tutorials/split-observations.md",
        "A Gaussian-process latent process" => "tutorials/gaussian-process.md",
        "Time-varying damping in an AR process" => "tutorials/time-varying-damping.md",
        "Partial pooling across groups" => "tutorials/hierarchy-stacked.md",
        "Coupled patch models" => "tutorials/patch-models.md",
        "Renewal modifiers: depletion and importation" => "getting-started/tutorials/renewal-modifiers.md",
    ],
    "API reference" => [
        "Public API" => "lib/public.md",
        "Internal API" => "lib/internals.md",
    ],
    "Benchmarks" => [
        "Overview" => "benchmarks.md",
        "Automatic differentiation" => "getting-started/tutorials/ad-backends.md",
    ],
]
