# Documentation navigation. `lib/public.md` and `lib/internals.md` are generated
# by make.jl from the module's documented bindings.
pages = [
    "Home" => "index.md",
    # Authored quickstart pages, distinct from the README-derived home page:
    # installation lives here, so the home page's ```julia blocks can all run.
    "Getting started" => [
        "Overview" => "getting-started/index.md",
        "Tutorials" => [
            "Automatic differentiation backends" => "getting-started/tutorials/ad-backends.md"
        ]
    ],
    "Overview" => "overview.md",
    "Composable design" => "design.md",
    "Case studies" => [
        "Overview" => "case-studies/index.md",
        "Renewal model with negative-binomial reporting" => "case-studies/renewal-negbin.md",
        "Reporting delays and day-of-week effects" => "case-studies/delays-dayofweek.md",
        "Real-time nowcasting: correcting right-truncation" => "case-studies/realtime-nowcast.md",
        "An SIR compartmental model" => "case-studies/sir-ode.md",
        "Declarative compartmental models with Catalyst" => "case-studies/catalyst-ode.md",
        "Multiple observation streams: cases, deaths, and strata" => "case-studies/split-observations.md",
        "Time-varying damping in an AR process" => "case-studies/time-varying-damping.md",
        "Partial pooling across groups" => "case-studies/hierarchy-stacked.md"
    ],
    "API reference" => [
        "Public API" => "lib/public.md",
        "Internal API" => "lib/internals.md"
    ],
    "Benchmarks" => "benchmarks.md"
]
