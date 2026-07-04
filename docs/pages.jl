# Documentation navigation. `lib/public.md` and `lib/internals.md` are generated
# by make.jl from the module's documented bindings.
pages = [
    "Home" => "index.md",
    "Composable design" => "design.md",
    "Case studies" => [
        "Overview" => "case-studies/index.md",
        "Renewal model with negative-binomial reporting" => "case-studies/renewal-negbin.md",
        "Reporting delays and day-of-week effects" => "case-studies/delays-dayofweek.md",
        "Real-time nowcasting: correcting right-truncation" => "case-studies/realtime-nowcast.md",
        "An SIR compartmental model" => "case-studies/sir-ode.md",
        "A declarative SEIR model with Catalyst" => "case-studies/catalyst-seir.md"
    ],
    "API reference" => [
        "Public API" => "lib/public.md",
        "Internal API" => "lib/internals.md"
    ]
]
