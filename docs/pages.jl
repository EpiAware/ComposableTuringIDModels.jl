# Documentation navigation. `lib/public.md` and `lib/internals.md` are generated
# by make.jl from the module's documented bindings.
pages = [
    "Home" => "index.md",
    "Composable design" => "design.md",
    "Case studies" => [
        "Overview" => "case-studies/index.md",
        "Renewal model with negative-binomial reporting" => "case-studies/renewal-negbin.md",
        "Reporting delays and day-of-week effects" => "case-studies/delays-dayofweek.md",
        "Joint cases and deaths from shared infections" => "case-studies/cases-deaths.md",
        "An SIR compartmental model" => "case-studies/sir-ode.md"
    ],
    "API reference" => [
        "Public API" => "lib/public.md",
        "Internal API" => "lib/internals.md"
    ]
]
