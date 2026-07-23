# PACKAGE-OWNED — scaffold writes this once and never overwrites it.
#
# Package-specific configuration read by the managed `make.jl`.

const LIGHT_TUTORIALS = String[]

const HEAVY_TUTORIALS = String["getting-started/tutorials/ad-backends"]

const TUTORIALS_SUBDIR = joinpath("getting-started", "tutorials")

const TUTORIAL_STUBS = Pair{String, String}[
    "getting-started/tutorials/ad-backends" => "# [Automatic differentiation backends](@id ad-backends)"
]

const FORCE_STUB_TUTORIALS = String[]

const ORG_BRANDING = true

const LINKCHECK_IGNORE = Regex[]

const INDEX_REWRITES = Pair{String, String}[]

const README_EXECUTE = true

const INDEX_STRIP_SECTIONS = String[]

const BENCHMARK_PAGE = false

const HISTORY_SUITES = String[]

const HISTORY_COMMITS = 5

const HISTORY_REGRESSION_THRESHOLD = 1.1
