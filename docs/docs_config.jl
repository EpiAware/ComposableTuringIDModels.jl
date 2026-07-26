# PACKAGE-OWNED — scaffold writes this once and never overwrites it.
#
# Package-specific configuration read by the managed `make.jl`.

const LIGHT_TUTORIALS = String[]

# Names are relative to `TUTORIALS_SUBDIR`: `HEAVY_TUTORIALS` holds the Literate
# `.jl` source names, `TUTORIAL_STUBS` is keyed by the rendered `.md` names.
const HEAVY_TUTORIALS = String["ad-backends.jl"]

const TUTORIALS_SUBDIR = joinpath("getting-started", "tutorials")

const TUTORIAL_STUBS = Pair{String, String}[
    "ad-backends.md" => "# [Automatic differentiation backends](@id ad-backends)"
]

const FORCE_STUB_TUTORIALS = String[]

const ORG_BRANDING = true

const LINKCHECK_IGNORE = Regex[]

const INDEX_REWRITES = Pair{String, String}[]

# The README's ```julia fences are illustrative, and the managed build converts
# EVERY one of them to an executed `@example readme` block. Executing them would
# run `Pkg.add("ComposableTuringIDModels")` from the Installation section against
# the registry, add a 1000-draw NUTS fit to every docs build, and call
# `summarystats`, which the README never brings into scope.
const README_EXECUTE = false

const INDEX_STRIP_SECTIONS = String[]

# The package benchmarks and publishes a timeline to the `benchmarks` branch,
# so the generated page (heading + `docs/benchmarks.md` prose + history) is on.
const BENCHMARK_PAGE = true

const HISTORY_SUITES = String[]

const HISTORY_COMMITS = 5

const HISTORY_REGRESSION_THRESHOLD = 1.1
