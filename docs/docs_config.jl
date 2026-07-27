# PACKAGE-OWNED — scaffold writes this once and never overwrites it.
#
# Package-specific configuration read by the managed `make.jl`.

# `renewal-modifiers.jl` is light: Literate emits `@example` blocks that
# Documenter runs in-process, exactly as the case studies do, so the page also
# renders under `--skip-notebooks`. Its fit is deliberately small (n = 60, 250
# draws) to keep that cheap.
const LIGHT_TUTORIALS = String["renewal-modifiers.jl"]

# Names are relative to `TUTORIALS_SUBDIR`: `HEAVY_TUTORIALS` holds the Literate
# `.jl` source names, `TUTORIAL_STUBS` is keyed by the rendered `.md` names.
const HEAVY_TUTORIALS = String["ad-backends.jl"]

const TUTORIALS_SUBDIR = joinpath("getting-started", "tutorials")

const TUTORIAL_STUBS = Pair{String, String}[
    "ad-backends.md" => "# [Automatic differentiation backends](@id ad-backends)"
]

const FORCE_STUB_TUTORIALS = String[]

const ORG_BRANDING = true

# `CITATION.cff` is package-owned and seeded only by the kit's `scaffold`, never
# by `update`, while the managed "How to cite" README section links to it
# unconditionally. This change adds the file, but linkcheck resolves the link
# against `main`, so it stays a 404 until this branch lands. Drop this entry
# once the file is on `main`.
const LINKCHECK_IGNORE = Regex[
    r"blob/main/CITATION\.cff$"
]

const _REPO_BLOB = "https://github.com/EpiAware/ComposableTuringIDModels.jl/blob/main"

# `NOTICE` and `LICENSE` are repo-root files. The README's relative links are
# right on GitHub, but the generated `index.md` is served from the docs site,
# where they resolve to nothing — so point the docs copy at the files on `main`.
const INDEX_REWRITES = Pair{String, String}[
    "(NOTICE)" => "($(_REPO_BLOB)/NOTICE)",
    "(LICENSE)" => "($(_REPO_BLOB)/LICENSE)"
]

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
