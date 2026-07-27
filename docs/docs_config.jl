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

# Run the README's ```julia fences on the home page so every block shows its
# printed output (#218). The managed build converts EVERY fence into an
# executed `@example readme` block and has no per-fence opt-out, so the one
# fence that must not run, `Pkg.add("ComposableTuringIDModels")` under
# Installation, is handled by stripping that section from the generated index
# (see `INDEX_STRIP_SECTIONS` below). Executing it would hit the registry and
# install a second copy of the package into the docs environment. Installation
# stays in `README.md` for GitHub readers, and the docs site carries it on the
# Overview page.
const README_EXECUTE = true

const INDEX_STRIP_SECTIONS = String["Installation"]

# The package benchmarks and publishes a timeline to the `benchmarks` branch,
# so the generated page (heading + `docs/benchmarks.md` prose + history) is on.
const BENCHMARK_PAGE = true

const HISTORY_SUITES = String[]

const HISTORY_COMMITS = 5

const HISTORY_REGRESSION_THRESHOLD = 1.1
