# PACKAGE-OWNED — scaffold writes this once and never overwrites it.
#
# Package-specific configuration read by the managed `make.jl`.

# `renewal-modifiers.jl` is light: Literate emits `@example` blocks that
# Documenter runs in-process, exactly as the tutorials do, so the page also
# renders under `--skip-notebooks`. Its fit is deliberately small (n = 60, 250
# draws) to keep that cheap.
const LIGHT_TUTORIALS = String["renewal-modifiers.jl"]

# Names are relative to `TUTORIALS_SUBDIR`: `HEAVY_TUTORIALS` holds the Literate
# `.jl` source names, `TUTORIAL_STUBS` is keyed by the rendered `.md` names.
const HEAVY_TUTORIALS = String[]

const TUTORIALS_SUBDIR = joinpath("getting-started", "tutorials")

const TUTORIAL_STUBS = Pair{String, String}[]

## `ad-comparison.jl` benchmarks every (backend, scenario) pair in-process,
## serially. This PR's own backend/scenario expansion (4 -> 7 backends,
## 32 -> 43 scenarios) pushed that past 2h without finishing (#232). Stubbed
## here as an interim measure pending an artefact-based redesign that reads
## pre-computed per-backend CI benchmark output instead of running it live
## during the docs build -- tracked upstream in EpiAwarePackageTools.jl.
## Drop this entry once that lands and the page is cheap again.
const FORCE_STUB_TUTORIALS = String["ad-comparison.jl"]

const HEAVY_BENCHMARKS = String["ad-comparison.jl"]

const BENCHMARK_STUBS = Pair{String, String}[
    "ad-comparison.md" => "# [AD backend comparison](@id ad-comparison)\n\n## [Choosing a backend](@id ad-backends)",
]

const ORG_BRANDING = true

# `CITATION.cff` is package-owned and seeded only by the kit's `scaffold`, never
# by `update`, while the managed "How to cite" README section links to it
# unconditionally. This change adds the file, but linkcheck resolves the link
# against `main`, so it stays a 404 until this branch lands. Drop this entry
# once the file is on `main`.
#
# `dev/tutorials` is the same self-reference problem: the README now points at
# the renamed Tutorials section, but linkcheck resolves it against the *live*
# deployed site, which still serves the old `case-studies` paths until this
# branch merges and redeploys. Drop this entry once the file is on `main`.
const LINKCHECK_IGNORE = Regex[
    r"blob/main/CITATION\.cff$",
    r"dev/tutorials$",
]

const _REPO_URL = "https://github.com/EpiAware/ComposableTuringIDModels.jl"

const _REPO_BLOB = "$(_REPO_URL)/blob/main"

# `NOTICE` and `LICENSE` are repo-root files. The README's relative links are
# right on GitHub, but the generated `index.md` is served from the docs site,
# where they resolve to nothing — so point the docs copy at the files on `main`.
#
# Nothing else is rewritten. The home page runs the README's code as written,
# unseeded, so a build shows a genuine draw from the prior rather than a curated
# one: the latent process is a random walk on the log scale, so the simulated
# series and the fit that follows differ from build to build.
const INDEX_REWRITES = Pair{String, String}[
    "(NOTICE)" => "($(_REPO_BLOB)/NOTICE)",
    "(LICENSE)" => "($(_REPO_BLOB)/LICENSE)",
]

# Run the README's ```julia fences on the home page so every block shows its
# printed output (#218). The managed build converts EVERY fence into an
# executed `@example readme` block and has no per-fence opt-out, so every fence
# in `README.md` has to be runnable. It is: the README carries no installation
# snippet, which is the one thing that must not run during a docs build (it
# would hit the registry and install a second copy of the package into the docs
# environment). Installation lives on the authored Installation page instead,
# where a plain ```julia fence renders without executing — the same arrangement
# as the other EpiAware packages.
#
# COST: the last block fits `NUTS()` for 1000 draws, which dominates the build
# time of the home page (minutes, not seconds). The draw count lives in
# `README.md` and is deliberately not reduced here.
const README_EXECUTE = true

# Nothing to omit from the home page: with installation moved to the
# Installation page, every README section is safe to render and run.
const INDEX_STRIP_SECTIONS = String[]

# The package benchmarks and publishes a timeline to the `benchmarks` branch,
# so the generated page (heading + `docs/benchmarks.md` prose + history) is on.
const BENCHMARK_PAGE = true

const HISTORY_SUITES = String[]

const HISTORY_COMMITS = 5

const HISTORY_REGRESSION_THRESHOLD = 1.1
