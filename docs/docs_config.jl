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

const _REPO_URL = "https://github.com/EpiAware/ComposableTuringIDModels.jl"

const _REPO_BLOB = "$(_REPO_URL)/blob/main"

# The exact README lines the rewrites below key off, kept as constants so both
# sides of each `from => to` pair stay in step with `README.md`.
const _USING_LINE = "using ComposableTuringIDModels, Distributions, Turing"

const _PRIOR_LINE = "prior_model = as_turing_model(model, missing, n)"

const _EXAMPLES_HEADING = "## Getting started"

const _SEED_LINES = "\nusing Random # hide\nRandom.seed!(27) # hide"

# Seeding and hidden output are both invisible in the rendered code, so say so
# on the page rather than letting a reader assume a fresh run reproduces it.
const _OUTPUT_NOTE = """
!!! note "About the output on this page"
    The blocks below run when this page is built, at a fixed seed set behind
    the scenes.
    This model's latent process is a random walk on the log scale, so its prior
    over case counts is heavy-tailed: running the same code yourself gives a
    different simulated series and a different fit.
    The block that builds the Turing model prints nothing, because its value
    displays as one very long type signature."""

# `README.md` is written for its GitHub readers; every adjustment the docs home
# page needs is made here, on the generated copy. `build_index` applies these
# rewrites line by line, in prose as well as inside ```julia fences, so a
# replacement containing `\n` inserts lines.
#
#  1. `NOTICE`/`LICENSE` are repo-root files: the README's relative links are
#     right on GitHub but resolve to nothing on the docs site.
#  2. The home page's blocks are executed (see `README_EXECUTE`) and start by
#     drawing from the prior, so unseeded the page differs on every build. Over
#     seeds 1:60 the simulated series ran from all-zero to a peak of 1e36
#     cases, and the fit that follows was usually unconverged. Seed 27 gives a
#     curve rising from 6 to a peak of 74 cases and decaying back to 9 over 30
#     days, and a fit that mixes (rhat 1.01, min ess_bulk 327). Re-check both
#     before changing it: the fit degrades sharply as the counts grow, because
#     a Poisson likelihood on large counts pins the log-scale latent series
#     into a very ill-conditioned posterior. Seed 15 (peak 5210 cases) leaves
#     NUTS at rhat 2.02 and min ess_bulk 2.6 after 1000 draws.
#  3. `as_turing_model` returns a `DynamicPPL.Model`, which displays as one
#     very long type signature — noise on the home page, so a `nothing # hide`
#     line drops it. `@example` shows the value of the block's last expression,
#     and unlike the REPL a trailing `;` does not suppress that.
#
# The inserted code lines all end in `# hide`, so they run but stay out of the
# rendered snippet: the code shown on the home page is character-for-character
# the code in `README.md`.
const INDEX_REWRITES = Pair{String, String}[
    "(NOTICE)" => "($(_REPO_BLOB)/NOTICE)",
    "(LICENSE)" => "($(_REPO_BLOB)/LICENSE)",
    _EXAMPLES_HEADING => "$(_EXAMPLES_HEADING)\n\n$(_OUTPUT_NOTE)",
    _USING_LINE => _USING_LINE * _SEED_LINES,
    _PRIOR_LINE => "$(_PRIOR_LINE)\nnothing # hide"
]

# The three rewrites that insert content are exact-string matches, and
# `build_index` applies them with a bare `replace` that cannot tell a no-match
# from a match, nor one match from two. A README edit that reflows a keyed line
# silently drops what it inserts; a second occurrence silently inserts it twice
# (re-seeding the page part-way through). The page still builds either way, so
# check here. The `NOTICE`/`LICENSE` rewrites are left out: they insert
# nothing and are harmless no-ops if they stop matching. (Appending to a keyed
# line still matches, and the rewrite still lands where it should.)
let readme = read(joinpath(@__DIR__, "..", "README.md"), String)
    for key in (_EXAMPLES_HEADING, _USING_LINE, _PRIOR_LINE)
        matches = length(findall(key, readme))
        matches == 1 ||
            error("docs/docs_config.jl: expected exactly one occurrence of " *
                  "\"$key\" in README.md but found $(matches), so the " *
                  "generated home page would lose the matching entry of " *
                  "INDEX_REWRITES or apply it twice")
    end
end

# Run the README's ```julia fences on the home page so every block shows its
# printed output (#218). The managed build converts EVERY fence into an
# executed `@example readme` block and has no per-fence opt-out, so every fence
# in `README.md` has to be runnable. It is: the README carries no installation
# snippet, which is the one thing that must not run during a docs build (it
# would hit the registry and install a second copy of the package into the docs
# environment). Installation lives on the authored Getting started page
# instead, where a plain ```julia fence renders without executing — the same
# arrangement as the other EpiAware packages.
#
# COST: the last block fits `NUTS()` for 1000 draws, which dominates the build
# time of the home page (minutes, not seconds). The draw count lives in
# `README.md` and is deliberately not reduced here.
const README_EXECUTE = true

# Nothing to omit from the home page: with installation moved to the Getting
# started page, every README section is safe to render and run.
const INDEX_STRIP_SECTIONS = String[]

# The package benchmarks and publishes a timeline to the `benchmarks` branch,
# so the generated page (heading + `docs/benchmarks.md` prose + history) is on.
const BENCHMARK_PAGE = true

const HISTORY_SUITES = String[]

const HISTORY_COMMITS = 5

const HISTORY_REGRESSION_THRESHOLD = 1.1
