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

const _GITHUB_LINE = "  [GitHub repository]($(_REPO_URL))."

const _INSTALL_POINTER = "- Ready to install? See [Installation](@ref)."

# `README.md` is owned by its GitHub readers and is edited only there; every
# adjustment the docs home page needs is made here, on the generated copy.
# `build_index` applies these rewrites line by line to the README, inside
# ```julia fences as well as in prose, so a replacement containing `\n` inserts
# lines into a code block.
#
#  1. `NOTICE`/`LICENSE` are repo-root files. The README's relative links are
#     right on GitHub, but the generated `index.md` is served from the docs
#     site, where they resolve to nothing — so point the docs copy at `main`.
#  2. The home page's blocks are executed (see `README_EXECUTE`), and the first
#     thing they do is draw from the prior. Unseeded, the page is different on
#     every build: this model's latent process is a random walk on the log
#     scale, so its prior over case counts is heavy-tailed. Over seeds 1:60 the
#     simulated series ran from all-zero to a peak of 1e36 cases, and the fit
#     that follows was anywhere from 4 s to 240 s and usually unconverged.
#     Seeding pins the page to one reproducible draw. Seed 27 gives an epidemic
#     curve that rises from 6 to a peak of 74 cases and decays back to 9 over
#     30 days, and a fit that mixes (rhat 1.01, min ess_bulk 327 over the 34
#     parameters). Do not change it without re-checking both: the fit degrades
#     sharply as the simulated counts grow, because a Poisson likelihood on
#     large counts pins the log-scale latent series into a very ill-conditioned
#     posterior. Seed 15 (peak 5210 cases) leaves NUTS stuck after 1000 draws,
#     at rhat 2.02 and min ess_bulk 2.6.
#  3. `as_turing_model` returns a `DynamicPPL.Model`, whose display is six lines
#     of nested type parameters. That is noise on the home page, so the docs
#     copy ends the line with `;`. The block's job is to build the model; the
#     blocks that follow carry the output worth reading.
#  4. The Installation section is stripped from the home page (see
#     `INDEX_STRIP_SECTIONS`), so add a pointer to where it went.
const INDEX_REWRITES = Pair{String, String}[
    "(NOTICE)" => "($(_REPO_BLOB)/NOTICE)",
    "(LICENSE)" => "($(_REPO_BLOB)/LICENSE)",
    _USING_LINE => "$(_USING_LINE)\nusing Random\nRandom.seed!(27)",
    _PRIOR_LINE => "$(_PRIOR_LINE);",
    _GITHUB_LINE => "$(_GITHUB_LINE)\n$(_INSTALL_POINTER)"
]

# Run the README's ```julia fences on the home page so every block shows its
# printed output (#218). The managed build converts EVERY fence into an
# executed `@example readme` block and has no per-fence opt-out, so the one
# fence that must not run, `Pkg.add("ComposableTuringIDModels")` under
# Installation, is handled by stripping that section from the generated index
# (see `INDEX_STRIP_SECTIONS` below). Executing it would hit the registry and
# install a second copy of the package into the docs environment. Installation
# stays in `README.md` for GitHub readers, is linked from the home page (see
# `INDEX_REWRITES`), and the docs site carries it on the Overview page.
#
# COST: the last block fits `NUTS()` for 1000 draws, which dominates the cost
# of the home page. At the pinned seed it took 226-319 s on an Apple-silicon
# laptop, the spread reflecting other load on the machine rather than the fit.
# Unseeded it ranged from 4 s to 240 s across seeds, the fast end being draws
# so degenerate that the sampler got stuck rather than converged. The draw
# count lives in `README.md` and is deliberately not reduced here.
const README_EXECUTE = true

const INDEX_STRIP_SECTIONS = String["Installation"]

# The package benchmarks and publishes a timeline to the `benchmarks` branch,
# so the generated page (heading + `docs/benchmarks.md` prose + history) is on.
const BENCHMARK_PAGE = true

const HISTORY_SUITES = String[]

const HISTORY_COMMITS = 5

const HISTORY_REGRESSION_THRESHOLD = 1.1
