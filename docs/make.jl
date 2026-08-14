# MANAGED by EpiAwarePackageTools.scaffold — do not edit by hand.
#
# Thin entry point for the standard EpiAware documentation build. All build
# logic lives in `EpiAwarePackageTools.DocsBuild.build_docs`; this file only
# wires `pages.jl` (managed, regenerated on every sync unless the package
# forked it -- see EpiAwarePackageTools' `_apply_pages`) + the package-owned
# `docs_config.jl` into that call, so it can be re-applied on every `update`
# without losing package content.
#
# `build_docs`:
#   - runs the Literate tutorial pipeline (light in-process, heavy one per
#     subprocess) driven by `docs_config.jl`; under `--skip-notebooks` light
#     tutorials still render in-process and only heavy ones fall back to
#     heading stubs; any `FORCE_STUB_TUTORIALS` entry always stubs,
#     independent of that flag, while its heavy siblings run normally, and
#     any `TUTORIAL_ENVIRONMENTS` entry resolves against the environment it
#     names instead of the shared docs one, and `HEAVY_TUTORIAL_WORKERS`
#     runs that many heavy tutorials at once,
#   - runs the same pipeline over `docs/src/benchmarks/`, driven by
#     `HEAVY_BENCHMARKS`/`BENCHMARK_STUBS`, so a benchmark report gets its
#     own nav group rather than sitting under Tutorials,
#   - generates `src/index.md` from the README (badges stripped, any
#     `INDEX_STRIP_SECTIONS` removed, link rewrites applied),
#   - generates `src/release-notes.md` from the repo's published GitHub
#     releases, fetched at build time (a link-out when they cannot be read),
#   - generates `src/benchmarks/over-time.md` (skeleton + package-owned
#     `docs/benchmarks.md` prose + rendered performance history),
#   - generates the API pages from the module's documented bindings, and
#   - renders + deploys with DocumenterVitepress.
#
# Build it with `task docs` (or `julia --project=docs docs/make.jl`).

using Pkg: Pkg
Pkg.instantiate()

using EpiAwarePackageTools
using ComposableTuringIDModels

# `pages.jl` (nav tree, managed -- see EpiAwarePackageTools' `_apply_pages`)
# and `docs_config.jl` (tutorial lists, link rewrites, linkcheck ignores,
# package-owned) can both be absent, e.g. an adopter predating either. Guard
# the include so a re-applied managed `make.jl` still loads and falls back to
# defaults (#163); `_cfg` then defaults any key a missing or older config
# predates. The fallback warns because a silently-defaulted `pages` publishes
# a Home-only nav that a green docs run would hide (#188).
for _f in ("pages.jl", "docs_config.jl")
    if isfile(joinpath(@__DIR__, _f))
        include(joinpath(@__DIR__, _f))
    else
        @warn "docs/$(_f) not found; building with defaults " *
            "(a missing pages.jl leaves the site with a Home-only nav). " *
            "Run `scaffold`/`update` to (re)write it."
    end
end

# Read a package-owned config const, defaulting when a missing or older
# `docs_config.jl`/`pages.jl` predates it.
_cfg(sym, default) = isdefined(@__MODULE__, sym) ?
    getfield(@__MODULE__, sym) : default

build_docs(
    ComposableTuringIDModels;
    repo = "EpiAware/ComposableTuringIDModels.jl",
    authors = "Sam Abbott",
    deploy_url = "https://composableturingidmodels.epiaware.org",
    pages = _cfg(:pages, ["Home" => "index.md"]),
    skip_notebooks = "--skip-notebooks" in ARGS ||
        get(ENV, "SKIP_NOTEBOOKS", "false") == "true",
    tutorials_subdir = _cfg(
        :TUTORIALS_SUBDIR,
        joinpath("getting-started", "tutorials")
    ),
    light_tutorials = _cfg(:LIGHT_TUTORIALS, String[]),
    heavy_tutorials = _cfg(:HEAVY_TUTORIALS, String[]),
    tutorial_stubs = _cfg(:TUTORIAL_STUBS, Pair{String, String}[]),
    force_stub_tutorials = _cfg(:FORCE_STUB_TUTORIALS, String[]),
    # Heavy tutorials that resolve against their own environment rather than
    # the shared docs one, as `"file.jl" => "environment/dir"` pairs. Empty
    # for a `docs_config.jl` predating the const, so every tutorial builds
    # against `docs/` exactly as before.
    tutorial_environments = _cfg(
        :TUTORIAL_ENVIRONMENTS, Pair{String, String}[]
    ),
    # How many heavy tutorials execute at once, each still in its own
    # subprocess, with the thread budget divided between them. Defaults to 1
    # (one after another, as before) for a `docs_config.jl` predating the
    # const, and for every package that has not opted in.
    heavy_tutorial_workers = _cfg(:HEAVY_TUTORIAL_WORKERS, 1),
    # The `docs/src/benchmarks/` pipeline (e.g. `ad-comparison.jl`), same
    # convention as the tutorials pipeline above but its own nav group
    # (#299/#305). `_cfg` defaults both to empty for a `docs_config.jl` that
    # predates this pipeline, so a package that has not yet added the two
    # consts still builds -- just without that page rendered.
    heavy_benchmarks = _cfg(:HEAVY_BENCHMARKS, String[]),
    benchmark_stubs = _cfg(:BENCHMARK_STUBS, Pair{String, String}[]),
    linkcheck_ignore = _cfg(:LINKCHECK_IGNORE, Regex[]),
    index_rewrites = _cfg(:INDEX_REWRITES, Pair{String, String}[]),
    readme_execute = _cfg(:README_EXECUTE, true),
    index_strip_sections = _cfg(:INDEX_STRIP_SECTIONS, String[]),
    benchmark_page = _cfg(:BENCHMARK_PAGE, false),
    # Performance-history rendering (#193): headline suites + cap on the
    # summary/detail to the most-recent revisions. Both default to the whole
    # timeline when a package predates these config keys.
    history_suites = _cfg(:HISTORY_SUITES, String[]),
    history_commits = _cfg(:HISTORY_COMMITS, 5),
    # Regression cutoff: ratio (against the oldest shown revision) at or
    # above which a suite's `Status` flags "⚠ reg".
    history_regression_threshold = _cfg(:HISTORY_REGRESSION_THRESHOLD, 1.1),
    # Extra docstring-owning modules for a re-export the alias walk cannot
    # reach; re-exported API owners are auto-discovered, so most packages
    # leave this empty (#175).
    extra_modules = _cfg(:EXTRA_MODULES, Module[])
)
