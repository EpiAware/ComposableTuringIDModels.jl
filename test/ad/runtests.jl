#!/usr/bin/env julia
# MANAGED by EpiAwarePackageTools.scaffold — do not edit by hand.
#
# AD gradient test entry, organised as `@testitem`s and run with
# TestItemRunner. The AD items live in their own environment (Enzyme, Mooncake,
# etc. are not main-test deps) and in dedicated per-backend CI.
#
#   julia --project=test/ad test/ad/runtests.jl              # all backends
#   julia --project=test/ad test/ad/runtests.jl enzyme_reverse  # one tag
#
# Per-backend tags let the per-backend CI run a single backend so a transiently
# unstable backend only reds its own badge. With no argument every AD item runs.
#
# A tag that matches nothing is an error, not a pass. The kit adds backends to
# the managed `ad.yaml` matrix over time, but `scenarios.jl` and `ADFixtures`
# are package-owned seeds a sync cannot extend, so a package can end up with a
# CI job for a backend it has no test items for. Left to report green, that job
# also uploads an empty coverage flag behind a public badge (kit#415).

using TestItemRunner

if isempty(ARGS)
    TestItemRunner.run_tests(@__DIR__)
else
    selected = Symbol.(ARGS)
    # `filter` is called once per discovered item, so it doubles as the
    # matched-anything check without reaching into TestItemRunner internals.
    matched = Ref(false)
    TestItemRunner.run_tests(
        @__DIR__; filter = function (ti)
            hit = any(in(ti.tags), selected)
            hit && (matched[] = true)
            return hit
        end
    )
    matched[] || error(
        "no AD test item carries any of the requested tags: " *
            join(ARGS, ", ") * ". Add a test item for the backend in " *
            "test/ad/scenarios.jl and the matching entry to " *
            "ADFixtures.backends(), or drop the backend from the ad.yaml " *
            "matrix. An empty run must not report green (kit#415)."
    )
end
