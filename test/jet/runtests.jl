#!/usr/bin/env julia
# MANAGED by EpiAwarePackageTools.scaffold — do not edit by hand.
#
# JET static-analysis runner, run in this isolated environment so JET's
# JuliaSyntax pin does not clash with the main test deps.
#
#   julia --project=test/jet test/jet/runtests.jl
#
# A DynamicPPL `@model` package gets spurious JET reports for every `~`/`:=`
# line (the tilde macro hides the assignment from JET). Suppress those via a
# package-owned `test/jet/jet_config.jl` defining `JET_REPORT_FILTER`
# (`EpiAwarePackageTools.dynamicppl_model_filter` is the ready-made one).
# Without it the runner fails on any report (the strict default).

using JET
using EpiAwarePackageTools: dynamicppl_model_filter
using ComposableTuringIDModels

const _CONFIG = joinpath(@__DIR__, "jet_config.jl")
isfile(_CONFIG) && include(_CONFIG)

if @isdefined(JET_REPORT_FILTER)
    result = JET.report_package(
        ComposableTuringIDModels;
        target_modules = (ComposableTuringIDModels,)
    )
    kept = filter(JET_REPORT_FILTER, JET.get_reports(result))
    for r in kept
        @info "JET report (not filtered)" report = sprint(show, r)
    end
    isempty(kept) || error("JET found $(length(kept)) report(s)")
    println("JET: no reports survived the configured filter")
else
    JET.test_package(ComposableTuringIDModels; target_modules = (ComposableTuringIDModels,))
end
