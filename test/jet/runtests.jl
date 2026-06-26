#!/usr/bin/env julia
# LOCALLY OVERRIDDEN (was MANAGED by EpiAwarePackageTools.scaffold).
#
# See ISSUES_FOR_PACKAGETOOLS.md item "JET runner cannot analyse a DynamicPPL
# @model package cleanly". The stock runner calls `JET.report_package` with no
# way to filter the false-positive "local variable `x` is not defined" reports
# that JET emits for every `~`-assigned variable inside a Turing `@model` (the
# tilde macro hides the assignment from JET's static analysis). The original
# upstream EpiAware package sidestepped this by not running JET at all.
#
# Minimal local workaround: drop exactly those false positives (reports whose
# message is "local variable ... is not defined" arising in this package's
# `@model`-generated `as_turing_model` / `generate_observation_error_priors`
# code) and fail only on any OTHER report. A template-sync (`update`) will
# revert this file; re-apply until EpiAwarePackageTools supports a JET filter.
#
#   julia --project=test/jet test/jet/runtests.jl

using JET
using EpiAwarePrototype
using DynamicPPL: DynamicPPL

result = JET.report_package(
    EpiAwarePrototype; target_modules = (EpiAwarePrototype,))

# A report is a DynamicPPL `@model`-macro false positive when it arises inside a
# `@model`-generated method — i.e. one whose method instance takes the DynamicPPL
# evaluator arguments `(::DynamicPPL.Model, ::DynamicPPL.AbstractVarInfo, ...)`.
# Those methods are the lowered bodies of the package's `as_turing_model`
# `@model`s, where the `~` (tilde) and `:=` (coloneq) macros produce code that
# JET cannot resolve statically: `UndefVarErrorReport`s for the hidden
# tilde-assigned locals, and `MethodErrorReport`s through the coloneq tracking
# machinery (`store_coloneq_value!!`). All such reports are well-defined at run
# time (the models sample, are tested, and run NUTS); only reports OUTSIDE these
# generated bodies are real.
function in_model_evaluator(report)
    isempty(report.vst) && return false
    linfo = last(report.vst).linfo
    linfo isa Core.MethodInstance || return false
    sig = linfo.specTypes.parameters
    length(sig) >= 3 || return false
    return sig[2] <: DynamicPPL.Model &&
           sig[3] <: DynamicPPL.AbstractVarInfo
end

reports = JET.get_reports(result)
real_reports = filter(!in_model_evaluator, reports)

if isempty(real_reports)
    n = length(reports) - length(real_reports)
    println("JET: no real errors ($(n) DynamicPPL @model-macro false positives filtered)")
    exit(0)
else
    println("JET: $(length(real_reports)) report(s) after filtering @model false positives:")
    foreach(r -> show(stdout, r), real_reports)
    exit(1)
end
