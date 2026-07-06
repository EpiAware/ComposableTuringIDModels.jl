# PACKAGE-OWNED — scaffold writes this once and never overwrites it.
#
# Optional JET configuration for the isolated runner (test/jet/runtests.jl). If
# this file defines `JET_REPORT_FILTER` (a `report -> Bool` predicate; a report
# is KEPT when it returns `true`), the runner switches from `test_package` to
# `report_package` + filter and fails only on reports the predicate keeps.
#
# ComposableTuringIDModels's public surface is DynamicPPL `@model` functions
# (`as_turing_model` and friends), so JET emits a false `UndefVarErrorReport`
# for every `~`-assigned local (and `MethodErrorReport`s through the `:=`
# tracker) because the tilde macro hides the assignment from JET's static
# analysis. `dynamicppl_model_filter` (shipped by EpiAwarePackageTools) drops
# exactly those and keeps any genuine report.

const JET_REPORT_FILTER = dynamicppl_model_filter
