#!/usr/bin/env julia
# MANAGED by EpiAwarePackageTools.scaffold — do not edit by hand.
#
# Run selected AD scenarios against selected backends, for fast
# diagnosis/fix/prototyping. CI is driven by `runtests.jl` and the ad.yaml
# matrix.
#
#   julia --project=test/ad test/ad/run_selected.jl --backend enzyme \
#       --scenario AR
#
# `--backend` and `--scenario` are repeatable case-insensitive substring
# filters; omit one to select everything.

using ADFixtures, ADTypes, DifferentiationInterface
# Load the backends so DifferentiationInterface registers their extensions.
using ForwardDiff, ReverseDiff, Enzyme, Mooncake
using EpiAwarePackageTools

function main()
    backend_filters = String[]
    scenario_filters = String[]
    i = 1
    while i <= length(ARGS)
        a = ARGS[i]
        if a == "--backend"
            i += 1
            push!(backend_filters, ARGS[i])
        elseif a == "--scenario"
            i += 1
            push!(scenario_filters, ARGS[i])
        else
            error("unknown argument: $a (expected --backend/--scenario)")
        end
        i += 1
    end

    EpiAwarePackageTools.run_selected(
        ADFixtures; backends = backend_filters, scenarios = scenario_filters
    )
    return nothing
end

main()
