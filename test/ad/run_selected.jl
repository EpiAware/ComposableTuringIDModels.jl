#!/usr/bin/env julia
# PACKAGE-OWNED — selectively run a subset of AD scenarios against a subset of
# backends for fast diagnosis/fix iteration, without a full per-backend suite.
#
# Usage (run from the repo, in the AD env):
#   julia --project=test/ad test/ad/run_selected.jl \
#       --backend "Enzyme forward" \
#       --backend "ReverseDiff (tape)" \
#       --scenario "LatentDelay" \
#       --scenario "ARIMA"
#
# `--backend` and `--scenario` are repeatable, case-insensitive SUBSTRING
# filters. Omit a filter to select everything. For each selected scenario the
# script runs each selected backend's `gradient` against the scenario's
# ForwardDiff reference and prints PASS / MISMATCH / ERROR (with a snippet of
# the error). This is a diagnostic helper only; it is not part of the CI suite
# (`test/ad/runtests.jl` + the ad.yaml matrix drive CI).
#
#   julia --project=test/ad test/ad/run_selected.jl --backend enzyme --scenario AR
#     # enzyme reverse + forward over every scenario whose name contains "AR"

using ADFixtures, ADTypes, DifferentiationInterface
# `using` the backend packages so DifferentiationInterface registers each
# backend's extension (without this, `gradient` falls back to an
# unloaded-extension `_prepare_pullback_aux`/`_prepare_pushforward_aux` error).
using ForwardDiff, ReverseDiff, Enzyme, Mooncake

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

    matches(filters, s) = isempty(filters) ||
                          any(f -> occursin(lowercase(f), lowercase(s)), filters)

    scens = ADFixtures.scenarios(with_reference = true)
    backends = ADFixtures.backends()

    sel_scens = [s for s in scens if matches(scenario_filters, String(s.name))]
    sel_backends = [b for b in backends if matches(backend_filters, b.name)]

    isempty(sel_scens) && error("no scenarios match: $scenario_filters")
    isempty(sel_backends) && error("no backends match: $backend_filters")

    println("Scenarios (", length(sel_scens), "): ", join(
        [String(s.name)
         for s in sel_scens], "; "))
    println("Backends (", length(sel_backends), "): ", join([b.name for b in sel_backends], "; "))
    println()

    npass = 0
    for s in sel_scens
        ref = s.res1
        for b in sel_backends
            try
                g = DifferentiationInterface.gradient(s.f, b.backend, s.x)
                ok = length(g) == length(ref) &&
                     maximum(abs.(g .- ref)) <= max(1e-6, 5e-2 * maximum(abs.(ref)))
                if ok && all(isfinite, g)
                    npass += 1
                    println(rpad(String(s.name), 46), rpad(b.name, 22), "PASS")
                else
                    println(rpad(String(s.name), 46), rpad(b.name, 22),
                        "MISMATCH (finite=", all(isfinite, g), ")")
                end
            catch e
                msg = replace(sprint(showerror, e), '\n' => ' ')
                println(rpad(String(s.name), 46), rpad(b.name, 22),
                    "ERROR: ", msg[1:min(end, 90)])
            end
        end
    end

    println()
    println("done: ", npass, " PASS")
end

main()
