#!/usr/bin/env julia
# PACKAGE-OWNED — run selected AD scenarios against selected backends, for
# diagnosis. CI is driven by `runtests.jl` and the ad.yaml matrix.
#
#   julia --project=test/ad test/ad/run_selected.jl --backend enzyme --scenario AR
#
# `--backend` and `--scenario` are repeatable case-insensitive substring
# filters; omit one to select everything.

using ADFixtures, ADTypes, DifferentiationInterface
# Load the backends so DifferentiationInterface registers their extensions.
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
                # Same test as the harness's `check_broken`, so a PASS here
                # cannot disagree with CI.
                ok = g isa AbstractVector && ref !== nothing &&
                     isapprox(g, ref; rtol = 5e-2, atol = 1e-6)
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
                    "ERROR: ", first(msg, 90))
            end
        end
    end

    println()
    println("done: ", npass, " PASS")
end

main()
