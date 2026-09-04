#!/usr/bin/env julia
# PACKAGE-OWNED.
#
# Runs the named sampling smoke scenarios for one backend and prints a result
# line per scenario. Driven by the `@testitem`s in `test/ad/sampling.jl`.
#
#   julia --project=test/ad --threads=2 test/ad/sampling_child.jl \
#       "Enzyme reverse" "Renewal+NegativeBinomial posterior (MCMCSerial)"
#
# This runs in a child process because the failure it guards against is a
# segfault, which no in-process `try` can catch. Scenarios are run in the given
# order and each result is flushed as it lands, so the parent can tell which
# scenario a dead process was on from the results that never arrived.

using ADFixtures
# Load every backend package so DifferentiationInterface registers its
# extensions, as `run_selected.jl` does.
using ForwardDiff, ReverseDiff, Enzyme, Mooncake
using ComposableTuringIDModels: apply_method, NUTSampler
using Random: Random

# Prefix marking a result line, so the parent can pick it out of the sampler's
# own chatter on the same stream.
const MARKER = "SMOKE"

# Collapse an error into the single line a result carries.
function _one_line(err)
    text = sprint(showerror, err)
    return first(replace(text, r"\s+" => " "), 300)
end

function main()
    backend_name = ARGS[1]
    wanted = ARGS[2:end]
    adtype = only(
        filter(e -> e.name == backend_name, ADFixtures.backends())
    ).backend
    by_name = Dict(s.name => s for s in ADFixtures.sampling_scenarios())
    for name in wanted
        scen = by_name[name]
        # Seed the sampler so a backend's pass or fail does not depend on the
        # trajectory it happens to draw.
        Random.seed!(1)
        status = try
            apply_method(
                scen.model,
                NUTSampler(; adtype = adtype, scen.method_kwargs...);
                progress = false
            )
            "PASS"
        catch err
            "FAIL " * _one_line(err)
        end
        println(stdout, MARKER, "\t", name, "\t", status)
        flush(stdout)
    end
    return nothing
end

main()
