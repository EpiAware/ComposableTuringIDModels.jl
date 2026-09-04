# PACKAGE-OWNED — not part of the scaffold template.
#
# Per-backend NUTS sampling smoke tests, the sampling half of the qualification
# the AD backend comparison page describes.
#
# Everything about a scenario comes from the same `ADFixtures` registry the
# gradient items use: the models, the sampler budget, and the per-backend
# broken and skip bookkeeping. A sampling scenario carries its ensemble
# strategy in its name, so it can be declared broken or skipped without
# touching the gradient scenario of the same model.

@testsnippet SamplingSmoke begin
    using Test
    using ADFixtures

    const CHILD = joinpath(@__DIR__, "sampling_child.jl")

    # Wall-clock cap on a single scenario, and how often the child's progress
    # is checked against it. Seven per-backend jobs share one CI time limit, so
    # a hung sampler is reported as a failure rather than left to spend the
    # budget. The cap is per scenario, not per child, so a slow first-use
    # compile cannot leave a healthy second scenario with no time.
    const SCENARIO_TIMEOUT = 900
    const POLL_INTERVAL = 5

    # The status a scenario carries when it never ran.
    const SKIPPED = "skipped by the registry"

    # Run the child over `names` and return its output plus a description of
    # how it exited. Output goes to a file rather than a pipe so it survives a
    # process that dies mid-write.
    function run_child(backend_name, names)
        cmd = `$(Base.julia_cmd()) --project=$(@__DIR__) --threads=2
            --startup-file=no $CHILD $backend_name $names`
        path, io = mktemp()
        proc = run(pipeline(cmd; stdout = io, stderr = io), wait = false)
        # Each result the child flushes restarts the clock, so the cap applies
        # to whichever scenario is currently running.
        deadline = Ref(time() + SCENARIO_TIMEOUT)
        reported = Ref(0)
        timer = Timer(POLL_INTERVAL; interval = POLL_INTERVAL) do _
            done = length(child_results(read(path, String)))
            if done > reported[]
                reported[] = done
                deadline[] = time() + SCENARIO_TIMEOUT
            end
            time() > deadline[] && process_running(proc) &&
                kill(proc, Base.SIGKILL)
        end
        try
            wait(proc)
        finally
            close(timer)
            close(io)
        end
        how = proc.termsignal != 0 ?
            "killed by signal " * string(proc.termsignal) :
            "exit code " * string(proc.exitcode)
        out = read(path, String)
        rm(path; force = true)
        return out, how
    end

    # Pull the child's result lines out of its combined output.
    function child_results(out)
        results = Dict{String, String}()
        for line in eachline(IOBuffer(out))
            startswith(line, "SMOKE\t") || continue
            _, name, status = split(line, '\t'; limit = 3)
            results[name] = status
        end
        return results
    end

    # Sample every non-skipped scenario for `backend_name`, one result string
    # per scenario name. A child that dies takes the first scenario it had not
    # yet reported with it, so that one is recorded as crashed and the rest are
    # retried in a fresh child.
    function smoke_statuses(backend_name, all_names, skip)
        statuses = Dict{String, String}(n => SKIPPED for n in all_names if n in skip)
        pending = filter(n -> !(n in skip), all_names)
        while !isempty(pending)
            out, how = run_child(backend_name, pending)
            for (name, status) in child_results(out)
                name in pending && (statuses[name] = status)
            end
            remaining = filter(n -> !haskey(statuses, n), pending)
            isempty(remaining) && break
            statuses[first(remaining)] = "CRASH child process died (" * how * ")"
            pending = remaining[2:end]
        end
        return statuses
    end

    # Drive the sampling smoke scenarios for one backend.
    #
    # A scenario the registry skips is reported skipped, one it declares broken
    # records `@test_broken` when it fails and passes when it does not (the
    # shape `check_broken` uses, so over-listing is safe), and anything else is
    # a hard test. The assertion compares the status string so a failure prints
    # the sampler error, or the signal that killed the child, rather than a
    # bare `false`.
    function test_sampling_smoke(backend_name)
        skip = get(
            ADFixtures.backend_skip_scenarios(), backend_name, Set{String}()
        )
        broken = union(
            Set(ADFixtures.broken_scenario_names()),
            get(
                ADFixtures.backend_broken_scenarios(), backend_name,
                Set{String}()
            )
        )
        all_names = [s.name for s in ADFixtures.sampling_scenarios()]
        statuses = smoke_statuses(backend_name, all_names, skip)
        for name in all_names
            status = statuses[name]
            @testset "$name" begin
                if name in skip
                    @test_skip status == "PASS"
                elseif name in broken && status != "PASS"
                    @test_broken status == "PASS"
                else
                    @test status == "PASS"
                end
            end
        end
        return nothing
    end
end

@testitem "ForwardDiff sampling smoke" tags = [:ad, :smoke, :forwarddiff] setup = [SamplingSmoke] begin
    test_sampling_smoke("ForwardDiff")
end

@testitem "ReverseDiff compiled sampling smoke" tags = [:ad, :smoke, :reversediff_compiled] setup = [SamplingSmoke] begin
    test_sampling_smoke("ReverseDiff (compiled)")
end

@testitem "ReverseDiff sampling smoke" tags = [:ad, :smoke, :reversediff] setup = [SamplingSmoke] begin
    test_sampling_smoke("ReverseDiff (tape)")
end

@testitem "Enzyme reverse sampling smoke" tags = [:ad, :smoke, :enzyme, :enzyme_reverse] setup = [SamplingSmoke] begin
    test_sampling_smoke("Enzyme reverse")
end

@testitem "Enzyme forward sampling smoke" tags = [:ad, :smoke, :enzyme, :enzyme_forward] setup = [SamplingSmoke] begin
    test_sampling_smoke("Enzyme forward")
end

@testitem "Mooncake reverse sampling smoke" tags = [:ad, :smoke, :mooncake, :mooncake_reverse] setup = [SamplingSmoke] begin
    test_sampling_smoke("Mooncake reverse")
end

@testitem "Mooncake forward sampling smoke" tags = [:ad, :smoke, :mooncake, :mooncake_forward] setup = [SamplingSmoke] begin
    test_sampling_smoke("Mooncake forward")
end
