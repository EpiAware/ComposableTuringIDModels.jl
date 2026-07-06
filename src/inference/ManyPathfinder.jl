# `ManyPathfinder` variational pre-sampler and the `manypathfinder` driver.

@doc raw"
Variational pre-sampler that runs Pathfinder several times and keeps the best run.

## Fields

  - `ndraws`: draws per Pathfinder run.
  - `nruns`: number of Pathfinder runs.
  - `maxiters`: optimiser iterations per run.
  - `max_tries`: extra tries if all runs fail.
"
@kwdef struct ManyPathfinder <: AbstractIDOptMethod
    "Draws per Pathfinder run."
    ndraws::Int = 10
    "Number of Pathfinder runs."
    nruns::Int = 4
    "Optimiser iterations per run."
    maxiters::Int = 100
    "Extra tries if all runs fail."
    max_tries::Int = 100
end

function _apply_method(model::DynamicPPL.Model, method::ManyPathfinder,
        prev_result = nothing; kwargs...)
    return _apply_pathfinder(model, method, prev_result; kwargs...)
end

function _apply_pathfinder(model, method, prev_result; kwargs...)
    return manypathfinder(model, method.ndraws; nruns = method.nruns,
        maxiters = method.maxiters, max_tries = method.max_tries, kwargs...)
end

function _apply_pathfinder(model, method, prev_result::Vector{<:Real}; kwargs...)
    return manypathfinder(model, method.ndraws; init = prev_result,
        nruns = method.nruns, maxiters = method.maxiters,
        max_tries = method.max_tries, kwargs...)
end

@doc raw"
Run Pathfinder several times and return the run with the largest ELBO estimate.

# Arguments

  - `mdl`: the `DynamicPPL.Model` to fit.
  - `ndraws`: draws per Pathfinder run.

# Keyword Arguments

  - `nruns`: number of Pathfinder runs (default `4`).
  - `maxiters`: optimiser iterations per run (default `50`).
  - `max_tries`: extra tries if all runs fail (default `100`).
  - `kwargs...`: forwarded to `pathfinder`.

# Examples
```@example manypathfinder
using ComposableTuringIDModels, Distributions
m = as_turing_model(
    IDModel(
        DirectInfections(; Z = RandomWalk(), initialisation = Normal()),
        PoissonError()), fill(10, 10), 10)
nothing
```
"
function manypathfinder(mdl::DynamicPPL.Model, ndraws; nruns = 4, maxiters = 50,
        max_tries = 100, kwargs...)
    return _run_manypathfinder(mdl; nruns, ndraws, maxiters, kwargs...) |>
           pfs -> _continue_manypathfinder!(pfs, mdl; max_tries, nruns, kwargs...) |>
                  pfs -> _get_best_elbo_pathfinder(pfs)
end

function _run_manypathfinder(mdl::DynamicPPL.Model; nruns, kwargs...)
    @info "Running pathfinder $nruns times"
    pfs = Vector{Union{PathfinderResult, Symbol}}(undef, nruns)
    Threads.@threads for i in 1:nruns
        try
            pfs[i] = pathfinder(mdl; kwargs...)
        catch
            pfs[i] = :fail
        end
    end
    return pfs
end

function _continue_manypathfinder!(pfs, mdl::DynamicPPL.Model; max_tries, nruns,
        kwargs...)
    tryiter = 1
    if all(pfs .== :fail)
        @warn "All initial pathfinder runs failed, trying again for $max_tries tries."
    end
    while all(pfs .== :fail) && tryiter <= max_tries
        new_pf = try
            pathfinder(mdl; kwargs...)
        catch
            :fail
        end
        pfs = vcat(pfs, new_pf)
        tryiter += 1
    end
    if all(pfs .== :fail)
        throw(ErrorException("All pathfinder runs failed after $max_tries tries."))
    end
    return pfs
end

function _get_best_elbo_pathfinder(pfs)
    # Rank runs by their best ELBO across optimiser iterations, matching how
    # Pathfinder itself selects a fit. The last iteration's ELBO
    # (`elbo_estimates[end]`) is not necessarily the ELBO of the returned fit, so
    # ranking on it could prefer a worse run.
    elbos = map(pfs) do pf_res
        pf_res == :fail ? -Inf : maximum(e.value for e in pf_res.elbo_estimates)
    end
    _, choice_of_pf = findmax(elbos)
    return pfs[choice_of_pf]
end
