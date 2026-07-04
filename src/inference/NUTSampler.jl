# NUTS sampling method.

@doc raw"
NUTS sampling method for a `DynamicPPL.Model`.

## Fields

  - `target_acceptance`: target acceptance rate.
  - `adtype`: automatic-differentiation backend.
  - `mcmc_parallel`: MCMC parallelisation strategy.
  - `nchains`: number of chains.
  - `max_depth`: NUTS tree-depth limit.
  - `Δ_max`: divergence threshold.
  - `init_ϵ`: initial step size (`0.0` lets NUTS find one).
  - `ndraws`: total draws.
  - `metricT`: HMC metric type.
  - `nadapts`: adaptation steps (`-1` uses the Turing default).
"
@kwdef struct NUTSampler{A <: ADTypes.AbstractADType,
    E <: AbstractMCMC.AbstractMCMCEnsemble, M} <: AbstractEpiSamplingMethod
    "Target acceptance rate."
    target_acceptance::Float64 = 0.8
    "Automatic-differentiation backend."
    adtype::A = AutoForwardDiff()
    "MCMC parallelisation strategy."
    mcmc_parallel::E = MCMCSerial()
    "Number of chains."
    nchains::Int = 1
    "NUTS tree-depth limit."
    max_depth::Int = 10
    "Divergence threshold."
    Δ_max::Float64 = 1000.0
    "Initial step size (`0.0` lets NUTS find one)."
    init_ϵ::Float64 = 0.0
    "Total draws."
    ndraws::Int
    "HMC metric type."
    metricT::M = DiagEuclideanMetric
    "Adaptation steps (`-1` uses the Turing default)."
    nadapts::Int = -1
end

function _apply_method(model::DynamicPPL.Model, method::NUTSampler, prev_result = nothing;
        kwargs...)
    return _apply_nuts(model, method, prev_result; kwargs...)
end

function _apply_nuts(model, method, prev_result; kwargs...)
    nuts = Turing.NUTS(method.target_acceptance; adtype = method.adtype,
        max_depth = method.max_depth, Δ_max = method.Δ_max,
        init_ϵ = method.init_ϵ, metricT = method.metricT)
    # The AbstractMCMC keyword is `n_adapts`; the old `nadapts` was silently
    # swallowed by `kwargs...`, so `NUTSampler.nadapts` never took effect and the
    # default adaptation was always used. `nadapts == -1` is the sentinel for
    # "use the Turing default", so only forward an explicit value.
    adapt_kwargs = method.nadapts == -1 ? (;) : (; n_adapts = method.nadapts)
    return sample(model, nuts,
        method.mcmc_parallel, method.ndraws ÷ method.nchains, method.nchains;
        adapt_kwargs..., kwargs...)
end

function _apply_nuts(model, method, prev_result::PathfinderResult; kwargs...)
    # A Pathfinder pre-step has run; thread its result through as the NUTS
    # initialisation. The mechanism by which earlier EpiAware seeded NUTS from a
    # Pathfinder draw (`init_params = eachrow(draws_transformed.value)`) is gone
    # in current Turing (`initial_params` now requires an `AbstractInitStrategy`,
    # not a vector) and Pathfinder (no `draws_transformed.value` array). The
    # `pathfinder` integration already initialises its own optimisation from the
    # model, so we run NUTS with the default strategy here. Note the Pathfinder
    # result is currently *not* threaded into the returned solution — this path
    # runs the optimisation and then samples from the default init. Warm-starting
    # NUTS from the draw will be reinstated once the init-strategy API stabilises
    # (tracked as a follow-up).
    return _apply_nuts(model, method, nothing; kwargs...)
end
