@doc raw"
A **prototype** for composable probabilistic infectious disease modelling in
Julia.

`EpiAwarePrototype` builds epidemiological models from small, reusable
components — latent processes, infection processes, and observation models —
each turned into a `Turing`/`DynamicPPL` model by the single generic constructor
[`as_turing_model`](@ref). Components compose by sampling one another as
submodels, so a full model is assembled rather than hand-written.

This package is **ported and adapted** from the open-source, Apache-2.0 licensed
`EpiAware` package; see the `NOTICE` file for attribution. It is exploratory and
clearly labelled as a prototype.

# Examples
```@example
using EpiAwarePrototype, Distributions
data = EpiData([0.2, 0.3, 0.5], exp)
model = EpiAwareModel(RandomWalk(),
    DirectInfections(; data = data, initialisation_prior = Normal()),
    PoissonError())
rand(as_turing_model(model, missing, 20))
```
"
module EpiAwarePrototype

# This package does NOT blanket-reexport Distributions/Turing (as the upstream
# EpiAware also did not): users `using EpiAwarePrototype, Distributions, Turing`.
# Only the names the prototype itself uses or extends are imported below, which
# keeps the public surface to the package's own exports.

using DynamicPPL: DynamicPPL, @model, to_submodel, fix, condition, prefix
using Turing: Turing, filldist, arraydist, sample, MCMCSerial
using CensoredDistributions: double_interval_censored
using LinearAlgebra: dot
using LogExpFunctions: softmax, xexpy, log1pexp
using OrdinaryDiffEq: ODEProblem, ODEFunction, solve, remake, AutoVern7, Rodas5P
using Random: AbstractRNG, randexp

# Inference-layer dependencies.
using ADTypes: ADTypes, AutoForwardDiff
using AbstractMCMC: AbstractMCMC
using AdvancedHMC: DiagEuclideanMetric
using MCMCChains: Chains
using Pathfinder: pathfinder, PathfinderResult
using DataFramesMeta: DataFrame, @rename!
using Tables: rowtable

# Distributions names used (and, for many, extended) by the prototype, imported
# explicitly (not reexported).
using Distributions: Distributions, Distribution, Sampleable,
                     ContinuousUnivariateDistribution, ContinuousDistribution,
                     Normal, Poisson, NegativeBinomial, Gamma, truncated,
                     cdf, ccdf, logcdf, logccdf, invlogcdf, pdf, logpdf, quantile,
                     params, mean, var, std, mode, skewness, kurtosis
using Statistics: Statistics

# --- core architecture ---
export AbstractEpiAwareModel, as_turing_model

# --- utilities and distributions ---
# (double-interval censoring is provided by CensoredDistributions.jl, used
# internally by `EpiData` / `LatentDelay`; it is not re-exported here.)
export accumulate_scan, get_state, HalfNormal, SafePoisson, SafeNegativeBinomial,
       NegativeBinomialMeanClust, condition_model

# --- latent models ---
export IID, HierarchicalNormal, RandomWalk, AR, MA, Intercept, FixedIntercept,
       Null, DiffLatentModel

# --- latent modifiers / manipulators / combinations / broadcasting ---
export TransformLatentModel, PrefixLatentModel, RecordExpectedLatent,
       CombineLatentModels, ConcatLatentModels, BroadcastLatentModel,
       RepeatEach, RepeatBlock, broadcast_rule, broadcast_n, broadcast_dayofweek,
       broadcast_weekly, equal_dimensions, arma, arima

# --- infection models ---
export EpiData, DirectInfections, ExpGrowthRate, Renewal,
       R_to_r, r_to_R, expected_Rt

# --- ODE compartmental models ---
export SIRParams, SEIRParams, ODEProcess

# --- observation models ---
export PoissonError, NegativeBinomialError, LatentDelay,
       observation_error, generate_observation_error_priors

# --- observation modifiers / manipulators ---
export Ascertainment, ascertainment_dayofweek, Aggregate, PrefixObservationModel,
       RecordExpectedObs, TransformObservationModel, StackObservationModels

# --- composition ---
export EpiAwareModel

# --- inference orchestration ---
export EpiProblem, EpiMethod, NUTSampler, ManyPathfinder, DirectSample,
       manypathfinder, apply_method, EpiAwareObservables, generated_observables,
       spread_draws, get_param_array

# --- core architecture ---
include("base/base.jl")
include("base/accumulation.jl")
include("base/prettyprinting.jl")

# --- utilities and distributions ---
include("utils/accumulate_scan.jl")
include("utils/HalfNormal.jl")
include("utils/SafeInt.jl")
include("utils/SafePoisson.jl")
include("utils/SafeNegativeBinomial.jl")
include("utils/censored_pmf.jl")
include("utils/turing-methods.jl")

# --- latent models ---
include("latent_models/utils.jl")
include("latent_models/models/IID.jl")
include("latent_models/models/HierarchicalNormal.jl")
include("latent_models/models/RandomWalk.jl")
include("latent_models/models/AR.jl")
include("latent_models/models/MA.jl")
include("latent_models/models/Intercept.jl")
include("latent_models/models/Null.jl")
include("latent_models/modifiers/DiffLatentModel.jl")
include("latent_models/modifiers/TransformLatentModel.jl")
include("latent_models/modifiers/PrefixLatentModel.jl")
include("latent_models/modifiers/RecordExpectedLatent.jl")
include("latent_models/manipulators/CombineLatentModels.jl")
include("latent_models/manipulators/ConcatLatentModels.jl")
include("latent_models/manipulators/broadcast/LatentModel.jl")
include("latent_models/manipulators/broadcast/rules.jl")
include("latent_models/manipulators/broadcast/helpers.jl")
include("latent_models/combinations/arma.jl")
include("latent_models/combinations/arima.jl")

# --- infection models ---
include("infection_models/EpiData.jl")
include("infection_models/utils.jl")
include("infection_models/DirectInfections.jl")
include("infection_models/ExpGrowthRate.jl")
include("infection_models/RenewalSteps.jl")
include("infection_models/Renewal.jl")

# --- ODE compartmental models ---
include("ode/SIRParams.jl")
include("ode/SEIRParams.jl")
include("ode/ODEProcess.jl")

# --- observation models ---
include("observation_models/ObservationErrorModels/methods.jl")
include("observation_models/ObservationErrorModels/PoissonError.jl")
include("observation_models/ObservationErrorModels/NegativeBinomialError.jl")
include("observation_models/modifiers/LatentDelay.jl")
include("observation_models/modifiers/ascertainment/Ascertainment.jl")
include("observation_models/modifiers/ascertainment/helpers.jl")
include("observation_models/modifiers/Aggregate.jl")
include("observation_models/modifiers/PrefixObservationModel.jl")
include("observation_models/modifiers/RecordExpectedObs.jl")
include("observation_models/modifiers/TransformObservationModel.jl")
include("observation_models/StackObservationModels.jl")

# --- composition ---
include("compose.jl")

# --- inference orchestration ---
include("inference/types.jl")
include("inference/EpiProblem.jl")
include("inference/EpiAwareObservables.jl")
include("inference/apply_method.jl")
include("inference/NUTSampler.jl")
include("inference/DirectSample.jl")
include("inference/get_param_array.jl")
include("inference/ManyPathfinder.jl")
include("inference/post-inference.jl")

end
