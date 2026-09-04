@doc raw"
Composable probabilistic infectious disease modelling in Julia.

`ComposableTuringIDModels` builds epidemiological models from small, reusable
components — infection processes (each owning its own latent parameter process)
and observation models — each turned into a `Turing`/`DynamicPPL` model by the
single generic constructor [`as_turing_model`](@ref). Components compose by
sampling one another as submodels, so a full model is assembled rather than
hand-written.

This package is **ported and adapted** from the open-source, Apache-2.0 licensed
`EpiAware` package; see the `NOTICE` file for attribution.

# Examples
```@example
using ComposableTuringIDModels, Distributions
model = IDModel(
    DirectInfections(; Z = RandomWalk(), initialisation = Normal()),
    PoissonError())
rand(as_turing_model(model, missing, 20))
```
"
module ComposableTuringIDModels

# Distributions and Turing are deliberately not blanket-reexported, so users
# write `using ComposableTuringIDModels, Distributions, Turing`.
# Only the names the package itself uses or extends are imported below.

using DynamicPPL: DynamicPPL, @model, to_submodel, prefix
using Turing: Turing, filldist, predict
using FlexiChains: FlexiChains
using CensoredDistributions: double_interval_censored
# `Accessors` and `rewrap` rebuild a component through
# `ConstructionBase.constructorof`. No component here declares one. `rewrap`'s
# docstring says what a derived field needs instead.
using ConstructionBase: ConstructionBase
using LinearAlgebra: dot, cholesky, Symmetric, I, UniformScaling
using LogExpFunctions: softmax, xexpy, log1pexp
using OrdinaryDiffEq: ODEProblem, ODEFunction, solve, remake, AutoVern7, Rodas5P
using Random: AbstractRNG, randexp, default_rng

# KernelFunctions.jl supplies the covariance kernel types.
# The package adds only the 1-D spectral densities the Hilbert-space
# approximation needs.
using KernelFunctions: Kernel, SqExponentialKernel, Matern32Kernel,
    Matern52Kernel, with_lengthscale, kernelmatrix

# `ReportTriangle` accepts reports as a long-form table of rows.
using Tables: rowtable

# Distributions names used (and, for many, extended) by the package, imported
# explicitly (not reexported).
using Distributions: Distributions, Distribution, Sampleable,
    UnivariateDistribution,
    ContinuousUnivariateDistribution, ContinuousDistribution,
    Normal, Poisson, NegativeBinomial, Binomial, Gamma, truncated,
    cdf, ccdf, logcdf, logccdf, invlogcdf, pdf, logpdf, quantile,
    params, mean, var, std, mode, skewness, kurtosis,
    product_distribution
using Statistics: Statistics

# --- core architecture ---
export AbstractComposableModel, as_turing_model
export AbstractPriorModel, AbstractLatentModel, AbstractInfectionModel,
    AbstractObservationModel, AbstractObservationErrorModel
export implements_prior_interface,
    implements_infection_interface, implements_observation_interface
export as_turing_submodel

# --- utilities and distributions ---
# Double-interval censoring comes from CensoredDistributions.jl and is not
# re-exported here.
export accumulate_scan, get_state, HalfNormal, SafePoisson, SafeNegativeBinomial,
    NegativeBinomialMeanClust

# --- latent models ---
export IID, HierarchicalNormal, RandomWalk, AR, MA, Intercept,
    FixedIntercept,
    Null, DiffLatentModel, HilbertSpaceGP, ExactGP
# Covariance kernels are re-exported from KernelFunctions.jl.
# `spectral_density` adds the Hilbert-space weights those kernels need.
export SqExponentialKernel, Matern32Kernel, Matern52Kernel, spectral_density,
    standardised_index

# --- latent modifiers / manipulators / combinations / broadcasting ---
export TransformLatentModel, PrefixLatentModel, RecordExpectedLatent,
    CombineLatentModels, ConcatLatentModels, BroadcastLatentModel,
    RepeatEach, RepeatBlock, InWindow, broadcast_rule, broadcast_n,
    broadcast_dayofweek, broadcast_weekly, broadcast_window, equal_dimensions,
    arma, arima, Hierarchy,
    Stratify, Replicate

# --- infection models ---
export DirectInfections, ExpGrowthRate, Renewal, SeedingPath,
    RenewalStep, SusceptibleDepletion, ImportedCases,
    R_to_r, r_to_R, expected_Rt
export CombineInfections

# --- coupling between strata ---
export renewal_pressure, pairwise_gen_int, AbstractMixingModel, MixingStep,
    Gravity, gravity

# --- ODE compartmental models ---
export SIRParams, SEIRParams, ODEProcess, CatalystODEParams

# --- observation models ---
export PoissonError, NegativeBinomialError, NormalError, BinomialError, LatentDelay,
    UncertainDelay, observation_error, generate_observation_error_priors,
    define_y_t

# --- observation modifiers / manipulators ---
export Ascertainment, ascertainment_dayofweek, Aggregate, RightTruncate,
    ReportingCDF, ReportTriangle, ReportingTriangle, ReportingPMF,
    PrefixObservationModel, RecordExpectedObs, TransformObservationModel

# --- observation composition ---
export Split, StrataMap

# --- observation diagnostics ---
# What data a model needs, and how much of the series its delays consume.
export data_requirements, data_fits, observation_lead_in
export DataRequirements, StreamRequirement

# --- composition ---
export IDModel

# --- inference orchestration ---
export IDProblem, forecast

# --- extension points ---
# Names a component author implements against but rarely calls.
# Public but not exported, so they are documented and supported without crowding
# the namespace of a `using` call.
public ModelShape, across_shape, infection_strata,
    AbstractAccumulationStep, AbstractConstantRenewalStep,
    ConstantRenewalStep, AbstractRenewalModifier, modifier_init_state,
    apply_modifier, renewal_foi, renewal_init_state, renewal_init_window,
    MissingObservations,
    at, path_prior, prior_order, assert_prior_length,
    wrapped_models, observation_components, rewrap, swap

# --- core architecture ---
include("base/base.jl")
include("base/roles.jl")
include("base/shapes.jl")
include("base/interfaces.jl")
include("base/priors.jl")
include("base/prettyprinting.jl")

# --- accumulation steps ---
# Included early so every step type is defined before the model components that
# construct them.
include("steps/AbstractAccumulationStep.jl")
include("steps/accumulate_scan.jl")
include("steps/RWStep.jl")
include("steps/ARStep.jl")
include("steps/TVARStep.jl")
include("steps/MAStep.jl")
include("steps/LDStep.jl")
include("steps/TimeVaryingLDStep.jl")
include("steps/RenewalSteps.jl")
include("steps/RenewalStep.jl")
include("steps/MixingStep.jl")
include("steps/Gravity.jl")
include("steps/ImportedCases.jl")

# --- utilities and distributions ---
include("utils/HalfNormal.jl")
include("utils/SafeInt.jl")
include("utils/SafePoisson.jl")
include("utils/SafeNegativeBinomial.jl")
include("utils/censored_pmf.jl")

# --- latent models ---
include("latent_models/models/IID.jl")
include("latent_models/models/HierarchicalNormal.jl")
include("latent_models/models/RandomWalk.jl")
include("latent_models/models/AR.jl")
include("latent_models/models/MA.jl")
include("latent_models/models/HilbertSpaceGP.jl")
include("latent_models/models/ExactGP.jl")
include("latent_models/models/Intercept.jl")
include("latent_models/models/Null.jl")
include("latent_models/modifiers/DiffLatentModel.jl")
include("latent_models/modifiers/TransformLatentModel.jl")
include("latent_models/modifiers/PrefixLatentModel.jl")
include("latent_models/modifiers/RecordExpectedLatent.jl")
include("latent_models/manipulators/CombineLatentModels.jl")
include("latent_models/manipulators/ConcatLatentModels.jl")
include("latent_models/manipulators/Hierarchy.jl")
include("latent_models/manipulators/Stratify.jl")
include("latent_models/manipulators/Replicate.jl")
include("latent_models/manipulators/broadcast/LatentModel.jl")
include("latent_models/manipulators/broadcast/rules.jl")
include("latent_models/manipulators/broadcast/helpers.jl")
include("latent_models/combinations/arma.jl")
include("latent_models/combinations/arima.jl")

# --- infection models ---
include("infection_models/DirectInfections.jl")
include("infection_models/ExpGrowthRate.jl")
include("infection_models/SeedingPath.jl")
include("infection_models/Renewal.jl")
# `utils.jl` defines the `R_to_r(::Renewal)` method, so it follows `Renewal`.
include("infection_models/utils.jl")
include("infection_models/CombineInfections.jl")

# --- ODE compartmental models ---
include("ode/SIRParams.jl")
include("ode/SEIRParams.jl")
include("ode/ODEProcess.jl")
include("ode/CatalystODEParams.jl")

# --- observation models ---
include("observation_models/ObservationErrorModels/methods.jl")
include("observation_models/ObservationErrorModels/PoissonError.jl")
include("observation_models/ObservationErrorModels/NegativeBinomialError.jl")
include("observation_models/ObservationErrorModels/NormalError.jl")
include("observation_models/ObservationErrorModels/BinomialError.jl")
include("observation_models/modifiers/LatentDelay.jl")
include("observation_models/modifiers/ascertainment/Ascertainment.jl")
include("observation_models/modifiers/ascertainment/helpers.jl")
include("observation_models/modifiers/Aggregate.jl")
include("observation_models/modifiers/RightTruncate.jl")
include("observation_models/modifiers/ReportTriangle.jl")
include("observation_models/modifiers/PrefixObservationModel.jl")
include("observation_models/modifiers/RecordExpectedObs.jl")
include("observation_models/modifiers/TransformObservationModel.jl")
include("observation_models/Split.jl")

# --- composition ---
include("compose.jl")

# Included after the composition so it can dispatch on `IDModel` as well as on
# every observation component.
include("observation_models/traversal.jl")

# What data an assembled model needs, read off the traversal seam above.
include("observation_models/requirements.jl")

# --- inference orchestration ---
include("inference/IDProblem.jl")
include("inference/forecast.jl")

end
