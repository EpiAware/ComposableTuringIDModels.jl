# API reference

This page documents the full surface of `EpiAwarePrototype`: the exported public
components first, then an [Internal types and helpers](@ref) section covering the
unexported building blocks the public components are built from (and that the
public docstrings cross-reference).

## Module

```@docs
EpiAwarePrototype
```

## Core architecture

```@docs
as_turing_model
AbstractEpiAwareModel
```

## Latent models

```@docs
IID
HierarchicalNormal
RandomWalk
AR
MA
Intercept
FixedIntercept
Null
DiffLatentModel
```

## Latent modifiers, manipulators, and combinations

```@docs
TransformLatentModel
PrefixLatentModel
RecordExpectedLatent
CombineLatentModels
ConcatLatentModels
BroadcastLatentModel
RepeatEach
RepeatBlock
broadcast_rule
broadcast_n
broadcast_dayofweek
broadcast_weekly
equal_dimensions
arma
arima
```

## Infection models

```@docs
EpiData
DirectInfections
ExpGrowthRate
Renewal
R_to_r
r_to_R
expected_Rt
```

## ODE compartmental models

```@docs
SIRParams
SEIRParams
ODEProcess
```

## Observation models

```@docs
PoissonError
NegativeBinomialError
LatentDelay
observation_error
generate_observation_error_priors
```

## Observation modifiers and manipulators

```@docs
Ascertainment
ascertainment_dayofweek
Aggregate
PrefixObservationModel
RecordExpectedObs
TransformObservationModel
StackObservationModels
```

## Composition

```@docs
EpiAwareModel
```

## Inference orchestration

```@docs
EpiProblem
EpiMethod
NUTSampler
ManyPathfinder
DirectSample
manypathfinder
apply_method
EpiAwareObservables
generated_observables
spread_draws
get_param_array
```

## Utilities and distributions

```@docs
accumulate_scan
get_state
HalfNormal
SafePoisson
SafeNegativeBinomial
NegativeBinomialMeanClust
condition_model
```

## Internal types and helpers

These are unexported. They are documented here because the public docstrings
cross-reference them and they describe how the public components are built; they
are not part of the stable public API.

### Accumulation and scan steps

```@docs
EpiAwarePrototype.AbstractAccumulationStep
EpiAwarePrototype.RWStep
EpiAwarePrototype.ARStep
EpiAwarePrototype.MAStep
```

### Renewal steps

```@docs
EpiAwarePrototype.AbstractConstantRenewalStep
EpiAwarePrototype.ConstantRenewalStep
EpiAwarePrototype.ConstantRenewalWithPopulationStep
EpiAwarePrototype.neg_MGF
```

### Broadcasting

```@docs
EpiAwarePrototype.AbstractBroadcastRule
```

### Observation internals

```@docs
EpiAwarePrototype.AbstractObservationErrorModel
EpiAwarePrototype.LDStep
```

### Inference method supertypes

```@docs
EpiAwarePrototype.AbstractEpiMethod
EpiAwarePrototype.AbstractEpiSamplingMethod
EpiAwarePrototype.AbstractEpiOptMethod
```

### Distribution and utility helpers

```@docs
EpiAwarePrototype.SafeIntValued
EpiAwarePrototype._expand_dist
```
