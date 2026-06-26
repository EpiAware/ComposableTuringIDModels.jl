# API reference

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

## Infection models

```@docs
EpiData
DirectInfections
```

## Observation models

```@docs
PoissonError
NegativeBinomialError
LatentDelay
observation_error
generate_observation_error_priors
```

## Composition

```@docs
EpiAwareModel
```

## Utilities and distributions

```@docs
accumulate_scan
get_state
HalfNormal
SafePoisson
SafeNegativeBinomial
NegativeBinomialMeanClust
censored_pmf
condition_model
```
