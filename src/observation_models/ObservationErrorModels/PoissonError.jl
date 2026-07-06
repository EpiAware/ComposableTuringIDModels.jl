# Poisson observation-error model.

@doc raw"
A Poisson observation-error model.

# Examples
```jldoctest PoissonError; output = false
using ComposableTuringIDModels
poi = PoissonError()
mdl = as_turing_model(poi, missing, fill(10, 10))
rand(mdl)
nothing
# output
```
"
struct PoissonError <: AbstractObservationErrorModel end

observation_error(::PoissonError, Y_t) = SafePoisson(Y_t)
