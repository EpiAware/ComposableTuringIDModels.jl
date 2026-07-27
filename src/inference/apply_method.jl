# `apply_method` driver: condition a model and run an inference method.

@doc raw"
Condition a model by fixing some parameters and conditioning on others, then run
an inference `method`.

# Arguments

  - `idproblem`: the [`IDProblem`](@ref) (or a `DynamicPPL.Model`).
  - `method`: the inference method (a sampler, e.g. [`NUTSampler`](@ref)).
  - `data`: the data to condition on (with a `y_t` field).

# Keyword Arguments

  - `fix_parameters`: a `NamedTuple` of parameters to fix.
  - `condition_parameters`: a `NamedTuple` of parameters to condition on.
  - `kwargs...`: forwarded to the inference method.

# Examples
```@example apply_method
using ComposableTuringIDModels, Distributions
problem = IDProblem(
    infection = DirectInfections(; Z = RandomWalk(), initialisation = Normal()),
    observation_model = PoissonError(),
    tspan = (1, 20))
y = rand(as_turing_model(problem, (; y_t = missing)))
nothing
```
"
function apply_method(idproblem::IDProblem, method::AbstractIDMethod, data;
        fix_parameters::NamedTuple = NamedTuple(),
        condition_parameters::NamedTuple = NamedTuple(), kwargs...)
    model = as_turing_model(idproblem, data)
    cond_model = condition_model(model, fix_parameters, condition_parameters)
    return apply_method(cond_model, method, data; kwargs...)
end

# Apply a method to a model and wrap the solution as observables. Mirrors the
# upstream two- and three-argument `apply_method`: run the method
# (`_apply_method`) and return an [`IDObservables`](@ref) via
# [`generated_observables`](@ref).
function apply_method(model::DynamicPPL.Model, method::AbstractIDMethod, data;
        kwargs...)
    solution = _run_method(model, method; kwargs...)
    return generated_observables(model, data, solution)
end

function apply_method(model::DynamicPPL.Model, method::AbstractIDMethod; kwargs...)
    return apply_method(model, method, nothing; kwargs...)
end

# Run a method to its raw solution.
function _run_method(model::DynamicPPL.Model, method::AbstractIDMethod; kwargs...)
    _apply_method(model, method, nothing; kwargs...)
end
