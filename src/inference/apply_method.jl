# `apply_method` driver: condition a model and run an inference method.

@doc raw"
Condition a model by fixing some parameters and conditioning on others, then run
an inference `method`.

# Arguments

  - `epiproblem`: the [`IDProblem`](@ref) (or a `DynamicPPL.Model`).
  - `method`: the inference method (a sampler or an [`IDMethod`](@ref)).
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
function apply_method(epiproblem::IDProblem, method::AbstractIDMethod, data;
        fix_parameters::NamedTuple = NamedTuple(),
        condition_parameters::NamedTuple = NamedTuple(), kwargs...)
    model = as_turing_model(epiproblem, data)
    cond_model = condition_model(model, fix_parameters, condition_parameters)
    return apply_method(cond_model, method, data; kwargs...)
end

# Apply a method to a model and wrap the solution as observables. Mirrors the
# upstream two- and three-argument `apply_method`: run the method (`_apply_method`,
# or the `IDMethod` pre-sampler→sampler chain) and return an
# [`IDObservables`](@ref) via [`generated_observables`](@ref).
function apply_method(model::DynamicPPL.Model, method::AbstractIDMethod, data;
        kwargs...)
    solution = _run_method(model, method; kwargs...)
    return generated_observables(model, data, solution)
end

function apply_method(model::DynamicPPL.Model, method::AbstractIDMethod; kwargs...)
    return apply_method(model, method, nothing; kwargs...)
end

# Run a method to its raw solution. An `IDMethod` threads its pre-sampler steps
# into the sampler; a bare method goes straight to its `_apply_method`.
function _run_method(model::DynamicPPL.Model, method::IDMethod; kwargs...)
    prev_result = nothing
    for pre_sampler in method.pre_sampler_steps
        prev_result = _apply_method(model, pre_sampler, prev_result; kwargs...)
    end
    return _apply_method(model, method.sampler, prev_result; kwargs...)
end

function _run_method(model::DynamicPPL.Model, method::AbstractIDMethod; kwargs...)
    _apply_method(model, method, nothing; kwargs...)
end
