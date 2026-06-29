# `ODEProcess` infection model: solve a compartmental ODE and map the solution
# to a latent-infection series.

@doc raw"
An infection process defined by solving an ODE.

`ODEProcess` combines a parameter struct (`params`, e.g. [`SIRParams`](@ref) or
[`SEIRParams`](@ref), whose `as_turing_model` samples `(u0, p)`) with a `solver`,
extra `solver_options`, and a `sol2infs` link mapping the ODE solution to a
latent-infection series. Its `as_turing_model` samples the parameters, solves the
ODE, and returns the mapped infections.

# Arguments

  - `epi_model`: the [`ODEProcess`](@ref).
  - `Z_t`: an optional latent path (may be `nothing` for self-contained ODE
    models); only its length is used.

# Examples
```@example ODEProcess
using EpiAwarePrototype, OrdinaryDiffEq, Distributions, LogExpFunctions
sirparams = SIRParams(
    tspan = (0.0, 100.0),
    infectiousness = LogNormal(log(0.3), 0.05),
    recovery_rate = LogNormal(log(0.1), 0.05),
    initial_prop_infected = Beta(1, 99))
N = 1000.0
sir_process = ODEProcess(
    params = sirparams,
    sol2infs = sol -> softplus.(N .* sol[2, :]),
    solver_options = Dict(:saveat => 1.0))
as_turing_model(sir_process, nothing)()
```

## Fields

  - `params`: the ODE parameter model (an [`AbstractLatentModel`](@ref), e.g.
    [`SIRParams`](@ref) / [`SEIRParams`](@ref), whose `as_turing_model` samples
    `(u0, p)`).
  - `solver`: the ODE solver (default `AutoVern7(Rodas5P())`).
  - `sol2infs`: link mapping the ODE solution to an infection series.
  - `solver_options`: extra options passed to `solve` (a `Dict` or `NamedTuple`).
"
@kwdef struct ODEProcess{P <: AbstractLatentModel, S, F <: Function,
    D <: Union{Dict, NamedTuple}} <: AbstractInfectionModel
    "The ODE parameter model."
    params::P
    "The ODE solver."
    solver::S = AutoVern7(Rodas5P())
    "Link mapping the ODE solution to an infection series."
    sol2infs::F
    "Extra options passed to `solve`."
    solver_options::D = Dict(:saveat => 1.0)
end

# Sample the ODE parameters and solve, returning the solution object.
@model function _generate_ode_solution(epi_model::ODEProcess, n)
    prob = epi_model.params.prob
    solver = epi_model.solver
    solver_options = epi_model.solver_options
    params ~ to_submodel(as_turing_model(epi_model.params, n), false)
    u0, p = params
    _prob = remake(prob; u0 = u0, p = p)
    sol = solve(_prob, solver; solver_options...)
    return sol
end

@model function as_turing_model(epi_model::ODEProcess, Z_t)
    n = isnothing(Z_t) ? 0 : size(Z_t, 1)
    sol ~ to_submodel(_generate_ode_solution(epi_model, n), false)
    return epi_model.sol2infs(sol)
end
