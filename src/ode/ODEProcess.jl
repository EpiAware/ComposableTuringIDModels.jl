# `ODEProcess` infection model: solve a compartmental ODE and map the solution
# to a latent-infection series.

@doc raw"
Rebuild an `ODEProcess`'s `ODEProblem` with freshly sampled state `u0` and
parameters `p`.

This is the single seam through which `ODEProcess` re-instantiates its problem on
every sample. The default applies the plain-vector `remake(prob; u0, p)` that the
hand-coded [`SIRParams`](@ref) / [`SEIRParams`](@ref) problems expect.

Parameter models whose problem stores parameters differently — e.g. a
`ModelingToolkit`/`Catalyst`-built problem carrying a structured `MTKParameters`
object plus an initialization system — specialise this on their own type to inject
the remake form they need. The Catalyst extension uses the plain-vector
`build_initializeprob = false` path to bypass the init system (and stay reverse-mode
differentiable). Dispatching on the parameter model keeps `ODEProcess` itself
backend-agnostic.

# Arguments

  - `params`: the ODE parameter model (the `ODEProcess`'s `params` field).
  - `prob`: the template `ODEProblem` to rebuild.
  - `u0`: the sampled initial state vector.
  - `p`: the sampled parameter vector.
"
remake_ode_problem(params::AbstractLatentModel, prob, u0, p) = remake(prob; u0 = u0, p = p)

@doc raw"
An infection process defined by solving an ODE.

`ODEProcess` combines a parameter struct (`params`, e.g. [`SIRParams`](@ref) or
[`SEIRParams`](@ref), whose `as_turing_model` samples `(u0, p)`) with a `solver`,
extra `solver_options`, and a `sol2infs` link mapping the ODE solution to a
latent-infection series. The compartmental dynamics are fully determined by the
sampled ODE parameters, so the model carries no separate latent ``R_t`` process:
its `as_turing_model` samples the parameters, solves the ODE, and returns
`(; I_t, Z_t)` with `Z_t = nothing` (no exposable latent path).

# Arguments

  - `epi_model`: the [`ODEProcess`](@ref).
  - `n`: the requested series length; passed through to the parameter model
    (the ODE dimension is fixed, so `n` is otherwise unused — `nothing` is also
    accepted).

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
    _prob = remake_ode_problem(epi_model.params, prob, u0, p)
    sol = solve(_prob, solver; solver_options...)
    return sol
end

@model function as_turing_model(epi_model::ODEProcess, n)
    n_steps = isnothing(n) ? 0 : n
    sol ~ to_submodel(_generate_ode_solution(epi_model, n_steps), false)
    I_t = epi_model.sol2infs(sol)
    return (; I_t, Z_t = nothing)
end
