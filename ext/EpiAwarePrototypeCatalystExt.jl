# Catalyst.jl-backed ODE infection models (opt-in package extension).
#
# Loads only when `Catalyst` and `ModelingToolkit` are present alongside
# `EpiAwarePrototype`. It adds a *declarative*, model-agnostic ODE parameter
# component: you hand it ANY Catalyst `@reaction_network` together with priors
# for its initial conditions and rate parameters, and Catalyst + ModelingToolkit
# generate the ODE system and a symbolic Jacobian for you, kept consistent by
# construction. There is no hand-written vector field or Jacobian, and nothing in
# this file is specialised to a particular compartmental model — SIR, SEIR, a
# vaccinated class, a second strain, etc. are all just different reaction
# networks passed to the same `CatalystODEParams`.
#
# The hand-coded `SIRParams` / `SEIRParams` path stays the zero-latency DEFAULT;
# this extension is purely additive and only pulls in the heavy symbolic stack
# when a user explicitly `using Catalyst`.
#
# Two Catalyst gotchas shape the code below (see the package's Catalyst case
# study and issue #46), handled GENERICALLY so they hold for any network:
#
#   1. Catalyst SORTS species and parameters when it builds the problem, so the
#      problem's internal state/parameter layout is generally NOT the order you
#      wrote the network in. We therefore resolve every species' and parameter's
#      stored position once, via symbolic indexing (`variable_index` /
#      `parameter_index`), and assemble the sampled `u0` / `p` vectors into those
#      positions — never assuming a positional layout.
#   2. A `ModelingToolkit`/`Catalyst` `ODEProblem` stores parameters as a
#      structured `MTKParameters` object and runs an initialization system, so
#      the plain-vector `remake(prob; u0 = vec, p = vec)` that the hand-coded
#      problems accept throws inside init reconstruction. We rebuild with
#      `build_initializeprob = false`, which bypasses the init system and is the
#      reverse-mode-AD-safe path, passing plain vectors already arranged into the
#      stored layout by (1).

module EpiAwarePrototypeCatalystExt

# `as_turing_model` and `remake_ode_problem` are extended below as
# `EpiAwarePrototype.<name>` (qualified), so they are not imported as bare names.
using EpiAwarePrototype: EpiAwarePrototype, AbstractLatentModel
using Catalyst: Catalyst, ReactionSystem
# `variable_index` / `parameter_index` (state/parameter index lookup) and
# `unknowns` / `parameters` (the system's symbolic species / rate handles) are
# owned by SymbolicIndexingInterface / ModelingToolkit; we call them qualified as
# `ModelingToolkit.<name>` rather than importing the bare names (keeps the
# extension's explicit imports public and owner-correct for ExplicitImports).
using ModelingToolkit: ModelingToolkit
using OrdinaryDiffEq: ODEProblem, remake
using DynamicPPL: DynamicPPL, @model, NamedDist
using Distributions: Distribution

# A specification for one species or parameter: its symbolic handle and the
# `spec` that gives its (initial) value — either a `Distribution` (sampled, with a
# flat Turing variable named after the symbol) or a plain `Real` (a fixed
# constant, NOT sampled — so fixed compartments / rates don't introduce discrete
# `Dirac` variables that break NUTS gradient checks). `name` is the flat Symbol;
# `index` is the symbol's position in the problem's stored (sorted)
# state / tunable-parameter vector.
struct SymbolSpec{S, D <: Union{Distribution, Real}}
    symbol::S
    name::Symbol
    index::Int
    spec::D
end

@doc raw"
Declarative, model-agnostic ODE parameter component built from **any** `Catalyst`
`ReactionSystem`, usable as the latent component of an
[`ODEProcess`](@ref EpiAwarePrototype.ODEProcess) in place of the hand-coded
[`SIRParams`](@ref EpiAwarePrototype.SIRParams) /
[`SEIRParams`](@ref EpiAwarePrototype.SEIRParams).

You declare a reaction network and give priors for its initial conditions and
rate parameters; Catalyst + ModelingToolkit generate the ODE system **and a
symbolic Jacobian** (`jac = true`), so there is no hand-written vector field or
Jacobian to keep in sync, and nothing here is specialised to a particular
compartmental model. Construct it for an SIR network, an SEIR network, or any
other network the same way.

This type lives in the optional `Catalyst` extension and is available once
`Catalyst` is loaded (`using EpiAwarePrototype, Catalyst`). The heavy symbolic
stack stays out of the default install; the hand-coded models remain the
zero-latency default.

Its `as_turing_model` samples each supplied prior into a flat Turing variable
named after the species / parameter symbol (e.g. `β`, `S`), then assembles the
initial-state vector `u0` and parameter vector `p` into the problem's stored
(Catalyst-sorted) layout, returning `(u0, p)` exactly like the hand-coded
parameter models — so it drops straight into [`ODEProcess`](@ref
EpiAwarePrototype.ODEProcess) and the prefix-off `to_submodel` flow.

# Arguments

  - `rn`: the `Catalyst` `ReactionSystem` (e.g. from `@reaction_network`).

# Keyword Arguments

  - `tspan`: the ODE solution time span.
  - `u0_priors`: the initial conditions, as symbolic-handle ⇒ spec pairs
    (`[rn.S => Beta(...), rn.R => 0.0, ...]`). Each spec is either a
    `Distribution` (sampled, as a flat variable named after the species) or a
    plain `Real` (a fixed initial value, not sampled). Every species of `rn` must
    appear.
  - `p_priors`: the rate parameters, as symbolic-handle ⇒ spec pairs
    (`[rn.β => LogNormal(...), ...]`), each a `Distribution` (sampled) or a fixed
    `Real`. Every parameter of `rn` must appear.

# Fields

  - `prob`: the `ODEProblem` built from `rn` (auto symbolic Jacobian).
  - `u0_specs`: per-species specs carrying each species' stored state index.
  - `p_specs`: per-parameter specs carrying each parameter's stored slot.
  - `species_index`: a `name ⇒ stored-state-index` map, exposed so a `sol2infs`
    link can pull a compartment out of the solution by name without hard-coding
    Catalyst's ordering (e.g. `sol[params.species_index[:I], :]`).

# Examples
```julia
using EpiAwarePrototype, Catalyst, OrdinaryDiffEq, Distributions
sir = @reaction_network begin
    β, S + I --> 2I
    γ, I --> R
end
params = CatalystODEParams(sir;
    tspan = (0.0, 30.0),
    u0_priors = [sir.S => Beta(99, 1), sir.I => Beta(1, 99), sir.R => 0.0],
    p_priors = [sir.β => LogNormal(log(0.3), 0.1), sir.γ => LogNormal(log(0.1), 0.1)])
rand(as_turing_model(params, nothing))
```
"
struct CatalystODEParams{P, U, R} <: AbstractLatentModel
    "The `ODEProblem` built from the reaction network (auto symbolic Jacobian)."
    prob::P
    "Per-species specs (symbol, flat name, stored state index, prior-or-fixed)."
    u0_specs::U
    "Per-parameter specs (symbol, flat name, stored tunable slot, prior-or-fixed)."
    p_specs::R
    "Map from species name to its stored state index, for `sol2infs` links."
    species_index::Dict{Symbol, Int}
end

# Build the per-symbol specs, resolving each symbol's stored index via symbolic
# indexing, and return them SORTED into the problem's stored layout (ascending
# index). `index_of` differs for species (state index) vs parameters
# (tunable-parameter slot), so it is passed in. Each `spec` is a `Distribution`
# (sampled) or a `Real` (fixed). Sorting lets `as_turing_model` build the `u0`/`p`
# vectors by a single ordered comprehension, so their element type promotes
# concretely (e.g. to `ForwardDiff.Dual`) instead of staying an abstract
# `Vector{Real}` — the abstract eltype breaks the stiff solver's linear-solve.
function _specs(prob, pairs, index_of)
    specs = SymbolSpec[]
    for (sym, spec) in pairs
        v = Catalyst.value(sym)
        name = Symbol(replace(string(Symbol(v)), "(t)" => ""))  # `S(t)` -> `S`
        push!(specs, SymbolSpec(v, name, index_of(prob, v), spec))
    end
    sort!(specs; by = s -> s.index)
    return Tuple(specs)
end

function CatalystODEParams(rn::ReactionSystem; tspan, u0_priors, p_priors)
    species = ModelingToolkit.unknowns(rn)
    rates = ModelingToolkit.parameters(rn)
    length(u0_priors) == length(species) || throw(ArgumentError(
        "u0_priors must give a prior for every species ($(length(species)) of them)"))
    length(p_priors) == length(rates) || throw(ArgumentError(
        "p_priors must give a prior for every parameter ($(length(rates)) of them)"))

    # Build the problem with placeholder symbolic maps. `jac = true` makes
    # ModelingToolkit emit the symbolic Jacobian the stiff/auto solver wants.
    u0_map = [Catalyst.value(sym) => 0.0 for (sym, _) in u0_priors]
    p_map = [Catalyst.value(sym) => 0.0 for (sym, _) in p_priors]
    prob = ODEProblem(rn, u0_map, tspan, p_map; jac = true)

    u0_specs = _specs(prob, u0_priors, ModelingToolkit.variable_index)
    p_specs = _specs(prob, p_priors, (p, s) -> ModelingToolkit.parameter_index(p, s).idx)
    species_index = Dict(s.name => s.index for s in u0_specs)
    return CatalystODEParams(prob, u0_specs, p_specs, species_index)
end

# Sample every distribution-valued spec into a flat, symbol-named Turing variable
# (`β`, `S`, ...) — fixed `Real` specs are used as constants, NOT sampled — and
# assemble the initial-state vector `u0` and parameter vector `p` in the problem's
# stored (Catalyst-sorted) layout, returning `(u0, p)`. The plain-vector remake
# below consumes them in exactly that layout.
#
# The specs are pre-sorted into stored-index order (see `_specs`), so we sample /
# read them in order, collect the scalars, and `reduce(vcat, …)` to a vector. The
# `vcat` promotes the element type concretely to whatever the sampler feeds in
# (e.g. a `ForwardDiff.Dual{…,Float64}`); building an abstract `Vector{Real}`
# instead would propagate an abstract eltype into the ODE solve and break the
# stiff solver's linear-solve cache.
@model function EpiAwarePrototype.as_turing_model(params::CatalystODEParams, n)
    u0_list = Vector{Any}(undef, length(params.u0_specs))
    p_list = Vector{Any}(undef, length(params.p_specs))
    for (k, s) in enumerate(params.u0_specs)
        if s.spec isa Distribution
            x ~ NamedDist(s.spec, DynamicPPL.VarName{s.name}())
            u0_list[k] = x
        else
            u0_list[k] = s.spec
        end
    end
    for (k, s) in enumerate(params.p_specs)
        if s.spec isa Distribution
            x ~ NamedDist(s.spec, DynamicPPL.VarName{s.name}())
            p_list[k] = x
        else
            p_list[k] = s.spec
        end
    end
    u0 = reduce(vcat, u0_list)
    p = reduce(vcat, p_list)
    return (u0, p)
end

# Catalyst/MTK `remake` hook: bypass the initialization system with plain
# vectors. `build_initializeprob = false` is the reverse-mode-AD-safe path; the
# vectors are already arranged into the problem's stored layout above.
function EpiAwarePrototype.remake_ode_problem(::CatalystODEParams, prob, u0, p)
    return remake(prob; u0 = u0, p = p, build_initializeprob = false)
end

end # module EpiAwarePrototypeCatalystExt
