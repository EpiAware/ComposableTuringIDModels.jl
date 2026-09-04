# Catalyst.jl-backed ODE infection models (opt-in package extension).
#
# Loads only when `Catalyst` and `ModelingToolkit` are present alongside
# `ComposableTuringIDModels`. The public, model-agnostic parameter component itself —
# `CatalystODEParams` — is defined in `src/` (so it is exported and `@ref`-able);
# this extension supplies only the parts that genuinely need the symbolic stack:
# the `ReactionSystem` constructor, the `as_turing_model` sampling method, and the
# `remake_ode_problem` / `ode_solution_view` hooks.
#
# Nothing here is specialised to a particular compartmental model. You hand
# `CatalystODEParams` ANY Catalyst `@reaction_network` plus priors for its initial
# conditions and rate parameters; Catalyst + ModelingToolkit generate the ODE
# system and a symbolic Jacobian, kept consistent by construction. SIR, SEIR, a
# vaccinated class, a second strain, etc. are all just different reaction networks
# passed to the same type.
#
# The hand-coded `SIRParams` path stays the zero-latency DEFAULT for SIR;
# this extension is purely additive and only pulls in the heavy symbolic stack
# when a user explicitly `using Catalyst`.
#
# Three things here exist so that reverse-mode AD works through the solve.
# Catalyst sorts species and parameters internally, so all three lean on one
# ordering resolved once at construction rather than on a symbolic lookup per
# sample.
#
#   1. Sampling returns PLAIN VECTORS in the problem's own order, and `remake`
#      places them with `build_initializeprob = false`. A symbolic
#      `symbol => value` remake instead reconstructs the problem's `MTKParameters`
#      and runs its initialization system, which drags ModelingToolkit's symbolic
#      caches onto the AD tape.
#   2. The problem is re-specialised `FullSpecialize`. An MTK-built problem
#      defaults to `AutoDespecialize`, which wraps the generated Jacobian in a
#      `Float64`-typed `FunctionWrapper`; the stiff half of `AutoVern7(Rodas5P())`
#      then calls it with `ForwardDiff.Dual` entries and no wrapper matches.
#   3. A solution is handed to `sol2infs` through a view that resolves the
#      network's own symbolic handles positionally, so `sol[rn.I, :]` keeps
#      working. A reverse-mode solve returns a solution without the symbolic
#      index provider that would otherwise serve that lookup.

module ComposableTuringIDModelsCatalystExt

# `CatalystODEParams` (the public struct), `as_turing_model`, `remake_ode_problem`
# and `ode_solution_view` are extended below as `ComposableTuringIDModels.<name>`
# (qualified), so they are not imported as bare names.
using ComposableTuringIDModels: ComposableTuringIDModels, CatalystODEParams
using Catalyst: Catalyst, ReactionSystem
# `unknowns` / `parameters` (a system's symbolic species / rate handles) and
# `variable_symbols` (the compiled problem's own state ordering) are owned by
# ModelingToolkit; we call them qualified as `ModelingToolkit.<name>` rather than
# importing the bare names (keeps the extension's explicit imports public and
# owner-correct for ExplicitImports).
using ModelingToolkit: ModelingToolkit
using OrdinaryDiffEq: ODEProblem, ODEFunction, remake
using SciMLBase: FullSpecialize
using DynamicPPL: DynamicPPL, @model, NamedDist
using Distributions: Distribution

# A specification for one species or parameter: its symbolic handle, the flat
# Symbol `name` that names the sampled Turing variable, and the `spec` giving its
# value — either a `Distribution` (sampled) or a plain `Real` (a fixed constant,
# NOT sampled, so fixed compartments / rates don't introduce discrete `Dirac`
# variables that break NUTS gradient checks).
struct SymbolSpec{S, D <: Union{Distribution, Real}}
    symbol::S
    name::Symbol
    spec::D
end

# The flat variable name a symbolic handle carries (`S(t)` -> `:S`).
_flat_name(sym) = Symbol(replace(string(Symbol(Catalyst.value(sym))), "(t)" => ""))

_spec(sym, spec) = SymbolSpec(Catalyst.value(sym), _flat_name(sym), spec)

# For each slot of the compiled problem's vector, the index of the spec that
# fills it. `symbols` is the problem's own ordering, `specs` the user's.
function _slot_sources(specs, symbols, what)
    spec_names = [s.name for s in specs]
    return Tuple(
        map(symbols) do sym
            k = findfirst(==(_flat_name(sym)), spec_names)
            isnothing(k) && throw(
                ArgumentError(
                    "the compiled problem has a $what named $(_flat_name(sym)) " *
                        "that no supplied prior covers"
                )
            )
            return k
        end
    )
end

function ComposableTuringIDModels.CatalystODEParams(
        rn::ReactionSystem; tspan, u0_priors, p_priors
    )
    species = ModelingToolkit.unknowns(rn)
    rates = ModelingToolkit.parameters(rn)
    length(u0_priors) == length(species) || throw(
        ArgumentError(
            "u0_priors must give a prior for every species ($(length(species)) of them)"
        )
    )
    length(p_priors) == length(rates) || throw(
        ArgumentError(
            "p_priors must give a prior for every parameter ($(length(rates)) of them)"
        )
    )

    # Build the problem with placeholder symbolic maps. `jac = true` makes
    # ModelingToolkit emit the symbolic Jacobian the stiff/auto solver wants.
    u0_map = [Catalyst.value(sym) => 0.0 for (sym, _) in u0_priors]
    p_map = [Catalyst.value(sym) => 0.0 for (sym, _) in p_priors]
    prob = ODEProblem(rn, u0_map, tspan, p_map; jac = true)
    # `ODEProblem{true, FullSpecialize}(rn, ...)` does not work here: the
    # constructor reads the `ReactionSystem` as a vector field. Re-specialising
    # the built problem's `ODEFunction` does.
    prob = remake(prob; f = ODEFunction{true, FullSpecialize}(prob.f))

    u0_specs = Tuple(_spec(sym, spec) for (sym, spec) in u0_priors)
    p_specs = Tuple(_spec(sym, spec) for (sym, spec) in p_priors)
    # The compiled problem's own ordering, which is not the order the priors were
    # written in. Rates come from the built system rather than from
    # `parameter_symbols(prob)`, which also lists the `Initial(...)` entries the
    # problem carries for its initialization system and which are no part of the
    # plain parameter vector.
    state_symbols = ModelingToolkit.variable_symbols(prob)
    rate_symbols = ModelingToolkit.parameters(prob.f.sys)
    layout = (
        u0 = _slot_sources(u0_specs, state_symbols, "species"),
        p = _slot_sources(p_specs, rate_symbols, "parameter"),
        slots = NamedTuple{Tuple(_flat_name(sym) for sym in state_symbols)}(
            Tuple(eachindex(state_symbols))
        ),
    )
    return CatalystODEParams(prob, u0_specs, p_specs, layout)
end

# Reorder sampled values from the order the user wrote the priors in into the
# order the compiled problem stores them. Fixed `Real` specs sit alongside
# sampled `Dual`s or tracked values, so `promote` is what keeps the result a
# concretely typed vector the solver can work in rather than a `Vector{Real}`.
_in_problem_order(values, sources) = collect(promote(map(k -> values[k], sources)...))

# Sample every distribution-valued spec into a flat, symbol-named Turing variable
# (`β`, `S`, ...) — fixed `Real` specs are used as constants, NOT sampled — and
# return the initial state and parameters as plain vectors in the compiled
# problem's own order.
@model function ComposableTuringIDModels.as_turing_model(
        params::CatalystODEParams, n::Union{Int, Nothing}
    )
    u0_values = Vector{Any}(undef, length(params.u0_specs))
    p_values = Vector{Any}(undef, length(params.p_specs))
    for (k, s) in enumerate(params.u0_specs)
        if s.spec isa Distribution
            x ~ NamedDist(s.spec, DynamicPPL.VarName{s.name}())
            u0_values[k] = x
        else
            u0_values[k] = s.spec
        end
    end
    for (k, s) in enumerate(params.p_specs)
        if s.spec isa Distribution
            x ~ NamedDist(s.spec, DynamicPPL.VarName{s.name}())
            p_values[k] = x
        else
            p_values[k] = s.spec
        end
    end
    u0 = _in_problem_order(u0_values, params.layout.u0)
    p = _in_problem_order(p_values, params.layout.p)
    return (u0, p)
end

# Catalyst/MTK `remake` hook: place the sampled vectors and bypass the
# initialization system.
function ComposableTuringIDModels.remake_ode_problem(::CatalystODEParams, prob, u0, p)
    return remake(prob; u0 = u0, p = p, build_initializeprob = false)
end

# A solution presented to `sol2infs`, with the network's symbolic handles
# resolved through the ordering stored at construction. Fields are read with
# `getfield` so `getproperty` can forward everything else (`sol.t`, `sol.u`) to
# the solution itself.
struct SymbolIndexedSolution{S, L}
    sol::S
    slots::L
end

Base.getproperty(v::SymbolIndexedSolution, name::Symbol) =
    getproperty(getfield(v, :sol), name)

# Anything Julia can already index an array with passes through; anything else is
# taken for a symbolic handle and resolved by name.
_resolve(::SymbolIndexedSolution, i::Union{Integer, Colon, AbstractRange}) = i
_resolve(::SymbolIndexedSolution, i::AbstractVector{<:Integer}) = i
_resolve(v::SymbolIndexedSolution, i) = getproperty(getfield(v, :slots), _flat_name(i))

function Base.getindex(v::SymbolIndexedSolution, i, rest...)
    return getindex(getfield(v, :sol), _resolve(v, i), rest...)
end

function ComposableTuringIDModels.ode_solution_view(params::CatalystODEParams, sol)
    return SymbolIndexedSolution(sol, params.layout.slots)
end

end # module ComposableTuringIDModelsCatalystExt
