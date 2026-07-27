# Catalyst.jl-backed ODE infection models (opt-in package extension).
#
# Loads only when `Catalyst` and `ModelingToolkit` are present alongside
# `ComposableTuringIDModels`. The public, model-agnostic parameter component itself —
# `CatalystODEParams` — is defined in `src/` (so it is exported and `@ref`-able);
# this extension supplies only the parts that genuinely need the symbolic stack:
# the `ReactionSystem` constructor, the `as_turing_model` sampling method, and the
# `remake_ode_problem` hook.
#
# Nothing here is specialised to a particular compartmental model. You hand
# `CatalystODEParams` ANY Catalyst `@reaction_network` plus priors for its initial
# conditions and rate parameters; Catalyst + ModelingToolkit generate the ODE
# system and a symbolic Jacobian, kept consistent by construction. SIR, SEIR, a
# vaccinated class, a second strain, etc. are all just different reaction networks
# passed to the same type.
#
# The hand-coded `SIRParams` / `SEIRParams` path stays the zero-latency DEFAULT;
# this extension is purely additive and only pulls in the heavy symbolic stack
# when a user explicitly `using Catalyst`.
#
# Sampling and problem rebuilding are fully SYMBOLIC. `as_turing_model` returns
# `symbol => value` maps and `remake` places each value into the problem by name,
# so we never resolve or sort stored indices and never assume a positional layout
# (Catalyst sorts species/parameters internally, so that order is not the order
# you wrote the network in). The solution is likewise indexed symbolically with
# the network's own handles, e.g. `sol2infs = sol -> sol[rn.I, :]`.
#
# The MTK/Catalyst `ODEProblem` stores parameters as a structured `MTKParameters`
# object and runs an initialization system, so we rebuild with
# `build_initializeprob = false`, which bypasses that init system and is the
# reverse-mode-AD-safe path.

module ComposableTuringIDModelsCatalystExt

# `CatalystODEParams` (the public struct), `as_turing_model` and
# `remake_ode_problem` are extended below as `ComposableTuringIDModels.<name>`
# (qualified), so they are not imported as bare names.
using ComposableTuringIDModels: ComposableTuringIDModels, CatalystODEParams
using Catalyst: Catalyst, ReactionSystem
# `unknowns` / `parameters` (the system's symbolic species / rate handles) are
# owned by ModelingToolkit; we call them qualified as `ModelingToolkit.<name>`
# rather than importing the bare names (keeps the extension's explicit imports
# public and owner-correct for ExplicitImports).
using ModelingToolkit: ModelingToolkit
using OrdinaryDiffEq: ODEProblem, remake
using DynamicPPL: DynamicPPL, @model, NamedDist
using Distributions: Distribution

# A specification for one species or parameter: its symbolic handle (used to build
# the symbolic `symbol => value` remake map), the flat Symbol `name` that names the
# sampled Turing variable, and the `spec` giving its value — either a
# `Distribution` (sampled) or a plain `Real` (a fixed constant, NOT sampled, so
# fixed compartments / rates don't introduce discrete `Dirac` variables that break
# NUTS gradient checks).
struct SymbolSpec{S, D <: Union{Distribution, Real}}
    symbol::S
    name::Symbol
    spec::D
end

# Build a `SymbolSpec` from a `symbol => spec` pair, deriving the flat variable
# name from the symbol (`S(t)` -> `S`).
function _spec(sym, spec)
    v = Catalyst.value(sym)
    name = Symbol(replace(string(Symbol(v)), "(t)" => ""))
    return SymbolSpec(v, name, spec)
end

function ComposableTuringIDModels.CatalystODEParams(
        rn::ReactionSystem; tspan, u0_priors, p_priors)
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

    u0_specs = Tuple(_spec(sym, spec) for (sym, spec) in u0_priors)
    p_specs = Tuple(_spec(sym, spec) for (sym, spec) in p_priors)
    return CatalystODEParams(prob, u0_specs, p_specs)
end

# Sample every distribution-valued spec into a flat, symbol-named Turing variable
# (`β`, `S`, ...) — fixed `Real` specs are used as constants, NOT sampled — and
# return the initial state and parameters as SYMBOLIC `symbol => value` maps. The
# symbolic `remake` below places each value into the problem by name, so the
# stored (Catalyst-sorted) layout is never assumed and the ForwardDiff `Dual`
# values propagate through `remake` unchanged.
@model function ComposableTuringIDModels.as_turing_model(params::CatalystODEParams, n)
    u0 = Vector{Pair}(undef, length(params.u0_specs))
    p = Vector{Pair}(undef, length(params.p_specs))
    for (k, s) in enumerate(params.u0_specs)
        if s.spec isa Distribution
            x ~ NamedDist(s.spec, DynamicPPL.VarName{s.name}())
            u0[k] = s.symbol => x
        else
            u0[k] = s.symbol => s.spec
        end
    end
    for (k, s) in enumerate(params.p_specs)
        if s.spec isa Distribution
            x ~ NamedDist(s.spec, DynamicPPL.VarName{s.name}())
            p[k] = s.symbol => x
        else
            p[k] = s.symbol => s.spec
        end
    end
    return (u0, p)
end

# Catalyst/MTK `remake` hook: place the sampled values by symbolic name and bypass
# the initialization system. `build_initializeprob = false` is the
# reverse-mode-AD-safe path; the symbolic maps keep ordering out of the picture and
# let promoted `Dual` values flow straight through `remake`.
function ComposableTuringIDModels.remake_ode_problem(::CatalystODEParams, prob, u0, p)
    return remake(prob; u0 = u0, p = p, build_initializeprob = false)
end

end # module ComposableTuringIDModelsCatalystExt
