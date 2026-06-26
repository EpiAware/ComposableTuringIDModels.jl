# PACKAGE-OWNED — scaffold writes this once and never overwrites it.
#
# Minimal AD-fixture registry implementing the EpiAwarePackageTools `ADRegistry`
# contract: scenarios (each with a ForwardDiff reference), a backend list, and
# broken/skip bookkeeping. The shared harness (driven from `test/ad/setup.jl`)
# consumes these. Replace the placeholder scenario with the package's own
# differentiable log densities and add the backends it supports.
module ADFixtures

using ADTypes: AutoForwardDiff
using DifferentiationInterface: DifferentiationInterface, Constant
import DifferentiationInterfaceTest as DIT
import ForwardDiff
using EpiAwarePrototype

export scenarios, backends, broken_scenario_names,
       backend_broken_scenarios, backend_skip_scenarios

# ForwardDiff reference gradient for a scenario function.
function _reference(f, θ, contexts)
    return DifferentiationInterface.gradient(
        f, AutoForwardDiff(), θ, contexts...)
end

"""
    scenarios(; with_reference = false, category = :marginal)

The AD gradient scenarios. Each is a `DIT.Scenario{:gradient, :out}` whose
`res1` carries a ForwardDiff reference when `with_reference = true`. Replace the
placeholder with the package's own differentiable log densities; group them by
`category` if the package distinguishes scenario groups (e.g. `:marginal` vs
`:latent`).
"""
function scenarios(; with_reference::Bool = false, category::Symbol = :marginal)
    out = DIT.Scenario{:gradient, :out}[]
    # Placeholder: a plain differentiable function so the harness runs out of
    # the box. Swap for a real log density, e.g.
    #   f = (θ, obs) -> sum(x -> logpdf(SomeDist(θ...), x), obs)
    θ = [1.0, 2.0]
    f = θ -> sum(abs2, θ)
    push!(out,
        DIT.Scenario{:gradient, :out}(f, θ; name = "placeholder sum_squares",
            res1 = with_reference ? _reference(f, θ, ()) : nothing))
    return out
end

"""
    backends()

The AD backends to test, as `(; name, backend)` named tuples. Add the package's
supported backends (ReverseDiff, Mooncake, Enzyme, ...).
"""
backends() = [(name = "ForwardDiff", backend = AutoForwardDiff())]

"Scenario names broken on every backend."
broken_scenario_names() = String[]

"Per-backend broken scenario names (`Dict{String, Set{String}}`)."
backend_broken_scenarios() = Dict{String, Set{String}}()

"Per-backend scenario names too unstable to run at all."
backend_skip_scenarios() = Dict{String, Set{String}}()

end # module ADFixtures
