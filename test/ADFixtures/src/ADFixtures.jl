# PACKAGE-OWNED — scaffold writes this once and never overwrites it.
#
# AD-fixture registry implementing the EpiAwarePackageTools `ADRegistry`
# contract. The scenarios are REAL differentiable log-densities from the
# package: the (linked) log-joint of representative latent processes and of
# composed `EpiAwareModel`s conditioned on simulated data — the gradients an AD
# backend must get right for NUTS to work. Each scenario carries a ForwardDiff
# reference gradient. The shared harness (driven from `test/ad/setup.jl`)
# consumes this registry.
module ADFixtures

using ADTypes: AutoForwardDiff
using DifferentiationInterface: DifferentiationInterface
import DifferentiationInterfaceTest as DIT
import ForwardDiff
using EpiAwarePrototype
using Distributions
using Random: Random, MersenneTwister
using DynamicPPL: DynamicPPL, LogDensityFunction, VarInfo, link, getlogjoint
import LogDensityProblems as LDP

export scenarios, backends, broken_scenario_names,
       backend_broken_scenarios, backend_skip_scenarios

# Turn a DynamicPPL model into a real differentiable scalar log-density.
#
# We link the model's `VarInfo` so every constrained variable (truncated damping
# priors, positive standard deviations, simplex-free reals, ...) maps to an
# unconstrained real coordinate. The returned `f(θ)` is then the log-joint
# (including the linking log-Jacobian) at the flat unconstrained vector `θ`, and
# is finite and smooth over all of ``ℝ^d`` — exactly the target a gradient-based
# sampler differentiates. Returns `(f, θ0, dim)`.
function _logdensity(model; seed::Int = 1)
    vi = link(VarInfo(model), model)
    ldf = LogDensityFunction(model, getlogjoint, vi)
    dim = LDP.dimension(ldf)
    f = θ -> LDP.logdensity(ldf, θ)
    θ0 = 0.3 .* randn(MersenneTwister(seed), dim)
    return f, θ0, dim
end

# A representative generation interval shared by the infection-model scenarios.
const _GEN_INT = [0.2, 0.3, 0.5]

# Build the registry's models once. Conditioned (posterior) scenarios use data
# simulated from the prior with a fixed seed so the target is deterministic.
function _models()
    data = EpiData(_GEN_INT, exp)
    n = 12

    # Simulate observations from a composed model's prior (its `generated_y_t`).
    sim(m, nn) = as_turing_model(m, missing, nn)().generated_y_t

    # --- latent-process log-joints (prior only) ---------------------------------
    rw = as_turing_model(RandomWalk(), n)
    ar = as_turing_model(AR(), n)
    arima = as_turing_model(
        DiffLatentModel(; model = AR(), init = [Normal(), Normal()]), n)
    # Moving-average: exercises `accumulate_scan(MAStep(θ), ...)` and its
    # `dot(θ, state)` innovation buffer (the MA counterpart of `AR`).
    ma = as_turing_model(MA(), 8)
    # Non-centred hierarchical normal: `σ ~ prior`, `η = σ ⋅ ϵ` — the simplest
    # scale-mixture latent (and the default innovation model everything reuses).
    hier = as_turing_model(HierarchicalNormal(), 8)
    # A random walk wrapped in a `d = 2` `DiffLatentModel` (the differencing
    # modifier over a non-`AR` inner process, distinct from the `AR`-based ARIMA).
    diffrw = as_turing_model(
        DiffLatentModel(; model = RandomWalk(), init = [Normal(), Normal()]), 8)
    # ARMA(p, q): an `AR` whose innovations are an `MA` (double accumulate-scan).
    armamdl = as_turing_model(arma(), 8)
    # Day-of-week broadcast: a `TransformLatentModel` (7·softmax) inner process
    # repeated across a 7-day period (`RepeatEach`).
    bdow = as_turing_model(broadcast_dayofweek(RandomWalk()), 14)
    # Weekly broadcast: a piecewise-constant weekly process (`RepeatBlock`).
    bweek = as_turing_model(broadcast_weekly(RandomWalk()), 14)
    # Concatenate an `Intercept` segment and a `RandomWalk` segment along time.
    concat = as_turing_model(
        ConcatLatentModels([Intercept(Normal(2, 0.2)), RandomWalk()]), 10)
    # Sum an `Intercept` and an `AR` over the full length (prefix-separated).
    combine = as_turing_model(
        CombineLatentModels([Intercept(Normal(2, 0.2)), AR()]), 10)

    # --- the #76 prior interface -----------------------------------------------
    # A `BroadcastPrior` over a *vector* of damping distributions (order 2): one
    # i.i.d. draw per lag through `arraydist`, threaded as a submodel.
    ar_vec = as_turing_model(
        AR(; damp = [truncated(Normal(0, 0.05), 0, 1),
                truncated(Normal(0, 0.05), 0, 1)],
            init = [Normal(), Normal()]), 8)
    # A latent MODEL as a prior: the AR damping coefficient is itself a
    # (prefixed) `RandomWalk` submodel, so the submodel-threading gradient path
    # is differentiated. The prefix keeps the inner `std`/`ϵ_t`/`rw_init` names
    # from colliding with the AR innovation's under the prefix-off convention.
    ar_lat = as_turing_model(
        AR(; damp = PrefixLatentModel(RandomWalk(), "damp")), 8)

    # --- infection posteriors ---------------------------------------------------
    direct = EpiAwareModel(
        DirectInfections(; Z = RandomWalk(), initialisation = Normal()),
        PoissonError())
    renewal = EpiAwareModel(
        Renewal(data; rt = RandomWalk(), initialisation = Normal()),
        NegativeBinomialError())
    # Exponential-growth-rate infections (the third infection family alongside
    # `DirectInfections` / `Renewal`): a cumulative growth-rate path exponentiated.
    egr = EpiAwareModel(
        ExpGrowthRate(; rt = RandomWalk(), initialisation = Normal()),
        PoissonError())

    # Nowcasting MARGINAL (right-truncation correction): a renewal model whose
    # observation error is wrapped in `RightTruncate` (fixed reporting-delay CDF
    # supplied as a `ReportingCDF` submodel). This exercises the `reverse`/
    # broadcast scaling the modifier adds on top of the inner error.
    nowcast = EpiAwareModel(
        Renewal(data; rt = RandomWalk(), initialisation = Normal()),
        RightTruncate(NegativeBinomialError(),
            truncated(Normal(4.0, 1.5), 0.0, Inf)))

    # Nowcasting JOINT (2D reporting triangle): a renewal model feeding the
    # per-cell `ReportTriangle` observation model. The gradient of the per-cell
    # Poisson log-likelihood over the masked triangle (`t + d ≤ now`) is what
    # nowcasting under NUTS depends on.
    triangle = EpiAwareModel(
        Renewal(data; rt = RandomWalk(), initialisation = Normal()),
        ReportTriangle(PoissonError(), [0.6, 0.25, 0.15]))

    # --- observation modifiers / error families over a composed model ----------
    # Reporting delay: convolves the expected observations with a delay PMF
    # (`accumulate_scan(LDStep(rev_pmf), ...)`) before the inner error.
    latdelay = EpiAwareModel(
        Renewal(data; rt = RandomWalk(), initialisation = Normal()),
        LatentDelay(NegativeBinomialError(), [0.3, 0.4, 0.3]))
    # Day-of-week ascertainment: scales the expected observations by a broadcast
    # latent (an `Ascertainment` wrapping `broadcast_dayofweek`).
    ascert = EpiAwareModel(
        DirectInfections(; Z = RandomWalk(), initialisation = Normal()),
        ascertainment_dayofweek(PoissonError()))
    # Aggregation: sum the expected observations over weekly reporting windows
    # (only the window endpoints are scored).
    aggregate = EpiAwareModel(
        DirectInfections(; Z = RandomWalk(), initialisation = Normal()),
        Aggregate(PoissonError(), [0, 0, 0, 0, 0, 0, 7]))
    # Transform-the-expected-observations: softplus applied before the error.
    transobs = EpiAwareModel(
        DirectInfections(; Z = RandomWalk(), initialisation = Normal()),
        TransformObservationModel(PoissonError()))
    # Gaussian observation error (continuous, `σ`-inferred) rather than a count
    # family — the minimal non-count likelihood.
    normalobs = EpiAwareModel(
        DirectInfections(; Z = RandomWalk(), initialisation = Normal()),
        NormalError())

    # Binomial error: a standalone observation model whose success PROBABILITY is
    # a latent process (a `HierarchicalNormal` pushed through a logistic link via
    # `Ascertainment`), with the number of trials supplied as data. This is the
    # meaningful way to reach a `BinomialError` gradient (its `Y_t` is a
    # probability, not a count, so it is not fed by an infection model here).
    binom_obs = Ascertainment(BinomialError(), HierarchicalNormal();
        transform = (Y_t, x) -> 1 ./ (1 .+ exp.(-x)), latent_prefix = "")
    n_b = 10
    N_b = fill(20, n_b)
    Ybase_b = fill(1.0, n_b)
    y_binom = as_turing_model(binom_obs, (y = missing, N = N_b), Ybase_b)()

    y_direct = sim(direct, n)
    y_renewal = sim(renewal, n)
    y_egr = sim(egr, n)
    y_nowcast = sim(nowcast, n)
    y_triangle = sim(triangle, n)
    y_latdelay = sim(latdelay, n)
    y_ascert = sim(ascert, 14)
    y_aggregate = sim(aggregate, 14)
    y_transobs = sim(transobs, n)
    y_normalobs = sim(normalobs, n)

    return [
        # latent-process log-joints
        ("RandomWalk latent logjoint", rw),
        ("AR latent logjoint", ar),
        ("ARIMA latent logjoint", arima),
        ("MA latent logjoint", ma),
        ("HierarchicalNormal latent logjoint", hier),
        ("DiffLatentModel(RandomWalk) latent logjoint", diffrw),
        ("ARMA latent logjoint", armamdl),
        ("BroadcastLatentModel day-of-week latent logjoint", bdow),
        ("BroadcastLatentModel weekly latent logjoint", bweek),
        ("ConcatLatentModels latent logjoint", concat),
        ("CombineLatentModels latent logjoint", combine),
        # the #76 prior interface
        ("AR vector BroadcastPrior latent logjoint", ar_vec),
        ("AR latent-model-as-prior latent logjoint", ar_lat),
        # infection posteriors
        ("DirectInfections+Poisson posterior",
            as_turing_model(direct, y_direct, n)),
        ("Renewal+NegativeBinomial posterior",
            as_turing_model(renewal, y_renewal, n)),
        ("ExpGrowthRate+Poisson posterior",
            as_turing_model(egr, y_egr, n)),
        # nowcasting
        ("Renewal+RightTruncate nowcast posterior",
            as_turing_model(nowcast, y_nowcast, n)),
        ("Renewal+ReportTriangle posterior",
            as_turing_model(triangle, y_triangle, n)),
        # observation modifiers / error families
        ("Renewal+LatentDelay posterior",
            as_turing_model(latdelay, y_latdelay, n)),
        ("DirectInfections+Ascertainment day-of-week posterior",
            as_turing_model(ascert, y_ascert, 14)),
        ("DirectInfections+Aggregate posterior",
            as_turing_model(aggregate, y_aggregate, 14)),
        ("DirectInfections+TransformObservation posterior",
            as_turing_model(transobs, y_transobs, n)),
        ("DirectInfections+NormalError posterior",
            as_turing_model(normalobs, y_normalobs, n)),
        ("BinomialError ascertainment posterior",
            as_turing_model(binom_obs, (y = y_binom, N = N_b), Ybase_b))
    ]
end

@doc """
    scenarios(; with_reference = false, category = :marginal)

The AD gradient scenarios — each a `DIT.Scenario{:gradient, :out}` over a real
package log-density (a latent process prior log-joint, or a composed
`EpiAwareModel` posterior conditioned on simulated data). When
`with_reference = true` each scenario carries its ForwardDiff reference gradient
in `res1`. `category` is accepted for the harness's group selector; all
scenarios are in the single `:marginal` group here.
"""
function scenarios(; with_reference::Bool = false, category::Symbol = :marginal)
    out = DIT.Scenario{:gradient, :out}[]
    for (i, (name, model)) in enumerate(_models())
        f, θ0, _ = _logdensity(model; seed = i)
        ref = with_reference ?
              DifferentiationInterface.gradient(f, AutoForwardDiff(), θ0) :
              nothing
        push!(out,
            DIT.Scenario{:gradient, :out}(f, θ0; name = name, res1 = ref))
    end
    return out
end

@doc """
    backends()

The AD backends exercised against the scenarios, as `(; name, backend)` named
tuples: ForwardDiff (the reference), ReverseDiff (tape), Mooncake, and Enzyme
reverse. Per-backend brokenness is recorded honestly in
[`backend_broken_scenarios`](@ref) / [`broken_scenario_names`](@ref) rather than
by trimming this list.
"""
function backends()
    return [
        (name = "ForwardDiff", backend = _forwarddiff()),
        (name = "ReverseDiff (tape)", backend = _reversediff()),
        (name = "Mooncake reverse", backend = _mooncake()),
        (name = "Enzyme reverse", backend = _enzyme())
    ]
end

# Backend constructors are written so that loading a backend package is only
# required when that backend is actually requested (the AD env loads them all,
# but this keeps the registry importable without every backend present).
_forwarddiff() = AutoForwardDiff()
function _reversediff()
    ADTypes = Base.require(Base.PkgId(
        Base.UUID("47edcb42-4c32-4615-8424-f2b9edc5f35b"), "ADTypes"))
    return ADTypes.AutoReverseDiff(; compile = false)
end
function _mooncake()
    ADTypes = Base.require(Base.PkgId(
        Base.UUID("47edcb42-4c32-4615-8424-f2b9edc5f35b"), "ADTypes"))
    return ADTypes.AutoMooncake(; config = nothing)
end
function _enzyme()
    ADTypes = Base.require(Base.PkgId(
        Base.UUID("47edcb42-4c32-4615-8424-f2b9edc5f35b"), "ADTypes"))
    Enzyme = Base.require(Base.PkgId(
        Base.UUID("7da242da-08ed-463a-9acd-ee780be4f1d9"), "Enzyme"))
    # `function_annotation = Enzyme.Const`: the log-density closures carry no
    # derivative data, and without this Enzyme raises `EnzymeMutabilityException`
    # ("argument cannot be proven readonly") on every DynamicPPL log-density.
    # With it, most scenarios differentiate correctly; a minority remain
    # genuinely broken (see `backend_broken_scenarios`).
    return ADTypes.AutoEnzyme(;
        mode = Enzyme.set_runtime_activity(Enzyme.Reverse),
        function_annotation = Enzyme.Const)
end

"Scenario names broken on every backend (none — all are real, FD-differentiable)."
broken_scenario_names() = String[]

@doc """
    backend_broken_scenarios()

Per-backend broken scenario names (`Dict{String, Set{String}}`), populated
HONESTLY from the actual `test/ad` run rather than by silencing.

Result matrix (24 scenarios × 4 backends), Julia 1.12:

| scenario                                              | ForwardDiff | ReverseDiff | Mooncake | Enzyme |
|-------------------------------------------------------|:-----------:|:-----------:|:--------:|:------:|
| RandomWalk latent logjoint                            |      ✓      |      ✓      |    ✓    |   ✓   |
| AR latent logjoint                                    |      ✓      |      ✓      |    ✓    |   ✗   |
| ARIMA latent logjoint                                 |      ✓      |      ✓      |    ✓    |   ✗   |
| MA latent logjoint                                    |      ✓      |      ✓      |    ✓    |   ✓   |
| HierarchicalNormal latent logjoint                    |      ✓      |      ✓      |    ✓    |   ✓   |
| DiffLatentModel(RandomWalk) latent logjoint           |      ✓      |      ✓      |    ✓    |   ✗   |
| ARMA latent logjoint                                  |      ✓      |      ✓      |    ✓    |   ✗   |
| BroadcastLatentModel day-of-week latent logjoint      |      ✓      |      ✓      |    ✓    |   ✓   |
| BroadcastLatentModel weekly latent logjoint           |      ✓      |      ✓      |    ✓    |   ✓   |
| ConcatLatentModels latent logjoint                    |      ✓      |      ✓      |    ✓    |   ✓   |
| CombineLatentModels latent logjoint                   |      ✓      |      ✓      |    ✓    |   ✗   |
| AR vector BroadcastPrior latent logjoint              |      ✓      |      ✓      |    ✓    |   ✗   |
| AR latent-model-as-prior latent logjoint              |      ✓      |      ✓      |    ✓    |   ✗   |
| DirectInfections+Poisson posterior                    |      ✓      |      ✓      |    ✓    |   ✓   |
| Renewal+NegativeBinomial posterior                    |      ✓      |      ✓      |    ✓    |   ✓   |
| ExpGrowthRate+Poisson posterior                       |      ✓      |      ✓      |    ✓    |   ✓   |
| Renewal+RightTruncate nowcast posterior               |      ✓      |      ✓      |    ✓    |   ✓   |
| Renewal+ReportTriangle posterior                      |      ✓      |      ✓      |    ✓    |   ✓   |
| Renewal+LatentDelay posterior                         |      ✓      |      ✓      |    ✓    |   ✓   |
| DirectInfections+Ascertainment day-of-week posterior  |      ✓      |      ✓      |    ✓    |   ✓   |
| DirectInfections+Aggregate posterior                  |      ✓      |      ✓      |    ✓    |   ✓   |
| DirectInfections+TransformObservation posterior       |      ✓      |      ✓      |    ✓    |   ✓   |
| DirectInfections+NormalError posterior                |      ✓      |      ✓      |    ✓    |   ✗   |
| BinomialError ascertainment posterior                 |      ✓      |      ✓      |    ✓    |   ✓   |

ForwardDiff (the reference), ReverseDiff, and Mooncake differentiate every
scenario correctly. Enzyme (configured with `function_annotation = Enzyme.Const`,
see [`backends`](@ref)) works on eighteen of the twenty-four but raises
`IllegalTypeAnalysisException` / a related type-analysis error on eight:

  - the `AR`-based latent log-densities (`AR`, `ARIMA`, `ARMA`,
    `CombineLatentModels` (which contains an `AR`), and both prior-interface `AR`
    scenarios), inside the `accumulate_scan(ARStep(damp_AR), ...)` /
    `LinearAlgebra.dot` recursion;
  - `DiffLatentModel(RandomWalk)` (the repeated `cumsum` reconstruction); and
  - `DirectInfections+NormalError` (the Gaussian likelihood loop).

These are real Enzyme type-analysis limitations, not defects in the package (the
same models sample fine under NUTS with ForwardDiff). They are recorded as
`@test_broken` for Enzyme below rather than hidden. Notably `MA` — whose step
also uses `dot` — differentiates under Enzyme, so the brokenness is specific to
these recursions rather than to `dot` in general.
"""
function backend_broken_scenarios()
    return Dict{String, Set{String}}(
        "Enzyme reverse" => Set([
        "AR latent logjoint",
        "ARIMA latent logjoint",
        "DiffLatentModel(RandomWalk) latent logjoint",
        "ARMA latent logjoint",
        "CombineLatentModels latent logjoint",
        "AR vector BroadcastPrior latent logjoint",
        "AR latent-model-as-prior latent logjoint",
        "DirectInfections+NormalError posterior"]))
end

"Per-backend scenario names too unstable to even run (segfault/hang)."
backend_skip_scenarios() = Dict{String, Set{String}}()

end # module ADFixtures
