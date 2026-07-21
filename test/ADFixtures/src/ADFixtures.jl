# PACKAGE-OWNED — scaffold writes this once and never overwrites it.
#
# AD-fixture registry implementing the EpiAwarePackageTools `ADRegistry`
# contract. The scenarios are REAL differentiable log-densities from the
# package: the (linked) log-joint of representative latent processes and of
# composed `IDModel`s conditioned on simulated data — the gradients an AD
# backend must get right for NUTS to work. Each scenario carries a ForwardDiff
# reference gradient. The shared harness (driven from `test/ad/setup.jl`)
# consumes this registry.
module ADFixtures

using ADTypes: AutoForwardDiff
using DifferentiationInterface: DifferentiationInterface
import DifferentiationInterfaceTest as DIT
import ForwardDiff
using ComposableTuringIDModels
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
    gen_int = _GEN_INT
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
    # A *vector* of damping distributions (order 2): one i.i.d. draw per lag
<<<<<<< HEAD
    # (identical priors, so the `filldist` branch of the seam), threaded as a
    # submodel via `as_turing_submodel`.
||||||| 25949fc
=======
    # through `arraydist`, threaded as a submodel via `as_turing_submodel`.
>>>>>>> origin/main
    ar_vec = as_turing_model(
        AR(; damp = [truncated(Normal(0, 0.05), 0, 1),
                truncated(Normal(0, 0.05), 0, 1)],
            init = [Normal(), Normal()]), 8)
    # A process as the damping prior: the bare `AR(damp = RandomWalk())` form —
    # now a genuinely TIME-VARYING coefficient path (issue #80 for the threading).
    # The AR damping coefficient is a length-(n-1) `RandomWalk` submodel mapped
    # through `tanh`, so the submodel-threading gradient path is differentiated.
    # The prior slot prefixes the latent-model prior (the `damp_AR` namespace) via
    # `as_turing_submodel(...; prefix = true)`, keeping the inner
    # `std`/`ϵ_t`/`rw_init` names from colliding with the AR innovation's — so this
    # linked log-density both evaluates and differentiates without a manual prefix.
    ar_lat = as_turing_model(AR(; damp = RandomWalk()), 8)

    # --- infection posteriors ---------------------------------------------------
    direct = IDModel(
        DirectInfections(; Z = RandomWalk(), initialisation = Normal()),
        PoissonError())
    renewal = IDModel(
        Renewal(; generation_time = gen_int, rt = RandomWalk(), initialisation = Normal()),
        NegativeBinomialError())
    # Exponential-growth-rate infections (the third infection family alongside
    # `DirectInfections` / `Renewal`): a cumulative growth-rate path exponentiated.
    egr = IDModel(
        ExpGrowthRate(; rt = RandomWalk(), initialisation = Normal()),
        PoissonError())

    # Nowcasting MARGINAL (right-truncation correction): a renewal model whose
    # observation error is wrapped in `RightTruncate` (fixed reporting-delay CDF
    # supplied as a `ReportingCDF` submodel). This exercises the `reverse`/
    # broadcast scaling the modifier adds on top of the inner error.
    nowcast = IDModel(
        Renewal(; generation_time = gen_int, rt = RandomWalk(), initialisation = Normal()),
        RightTruncate(NegativeBinomialError(),
            truncated(Normal(4.0, 1.5), 0.0, Inf)))

    # Nowcasting JOINT (2D reporting triangle): a renewal model feeding the
    # per-cell `ReportTriangle` observation model. The gradient of the per-cell
    # Poisson log-likelihood over the masked triangle (`t + d ≤ now`) is what
    # nowcasting under NUTS depends on.
    triangle = IDModel(
        Renewal(; generation_time = gen_int, rt = RandomWalk(), initialisation = Normal()),
        ReportTriangle(PoissonError(), [0.6, 0.25, 0.15]))

    # --- observation modifiers / error families over a composed model ----------
    # Reporting delay: convolves the expected observations with a delay PMF
    # (`accumulate_scan(LDStep(rev_pmf), ...)`) before the inner error.
    latdelay = IDModel(
        Renewal(; generation_time = gen_int, rt = RandomWalk(), initialisation = Normal()),
        LatentDelay(NegativeBinomialError(), [0.3, 0.4, 0.3]))
    # Uncertain reporting delay: the delay distribution's parameters are prior
    # slots (a `LogNormal` whose meanlog/sdlog carry priors), sampled through the
    # priors seam and rediscretised into a PMF per draw before the same
    # convolution. The gradient must flow through the discretisation
    # (`_discretised_pmf` / `double_interval_censored` `pdf`) — the AD-sensitive
    # part of an inferred delay.
    udelay = IDModel(
        Renewal(; generation_time = gen_int, rt = RandomWalk(), initialisation = Normal()),
        LatentDelay(NegativeBinomialError(),
            UncertainDelay(LogNormal,
                [Normal(1.0, 0.3), truncated(Normal(0.4, 0.2), 0, Inf)];
                D = 6.0)))
    # Time-varying reporting delay: the delay distribution's meanlog is a latent
    # process (a `RandomWalk`), so the delay — and its discretised pmf — varies
    # with time. Each per-time pmf is built through the priors seam and the
    # time-indexed convolution (`TimeVaryingLDStep`) is driven by a reversed kernel
    # per step. The gradient must flow through the per-time discretisation
    # (`_discretised_pmf`) and the process submodel threading — the time-varying
    # counterpart of `udelay`.
    tvdelay = IDModel(
        Renewal(; generation_time = gen_int, rt = RandomWalk(), initialisation = Normal()),
        LatentDelay(NegativeBinomialError(),
            UncertainDelay(LogNormal,
                [RandomWalk(), truncated(Normal(0.4, 0.2), 0, Inf)];
                D = 6.0)))
    # Uncertain generation interval: the renewal generation interval is itself
    # inferred — a `LogNormal` whose meanlog/sdlog carry priors, sampled through
    # the priors seam and discretised into a pmf per draw (lag-0 bin dropped,
    # renormalised) before the renewal step is built. The gradient must flow
    # through the discretisation (`_discretised_pmf`) and the renewal recursion
    # built from the sampled interval — the AD-sensitive part of an inferred
    # generation interval, the renewal counterpart of `udelay`.
    ugen = IDModel(
        Renewal(;
            generation_time = UncertainDelay(LogNormal,
                [Normal(0.7, 0.3), truncated(Normal(0.4, 0.2), 0, Inf)];
                D = 6.0),
            rt = RandomWalk(), initialisation = Normal()),
        PoissonError())
    # Day-of-week ascertainment: scales the expected observations by a broadcast
    # latent (an `Ascertainment` wrapping `broadcast_dayofweek`).
    ascert = IDModel(
        DirectInfections(; Z = RandomWalk(), initialisation = Normal()),
        ascertainment_dayofweek(PoissonError()))
    # Aggregation: sum the expected observations over weekly reporting windows
    # (only the window endpoints are scored).
    aggregate = IDModel(
        DirectInfections(; Z = RandomWalk(), initialisation = Normal()),
        Aggregate(PoissonError(), [0, 0, 0, 0, 0, 0, 7]))
    # Transform-the-expected-observations: softplus applied before the error.
    transobs = IDModel(
        DirectInfections(; Z = RandomWalk(), initialisation = Normal()),
        TransformObservationModel(PoissonError()))
    # Gaussian observation error (continuous, `σ`-inferred) rather than a count
    # family — the minimal non-count likelihood.
    normalobs = IDModel(
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
    y_binom = as_turing_model(binom_obs, (y = missing, N = N_b), Ybase_b)().y_t

    # `Split` observation composition: a renewal model observed through two
    # streams, with `deaths` cascaded downstream of `cases` by sharing the case
    # delay and splitting after it. Exercises the per-stream prefixing and the
    # expected-series threading gradient path.
    split = IDModel(
        Renewal(; generation_time = gen_int, rt = RandomWalk(), initialisation = Normal()),
        LatentDelay(
            Split((
                cases = NegativeBinomialError(),
                deaths = LatentDelay(
                    Ascertainment(NegativeBinomialError(), FixedIntercept(log(0.1))),
                    [0.2, 0.3, 0.5]))),
            [0.4, 0.3, 0.2, 0.1]))

    y_direct = sim(direct, n)
    y_renewal = sim(renewal, n)
    y_egr = sim(egr, n)
    y_nowcast = sim(nowcast, n)
    y_triangle = sim(triangle, n)
    y_latdelay = sim(latdelay, n)
    y_udelay = sim(udelay, n)
    y_tvdelay = sim(tvdelay, n)
    y_ugen = sim(ugen, n)
    y_ascert = sim(ascert, 14)
    y_aggregate = sim(aggregate, 14)
    y_transobs = sim(transobs, n)
    y_normalobs = sim(normalobs, n)
    y_split = sim(split, n)

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
        ("AR vector-prior latent logjoint", ar_vec),
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
        ("Renewal+UncertainLatentDelay posterior",
            as_turing_model(udelay, y_udelay, n)),
        ("Renewal+TimeVaryingLatentDelay posterior",
            as_turing_model(tvdelay, y_tvdelay, n)),
        ("Renewal+UncertainGenInterval posterior",
            as_turing_model(ugen, y_ugen, n)),
        ("DirectInfections+Ascertainment day-of-week posterior",
            as_turing_model(ascert, y_ascert, 14)),
        ("DirectInfections+Aggregate posterior",
            as_turing_model(aggregate, y_aggregate, 14)),
        ("DirectInfections+TransformObservation posterior",
            as_turing_model(transobs, y_transobs, n)),
        ("DirectInfections+NormalError posterior",
            as_turing_model(normalobs, y_normalobs, n)),
        ("BinomialError ascertainment posterior",
            as_turing_model(binom_obs, (y = y_binom, N = N_b), Ybase_b)),
        # unified Split observation composition
        ("Renewal+Split cascade posterior",
            as_turing_model(split, y_split, n))
    ]
end

@doc """
    scenarios(; with_reference = false, category = :marginal)

The AD gradient scenarios — each a `DIT.Scenario{:gradient, :out}` over a real
package log-density (a latent process prior log-joint, or a composed
`IDModel` posterior conditioned on simulated data). When
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

Result matrix (28 scenarios × 4 backends), Julia 1.12:

| scenario                                              | ForwardDiff | ReverseDiff | Mooncake | Enzyme |
|-------------------------------------------------------|:-----------:|:-----------:|:--------:|:------:|
| RandomWalk latent logjoint                            |      ✓      |      ✓      |    ✓    |   ✓   |
| AR latent logjoint                                    |      ✓      |      ✓      |    ✓    |   ✗   |
| ARIMA latent logjoint                                 |      ✓      |      ✓      |    ✓    |   ✗   |
| MA latent logjoint                                    |      ✓      |      ✓      |    ✓    |   ✗   |
| HierarchicalNormal latent logjoint                    |      ✓      |      ✓      |    ✓    |   ✓   |
| DiffLatentModel(RandomWalk) latent logjoint           |      ✓      |      ✓      |    ✓    |   ✗   |
| ARMA latent logjoint                                  |      ✓      |      ✓      |    ✓    |   ✗   |
| BroadcastLatentModel day-of-week latent logjoint      |      ✓      |      ✓      |    ✓    |   ✓   |
| BroadcastLatentModel weekly latent logjoint           |      ✓      |      ✓      |    ✓    |   ✓   |
| ConcatLatentModels latent logjoint                    |      ✓      |      ✓      |    ✓    |   ✓   |
| CombineLatentModels latent logjoint                   |      ✓      |      ✓      |    ✓    |   ✗   |
| AR vector-prior latent logjoint                       |      ✓      |      ✓      |    ✓    |   ✗   |
| AR latent-model-as-prior latent logjoint              |      ✓      |      ✓      |    ✓    |   ✗   |
| DirectInfections+Poisson posterior                    |      ✓      |      ✓      |    ✓    |   ✓   |
| Renewal+NegativeBinomial posterior                    |      ✓      |      ✓      |    ✓    |   ✗   |
| ExpGrowthRate+Poisson posterior                       |      ✓      |      ✓      |    ✓    |   ✓   |
| Renewal+RightTruncate nowcast posterior               |      ✓      |      ✓      |    ✓    |   ✓   |
| Renewal+ReportTriangle posterior                      |      ✓      |      ✓      |    ✓    |   ✓   |
| Renewal+LatentDelay posterior                         |      ✓      |      ✓      |    ✓    |   ✓   |
<<<<<<< HEAD
| Renewal+UncertainLatentDelay posterior                |      ✓      |      ✓      |    ✓    |   ✓   |
| Renewal+TimeVaryingLatentDelay posterior              |      ✓      |      ✓      |    ✓    |   ✗   |
| Renewal+UncertainGenInterval posterior                |      ✓      |      ✓      |    ✓    |   ✓   |
||||||| 25949fc
| scenario                              | ForwardDiff | ReverseDiff | Mooncake | Enzyme |
|---------------------------------------|:-----------:|:-----------:|:--------:|:------:|
| RandomWalk latent logjoint            |      ✓      |      ✓      |    ✓    |   ✓   |
| AR latent logjoint                    |      ✓      |      ✓      |    ✓    |   ✗   |
| ARIMA latent logjoint                 |      ✓      |      ✓      |    ✓    |   ✗   |
| DirectInfections+Poisson posterior    |      ✓      |      ✓      |    ✓    |   ✓   |
| Renewal+NegativeBinomial posterior    |      ✓      |      ✓      |    ✓    |   ✓   |
| Renewal+RightTruncate nowcast posterior |    ✓      |      ✓      |    ✓    |   ✓   |
| Renewal+ReportTriangle posterior      |      ✓      |      ✓      |    ✓    |   ✓   |
| Renewal+Split cascade posterior       |      ✓      |      ✓      |    ✓    |   ✗   |
=======
| Renewal+UncertainLatentDelay posterior                |      ✓      |      ✓      |    ✓    |   ✗   |
| Renewal+TimeVaryingLatentDelay posterior              |      ✓      |      ✓      |    ✓    |   ✗   |
| Renewal+UncertainGenInterval posterior                |      ✓      |      ✓      |    ✓    |   ✗   |
>>>>>>> origin/main
| DirectInfections+Ascertainment day-of-week posterior  |      ✓      |      ✓      |    ✓    |   ✗   |
| DirectInfections+Aggregate posterior                  |      ✓      |      ✓      |    ✓    |   ✓   |
| DirectInfections+TransformObservation posterior       |      ✓      |      ✓      |    ✓    |   ✓   |
| DirectInfections+NormalError posterior                |      ✓      |      ✓      |    ✓    |   ✗   |
| BinomialError ascertainment posterior                 |      ✓      |      ✓      |    ✓    |   ✓   |
| Renewal+Split cascade posterior                       |      ✓      |      ✓      |    ✓    |   ✗   |

scenario correctly. Enzyme (configured with `function_annotation = Enzyme.Const`,
<<<<<<< HEAD
see [`backends`](@ref)) works on fifteen of the twenty-eight but raises
`IllegalTypeAnalysisException` / a related type-analysis or shadow error on
thirteen:

  - the `AR`-based latent log-densities (`AR`, `ARIMA`, `ARMA`,
    `CombineLatentModels` (which contains an `AR`), and both prior-interface `AR`
    scenarios), inside the `accumulate_scan(ARStep(damp_AR), ...)` /
    `LinearAlgebra.dot` recursion;
  - `DiffLatentModel(RandomWalk)` (the repeated `cumsum` reconstruction);
  - `DirectInfections+NormalError` (the Gaussian likelihood loop);
  - `DirectInfections+Ascertainment day-of-week` — an `EnzymeNoShadowError`
    through the `PrefixLatentModel`-wrapped day-of-week `BroadcastLatentModel`
    submodel threading;
  - `Renewal+TimeVaryingLatentDelay` — the process-parameter (time-varying) delay
    collects its per-parameter submodel returns in a `Vector{Any}` and splats a
    runtime-length generator into the delay constructor, so Enzyme boxes the
    elements in a `Base.RefValue` the `_at` accessor cannot index; and
  - the deeply-nested `Renewal+Split` cascade, through its per-stream submodel
    threading.

The uncertain-delay / uncertain-generation-interval scenarios
(`Renewal+UncertainLatentDelay`, `Renewal+UncertainGenInterval`) previously
raised an `EnzymeNoShadowError` here. That was **not** the prefix seam: a
heterogeneous prior vector (e.g. a `LogNormal` delay's `[meanlog, sdlog]`
priors) was drawn through DynamicPPL's `arraydist`, which routes a vector of
univariate distributions through the **deprecated** `Distributions.Product`
constructor. Its `Base.depwarn`/`invokelatest` world-age machinery has no shadow
under Enzyme reverse. Building the identical product with `product_distribution`
instead (see
`ComposableTuringIDModels`'s prior seam) removes the deprecation path, and both
scenarios now differentiate correctly under Enzyme (#97).
||||||| 25949fc
ForwardDiff (the reference), ReverseDiff, and Mooncake differentiate every
scenario correctly — including both nowcasting models (the `RightTruncate`
marginal and the `ReportTriangle` joint triangle) and the unified `Split`
observation composition. Enzyme works on five of the eight once configured with
`function_annotation = Enzyme.Const` (see [`backends`](@ref)). The two AR-based
latent log-densities raise `IllegalTypeAnalysisException` inside the
`accumulate_scan(ARStep(damp_AR), ...)` / `LinearAlgebra.dot` recursion, and the
deeply-nested `Split` composition raises the same exception through its
per-stream submodel threading — real Enzyme type-analysis limitations, not
defects in the package (all three sample fine under NUTS with ForwardDiff or
Mooncake). They are recorded as `@test_broken` for Enzyme below rather than
hidden.
=======
see [`backends`](@ref)) works on sixteen of the twenty-eight but raises
`IllegalTypeAnalysisException` / a related type-analysis or shadow error on
twelve:

  - the `AR`-based latent log-densities (`AR`, `ARIMA`, `ARMA`,
    `CombineLatentModels` (which contains an `AR`), and both prior-interface `AR`
    scenarios), inside the `accumulate_scan(ARStep(damp_AR), ...)` /
    `LinearAlgebra.dot` recursion;
  - `DiffLatentModel(RandomWalk)` (the repeated `cumsum` reconstruction);
  - `DirectInfections+NormalError` (the Gaussian likelihood loop);
  - `DirectInfections+Ascertainment day-of-week` — an `EnzymeNoShadowError`
    through the `PrefixLatentModel`-wrapped day-of-week `BroadcastLatentModel`
    submodel threading, surfaced by the #76 prefix-on prior collapse;
  - `Renewal+TimeVaryingLatentDelay` — the process-parameter (time-varying) delay
    threads a `RandomWalk` submodel through the priors seam per delay parameter,
    where Enzyme wraps the constant `sdlog` scalar in a `Base.RefValue` its
    type-analysis cannot resolve against the `_at` accessor; and
  - the deeply-nested `Renewal+Split` cascade, through its per-stream submodel
    threading.
>>>>>>> origin/main

These are real Enzyme type-analysis limitations, not defects in the package (the
same models sample fine under NUTS with ForwardDiff or Mooncake). They are
recorded as `@test_broken` for Enzyme below rather than hidden. Notably `MA` —
whose step also uses `dot` — differentiates under Enzyme, so the brokenness is
specific to these recursions rather than to `dot` in general.
"""
function backend_broken_scenarios()
    return Dict{String, Set{String}}(
        "Enzyme reverse" => Set([
        "AR latent logjoint",
        "ARIMA latent logjoint",
        # Plain `MA` threads its `HierarchicalNormal` innovation (and vector-`θ`)
        # submodels through the prior seam, the same `EnzymeNoShadowError`
        # submodel-threading limit as its `AR`/`ARIMA`/`ARMA` siblings above.
        "MA latent logjoint",
        "DiffLatentModel(RandomWalk) latent logjoint",
        "ARMA latent logjoint",
        "CombineLatentModels latent logjoint",
        "AR vector-prior latent logjoint",
        "AR latent-model-as-prior latent logjoint",
        "DirectInfections+NormalError posterior",
        "Renewal+NegativeBinomial posterior",
        # Enzyme type-analysis brokenness tracked in #97.
        "Renewal+Split cascade posterior",
        # `EnzymeNoShadowError` through the `Ascertainment` +
        # `broadcast_dayofweek` submodel threading after the #76 prior collapse
        # (prefix-on prior slots): Enzyme cannot find a shadow for the
        # `PrefixLatentModel`-wrapped day-of-week `BroadcastLatentModel`. Enzyme
        # only; ForwardDiff/ReverseDiff/Mooncake differentiate it correctly.
        # Tracked in #97.
        "DirectInfections+Ascertainment day-of-week posterior",
<<<<<<< HEAD
        # The process-parameter (time-varying) reporting delay collects its
        # per-parameter submodel returns in a `Vector{Any}` and splats a
        # runtime-length generator into the delay constructor; Enzyme boxes the
        # elements in a `Base.RefValue` and the generator splat reaches `_at` on
        # the box (`MethodError: no method matching _at(::RefValue, ::Int)`). A
        # distinct cause from the `arraydist` shadow error fixed for the two
        # uncertain scenarios; Enzyme only, correct under the other three
        # backends. Tracked in #97.
        "Renewal+TimeVaryingLatentDelay posterior"]))
||||||| 25949fc
        "AR latent logjoint", "ARIMA latent logjoint",
        "Renewal+Split cascade posterior"]))
=======
        # The process-parameter (time-varying) reporting delay threads a
        # `RandomWalk` submodel per delay parameter; Enzyme wraps the constant
        # `sdlog` scalar in a `Base.RefValue` its type-analysis cannot resolve
        # against `_at` (`MethodError: no method matching _at(::RefValue, ::Int)`).
        # Enzyme only; ForwardDiff/ReverseDiff/Mooncake differentiate it correctly.
        # Tracked in #97.
        "Renewal+TimeVaryingLatentDelay posterior",
        # An inferred (uncertain) reporting delay draws its distribution
        # parameters through the prefix-on prior-slot seam
        # (`as_turing_submodel(delay; prefix = true)`) — the same
        # `EnzymeNoShadowError` submodel-threading limitation as the day-of-week
        # ascertainment above. Enzyme only; ForwardDiff/ReverseDiff/Mooncake
        # differentiate it correctly. Tracked in #97.
        "Renewal+UncertainLatentDelay posterior",
        # An inferred (uncertain) generation interval threads its `UncertainDelay`
        # parameters through the same prefix-on prior-slot seam — the same
        # `EnzymeNoShadowError` limit. Enzyme only; ForwardDiff/ReverseDiff/
        # Mooncake differentiate it correctly. Tracked in #97.
        "Renewal+UncertainGenInterval posterior"]))
>>>>>>> origin/main
end

"Per-backend scenario names too unstable to even run (segfault/hang)."
backend_skip_scenarios() = Dict{String, Set{String}}()

end # module ADFixtures
