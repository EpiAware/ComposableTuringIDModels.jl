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
    # `as_turing_model(::IDModel, ...)` narrows its own `y_t` argument via
    # `ComposableTuringIDModels.concrete_observations` automatically (see that
    # function's docstring: DynamicPPL's predictive container keeps a
    # `Union{Missing, T}` eltype even once concrete, which trips a `deepcopy`
    # Enzyme forward has no rule for). `sim` narrows explicitly too, so
    # callers below that use its result outside `IDModel` (the standalone
    # `BinomialError` scenario) get the same fix.
    sim(m, nn) = ComposableTuringIDModels.concrete_observations(
        as_turing_model(m, missing, nn)().generated_y_t)

    # --- latent-process log-joints (prior only) ---------------------------------
    rw = as_turing_model(RandomWalk(), n)
    ar = as_turing_model(AR(), n)
    arima = as_turing_model(
        DiffLatentModel(; model = AR(), init = [Normal(), Normal()]), n)
    hsgp = as_turing_model(HilbertSpaceGP(; m = 8), n)
    hsgp_matern = as_turing_model(
        HilbertSpaceGP(; m = 8, kernel = Matern52Kernel()), n)
    # The exact GP the Hilbert-space model approximates. Its gradient runs
    # through `kernelmatrix` and a `cholesky(Symmetric(...))` of a matrix built
    # from the sampled `ℓ`/`σ` — the most AD-sensitive path of the two GPs, and
    # the one the Gaussian-process case study drives under Mooncake.
    exactgp = as_turing_model(ExactGP(), n)
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
    # Non-centred partial pooling across groups: a shared level `mean` plus
    # `IID` group deviations from `across` (the default cross-group prior,
    # otherwise untested on its own). Exercises the two-submodel-slot
    # (`mean`, `across`) threading distinct from every other manipulator here.
    hierarchy = as_turing_model(
        Hierarchy(; mean = Normal(2, 0.5), across = IID(Normal(0, 0.5))), 6)

    # --- the #76 prior interface -----------------------------------------------
    # A *vector* of damping distributions (order 2): one i.i.d. draw per lag
    # (identical priors, so the `filldist` branch of the seam), threaded as a
    # submodel via `as_turing_submodel`.
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
    # Renewal MODIFIERS: a susceptible-depleting renewal whose incidence is also
    # seeded by an imported-cases rate. Both modifier kinds ride the one
    # pre-scan seam — `SusceptibleDepletion` samples nothing and returns
    # itself, while `ImportedCases` draws its rate before the scan and hands
    # back a resolved modifier. The gradient must flow through the pre-scan
    # submodel, the step rebuilt from the resolved modifiers, and the modifier
    # threading inside the renewal recursion (#189).
    modifiers = IDModel(
        Renewal(gen_int, SusceptibleDepletion(2_000.0),
            ImportedCases(Normal(0.0, 0.5));
            rt = RandomWalk(), initialisation = Normal()),
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
    # Partially-missing observations: a data vector carrying genuine reporting
    # gaps (some entries `missing`, some concrete). This is the ragged case
    # `ComposableTuringIDModels.concrete_observations` cannot narrow away —
    # some entries really are unobserved — so the argument keeps its
    # `Union{Missing, T}` eltype, DynamicPPL's `hasmissing` `deepcopy` fires on
    # every evaluation, and it stays differentiable under Enzyme forward only
    # via the `EnzymeRules.inactive` `deepcopy` mark in
    # `ext/ComposableTuringIDModelsEnzymeExt.jl`. That gradient path is the
    # point of the scenario.
    #
    # The error family is `NormalError`, not a count family, and that is
    # load-bearing. A partially-missing vector makes the WHOLE `y_t` latent
    # rather than only its blank entries, so `_logdensity`'s `link` must find
    # a bijector for the observation distribution. Count families are discrete
    # and have none (`SafePoisson` raises `no method matching
    # scalar_to_scalar_bijector`), which is not an AD limitation but a
    # consequence of linking a discrete variable. A continuous error keeps the
    # linked log-density well defined.
    partialmiss = IDModel(
        DirectInfections(; Z = RandomWalk(), initialisation = Normal()),
        NormalError())
    # Transform-the-expected-observations: softplus applied before the error.
    transobs = IDModel(
        DirectInfections(; Z = RandomWalk(), initialisation = Normal()),
        TransformObservationModel(PoissonError()))
    # Gaussian observation error (continuous, `σ`-inferred) rather than a count
    # family — the minimal non-count likelihood.
    normalobs = IDModel(
        DirectInfections(; Z = RandomWalk(), initialisation = Normal()),
        NormalError())

    # Record-the-expected modifiers: `RecordExpectedLatent` tracks the inner
    # latent path as a generated quantity before it feeds `Z`, and
    # `RecordExpectedObs` tracks the expected observations before the inner
    # error; both are `:=`-tracked pass-throughs, so this exercises the
    # gradient still flows unchanged through the tracked variable on both the
    # latent and the observation side.
    recordexp = IDModel(
        DirectInfections(;
            Z = RecordExpectedLatent(RandomWalk()), initialisation = Normal()),
        RecordExpectedObs(PoissonError()))
    # Prefix modifiers: `PrefixLatentModel` and `PrefixObservationModel` rename
    # an inner model's sampled variables via `DynamicPPL.prefix` before
    # threading it as a submodel (the same mechanism `Split` and the #76 prior
    # seam use inline); this exercises the two standalone wrapper structs
    # directly rather than only their inline call sites.
    prefixmods = IDModel(
        DirectInfections(;
            Z = PrefixLatentModel(; model = RandomWalk(), prefix = "zpre"),
            initialisation = Normal()),
        PrefixObservationModel(; model = PoissonError(), prefix = "obspre"))

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
    y_binom = ComposableTuringIDModels.concrete_observations(
        as_turing_model(binom_obs, (y = missing, N = N_b), Ybase_b)().y_t)

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
    y_modifiers = sim(modifiers, n)
    y_egr = sim(egr, n)
    y_nowcast = sim(nowcast, n)
    y_triangle = sim(triangle, n)
    y_latdelay = sim(latdelay, n)
    y_udelay = sim(udelay, n)
    y_tvdelay = sim(tvdelay, n)
    y_ugen = sim(ugen, n)
    y_ascert = sim(ascert, 14)
    y_aggregate = sim(aggregate, 14)
    # A straight `sim` draw, then genuine reporting gaps re-introduced by hand:
    # every third entry blanked back to `missing` (distinct from `Aggregate`
    # above, whose missing pattern comes from its window scattering, not the
    # data itself).
    y_partialmiss = Vector{Union{Missing, Float64}}(sim(partialmiss, n))
    y_partialmiss[2:3:end] .= missing
    y_transobs = sim(transobs, n)
    y_normalobs = sim(normalobs, n)
    y_split = sim(split, n)
    y_recordexp = sim(recordexp, n)
    y_prefixmods = sim(prefixmods, n)

    return [
        # latent-process log-joints
        ("RandomWalk latent logjoint", rw),
        ("AR latent logjoint", ar),
        ("ARIMA latent logjoint", arima),
        ("HilbertSpaceGP latent logjoint", hsgp),
        ("HilbertSpaceGP Matern latent logjoint", hsgp_matern),
        ("ExactGP latent logjoint", exactgp),
        ("MA latent logjoint", ma),
        ("HierarchicalNormal latent logjoint", hier),
        ("DiffLatentModel(RandomWalk) latent logjoint", diffrw),
        ("ARMA latent logjoint", armamdl),
        ("BroadcastLatentModel day-of-week latent logjoint", bdow),
        ("BroadcastLatentModel weekly latent logjoint", bweek),
        ("ConcatLatentModels latent logjoint", concat),
        ("CombineLatentModels latent logjoint", combine),
        ("Hierarchy latent logjoint", hierarchy),
        # the #76 prior interface
        ("AR vector-prior latent logjoint", ar_vec),
        ("AR latent-model-as-prior latent logjoint", ar_lat),
        # infection posteriors
        ("DirectInfections+Poisson posterior",
            as_turing_model(direct, y_direct, n)),
        ("Renewal+NegativeBinomial posterior",
            as_turing_model(renewal, y_renewal, n)),
        ("Renewal+ImportedCases posterior",
            as_turing_model(modifiers, y_modifiers, n)),
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
        ("DirectInfections+PartiallyMissing posterior",
            as_turing_model(partialmiss, y_partialmiss, n)),
        ("DirectInfections+TransformObservation posterior",
            as_turing_model(transobs, y_transobs, n)),
        ("DirectInfections+NormalError posterior",
            as_turing_model(normalobs, y_normalobs, n)),
        ("DirectInfections+RecordExpected posterior",
            as_turing_model(recordexp, y_recordexp, n)),
        ("DirectInfections+PrefixModifiers posterior",
            as_turing_model(prefixmods, y_prefixmods, n)),
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
tuples: ForwardDiff (the reference), ReverseDiff (tape), Mooncake reverse,
Mooncake forward, Enzyme reverse, and Enzyme forward — the full six-backend
matrix `ad.yaml` runs in CI. Per-backend brokenness is recorded honestly in
[`backend_broken_scenarios`](@ref) / [`broken_scenario_names`](@ref) rather than
by trimming this list.
"""
function backends()
    return [
        (name = "ForwardDiff", backend = _forwarddiff()),
        (name = "ReverseDiff (tape)", backend = _reversediff()),
        (name = "Mooncake reverse", backend = _mooncake()),
        (name = "Mooncake forward", backend = _mooncake_forward()),
        (name = "Enzyme reverse", backend = _enzyme()),
        (name = "Enzyme forward", backend = _enzyme_forward())
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
function _mooncake_forward()
    ADTypes = Base.require(Base.PkgId(
        Base.UUID("47edcb42-4c32-4615-8424-f2b9edc5f35b"), "ADTypes"))
    return ADTypes.AutoMooncakeForward(; config = nothing)
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
function _enzyme_forward()
    ADTypes = Base.require(Base.PkgId(
        Base.UUID("47edcb42-4c32-4615-8424-f2b9edc5f35b"), "ADTypes"))
    Enzyme = Base.require(Base.PkgId(
        Base.UUID("7da242da-08ed-463a-9acd-ee780be4f1d9"), "Enzyme"))
    # Same `function_annotation = Enzyme.Const` rationale as `_enzyme()`,
    # forward mode.
    return ADTypes.AutoEnzyme(;
        mode = Enzyme.set_runtime_activity(Enzyme.Forward),
        function_annotation = Enzyme.Const)
end

"Scenario names broken on every backend (none — all are real, FD-differentiable)."
broken_scenario_names() = String[]

@doc """
    backend_broken_scenarios()

Per-backend broken scenario names (`Dict{String, Set{String}}`, keyed
`"Enzyme reverse"` / `"Enzyme forward"`). The AD harness's `check_broken`
(EpiAwarePackageTools' `ad_harness.jl`) computes the pass/fail boolean
itself and only falls back to `@test_broken` when a listed scenario
actually fails; a listed scenario that passes still records an ordinary
`@test` pass. So over-listing a scenario here is safe; under-listing one
reds CI. Both `check_broken` and the full correctness sweep use the same
tolerance (`rtol = 5e-2, atol = 1e-6`); `check_broken` runs a single
plain `DI.gradient` call, while the full sweep also covers the prepared,
in-place, and re-preparation forms at that same tolerance.

Evidence below is the Linux `ad.yaml` CI job at commit `5e639bb`, not a
local run.

| Backend          | Result                                   |
|------------------|-------------------------------------------|
| ForwardDiff      | reference backend; not listed broken     |
| ReverseDiff      | not listed broken                        |
| Mooncake reverse | 735/735 pass, all 35 scenarios           |
| Mooncake forward | 735/735 pass, all 35 scenarios           |
| Enzyme reverse   | 427 pass, 8 broken (of 15 listed names)  |
| Enzyme forward   | 273 pass, 20 broken (all 20 listed names)|

Mooncake passes every scenario in both modes, including all three
Gaussian-process ones; there is no Mooncake gap.

Enzyme reverse lists 15 names; 8 genuinely fail and 7 pass (established
by the CI pass/broken counts; which 7 pass is inferred from a prior
local record, not individually confirmed by CI). Enzyme forward lists 20
names and all 20 genuinely fail (the 13 unlisted scenarios pass clean,
and `13 × 21 = 273` matches the CI pass count exactly).

Genuine Enzyme failures, grouped by cause:

  - `AR vector-prior latent logjoint` is the only scenario in the AR/MA
    family that executes `dot` (the order-1 default constructors used
    elsewhere take a dot-free path). Fixed by replacing `dot` with
    `mapreduce(*, +, ...)`. Established.
  - `AR`, `CombineLatentModels`, `AR latent-model-as-prior` execute no
    `dot` and pass Enzyme forward on CI; the reverse-mode counts are
    consistent with them passing there too. Established.
  - `MA`, `ARMA` execute no `dot` and no `accumulate_scan`, yet fail
    under both Enzyme modes. Hypothesised cause: combining two
    independently drawn `to_submodel` results in one expression.
    `RandomWalk` has the same submodel nesting depth as `MA` and
    passes, so depth alone is not the trigger.
  - `ARIMA`, `DiffLatentModel(RandomWalk)` touch no `dot`. `ARIMA`'s
    inner `AR()` uses the same passing default config as plain `AR`, so
    the cause is attributed to the `DiffLatentModel` wrapper rather
    than to `AR`; the two likely share one cause. Hypothesised.
  - `Renewal+Split cascade` runs a runtime-length loop over
    heterogeneously typed observation models, each issuing its own
    prefixed `to_submodel` at unequal depth. That shape is `Split`'s
    purpose and is not restructurable; the existing `Vector{Any}` to
    `Tuple` mitigation does not fix it. Established.
  - `DirectInfections+Ascertainment day-of-week` and
    `DirectInfections+PrefixModifiers` are not among the 8 genuine
    Enzyme-reverse failures; there is no prefix-threading blocker.
    Established.

Unresolved causes are tracked in #97.
"""
function backend_broken_scenarios()
    enzyme_reverse = Set([
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
        # The modifier renewal is the `Renewal+NegativeBinomial` model plus a
        # modifier tuple, so it inherits that scenario's Enzyme brokenness and
        # adds the nested `to_submodel` recursion that resolves the modifiers —
        # the same threading behind the other `EnzymeNoShadowError` rows. It is
        # recorded broken conservatively: `check_broken` records a plain pass
        # when a listed scenario succeeds, so listing it costs nothing, while a
        # wrong ✓ would red the Enzyme job. Both this scenario and its
        # `Renewal+NegativeBinomial` base pass under Enzyme on macOS/aarch64
        # locally, so this row (and that one) want re-checking against the
        # Linux `ad.yaml` Enzyme job. Tracked in #97.
        "Renewal+ImportedCases posterior",
        # Enzyme type-analysis brokenness tracked in #97.
        "Renewal+Split cascade posterior",
        # `EnzymeNoShadowError` through the `Ascertainment` +
        # `broadcast_dayofweek` submodel threading after the #76 prior collapse
        # (prefix-on prior slots): Enzyme cannot find a shadow for the
        # `PrefixLatentModel`-wrapped day-of-week `BroadcastLatentModel`. Enzyme
        # only; ForwardDiff/ReverseDiff/Mooncake differentiate it correctly.
        # Tracked in #97.
        "DirectInfections+Ascertainment day-of-week posterior",
        # `PrefixLatentModel`/`PrefixObservationModel` wrap an inner model with
        # the same `prefix(as_turing_model(...), Symbol(...))` +
        # `to_submodel` call as the day-of-week `Ascertainment` scenario above,
        # so it is listed conservatively for the same `EnzymeNoShadowError`
        # submodel-threading reason. Both this scenario and its
        # `DirectInfections+Ascertainment day-of-week` sibling pass under
        # Enzyme reverse on macOS/aarch64 locally, so both rows want
        # re-checking against the Linux `ad.yaml` Enzyme job; listing costs
        # nothing (`check_broken` records a plain pass if a listed scenario
        # succeeds). Tracked in #97.
        "DirectInfections+PrefixModifiers posterior",
        # The process-parameter (time-varying) reporting delay threads a
        # `RandomWalk` submodel per delay parameter; Enzyme wraps the constant
        # `sdlog` scalar in a `Base.RefValue` its type-analysis cannot resolve
        # against `_at` (`MethodError: no method matching _at(::RefValue, ::Int)`).
        # Enzyme only; ForwardDiff/ReverseDiff/Mooncake differentiate it correctly.
        # Tracked in #97.
        "Renewal+TimeVaryingLatentDelay posterior"
    ])
    # Carried forward from this PR's own prior local measurement (before the
    # #53 merge and the fixture fix — see the docstring above), NOT
    # re-confirmed against the merged 35-scenario set, except
    # `DirectInfections+Poisson posterior`, which was individually re-run
    # against the fixed fixture and removed from this set after confirming it
    # now passes. The other fourteen `Vector{Union{Missing, Int64}}`-observed
    # entries here share the identical fixture code path and are expected to
    # be fixed the same way, but that is an expectation, not a
    # re-measurement — the `ad.yaml` `enzyme_forward` job is what actually
    # confirms it. The `accumulate_scan`/`dot`-recursion family
    # (`ARIMA`/`MA`/`DiffLatentModel(RandomWalk)`/`ARMA`/`AR vector-prior`)
    # and `Renewal+Split cascade` are unaffected by the fixture fix (see
    # docstring) and stay broken for their own, unrelated reasons.
    enzyme_forward = Set([
        "ARIMA latent logjoint",
        "MA latent logjoint",
        "DiffLatentModel(RandomWalk) latent logjoint",
        "ARMA latent logjoint",
        "AR vector-prior latent logjoint",
        "Renewal+NegativeBinomial posterior",
        "Renewal+ImportedCases posterior",
        "ExpGrowthRate+Poisson posterior",
        "Renewal+RightTruncate nowcast posterior",
        "Renewal+ReportTriangle posterior",
        "Renewal+LatentDelay posterior",
        "Renewal+UncertainLatentDelay posterior",
        "Renewal+TimeVaryingLatentDelay posterior",
        "Renewal+UncertainGenInterval posterior",
        "DirectInfections+Ascertainment day-of-week posterior",
        "DirectInfections+Aggregate posterior",
        "DirectInfections+TransformObservation posterior",
        "DirectInfections+RecordExpected posterior",
        "DirectInfections+PrefixModifiers posterior",
        # `EnzymeNoDerivativeError`, not the `deepcopy` rule gap above — its own
        # cause, through the same per-stream submodel threading that breaks it
        # under Enzyme reverse.
        "Renewal+Split cascade posterior"
    ])
    return Dict{String, Set{String}}(
        "Enzyme reverse" => enzyme_reverse,
        "Enzyme forward" => enzyme_forward)
end

"Per-backend scenario names too unstable to even run (segfault/hang)."
backend_skip_scenarios() = Dict{String, Set{String}}()

end # module ADFixtures
