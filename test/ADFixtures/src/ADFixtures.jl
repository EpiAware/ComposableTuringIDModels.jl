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
    # Some renewal/delay processes generate genuine leading `missing`
    # (warm-up) entries. We do not infer missing events here, so `sim`
    # coalesces them into a concrete observed vector (`hasmissing == false`)
    # at inference/conditioning time. That keeps the AD gradient path on a
    # concrete `Vector{Int}`/`Vector{Float64}` instead of a
    # `Union{Missing, T}` vector, which Enzyme forward cannot compile (its
    # `deepcopy` of a `Union{Missing, Int}` observed argument has no working
    # rule). Only plain arrays with a `Missing` eltype are touched; non-array
    # outputs (e.g. `ReportingTriangle`, `Split` stream tuples) pass through.
    _concretize(x::AbstractArray) = eltype(x) >: Missing ?
                                    coalesce.(x, zero(nonmissingtype(eltype(x)))) : x
    # Split-style streams are returned as a NamedTuple of observed vectors;
    # concretize each stream the same way as a plain array.
    _concretize(x::NamedTuple) = map(_concretize, x)
    _concretize(x) = x
    sim(m, nn) = _concretize(
        ComposableTuringIDModels.concrete_observations(
        as_turing_model(m, missing, nn)().generated_y_t))

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
    # Partially-missing observations (NormalError). We do not infer missing
    # events, so the observed data is kept concrete at inference time (the
    # same `sim` narrowing applied elsewhere) rather than threading a ragged
    # `Union{Missing, Float64}` vector — which DynamicPPL `deepcopy`s and
    # Enzyme cannot compile. A continuous error keeps the linked log-density
    # well defined (count families have no scalar bijector).
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
    y_partialmiss = sim(partialmiss, n)
    y_transobs = sim(transobs, n)
    y_normalobs = sim(normalobs, n)
    y_split = sim(split, n)
    y_recordexp = sim(recordexp, n)
    y_prefixmods = sim(prefixmods, n)

    # --- coupled / stratified 'patch' models (the patch-models tutorial) -----
    # Stratified latent seams (2-D `Dims{2}` shapes): a shared `RandomWalk` with
    # partially pooled per-stratum deviations (`Stratify`), and fully independent
    # per-stratum paths (`Replicate`).
    strat_lat = as_turing_model(
        Stratify(RandomWalk(), Hierarchy(; across = IID(Normal(0.0, 0.3)))),
        (3, 8))
    repl_lat = as_turing_model(Replicate(RandomWalk()), (3, 8))
    # A fixed diagonally-dominant off-diagonal coupling matrix and a `Gravity`
    # operator with inferred exponent priors, the two `mixing` choices.
    patch_K = [1.0 0.1 0.05
               0.1 1.0 0.1
               0.05 0.1 1.0]
    patch_pop = [1.0e5, 5.0e4, 2.0e5]
    patch_dist = [0.0 10.0 30.0
                  10.0 0.0 25.0
                  30.0 25.0 0.0]
    rt_process = Stratify(RandomWalk(),
        Hierarchy(; across = IID(Normal(0.0, 0.3))))
    # Stratified `Renewal` posteriors over a 3-stratum x 12-time panel:
    # independent patches, shared/partially-pooled `R_t`, fixed `K` mixing,
    # `Gravity` mixing, and stratified exogenous `ImportedCases`.
    strat_panel = (3, 12)
    patch_independent = IDModel(
        Renewal(; generation_time = gen_int, rt = Replicate(RandomWalk()),
            initialisation = Normal(log(50.0), 0.2)),
        PoissonError())
    patch_pooled = IDModel(
        Renewal(; generation_time = gen_int, rt = rt_process,
            initialisation = Normal(log(50.0), 0.2)),
        PoissonError())
    patch_mixk = IDModel(
        Renewal(; generation_time = gen_int, rt = rt_process,
            initialisation = Normal(log(50.0), 0.2), mixing = patch_K),
        PoissonError())
    patch_mixg = IDModel(
        Renewal(; generation_time = gen_int, rt = rt_process,
            initialisation = Normal(log(50.0), 0.2),
            mixing = Gravity(patch_pop, patch_dist; α = HalfNormal(1.0),
                β = HalfNormal(1.0), γ = HalfNormal(2.0))),
        PoissonError())
    patch_imports = IDModel(
        Renewal(
            gen_int, ImportedCases(
                Stratify(FixedIntercept(-2.0), IID(Normal(0.0, 0.3))));
            rt = rt_process, initialisation = Normal(log(50.0), 0.2)),
        PoissonError())
    y_independent = sim(patch_independent, strat_panel)
    y_pooled = sim(patch_pooled, strat_panel)
    y_mixk = sim(patch_mixk, strat_panel)
    y_mixg = sim(patch_mixg, strat_panel)
    y_imports = sim(patch_imports, strat_panel)

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
            as_turing_model(split, y_split, n)),
        # coupled / stratified 'patch' models (the patch-models tutorial)
        ("Stratify latent logjoint", strat_lat),
        ("Replicate latent logjoint", repl_lat),
        ("Renewal+IndependentPatches posterior",
            as_turing_model(patch_independent, y_independent, strat_panel)),
        ("Renewal+StratifiedRt posterior",
            as_turing_model(patch_pooled, y_pooled, strat_panel)),
        ("Renewal+FixedMixingK posterior",
            as_turing_model(patch_mixk, y_mixk, strat_panel)),
        ("Renewal+GravityMixing posterior",
            as_turing_model(patch_mixg, y_mixg, strat_panel)),
        ("Renewal+StratifiedImportedCases posterior",
            as_turing_model(patch_imports, y_imports, strat_panel))
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
tuples: ForwardDiff (the reference), ReverseDiff (compiled), ReverseDiff (tape),
Mooncake reverse, Mooncake forward, Enzyme reverse, and Enzyme forward — the
full seven-backend matrix `ad.yaml` runs in CI. Per-backend brokenness is
recorded honestly in [`backend_broken_scenarios`](@ref) /
[`broken_scenario_names`](@ref) rather than by trimming this list.
"""
function backends()
    return [
        (name = "ForwardDiff", backend = _forwarddiff()),
        (name = "ReverseDiff (compiled)", backend = _reversediff_compiled()),
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
function _reversediff_compiled()
    ADTypes = Base.require(Base.PkgId(
        Base.UUID("47edcb42-4c32-4615-8424-f2b9edc5f35b"), "ADTypes"))
    return ADTypes.AutoReverseDiff(; compile = true)
end
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

Evidence below is measured locally (`DI.gradient` vs the ForwardDiff
reference, `rtol = 5e-2, atol = 1e-6`); the harness's full sweep additionally
covers the prepared/in-place forms.

| Backend                | Result                          |
|------------------------|---------------------------------|
| ForwardDiff            | reference; not listed broken    |
| ReverseDiff (tape)     | not listed broken               |
| ReverseDiff (compiled) | all 43 listed (DI batch-sweep gap) |
| Enzyme reverse         | 39 pass / 4 broken             |
| Enzyme forward         | 39 pass / 4 broken             |

Both modes share the same 4 remaining broken scenarios: the stratified patch
posters (`Renewal+IndependentPatches`, `Renewal+StratifiedRt`, `Renewal+
StratifiedImportedCases`) and `Renewal+GravityMixing`. They fail with
`IllegalTypeAnalysisException` (a Union through the Stratify/Replicate renewal
path); Gravity additionally surfaces a `BigFloat`/`rewrite_union_returns_as_ref`
issue. Root cause not yet established.

Previously fixed along the way: the four AR scenarios (type-stabilised
scan/diff seams, 47411fb); the LatentDelay / UncertainLatentDelay /
TimeVaryingLatentDelay posters and `PartiallyMissing` / `Split cascade` (making
observed data concrete at inference time); and `Renewal+FixedMixingK` plus the
delay-convolution scenarios (dot-free reductions).

Unresolved causes are tracked in #97.
"""
function backend_broken_scenarios()
    # Enzyme reverse now has one remaining broken scenario.
    #
    # The five scenarios that shared a union-typed value through the prior
    # seam (`as_turing_submodel`'s `filldist`/`product_distribution` runtime
    # ternary) were fixed by 375eee3, which unconditionally returns
    # `product_distribution` and removes the `Union{Product, FillDist}`
    # return type that Enzyme's type analysis could not resolve.
    #
    # `Renewal+TimeVaryingLatentDelay posterior` was measured passing after
    # the `_at_all` fix (0e740ef) and is now removed from the broken list.
    # MA latent logjoint was fixed by the priors.jl filldist ternary removal
    # (375eee3) — the docstring identified the prior-seam union as the exact
    # cause. It now passes under Enzyme reverse (21 tests).
    #
    # The four AR-related scenarios (ARIMA, DiffLatentModel(RW), ARMA,
    # AR vector-prior) still fail with `IllegalTypeAnalysisException` under
    # Enzyme reverse despite the same fix. The common thread is the
    # `accumulate_scan(ARStep(...), ...)` / `DiffLatentModel._combine_diff`
    # path, which differs from the `MA`-side `accumulate_scan(MAStep(...), ...)`
    # in that it composes further accumulation steps (AR inside ARMA inside
    # DiffLatentModel). The remaining union type is likely introduced by
    # `vcat`/`cumsum` in `_combine_diff` or by the `%` branch in
    # `combine_correct` — the root cause is not yet established.
    enzyme_reverse = Set([
        # Stratified/mixing patch posters: Union through the renewal path
        # (Gravity also surfaces a `BigFloat`/`rewrite_union_returns_as_ref`
        # issue). Not yet root-caused.
        "Renewal+IndependentPatches posterior",
        "Renewal+StratifiedRt posterior",
        "Renewal+StratifiedImportedCases posterior",
        "Renewal+GravityMixing posterior"
    ])
    # Enzyme forward still has 12 broken scenarios.
    #
    # The four AR-related scenarios (ARIMA, DiffLatentModel(RW), ARMA,
    # AR vector-prior) fail with `IllegalTypeAnalysisException` (same cause
    # as Enzyme reverse — the `accumulate_scan(ARStep(...), ...)` path).
    #
    # The eight renewal/delay-convolution scenarios fail with
    # `EnzymeNoDerivativeError` from the BLAS `dot` call in
    # `renewal_foi`/`LDStep`/`TimeVaryingLDStep`. The `EnzymeRules.forward`
    # rule for `LinearAlgebra.dot` (375eee3) is loaded and works for simple
    # cases, but Enzyme inlines the BLAS `cblas_ddot64_` call before the
    # rule can intercept it, so the fallback BLAS replacement fires instead.
    # This is a known gap — either the rule needs a narrower type signature
    # that Enzyme can match before inlining, or the model code needs to avoid
    # `dot` on BLAS-visible vector lengths (the `mapreduce` alternative used
    # in `ARStep`/`MAStep` regressed ReverseDiff — see the revert in
    # f69939f — but a fresh approach may work).
    #
    # `Renewal+ImportedCases posterior` is also listed as broken: 12 of 21
    # differentiation tests fail (the `pairwise_gen_int` gravity-coupling
    # path likely hits the same `dot` BLAS gap, plus the dense `*` matrix
    # operations in `renewal_pressure` that Enzyme does support but the
    # combined gradient path triggers a `BoundsError` in).
    enzyme_forward = Set([
        # Stratified/mixing patch posters (same as reverse): Union through the
        # Stratify/Replicate/bernoulli renewal path; Gravity also BigFloat.
        "Renewal+IndependentPatches posterior",
        "Renewal+StratifiedRt posterior",
        "Renewal+StratifiedImportedCases posterior",
        "Renewal+GravityMixing posterior"
    ])
    # ReverseDiff compiled (`AutoReverseDiff(; compile = true)`) differentiates
    # the DynamicPPL log-densities correctly: a plain `DI.gradient` call and an
    # isolated per-scenario `DifferentiationInterfaceTest.test_differentiation`
    # run both pass all 43 scenarios, matching the ForwardDiff reference. It is
    # still listed broken below because this DifferentiationInterface/
    # DifferentiationInterfaceTest version cannot run the harness's batched
    # `test_differentiation` correctness sweep for `AutoReverseDiff{true}`: the
    # `PullbackFast` path it selects has no `_prepare_pullback_aux` method for
    # these closure types (only `PullbackSlow` is implemented). That is a
    # DifferentiationInterface integration gap, not a ReverseDiff or DynamicPPL
    # limitation — with ReverseDiff properly loaded the compiled mode works, as
    # it does in packages that drive it outside DI's batched sweep. Listing all
    # 43 keeps the harness green (over-listed entries still record an ordinary
    # `@test` pass via `check_broken`).
    reversediff_compiled = Set([
        "RandomWalk latent logjoint",
        "AR latent logjoint",
        "ARIMA latent logjoint",
        "HilbertSpaceGP latent logjoint",
        "HilbertSpaceGP Matern latent logjoint",
        "ExactGP latent logjoint",
        "MA latent logjoint",
        "HierarchicalNormal latent logjoint",
        "DiffLatentModel(RandomWalk) latent logjoint",
        "ARMA latent logjoint",
        "BroadcastLatentModel day-of-week latent logjoint",
        "BroadcastLatentModel weekly latent logjoint",
        "ConcatLatentModels latent logjoint",
        "CombineLatentModels latent logjoint",
        "Hierarchy latent logjoint",
        "AR vector-prior latent logjoint",
        "AR latent-model-as-prior latent logjoint",
        "DirectInfections+Poisson posterior",
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
        "DirectInfections+PartiallyMissing posterior",
        "DirectInfections+TransformObservation posterior",
        "DirectInfections+NormalError posterior",
        "DirectInfections+RecordExpected posterior",
        "DirectInfections+PrefixModifiers posterior",
        "BinomialError ascertainment posterior",
        "Renewal+Split cascade posterior",
        "Stratify latent logjoint",
        "Replicate latent logjoint",
        "Renewal+IndependentPatches posterior",
        "Renewal+StratifiedRt posterior",
        "Renewal+FixedMixingK posterior",
        "Renewal+GravityMixing posterior",
        "Renewal+StratifiedImportedCases posterior"
    ])
    return Dict{String, Set{String}}(
        "Enzyme reverse" => enzyme_reverse,
        "Enzyme forward" => enzyme_forward,
        "ReverseDiff (compiled)" => reversediff_compiled)
end

"Per-backend scenario names too unstable to even run (segfault/hang)."
backend_skip_scenarios() = Dict{String, Set{String}}()

end # module ADFixtures
