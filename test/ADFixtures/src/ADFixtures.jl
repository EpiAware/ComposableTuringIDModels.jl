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

    # Coalesce warm-up `missing` entries into a concrete observed vector. We do
    # not infer missing events, so conditioning on a concrete `Vector{Int}` /
    # `Vector{Float64}` keeps the gradient path off a `Union{Missing, T}`
    # eltype. Non-array outputs (`ReportingTriangle`) pass through.
    _coalesce_warmup(x::AbstractArray) = eltype(x) >: Missing ?
                                         coalesce.(x, zero(nonmissingtype(eltype(x)))) : x
    _coalesce_warmup(x::NamedTuple) = map(_coalesce_warmup, x)
    _coalesce_warmup(x) = x
    # Simulate observations from a composed model's prior. The seed is per
    # scenario so every backend process conditions on identical data.
    sim(m, nn, seed) = _coalesce_warmup(
        as_turing_model(m, missing, nn)(MersenneTwister(seed)).generated_y_t)

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

    # --- prior slots -----------------------------------------------------------
    # A vector of damping distributions: one i.i.d. draw per lag.
    ar_vec = as_turing_model(
        AR(; damp = [truncated(Normal(0, 0.05), 0, 1),
                truncated(Normal(0, 0.05), 0, 1)],
            init = [Normal(), Normal()]), 8)
    # A process as the damping prior: a time-varying coefficient path drawn as a
    # `RandomWalk` submodel and mapped through `tanh`.
    ar_lat = as_turing_model(AR(; damp = RandomWalk()), 8)

    # --- infection posteriors ---------------------------------------------------
    direct = IDModel(
        DirectInfections(; Z = RandomWalk(), initialisation = Normal()),
        PoissonError())
    renewal = IDModel(
        Renewal(; generation_time = gen_int, rt = RandomWalk(), initialisation = Normal()),
        NegativeBinomialError())
    # Both renewal modifier kinds at once: `SusceptibleDepletion` samples
    # nothing, `ImportedCases` draws its rate in the pre-scan seam.
    modifiers = IDModel(
        Renewal(gen_int, SusceptibleDepletion(2_000.0),
            ImportedCases(Normal(0.0, 0.5));
            rt = RandomWalk(), initialisation = Normal()),
        NegativeBinomialError())
    # Exponential-growth-rate infections: a cumulative growth-rate path
    # exponentiated.
    egr = IDModel(
        ExpGrowthRate(; rt = RandomWalk(), initialisation = Normal()),
        PoissonError())

    # Nowcasting, marginal: right-truncation correction from a reporting-delay
    # distribution, scaling the inner error.
    nowcast = IDModel(
        Renewal(; generation_time = gen_int, rt = RandomWalk(), initialisation = Normal()),
        RightTruncate(NegativeBinomialError(),
            truncated(Normal(4.0, 1.5), 0.0, Inf)))

    # Nowcasting, joint: a per-cell Poisson likelihood over the masked
    # reporting triangle (`t + d ≤ now`).
    triangle = IDModel(
        Renewal(; generation_time = gen_int, rt = RandomWalk(), initialisation = Normal()),
        ReportTriangle(PoissonError(), [0.6, 0.25, 0.15]))

    # --- observation modifiers / error families over a composed model ----------
    # Reporting delay: convolves the expected observations with a delay PMF
    # (`accumulate_scan(LDStep(rev_pmf), ...)`) before the inner error.
    latdelay = IDModel(
        Renewal(; generation_time = gen_int, rt = RandomWalk(), initialisation = Normal()),
        LatentDelay(NegativeBinomialError(), [0.3, 0.4, 0.3]))
    # Uncertain reporting delay: the delay parameters are prior slots, so the
    # pmf is rediscretised per draw. Differentiates through `_discretised_pmf`.
    udelay = IDModel(
        Renewal(; generation_time = gen_int, rt = RandomWalk(), initialisation = Normal()),
        LatentDelay(NegativeBinomialError(),
            UncertainDelay(LogNormal,
                [Normal(1.0, 0.3), truncated(Normal(0.4, 0.2), 0, Inf)];
                D = 6.0)))
    # Time-varying reporting delay: the meanlog is a `RandomWalk`, so the pmf is
    # rebuilt per time step and convolved through `TimeVaryingLDStep`.
    tvdelay = IDModel(
        Renewal(; generation_time = gen_int, rt = RandomWalk(), initialisation = Normal()),
        LatentDelay(NegativeBinomialError(),
            UncertainDelay(LogNormal,
                [RandomWalk(), truncated(Normal(0.4, 0.2), 0, Inf)];
                D = 6.0)))
    # Uncertain generation interval: inferred and rediscretised per draw (lag-0
    # bin dropped, renormalised) before the renewal step is built.
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
    # Partially-missing observations: only the blank entries are latent. The
    # error must be continuous, since a count family has no scalar bijector to
    # link them through.
    partialmiss = IDModel(
        DirectInfections(; Z = RandomWalk(), initialisation = Normal()),
        NormalError())
    # Transform-the-expected-observations: softplus applied before the error.
    transobs = IDModel(
        DirectInfections(; Z = RandomWalk(), initialisation = Normal()),
        TransformObservationModel(PoissonError()))
    # The same model fully observed: the minimal non-count likelihood.
    normalobs = IDModel(
        DirectInfections(; Z = RandomWalk(), initialisation = Normal()),
        NormalError())

    # `:=`-tracked pass-throughs on both the latent and the observation side.
    recordexp = IDModel(
        DirectInfections(;
            Z = RecordExpectedLatent(RandomWalk()), initialisation = Normal()),
        RecordExpectedObs(PoissonError()))
    # `PrefixLatentModel` / `PrefixObservationModel` as standalone wrappers,
    # rather than the prefixing `Split` and the prior seam do inline.
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
    y_binom = _coalesce_warmup(
        as_turing_model(binom_obs, (y = missing, N = N_b), Ybase_b)(
        MersenneTwister(124)).y_t)

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

    y_direct = sim(direct, n, 101)
    y_renewal = sim(renewal, n, 102)
    y_modifiers = sim(modifiers, n, 103)
    y_egr = sim(egr, n, 104)
    y_nowcast = sim(nowcast, n, 105)
    y_triangle = sim(triangle, n, 106)
    y_latdelay = sim(latdelay, n, 107)
    y_udelay = sim(udelay, n, 108)
    y_tvdelay = sim(tvdelay, n, 109)
    y_ugen = sim(ugen, n, 110)
    y_ascert = sim(ascert, 14, 111)
    y_aggregate = sim(aggregate, 14, 112)
    # Blank two entries so the scenario really is partially missing.
    y_partialmiss = Vector{Union{Missing, Float64}}(sim(partialmiss, n, 113))
    y_partialmiss[[3, 7]] .= missing
    y_transobs = sim(transobs, n, 114)
    y_normalobs = sim(normalobs, n, 115)
    y_split = sim(split, n, 116)
    y_recordexp = sim(recordexp, n, 117)
    y_prefixmods = sim(prefixmods, n, 118)

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
    y_independent = sim(patch_independent, strat_panel, 119)
    y_pooled = sim(patch_pooled, strat_panel, 120)
    y_mixk = sim(patch_mixk, strat_panel, 121)
    y_mixg = sim(patch_mixg, strat_panel, 122)
    y_imports = sim(patch_imports, strat_panel, 123)

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
        # prior slots
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
full seven-backend matrix `ad.yaml` runs in CI. Failures are recorded in
[`backend_broken_scenarios`](@ref) rather than by trimming this list.
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
    # derivative data, and without it Enzyme raises `EnzymeMutabilityException`
    # ("argument cannot be proven readonly") on every DynamicPPL log-density.
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

"Scenario names broken on every backend."
broken_scenario_names() = String[]

@doc """
    backend_broken_scenarios()

Per-backend broken scenario names (`Dict{String, Set{String}}`).

The harness's `check_broken` computes the pass/fail boolean itself and only
falls back to `@test_broken` when a listed scenario actually fails, so
over-listing is safe and under-listing reds CI. Listing every scenario for a
backend is not safe in a different sense: it empties the set the full
correctness sweep runs over, so that backend stops being tested.

Add an entry only with a measured failure.
"""
function backend_broken_scenarios()
    # DynamicPPL `deepcopy`s a model argument whose type admits `Missing`,
    # before the model body runs. Enzyme's generic reverse path differentiates
    # that copy without consulting the `inactive` rule the package ships for
    # it, and raises `FieldError: type Const has no field dval`. Forward mode
    # takes the rule, so it and the other five backends differentiate this.
    return Dict{String, Set{String}}(
        "Enzyme reverse" => Set(["DirectInfections+PartiallyMissing posterior"]))
end

"Per-backend scenario names too unstable to even run (segfault/hang)."
backend_skip_scenarios() = Dict{String, Set{String}}()

end # module ADFixtures
