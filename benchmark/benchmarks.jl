# PACKAGE-OWNED — scaffold writes this once and never overwrites it.
#
# Benchmark suite for ComposableTuringIDModels. Defines a BenchmarkTools `BenchmarkGroup`
# named `SUITE` that the managed `run.jl` / `compare.jl` consume.
#
# Groups:
#   "Model evaluation" — building + evaluating representative models: a prior
#       draw (`rand`) and the generated-quantities forward pass (`model()`) for
#       latent processes (AR, RandomWalk) and composed `IDModel`s
#       (DirectInfections+Poisson, Renewal+NegativeBinomial).
#   "Sampling"         — a short NUTS run on a composed model.
#   "AD gradients"     — gradient of a representative log-density across AD
#       backends, keyed `["AD gradients"][scenario][backend]` so `compare.jl`
#       folds it into a per-(scenario × backend) matrix.

using BenchmarkTools
using ComposableTuringIDModels
using Distributions
using Random: MersenneTwister
using Turing: NUTS, sample
using DynamicPPL: LogDensityFunction, VarInfo, link, getlogjoint
import LogDensityProblems as LDP
import DifferentiationInterface as DI
using ADTypes: AutoForwardDiff
import ForwardDiff, ReverseDiff, Mooncake, Enzyme
using ADTypes: AutoReverseDiff, AutoMooncake, AutoEnzyme

const SUITE = BenchmarkGroup()

# --- shared fixtures --------------------------------------------------------

const GEN_INT = [0.2, 0.3, 0.5]
const N = 30

# The representative models, each as a Turing model ready to evaluate. Posterior
# models are conditioned on data simulated from the prior with a fixed seed.
function _eval_models()
    direct = IDModel(
        DirectInfections(; Z = RandomWalk(), initialisation = Normal()),
        PoissonError())
    renewal = IDModel(
        Renewal(GEN_INT; rt = RandomWalk(), initialisation = Normal()),
        NegativeBinomialError())
    y_direct = as_turing_model(direct, missing, N)().generated_y_t
    y_renewal = as_turing_model(renewal, missing, N)().generated_y_t
    return [
        ("AR latent", as_turing_model(AR(), N)),
        ("RandomWalk latent", as_turing_model(RandomWalk(), N)),
        ("DirectInfections+Poisson", as_turing_model(direct, y_direct, N)),
        ("Renewal+NegativeBinomial", as_turing_model(renewal, y_renewal, N))
    ]
end

# --- Model evaluation -------------------------------------------------------

let eval_grp = SUITE["Model evaluation"] = BenchmarkGroup()
    for (name, model) in _eval_models()
        eval_grp[name] = BenchmarkGroup()
        # A prior draw (samples every random variable).
        eval_grp[name]["rand"] = @benchmarkable rand($model)
        # The forward pass returning the generated quantities.
        eval_grp[name]["forward"] = @benchmarkable $model()
    end
end

# --- Sampling ---------------------------------------------------------------

let samp_grp = SUITE["Sampling"] = BenchmarkGroup()
    model = IDModel(
        DirectInfections(; Z = RandomWalk(), initialisation = Normal()),
        PoissonError())
    y = as_turing_model(model, missing, N)().generated_y_t
    cond = as_turing_model(model, y, N)
    # A short NUTS run; `seconds` in run.jl caps wall time, and a low draw count
    # keeps this representative without dominating the suite.
    samp_grp["NUTS (DirectInfections+Poisson, 50 draws)"] = @benchmarkable sample(
        $cond, NUTS(), 50; progress = false)
end

# --- AD gradients -----------------------------------------------------------

# A real differentiable log-density: the linked log-joint over unconstrained ℝᵈ.
function _logdensity(model; seed::Int = 1)
    ldf = LogDensityFunction(model, getlogjoint, link(VarInfo(model), model))
    θ0 = 0.3 .* randn(MersenneTwister(seed), LDP.dimension(ldf))
    return (θ -> LDP.logdensity(ldf, θ)), θ0
end

# Backends. Enzyme needs `function_annotation = Const` (the closures carry no
# derivative data); the AR-based log-densities are excluded from Enzyme as they
# hit a known Enzyme type-analysis limitation (see test/ADFixtures).
const _ENZYME = AutoEnzyme(;
    mode = Enzyme.set_runtime_activity(Enzyme.Reverse),
    function_annotation = Enzyme.Const)

# (backend display name, backend, runs_on_AR-based-scenarios?)
const _AD_BACKENDS = [
    ("ForwardDiff", AutoForwardDiff(), true),
    ("ReverseDiff (tape)", AutoReverseDiff(; compile = false), true),
    ("Mooncake reverse", AutoMooncake(; config = nothing), true),
    ("Enzyme reverse", _ENZYME, false)
]

# (scenario name, model, is_AR_based?)
function _ad_scenarios()
    direct = IDModel(
        DirectInfections(; Z = RandomWalk(), initialisation = Normal()),
        PoissonError())
    y = as_turing_model(direct, missing, N)().generated_y_t
    return [
        ("AR latent logjoint", as_turing_model(AR(), N), true),
        ("DirectInfections+Poisson posterior",
            as_turing_model(direct, y, N), false)
    ]
end

let ad_grp = SUITE["AD gradients"] = BenchmarkGroup()
    for (i, (sname, model, ar_based)) in enumerate(_ad_scenarios())
        f, θ0 = _logdensity(model; seed = i)
        ad_grp[sname] = BenchmarkGroup()
        for (bname, backend, runs_ar) in _AD_BACKENDS
            (ar_based && !runs_ar) && continue
            # `prep` is built once; the benchmark times the gradient itself.
            prep = DI.prepare_gradient(f, backend, θ0)
            ad_grp[sname][bname] = @benchmarkable DI.gradient(
                $f, $prep, $backend, $θ0)
        end
    end
end
