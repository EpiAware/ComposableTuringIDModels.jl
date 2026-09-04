# Renewal infection process model.

@doc raw"
Model unobserved infections via a time-varying renewal process driven by an
internally generated (log) reproduction number.

```math
Z_t \sim \text{latent}, \qquad
\mathcal R_t = g(Z_t), \qquad
I_t = \mathcal R_t \sum_{i=1}^{n-1} I_{t-i} g_i
```

where the latent model `rt` supplies the (log) reproduction number ``Z_t``, ``g``
is `transformation`, ``g_i`` is the discrete generation interval, and the
pre-window infections decay at the growth rate implied by ``\mathcal R_1``. The
``R_t`` process is generated *inside* the model, so `as_turing_model` takes a
[`ModelShape`](@ref) `n` and returns the named tuple `(; I_t, Z_t, I_seed)` with
`Z_t` the (log) ``R_t`` path and `I_seed` the seeding window the scan started
from, which is not part of `I_t`.

The seeding window is deterministic unless `initialisation` is wrapped in a
[`SeedingPath`](@ref), which estimates the whole run-up instead of decaying a
single level.

Renewal is the one infection model that needs a generation interval, so it takes
one directly through the `generation_time` keyword, which dispatches on the value:
a discrete probability vector is used as-is, a continuous `Distribution` is
discretised internally (see the constructor), and a pmf-producing prior model
(e.g. an [`UncertainDelay`](@ref)) lets the generation interval itself be
**inferred** — its distribution's parameters carry priors and the interval is
rediscretised per draw through the [`as_turing_submodel`](@ref) seam.

`Renewal` is a step-composing helper: [`AbstractRenewalModifier`](@ref)s, given
positionally or through the `modifiers` keyword, are composed onto the renewal
[`RenewalStep`](@ref). Either form takes any generation time and the
discretisation keywords, so a modifier and a continuous generation time compose:
`Renewal(Gamma(2, 1.5), SusceptibleDepletion(N); D_gen = 15.0)`. Passing a
[`SusceptibleDepletion`](@ref)`(N)` gives a renewal process with a fixed
population ``N`` and susceptible depletion
```math
I_t = \frac{S_{t-1}}{N} \mathcal R_t \sum_{i=1}^{n-1} I_{t-i} g_i, \qquad
S_t = S_{t-1} - I_t.
```
Modifiers apply in the order given, and a modifier carrying priors (e.g. an
[`ImportedCases`](@ref) importation rate) draws them once before the scan, so
composing one takes no extra wiring:
`Renewal(gen_int, SusceptibleDepletion(N), ImportedCases(Normal(-1, 0.5)))`.

## Strata and coupling

Handing `rt` a [`Stratify`](@ref) (or [`Replicate`](@ref)) and calling
`as_turing_model(renewal, (n_strata, n_time))` runs one renewal recursion per
stratum. `R_t` is a `n_strata × n_time` matrix, each stratum has its own seed
(from a vector-valued `initialisation`, or the same scalar seed broadcast to
every stratum), and each stratum's incidence window advances on its own. The
`mixing` keyword couples the strata by transforming the incidence window each
stratum's force of infection is convolved against — see
[`renewal_pressure`](@ref) for the extension point. It defaults to `I`
(`LinearAlgebra.I`), leaving strata uncoupled: a `Dims{2}` shape with no mixing
is `n_strata` independent single-series renewal processes sharing one model.

The generation interval, when inferred, is drawn at the shape the renewal needs:
one pmf for `n::Int`, one pmf per stratum for `n::Dims{2}`. Because each
parameter of an [`UncertainDelay`](@ref) is itself a prior slot, giving one a
[`Hierarchy`](@ref) partially pools that parameter across strata.
A partially pooled generation interval therefore needs no new component, e.g.
`UncertainDelay(LogNormal, [Hierarchy(), σ_prior]; D = 14.0)` drawn at
`n_strata`.

## Fields

  - `gen_int`: the discrete generation interval vector (non-negative, sums to 1),
    or, for an inferred generation interval, the pmf-producing prior model.
  - `transformation`: the transformation between the unconstrained and
    constrained domains (default `exp`).
  - `rt`: the latent process model (an [`AbstractLatentModel`](@ref)) generating
    the (log) reproduction number. A length-`n` PATH slot: a bare `Distribution`
    here is auto-wrapped in an [`Intercept`](@ref), giving a constant path (one
    shared draw broadcast to length `n`); use [`IID`](@ref) for `n` independent
    draws, or [`Stratify`](@ref)/[`Replicate`](@ref) for a strata axis.
  - `initialisation`: prior for the unconstrained infections at ``t_0`` — a
    LEVEL, one value or one per stratum (a `Distribution` or prior model,
    sampled through [`as_turing_submodel`](@ref)). Wrap it in a
    [`SeedingPath`](@ref) to estimate the whole seeding window instead.
  - `recurrent_step`: the renewal accumulation step (an
    [`AbstractConstantRenewalStep`](@ref)), or `nothing` when the generation
    interval is inferred and the step is built per draw.
  - `mixing`: the coupling operator applied to the incidence window before
    `R_t` (default `I`, uncoupled); see [`renewal_pressure`](@ref).
  - `modifiers`: the [`AbstractRenewalModifier`](@ref)s composed onto the step.
    Held here as well as in `recurrent_step` so the step can be rebuilt per
    draw when the generation interval is inferred.

## Constructor

  - `Renewal(; generation_time, modifiers = (), rt, initialisation,
    transformation = exp, mixing = I, D_gen = nothing, Δd = 1.0)` — one keyword
    constructor that dispatches on `generation_time`:

      + a discrete probability **vector** (non-negative, sums to 1) is used
        directly as the generation interval;
      + a continuous **`Distribution`** is discretised via double-interval
        censoring (CensoredDistributions.jl), using `D_gen`/`Δd`, with the
        delay-0 bin dropped and the remainder renormalised; and
      + a pmf-producing **prior model** (an [`AbstractPriorModel`](@ref) such as
        an [`UncertainDelay`](@ref)) is held as-is and sampled per draw, giving an
        **inferred** generation interval (an uncertain discretised distribution
        used as the generation interval). Its fixed horizon keeps the interval
        length constant across draws; the lag-0 bin is dropped and the remainder
        renormalised per draw, exactly as for the fixed distribution.

  - `Renewal(generation_time, modifiers...; rt, initialisation, ...)` — the same
    constructor with the generation time and the modifiers given positionally.

# Examples

A fixed generation interval:

```@example Renewal
using ComposableTuringIDModels, Distributions
renewal = Renewal(; generation_time = [0.2, 0.3, 0.5], rt = RandomWalk(),
    initialisation = Normal())
rand(as_turing_model(renewal, 20))
```

An inferred generation interval — an uncertain discretised distribution used as
the generation interval, whose `LogNormal` parameters carry priors:

```@example Renewal
gen = UncertainDelay(
    LogNormal, [Normal(1.9, 0.2), truncated(Normal(0.5, 0.2), 0, Inf)]; D = 14.0)
renewal = Renewal(; generation_time = gen, rt = RandomWalk(),
    initialisation = Normal())
rand(as_turing_model(renewal, 20))

# With a fixed population and susceptible depletion.
depleting = Renewal([0.2, 0.3, 0.5], SusceptibleDepletion(1000.0);
    rt = RandomWalk(), initialisation = Normal())
rand(as_turing_model(depleting, 20))
```

A continuous generation time discretised by the constructor, carrying a
modifier:

```@example Renewal
discretised = Renewal(Gamma(2, 1.5), SusceptibleDepletion(1000.0);
    D_gen = 15.0, rt = RandomWalk(), initialisation = Normal(log(50), 0.2))
length(discretised.gen_int)
```

An estimated seeding window rather than a decaying one, through
[`SeedingPath`](@ref):

```@example Renewal
seeded = Renewal(; generation_time = [0.2, 0.3, 0.5], rt = RandomWalk(),
    initialisation = SeedingPath(RandomWalk(; init = Normal(log(50), 0.5))))
as_turing_model(seeded, 20)().I_seed
```

A stratified renewal: one recursion per stratum, sharing a random walk in
`R_t`-space with partially pooled per-stratum deviations, left uncoupled
(`mixing` defaults to `I`):

```@example Renewal
strat = Renewal(; generation_time = [0.2, 0.3, 0.5],
    rt = Stratify(RandomWalk(), Hierarchy()), initialisation = Normal())
size(as_turing_model(strat, (3, 20))().I_t)
```
"
struct Renewal{
        G, F <: Function, L <: PriorLike, S <: PriorLike, A, K, M <: Tuple,
    } <: AbstractInfectionModel
    "Discrete generation interval, or a pmf-producing prior model when inferred."
    gen_int::G
    "Transformation between unconstrained and constrained domains."
    transformation::F
    "Latent process model generating the (log) reproduction number."
    rt::L
    "Prior for the unconstrained infections at ``t_0``, or a `SeedingPath`."
    initialisation::S
    "The renewal accumulation step (`nothing` when the interval is inferred)."
    recurrent_step::A
    "The coupling operator applied to the incidence window before `R_t`."
    mixing::K
    "The renewal modifiers composed onto the step."
    modifiers::M

    function Renewal(
            gen_int::G, transformation::F, rt, initialisation::S,
            recurrent_step::A, mixing::K, modifiers::M
        ) where {G, F <: Function, S <: PriorLike, A, K, M <: Tuple}
        # `rt` is a length-`n` path slot, so a bare `Distribution` is wrapped in
        # an `Intercept` rather than left as a scalar.
        # The widening runs here rather than in the keyword constructor so the
        # positional form cannot bypass it.
        # `path_prior` is idempotent, so rebuilding a `Renewal` from its own
        # fields changes nothing.
        path = path_prior(rt)
        return new{G, F, typeof(path), S, A, K, M}(
            gen_int, transformation, path, initialisation, recurrent_step,
            mixing, modifiers
        )
    end
end

function Renewal(;
        generation_time, modifiers = (), rt = RandomWalk(),
        initialisation = Normal(), transformation::Function = exp,
        mixing = I, D_gen = nothing, Δd = 1.0
    )
    mods = _modifier_tuple(modifiers)
    gen_int, recurrent_step = _renewal_fields(
        generation_time, mixing, mods; D_gen = D_gen, Δd = Δd
    )
    return Renewal(
        gen_int, transformation, rt, initialisation,
        recurrent_step, mixing, mods
    )
end

# Positional modifier form, taking the generation time first and then the
# [`AbstractRenewalModifier`](@ref)s composed onto the renewal step.
# It delegates to the keyword form, so there is one construction path.
function Renewal(
        generation_time::Union{AbstractVector, Distribution, AbstractPriorModel},
        modifiers::AbstractRenewalModifier...;
        rt = RandomWalk(), mixing = I, initialisation = Normal(),
        transformation::Function = exp, D_gen = nothing, Δd = 1.0
    )
    return Renewal(;
        generation_time = generation_time, modifiers = modifiers, rt = rt,
        mixing = mixing, initialisation = initialisation,
        transformation = transformation, D_gen = D_gen, Δd = Δd
    )
end

# The modifier slot takes one modifier or a collection of them, and a tuple keeps
# the step's modifier types concrete.
# The keyword form has no signature to constrain it, so what it was given is
# checked here rather than failing as a `MethodError` from inside the step's seam.
_modifier_tuple(mod::AbstractRenewalModifier) = (mod,)

function _modifier_tuple(mods)
    tup = Tuple(mods)
    all(m -> m isa AbstractRenewalModifier, tup) || throw(
        ArgumentError(
            "`modifiers` takes `AbstractRenewalModifier`s, but was given " *
                "$(join(nameof.(typeof.(tup)), ", "))."
        )
    )
    return tup
end

# The step a core and a modifier tuple make.
# With no modifiers this is the bare core, so a modifier-free renewal keeps the
# step, and so the variable names, it has always had.
_bake_step(core, mods::Tuple) = RenewalStep(core, mods)
_bake_step(core, ::Tuple{}) = core

# A fixed generation interval bakes the discretised interval and its renewal
# step, with `mixing` and the modifiers folded in, at construction.
function _renewal_fields(
        generation_time, mixing, mods; D_gen = nothing, Δd = 1.0
    )
    gen_int = _renewal_gen_int(generation_time; D_gen = D_gen, Δd = Δd)
    return gen_int, _bake_step(_renewal_step(gen_int, mixing), mods)
end

# An inferred generation interval holds the pmf-producing prior model and builds
# the renewal step per draw inside the `@model`, so no interval or step is baked.
# `mixing` and the modifiers live on the struct and are folded in there instead.
function _renewal_fields(
        generation_time::AbstractPriorModel, mixing, mods;
        D_gen = nothing, Δd = 1.0
    )
    return generation_time, nothing
end

# `generation_time` as a discrete PMF: use it directly (must be a valid pmf).
function _renewal_gen_int(gen_int::AbstractVector; D_gen = nothing, Δd = 1.0)
    @assert all(gen_int .>= 0) "Generation interval must be non-negative"
    @assert sum(gen_int) ≈ 1 "Generation interval must sum to 1"
    return collect(gen_int)
end

# A continuous `generation_time` is discretised by double-interval censoring.
# The delay-0 bin is dropped, because a generation interval has no mass at lag 0,
# and the remainder renormalised.
function _renewal_gen_int(
        gen_distribution::ContinuousDistribution;
        D_gen = nothing, Δd = 1.0
    )
    return _discretised_pmf(gen_distribution; Δd = Δd, D = D_gen) |>
        p -> p[2:end] ./ sum(p[2:end])
end

# The shape a `gen_int` prior slot is drawn at.
# A single series takes no shape argument and a stratified renewal takes
# `n_strata`, giving one pmf per stratum.
# Returned as a tuple so it splats straight into the `as_turing_submodel` call.
_gen_int_shape(n::Int) = ()
_gen_int_shape(n::Dims{2}) = (n[1],)

# Drop the lag-0 bin and renormalise, which is the convention the
# fixed-distribution path applies at construction, applied here to a drawn pmf.
# `_stack_pmfs` assembles a vector of per-stratum pmfs into a `strata × lags`
# matrix before the matrix method renormalises each row.
_drop_lag_zero(p::AbstractVector) = p[2:end] ./ sum(p[2:end])
function _drop_lag_zero(ps::AbstractVector{<:AbstractVector})
    return _drop_lag_zero(_stack_pmfs(ps))
end
function _drop_lag_zero(P::AbstractMatrix)
    Q = P[:, 2:end]
    return Q ./ sum(Q; dims = 2)
end

_stack_pmfs(ps) = permutedims(reduce(hcat, ps))

# Initial renewal state from sampled I₀ and R₀, decaying at the implied rate.
# The generation interval and its step are passed in so the baked and per-draw
# paths share one initialiser.
function _make_renewal_init(step::AbstractConstantRenewalStep, gen_int, I₀, Rt₀)
    r_approx = _init_rate(Rt₀, gen_int)
    return renewal_init_state(step, I₀, r_approx, _n_lags(gen_int))
end

@model function as_turing_model(infection::Renewal, n::ModelShape)
    setup ~ to_submodel(_renewal_setup(infection, n), false)
    I_t = accumulate_scan(setup.scan_step, setup.init, setup.Rts)
    return (; I_t, Z_t = setup.Z_t, I_seed = setup.init.window)
end

# The scan runs on a `RecordingRenewalStep`, whose state keeps each step's
# expectation, so both series come off one scan. `accumulate_scan` is unrolled
# because the accumulated states it drops are what the second series is read
# from.
@model function with_expected_infections(infection::Renewal, n::ModelShape)
    setup ~ to_submodel(_renewal_setup(infection, n), false)
    step = RecordingRenewalStep(setup.scan_step)
    init = _recording_state(setup.init)
    states = accumulate(step, setup.Rts; init = init)
    I_t = get_state(step, init, states)
    exp_I_t = get_expected_state(step, init, states)
    return (; I_t, Z_t = setup.Z_t, I_seed = setup.init.window, exp_I_t)
end
