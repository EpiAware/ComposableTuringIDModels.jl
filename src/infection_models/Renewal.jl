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
[`ModelShape`](@ref) `n` and returns the named tuple `(; I_t, Z_t)` with `Z_t`
the (log) ``R_t`` path.

Renewal is the one infection model that needs a generation interval, so it takes
one directly through the `generation_time` keyword, which dispatches on the value:
a discrete probability vector is used as-is, a continuous `Distribution` is
discretised internally (see the constructor), and a pmf-producing prior model
(e.g. an [`UncertainDelay`](@ref)) lets the generation interval itself be
**inferred** — its distribution's parameters carry priors and the interval is
rediscretised per draw through the [`as_turing_submodel`](@ref) seam.

`Renewal` is a step-composing helper: positional [`AbstractRenewalModifier`](@ref)
arguments are composed onto the renewal [`RenewalStep`](@ref). Passing a
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
  - `initialisation`: prior for the unconstrained initial infections (a
    `Distribution` or prior model, sampled through [`as_turing_submodel`](@ref)).
  - `recurrent_step`: the renewal accumulation step (an
    [`AbstractConstantRenewalStep`](@ref)), or `nothing` when the generation
    interval is inferred and the step is built per draw.
  - `mixing`: the coupling operator applied to the incidence window before
    `R_t` (default `I`, uncoupled); see [`renewal_pressure`](@ref).

## Constructor

  - `Renewal(; generation_time, rt, initialisation, transformation = exp,
    mixing = I, D_gen = nothing, Δd = 1.0)` — one keyword constructor that
    dispatches on `generation_time`:

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

A stratified renewal: one recursion per stratum, sharing a random walk in
`R_t`-space with partially pooled per-stratum deviations, left uncoupled
(`mixing` defaults to `I`):

```@example Renewal
strat = Renewal(; generation_time = [0.2, 0.3, 0.5],
    rt = Stratify(RandomWalk(), Hierarchy()), initialisation = Normal())
size(as_turing_model(strat, (3, 20))().I_t)
```
"
struct Renewal{G, F <: Function, L <: PriorLike, S <: PriorLike, A, K} <:
    AbstractInfectionModel
    "Discrete generation interval, or a pmf-producing prior model when inferred."
    gen_int::G
    "Transformation between unconstrained and constrained domains."
    transformation::F
    "Latent process model generating the (log) reproduction number."
    rt::L
    "Prior for the unconstrained initial infections."
    initialisation::S
    "The renewal accumulation step (`nothing` when the interval is inferred)."
    recurrent_step::A
    "The coupling operator applied to the incidence window before `R_t`."
    mixing::K
end

function Renewal(;
        generation_time, rt = RandomWalk(),
        initialisation = Normal(), transformation::Function = exp,
        mixing = I, D_gen = nothing, Δd = 1.0
    )
    gen_int, recurrent_step = _renewal_fields(
        generation_time, mixing; D_gen = D_gen, Δd = Δd
    )
    return Renewal(
        gen_int, transformation, path_prior(rt), initialisation,
        recurrent_step, mixing
    )
end

# Positional modifier constructor: a discrete generation interval (a non-negative
# pmf that sums to 1) with positional [`AbstractRenewalModifier`](@ref)s composed
# onto the renewal step — e.g.
# `Renewal([0.2, 0.3, 0.5], SusceptibleDepletion(1000.0); rt = RandomWalk())`.
function Renewal(
        gen_int::AbstractVector,
        modifiers::AbstractRenewalModifier...;
        rt = RandomWalk(), mixing = I,
        initialisation = Normal(), transformation::Function = exp
    )
    @assert all(gen_int .>= 0) "Generation interval must be non-negative"
    @assert sum(gen_int) ≈ 1 "Generation interval must sum to 1"
    core = _renewal_step(gen_int, mixing)
    recurrent_step = RenewalStep(core, modifiers)
    return Renewal(
        gen_int, transformation, path_prior(rt), initialisation,
        recurrent_step, mixing
    )
end

# Fixed generation interval (a pmf vector or a continuous distribution): bake the
# discretised interval and its renewal step (with `mixing` folded in) at
# construction, exactly as before.
function _renewal_fields(generation_time, mixing; D_gen = nothing, Δd = 1.0)
    gen_int = _renewal_gen_int(generation_time; D_gen = D_gen, Δd = Δd)
    return gen_int, _renewal_step(gen_int, mixing)
end

# Inferred generation interval: hold the pmf-producing prior model and build the
# renewal step per draw inside the `@model`, so no interval or step is baked.
# `mixing` is still folded in there (it lives on the struct, not the step).
function _renewal_fields(
        generation_time::AbstractPriorModel, mixing;
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

# `generation_time` as a continuous distribution: discretise via double-interval
# censoring, drop the delay-0 bin (a generation interval has no mass at lag 0) and
# renormalise.
function _renewal_gen_int(
        gen_distribution::ContinuousDistribution;
        D_gen = nothing, Δd = 1.0
    )
    return _discretised_pmf(gen_distribution; Δd = Δd, D = D_gen) |>
        p -> p[2:end] ./ sum(p[2:end])
end

# The generation-interval shape a `gen_int` prior slot is drawn at: no shape
# argument (giving one pmf) for a single series, `n_strata` (giving one pmf per
# stratum) for a stratified renewal. Returned as a tuple so it splats straight
# into the `as_turing_submodel` call — `()` splats to no extra argument.
_gen_int_shape(n::Int) = ()
_gen_int_shape(n::Dims{2}) = (n[1],)

# Drop the lag-0 bin and renormalise, the generation-interval convention the
# fixed-distribution path applies at construction, applied here to a drawn pmf
# (or one pmf per stratum). `_stack_pmfs` assembles a vector of per-stratum pmfs
# into a `strata × lags` matrix before the matrix method renormalises each row.
_drop_lag_zero(p::AbstractVector) = p[2:end] ./ sum(p[2:end])
function _drop_lag_zero(ps::AbstractVector{<:AbstractVector})
    return _drop_lag_zero(_stack_pmfs(ps))
end
function _drop_lag_zero(P::AbstractMatrix)
    Q = P[:, 2:end]
    return Q ./ sum(Q; dims = 2)
end

_stack_pmfs(ps) = permutedims(reduce(hcat, ps))

# Initial renewal state from sampled I₀ and R₀, decaying at the implied rate. The
# generation interval and its step are passed in so the fixed (baked) and
# inferred (per-draw) paths share one initialiser; `_init_rate`/`_n_lags` widen
# it to a per-stratum seed, `R_t` and/or generation interval.
function _make_renewal_init(step::AbstractConstantRenewalStep, gen_int, I₀, Rt₀)
    r_approx = _init_rate(Rt₀, gen_int)
    return renewal_init_state(step, I₀, r_approx, _n_lags(gen_int))
end

@model function as_turing_model(infection::Renewal, n::ModelShape)
    Z_t ~ as_turing_submodel(infection.rt, n)
    init_incidence ~ as_turing_submodel(
        infection.initialisation, _n_strata(n); prefix = true
    )
    I₀ = infection.transformation.(_seed(init_incidence, n))
    Rt = infection.transformation.(Z_t)

    # The generation interval is either fixed (its baked renewal step used
    # directly) or inferred: a pmf-producing prior model sampled through the
    # single seam at the shape the renewal needs, with the lag-0 bin dropped
    # and the remainder renormalised per draw. The step is rebuilt per draw so
    # the gradient flows through the discretisation, and `mixing` is folded in
    # as the constructors fold it into the baked step.
    if infection.gen_int isa AbstractPriorModel
        gen ~ as_turing_submodel(
            infection.gen_int, _gen_int_shape(n)...; prefix = true
        )
        gen_int = _drop_lag_zero(gen)
        step = _renewal_step(gen_int, infection.mixing)
    else
        gen_int = infection.gen_int
        step = infection.recurrent_step
    end

    # Resolve the step before scanning: any modifier carrying priors (e.g.
    # `ImportedCases`) draws them here and hands back the modifier the scan
    # uses, and a drawn `mixing` model (e.g. `Gravity`) draws its parameters and
    # hands back a resolved core, while a purely deterministic step returns
    # itself. One seam, no branch on what the step's modifiers or mixing are.
    scan_step ~ as_turing_submodel(step, n)

    Rts = _steps(Rt)
    init = _make_renewal_init(scan_step, gen_int, I₀, first(Rts))
    I_t = accumulate_scan(scan_step, init, Rts)
    return (; I_t, Z_t)
end
