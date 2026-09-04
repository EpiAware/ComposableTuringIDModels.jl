# The centred stochastic renewal process.
#
# A centred draw needs `I_t ~ dist(ι_t, …)` with `ι_t` known during the
# recursion, and a scan step is a deterministic function. So this is an
# infection model rather than a renewal modifier. It runs the recursion as a
# `@model` loop over the propose/commit halves of a step, drawing in between.

@doc raw"
A renewal process whose infections are drawn around the renewal expectation
rather than fixed at it, **centred**.

```math
Z_t \sim \text{latent}, \qquad
\mathcal R_t = g(Z_t), \qquad
\iota_t = \mathcal R_t \sum_{i=1}^{n-1} I_{t-i} g_i, \qquad
I_t \sim F(\iota_t, c_t \iota_t)
```

``F`` and the coefficient of variation ``c_t`` come from the `noise` specification ([`InfectionNoise`](@ref)), matching a negative binomial's first two moments.
Everything else is [`Renewal`](@ref)'s, and the two take the same arguments.

Infections are the sampled parameter, so the likelihood informs them directly.
The same noise applied through the renewal modifier seam can only be non-centred, which is the worse parameterisation when the data are informative about the latent path.
Reach for this first, and for the modifier only when a non-centred draw is wanted or when the renewal is stratified.

The default `LogNormal` noise family keeps infections positive, which is where the centred draw's support constraint comes from.
A `Normal` family has no such bound.

Single-series only.
A stratified renewal needs one draw per stratum per step, which this does not provide.

## Fields

  - `gen_int`: the discrete generation interval, or the pmf-producing prior
    model when it is inferred.
  - `transformation`: the transformation between the unconstrained and
    constrained domains (default `exp`).
  - `rt`: the latent process model generating the (log) reproduction number.
  - `initialisation`: prior for the unconstrained infections at ``t_0``, or a
    [`SeedingPath`](@ref).
  - `recurrent_step`: the renewal accumulation step (`nothing` when the
    generation interval is inferred).
  - `mixing`: the coupling operator, which is `I` because this is single-series.
  - `modifiers`: the renewal modifiers composed onto the step.
  - `noise`: the [`InfectionNoise`](@ref) specification.

## Constructors

Both mirror [`Renewal`](@ref)'s, with a `noise` keyword added:

  - `StochasticRenewal(; generation_time, rt, initialisation, noise, …)`
  - `StochasticRenewal(generation_time, modifiers...; noise, …)`

# Examples

Infections are stochastic, so two draws differ even with ``R_t`` and the
seeding held fixed:

```@example StochasticRenewal
using ComposableTuringIDModels, Distributions, Random
sr = StochasticRenewal(
    [0.2, 0.3, 0.5]; rt = FixedIntercept(0.0),
    initialisation = Normal(log(500.0), 0.0)
)
a = as_turing_model(sr, 12)(Xoshiro(1)).I_t
b = as_turing_model(sr, 12)(Xoshiro(2)).I_t
all(>(0), a), isapprox(a, b)
```

Infections carry their own variable in the chain:

```@example StochasticRenewal
keys(rand(Xoshiro(1), as_turing_model(sr, 12)))
```
"
struct StochasticRenewal{
        G, F <: Function, L <: PriorLike, S <: PriorLike, A, K, M <: Tuple, N,
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
    "The infection-noise specification."
    noise::N

    # `rt` is widened here rather than in the keyword constructor, as
    # `Renewal` widens its own, so no construction path can bypass it.
    function StochasticRenewal(
            gen_int::G, transformation::F, rt, initialisation::S,
            recurrent_step::A, mixing::K, modifiers::M, noise::N
        ) where {G, F <: Function, S <: PriorLike, A, K, M <: Tuple, N}
        path = path_prior(rt)
        return new{G, F, typeof(path), S, A, K, M, N}(
            gen_int, transformation, path, initialisation, recurrent_step,
            mixing, modifiers, noise
        )
    end
end

function StochasticRenewal(;
        generation_time, modifiers = (), rt = RandomWalk(),
        initialisation = Normal(), transformation::Function = exp,
        mixing = I, D_gen = nothing, Δd = 1.0, noise = InfectionNoise()
    )
    mods = _modifier_tuple(modifiers)
    gen_int, recurrent_step = _renewal_fields(
        generation_time, mixing, mods; D_gen = D_gen, Δd = Δd
    )
    return StochasticRenewal(
        gen_int, transformation, rt, initialisation, recurrent_step,
        mixing, mods, noise
    )
end

# Positional modifier form, delegating to the keyword one exactly as
# `Renewal`'s does, so there is one construction path.
function StochasticRenewal(
        generation_time::Union{AbstractVector, Distribution, AbstractPriorModel},
        modifiers::AbstractRenewalModifier...;
        rt = RandomWalk(), mixing = I, initialisation = Normal(),
        transformation::Function = exp, D_gen = nothing, Δd = 1.0,
        noise = InfectionNoise()
    )
    return StochasticRenewal(;
        generation_time = generation_time, modifiers = modifiers, rt = rt,
        mixing = mixing, initialisation = initialisation,
        transformation = transformation, D_gen = D_gen, Δd = Δd, noise = noise
    )
end

@model function as_turing_model(infection::StochasticRenewal, n::Int)
    setup ~ to_submodel(_renewal_setup(infection, n), false)
    ξ ~ to_submodel(_noise_overdispersion(infection.noise, n), false)

    scan_step, init, Rts = setup.scan_step, setup.init, setup.Rts

    # The recursion, written out because each step's distribution depends on
    # the step before it. `I_t` is typed from the state it starts in, so the
    # container holds whatever the AD backend in use puts in it.
    state = init
    I_t = Vector{typeof(init.val * first(Rts))}(undef, n)
    for t in 1:n
        ι, _, substates = _propose(scan_step, state, Rts[t])
        I_t[t] ~ _noise_dist(infection.noise, ι, at(ξ, t))
        state = _commit(scan_step, state, I_t[t], substates)
    end
    return (; I_t, Z_t = setup.Z_t, I_seed = init.window)
end

# The recursion again, keeping each step's expectation. The loop is written
# twice rather than always building the series and discarding it, so a model
# nobody has wrapped in a `RecordExpectedInfections` allocates nothing for a
# quantity it never reports.
@model function with_expected_infections(infection::StochasticRenewal, n::Int)
    setup ~ to_submodel(_renewal_setup(infection, n), false)
    ξ ~ to_submodel(_noise_overdispersion(infection.noise, n), false)

    scan_step, init, Rts = setup.scan_step, setup.init, setup.Rts

    state = init
    T = typeof(init.val * first(Rts))
    I_t = Vector{T}(undef, n)
    exp_I_t = Vector{T}(undef, n)
    for t in 1:n
        ι, _, substates = _propose(scan_step, state, Rts[t])
        exp_I_t[t] = ι
        I_t[t] ~ _noise_dist(infection.noise, ι, at(ξ, t))
        state = _commit(scan_step, state, I_t[t], substates)
    end
    return (; I_t, Z_t = setup.Z_t, I_seed = init.window, exp_I_t)
end

# One draw per step, and a stratified renewal would need one per stratum per
# step. Say so rather than seeding a strata axis the loop cannot fill.
function as_turing_model(infection::StochasticRenewal, n::Dims{2})
    return error(
        "`StochasticRenewal` runs one series. A stratified renewal needs " *
            "one draw per stratum at every step, and what that draw should " *
            "be is an open question: " *
            "https://github.com/EpiAware/ComposableTuringIDModels.jl/issues/355" *
            ". Use `Renewal` with an `InfectionNoise` modifier for a " *
            "stratified stochastic renewal, at the cost of a non-centred " *
            "parameterisation."
    )
end
