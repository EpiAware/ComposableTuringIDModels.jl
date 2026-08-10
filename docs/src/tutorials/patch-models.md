# [Coupled patch models](@id tutorial-patches)

A spatial or age-stratified epidemic runs several patches over the same time
axis.
Patches can be independent, share a process with partial pooling, or exchange
infection pressure with each other.
[`Stratify`](@ref) puts the patch axis on a single process.
[`Renewal`](@ref)'s `mixing` slot couples the patches that axis produces.
[`CombineInfections`](@ref) draws several genuinely separate processes
instead.
One model over an axis is couplable.
Several independent models are not.

Every pattern below shares the same time axis and generation interval.

```@example patches
using ComposableTuringIDModels, Distributions, Random
Random.seed!(2026)

n_strata, n_time = 3, 25
gen_int = [0.2, 0.3, 0.3, 0.2]
```

## Independent patches

[`Replicate`](@ref) draws one process per patch, each independent over the
time axis.
It stacks the results into the `strata × time` matrix a renewal process needs.

```@example patches
independent = Renewal(gen_int; rt = Replicate(RandomWalk()),
    initialisation = Normal(log(50.0), 0.2))
size(as_turing_model(independent, (n_strata, n_time))().I_t)
```

Nothing here shares information across patches.
There is no common driver and no infection pressure moving between them.

## A shared R_t, partially pooled

[`Stratify`](@ref) draws one shared process over time and a per-patch
deviation.
It combines the two with `+` by default, a multiplicative effect once
exponentiated.

```@example patches
rt_process = Stratify(RandomWalk(),
    Hierarchy(; across = IID(Normal(0.0, 0.3))))
pooled = Renewal(gen_int; rt = rt_process,
    initialisation = Normal(log(50.0), 0.2))
size(as_turing_model(pooled, (n_strata, n_time))().I_t)
```

Every patch's `R_t` moves with the shared random walk.
The [`Hierarchy`](@ref) shrinks each patch's deviation toward the shared
level.
An `IID` `across` slot gives no pooling.
`FixedIntercept(0.0)` collapses every patch onto the identical `R_t`.

## Cross-patch imports: exogenous and endogenous

Infections can arrive at a patch from outside the modelled system.
They can also arrive from another patch inside it.
[`ImportedCases`](@ref) with a stratified rate covers the first case.
Its rate is drawn before the scan runs.
It never depletes a susceptible pool.
It never comes from another patch's own transmission chain.

```@example patches
exogenous = Renewal(gen_int,
    ImportedCases(Stratify(FixedIntercept(-2.0), IID(Normal(0.0, 0.3))));
    rt = rt_process, initialisation = Normal(log(50.0), 0.2))
nothing # hide
```

An off-diagonal `mixing` matrix covers the second case.
Pressure generated inside one patch's own renewal recursion reaches another
patch directly.

```@example patches
K = [1.0 0.1 0.05
     0.1 1.0 0.1
     0.05 0.1 1.0]
endogenous = Renewal(gen_int; rt = rt_process,
    initialisation = Normal(log(50.0), 0.2), mixing = K)
size(as_turing_model(endogenous, (n_strata, n_time))().I_t)
```

## A fixed mixing matrix

Every renewal step convolves each patch's own incidence history with the
generation interval.
`K` then mixes the resulting pressures across patches.
The mixed pressure scales that patch's `R_t`:

```math
\Lambda_{g,t} = \mathcal R_{g,t} \sum_h K_{gh} \sum_i g_i I_{h,t-i}
```

`K`'s diagonal is each patch's own weight.
Its off-diagonal is imported pressure from other patches.
`K = I`, [`Renewal`](@ref)'s default, is the uncoupled case above.
The dispatch point is [`renewal_pressure`](@ref), a plain function.
A new mixing structure needs only a new method.

## A gravity model with inferred exponents

A fixed `K` is one choice.
[`Gravity`](@ref) instead infers the mixing weights from population sizes and
pairwise distances.

```@example patches
pop = [1.0, 0.5, 2.0]   # hundreds of thousands
dist = [0.0 10.0 30.0
        10.0 0.0 25.0
        30.0 25.0 0.0]

movement = Gravity(pop, dist; α = HalfNormal(1.0), β = HalfNormal(1.0),
    γ = HalfNormal(2.0), normalise = true)
gravity_model = Renewal(gen_int; rt = rt_process,
    initialisation = Normal(log(50.0), 0.2), mixing = movement)
size(as_turing_model(gravity_model, (n_strata, n_time))().I_t)
```

[`gravity`](@ref) is the plain function behind it.
`K[g, h]` for `g ≠ h` is the pressure patch `h` puts on patch `g`, relative to
`g`'s own weight.
That own weight is the `within` keyword, default `1.0`.
`Gravity` draws priors on `α`, `β` and `γ`, then calls the same function
inside the model, so the fixed and inferred paths cannot drift.

The units of `pop` matter here.
`gravity` is unnormalised, so `pop^α` carries its scale straight into `K`.
Raw counts of `1e5` give off-diagonal weights around `1e7` against a `within`
diagonal of `1.0`, and the recursion overflows.
Scaling `pop` keeps both sides comparable, which is why it is in hundreds of
thousands above.

`normalise = true` then makes each row sum to one, so the operator conserves
the total force of infection.
It is not a substitute for scaling.
Normalising raw counts bounds the total but leaves the diagonal at around
`1e-7`, which says a patch infects itself essentially never.

## Many-to-one observation: age bands into one stream

A hospitalisation stream rarely reports by age even when transmission does.
[`Split`](@ref) with a weight matrix aggregates infection strata into fewer
observation streams.
Here three age bands feed one stream, weighted by each band's
hospitalisation rate.

```@example patches
W = [0.05 0.15 0.6]                     # young, middle, elderly
hosp = Split(NegativeBinomialError(cluster_factor = HalfNormal(0.1)), W)
hosp_model = IDModel(pooled, hosp)

Ymiss = Matrix{Union{Missing, Float64}}(missing, 1, n_time)
sim = as_turing_model(hosp_model, Ymiss)()
sim.generated_y_t
```

The `1 × n_time` data matrix carries no strata count.
`IDModel`'s two-argument form reads the infection strata straight off the
observation model's weight matrix.
The same three-patch process built earlier assembles itself from the data's
shape alone.
