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
The first sections build and inspect each pattern; the last two simulate a
known truth, fit it under NUTS, and check what actually comes back out.

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
pop = [1.0e5, 5.0e4, 2.0e5]
dist = [0.0 10.0 30.0
        10.0 0.0 25.0
        30.0 25.0 0.0]

movement = Gravity(pop, dist; α = HalfNormal(1.0), β = HalfNormal(1.0),
    γ = HalfNormal(2.0))
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

`pop` is in whatever units you have.
`gravity` works in units of the mean population and returns rows that sum to
one, so `K` says only where each patch's force comes from.
`R_t` carries its size.
Without that split the two are confounded, since any overall scale on `K` is
the same as scaling `R_t`.

`within` sets how much of a patch's force is its own, relative to a typical
pairwise term.
It has to be non-zero for `α` to be identifiable.

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

## Does the coupling actually recover? Simulate and fit

Every pattern above only ever ran forward: build a model, simulate or check a
shape, move on.
None of it was fitted, so none of it answers whether a coupled, partially
pooled patch process is actually **identifiable** from data.
This section closes that gap for the fixed-`K` pattern above: simulate from a
known truth, fit under NUTS, and check whether the three patches' `R_t`
paths come back out.

The model reuses the shared, partially pooled `R_t` process and the `K` from
above, with two changes purely so the fit finishes in reasonable time for a
docs build: a tighter `RandomWalk` innovation (so `R_t` drifts slowly rather
than wandering) and a narrower prior on the initial level (so simulated case
counts stay in the hundreds rather than the millions).

```@example patches
using Turing, Statistics
Random.seed!(20260817)

n_time = 24
K_fit = [1.0 0.15 0.05
         0.15 1.0 0.1
         0.05 0.1 1.0]
rt_fit = Stratify(
    RandomWalk(init = Normal(log(1.15), 0.05),
        ϵ_t = HierarchicalNormal(std = HalfNormal(0.03))),
    Hierarchy(; mean = Normal(0.0, 0.1), across = IID(Normal(0.0, 0.2))))
coupled = Renewal(gen_int; rt = rt_fit,
    initialisation = Normal(log(50.0), 0.2), mixing = K_fit)
coupled_model = IDModel(
    coupled, Split(NegativeBinomialError(cluster_factor = HalfNormal(0.1))))

Ymiss = Matrix{Union{Missing, Float64}}(missing, n_strata, n_time)
sim = as_turing_model(coupled_model, Ymiss)()
Rt_true = exp.(sim.Z_t)
(Rt_range = extrema(Rt_true),
    y_range = extrema(reduce(vcat, collect(values(sim.generated_y_t)))))
```

`R_t` stays between roughly 1.1 and 1.36 across the three patches: sustained
but modest growth, giving case counts that reach a few thousand by day 24
rather than the runaway millions a wider prior produces.
Conditioning on the simulated counts and sampling recovers the posterior:

```@example patches
Ydata = Float64.(reduce(vcat,
    [permutedims(sim.generated_y_t[Symbol(:group, g)]) for g in 1:n_strata]))
posterior = as_turing_model(coupled_model, Ydata)
chain = sample(
    posterior, NUTS(0.85; adtype = Turing.AutoForwardDiff()), 300;
    progress = false)
nothing # hide
```

`Z_t` is a generated quantity per patch, exactly as in the [Renewal
tutorial](@ref tutorial-renewal), so the fitted `R_t` path is recovered per
draw with [`Turing.returned`](https://turinglang.org/) and stacked into
credible bands:

```@example patches
using Turing: returned
draws = vec(returned(posterior, chain))
Rt_stack = cat((exp.(d.Z_t) for d in draws)...; dims = 3)
Rt_med = dropdims(median(Rt_stack; dims = 3); dims = 3)
Rt_lo = dropdims(mapslices(x -> quantile(x, 0.1), Rt_stack; dims = 3); dims = 3)
Rt_hi = dropdims(mapslices(x -> quantile(x, 0.9), Rt_stack; dims = 3); dims = 3)
(rmse = round(sqrt(mean((Rt_true .- Rt_med) .^ 2)); digits = 4),
    coverage80 = round(
        mean((Rt_true .>= Rt_lo) .& (Rt_true .<= Rt_hi)); digits = 3),
    true_range = extrema(Rt_true), posterior_median_range = extrema(Rt_med))
```

```@example patches
using CairoMakie
fig = Figure(; size = (760, 640))
for g in 1:n_strata
    ax = Axis(fig[g, 1]; ylabel = "Patch $g  Rₜ",
        xlabel = g == n_strata ? "Day" : "")
    band!(ax, 1:n_time, Rt_lo[g, :], Rt_hi[g, :]; color = (:purple, 0.25))
    lines!(ax, 1:n_time, Rt_med[g, :]; color = :purple, linewidth = 2,
        label = "posterior median")
    lines!(ax, 1:n_time, Rt_true[g, :]; color = :black, linestyle = :dash,
        linewidth = 2, label = "simulated truth")
    hlines!(ax, [1.0]; color = :grey, linestyle = :dot)
    g == 1 && axislegend(ax; position = :lt)
end
fig
```

The posterior median tracks each patch's true `R_t` path closely, and the
80% credible interval covers the truth in 76% of the 72 patch-day cells,
close to the nominal rate.
A coupled, partially pooled renewal process — `Stratify`, `Hierarchy`, and an
off-diagonal `mixing` matrix composed together — is identifiable from data at
this signal-to-noise level: the machinery in the sections above is not just
constructible, it fits.

## Does the gravity model recover its own coupling?

The `K` above was fixed and known.
[`Gravity`](@ref) instead *infers* `K`'s shape from three exponents, `α`,
`β`, and `γ`, so the natural next question is whether those exponents
themselves come back out of a fit — a stronger claim than merely running
without error.

The same `pop` and `dist` from the gravity example above are reused, with
narrower priors on the exponents so the truth is a value the small panel
below has a realistic chance of pinning down, and the same tight `R_t`
process as the coupled fit above.

```@example patches
movement_fit = Gravity(pop, dist; α = HalfNormal(0.5), β = HalfNormal(0.5),
    γ = HalfNormal(1.0))
grav_rt = Stratify(
    RandomWalk(init = Normal(log(1.08), 0.03),
        ϵ_t = HierarchicalNormal(std = HalfNormal(0.02))),
    Hierarchy(; mean = Normal(0.0, 0.05), across = IID(Normal(0.0, 0.15))))
gravity_fit_model = Renewal(gen_int; rt = grav_rt,
    initialisation = Normal(log(50.0), 0.2), mixing = movement_fit)
grav_model = IDModel(
    gravity_fit_model,
    Split(NegativeBinomialError(cluster_factor = HalfNormal(0.1))))
nothing # hide
```

Unlike the fixed-`K` fit, the truth here is itself a set of *sampled*
parameters, so a single draw from the prior — rather than the bare
`as_turing_model(...)()` call used above — is what fixes a reproducible
truth alongside its exponent values:

```@example patches
Random.seed!(20260817)
n_time_g = 20
Ymiss_g = Matrix{Union{Missing, Float64}}(missing, n_strata, n_time_g)
grav_prior = as_turing_model(grav_model, Ymiss_g)
truth_chain = sample(grav_prior, Prior(), 1; progress = false)
truth = only(vec(returned(grav_prior, truth_chain)))
Ydata_g = Float64.(reduce(vcat,
    [permutedims(truth.generated_y_t[Symbol(:group, g)]) for g in 1:n_strata]))
(y_range = extrema(reduce(vcat, collect(values(truth.generated_y_t)))),
    Rt_range = extrema(exp.(truth.Z_t)))
```

Fitting proceeds exactly as before, and the exponents' variables are
namespaced under `core.mixing`, as [`MixingStep`](@ref) documents:

```@example patches
grav_posterior = as_turing_model(grav_model, Ydata_g)
grav_chain = sample(
    grav_posterior, NUTS(0.85; adtype = Turing.AutoForwardDiff()), 300;
    progress = false)

recovery = map((:α, :β, :γ)) do sym
    label = "Parameter(core.mixing.$(sym))"
    vn = only(filter(v -> string(v) == label, collect(keys(grav_chain))))
    post = vec(grav_chain[vn])
    true_val = only(vec(truth_chain[vn]))
    (parameter = sym, true_value = round(true_val; digits = 3),
        posterior_mean = round(mean(post); digits = 3),
        q10 = round(quantile(post, 0.1); digits = 3),
        q90 = round(quantile(post, 0.9); digits = 3))
end
recovery
```

`γ`, the distance exponent, recovers well: the truth sits inside its 80%
interval and close to the posterior mean.
`β` is roughly captured too, if with a wide interval.
`α` is not: the posterior mean overshoots the truth and the 80% interval
misses it entirely.
With only three patches there are only three off-diagonal entries of `K` to
learn from, and `α` (destination population) and `β` (origin population)
trade off against each other, and against `γ`, over that little data — the
gravity form is *constructible* and *runs* cleanly (no overflow, no
divergences beyond the odd one in 300 draws), but three strata is not enough
data to identify all three exponents separately.
A real application would need more strata, an informative prior on one
exponent, or holding one fixed (`β = 1`, the standard gravity form) to
identify the rest cleanly.
