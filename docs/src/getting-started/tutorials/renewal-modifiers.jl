md"""
# [Renewal modifiers](@id renewal-modifiers)

A [`Renewal`](@ref) process can be extended by *modifiers* that change how the next incidence is formed.
This page takes one delayed renewal process and composes three of them onto it.
[`SusceptibleDepletion`](@ref) bounds the epidemic by a finite population.
[`ImportedCases`](@ref) seeds infections from outside it.
[`InfectionNoise`](@ref) gives infections a distribution of their own rather than fixing them at the renewal expectation.
Each is added against the same baseline, so its contribution is visible on its own.

Everything else is held fixed.
The generation interval, the reproduction number, the reporting delay and the initial incidence are the same in every run, and only the modifier list changes.
The last two sections fit models carrying modifiers back to simulated data.
"""

using ComposableTuringIDModels, Distributions, Random, Turing, Mooncake
using CairoMakie
using ADTypes: AutoMooncake
using DynamicPPL: fix
using Statistics: quantile

Random.seed!(20260727)
nothing # hide

md"""
## A delayed renewal process

The infection process is a renewal equation with a discretised generation
interval and a constant ``R_t = 1.3``.
Fixing ``\log R_t`` with a [`FixedIntercept`](@ref) keeps the comparison clean:
any difference between the runs below is the modifiers' doing, not a different
draw of the reproduction number.
"""

gen_int = [0.2, 0.3, 0.3, 0.2]
n = 60
pop_size = 2_000.0

function renewal(modifiers...)
    return Renewal(
        gen_int, modifiers...;
        rt = FixedIntercept(log(1.3)), initialisation = Normal()
    )
end
nothing # hide

md"""
Reported cases are the infections convolved with a reporting delay and observed
with negative-binomial noise — a [`LatentDelay`](@ref) wrapped around a
[`NegativeBinomialError`](@ref).
"""

delay = [0.1, 0.4, 0.3, 0.2]
obs = LatentDelay(
    NegativeBinomialError(cluster_factor = HalfNormal(0.1)), delay
)

md"""
## Three models

The three models differ only in the modifiers composed onto the renewal step.
Modifiers apply in the order given, so importation here is added *after*
depletion: imported infections are not scaled by the susceptible fraction and
do not themselves deplete the pool.
"""

models = (
    plain = IDModel(renewal(), obs),
    depleting = IDModel(renewal(SusceptibleDepletion(pop_size)), obs),
    seeded = IDModel(
        renewal(
            SusceptibleDepletion(pop_size),
            ImportedCases(FixedIntercept(log(2.0)))
        ),
        obs
    ),
)
nothing # hide

md"""
Simulating from each model with `missing` observations and the same fixed
initial incidence returns its infections and reported cases.
"""

simulate(model) = fix(
    as_turing_model(model, missing, n),
    (init_incidence = log(1.0),)
)()

sims = map(simulate, models)
nothing # hide

md"""
## What each modifier does

Both panels use a log scale; the reported counts in the lower panel are floored
at one so that zero-report days stay visible.
The reporting delay leaves the first `length(delay) - 1` days without a reported
count, so the lower panel starts on day `length(delay)`.
`n` is the number of reports, and the infection process is run over the delay's lead-in on top of it, so the two panels share a day axis once that is added.
"""

lead = observation_lead_in(obs)
inf_days = 1:(n + lead)
obs_days = (lead + 1):(n + lead)

fig = Figure(; size = (760, 560))
ax1 = Axis(fig[1, 1]; ylabel = "Infections Iₜ", yscale = log10)
ax2 = Axis(fig[2, 1]; xlabel = "Day", ylabel = "Reported cases", yscale = log10)
colours = (plain = :grey30, depleting = :purple, seeded = :teal)
labels = (
    plain = "renewal", depleting = "+ susceptible depletion",
    seeded = "+ importation",
)
for k in keys(sims)
    lines!(
        ax1, inf_days, sims[k].I_t; color = colours[k], linewidth = 2,
        label = labels[k]
    )
    scatter!(
        ax2, obs_days, max.(sims[k].generated_y_t, 1);
        color = colours[k], markersize = 7
    )
end
axislegend(ax1; position = :lt)
fig

md"""
The plain renewal grows exponentially without bound.
Susceptible depletion turns that growth over: as susceptibles are used up the
effective reproduction number falls below one, so incidence peaks and declines.
Importation keeps seeding the epidemic instead of leaving it to grow from the
initial incidence alone, which front-loads it, so it peaks earlier and higher
and spends the susceptible pool sooner.
"""

function summarise(sim)
    return (;
        peak_day = argmax(sim.I_t),
        peak = round(maximum(sim.I_t); digits = 1),
        final = round(last(sim.I_t); digits = 1),
    )
end
map(summarise, (depleting = sims.depleting, seeded = sims.seeded))

md"""
Late incidence in the seeded run is therefore the lower of the two, its epidemic
already past.

Importation also stops a process from dying out altogether, which is the other
thing it is for.
With ``R_t = 0.8`` and no susceptible depletion the plain process decays away,
while a seeded one levels off where importation and decay balance — at
``\iota / (1 - R_t) = 10`` infections per day here.
"""

function subcritical(modifiers...)
    return fix(
        as_turing_model(
            IDModel(
                Renewal(
                    gen_int, modifiers...; rt = FixedIntercept(log(0.8)),
                    initialisation = Normal()
                ),
                obs
            ),
            missing, n
        ),
        (init_incidence = log(50.0),)
    )()
end

sub_plain = subcritical()
sub_seeded = subcritical(ImportedCases(FixedIntercept(log(2.0))))

fig2 = Figure(; size = (760, 300))
ax3 = Axis(
    fig2[1, 1]; xlabel = "Day", ylabel = "Infections Iₜ",
    yscale = log10
)
lines!(
    ax3, inf_days, sub_plain.I_t; color = :grey30, linewidth = 2,
    label = "Rₜ = 0.8"
)
lines!(
    ax3, inf_days, sub_seeded.I_t; color = :teal, linewidth = 2,
    label = "Rₜ = 0.8 + importation"
)
axislegend(ax3; position = :rt)
fig2

md"""
## Stochastic infections

Both modifiers so far are deterministic.
Each transforms the incidence the renewal equation implies, and the result is still one number per day.
[`InfectionNoise`](@ref) instead gives infections a distribution of their own, matched to the first two moments of a negative binomial about the renewal expectation ``\iota_t``.
It composes positionally like the others, so a noisy, depleting, seeded renewal is one more argument.
"""

noise = InfectionNoise(; overdispersion = 0.2)
noisy = IDModel(
    renewal(
        SusceptibleDepletion(pop_size),
        ImportedCases(FixedIntercept(log(2.0))),
        noise
    ),
    obs
)

paths = [simulate(noisy).I_t for _ in 1:20]
nothing # hide

md"""
Each day's draw is centred on the incidence the renewal equation implies from the history it has reached.
Noise moves that day off the deterministic path, and the recursion carries the displacement forward, so the paths separate rather than scattering about a common line.
"""

fig3 = Figure(; size = (760, 300))
ax4 = Axis(
    fig3[1, 1]; xlabel = "Day", ylabel = "Infections Iₜ", yscale = log10
)
for path in paths
    lines!(ax4, inf_days, path; color = (:teal, 0.3), linewidth = 1)
end
lines!(
    ax4, inf_days, sims.seeded.I_t; color = :black, linewidth = 2,
    label = "deterministic"
)
axislegend(ax4; position = :lt)
fig3

md"""
## Centred and non-centred

A modifier is a deterministic function of the incidence reaching it, so as a modifier the noise has to be non-centred.
The sampled parameters are standard normals and the location and scale are applied inside the scan.

[`StochasticRenewal`](@ref) runs the same specification centred.
It draws `I_t` with the expectation in hand, so infections are themselves the sampled parameters and the likelihood informs them directly.
It takes [`Renewal`](@ref)'s arguments with a `noise` keyword added, and the deterministic modifiers compose onto it unchanged.
Reach for it first, and for the modifier when a non-centred draw is wanted or when the renewal is stratified.
"""

stochastic = IDModel(
    StochasticRenewal(
        gen_int,
        SusceptibleDepletion(pop_size),
        ImportedCases(FixedIntercept(log(2.0)));
        rt = FixedIntercept(log(1.3)), initialisation = Normal(),
        noise = noise
    ),
    obs
)
sim_stochastic = simulate(stochastic)
nothing # hide

md"""
## Fitting a model with modifiers

Modifiers are part of the model, so a model carrying them is fitted exactly like
any other.
An [`ImportedCases`](@ref) rate is a prior slot, so it is *estimated* rather than
assumed: here it is a single unknown constant on the log scale, drawn once
before the renewal recursion runs.
We fit that model to the reported cases simulated above.
The reproduction number is still held at its simulated value, so the importation
rate is what is being learned; with ``R_t`` free as well the two compete to
explain the same growth and the fit is much less sharp.
Every simulated report is conditioned on.
The model runs the infection process over the delay's lead-in itself, so nothing is dropped from the head.
"""

model = IDModel(
    renewal(
        SusceptibleDepletion(pop_size),
        ImportedCases(Normal(0.0, 1.0))
    ),
    obs
)
y_obs = sims.seeded.generated_y_t
chain = sample(
    as_turing_model(model, y_obs, n),
    NUTS(0.95; adtype = AutoMooncake(; config = nothing)),
    MCMCThreads(), 250, 2; progress = false
)
nothing # hide

md"""
The rate is namespaced by the modifier's position on the renewal step, so
several modifiers carrying priors can never collide.
Exponentiating the draws puts it back on the scale of infections per day, where
the true value was two.
"""

import_draws = exp.(vec(chain[@varname(modifier_2.import_rates)]))
(
    posterior = round.(quantile(import_draws, [0.05, 0.5, 0.95]), digits = 2),
    truth = 2.0,
)

md"""
The posterior predictive then tracks the simulated series.
"""

pred = predict(as_turing_model(model, fill(missing, n), n), chain)
y_draws(i) = Float64.(vec(pred[@varname(y_t[i])]))
quantiles(i) = quantile(y_draws(i), [0.05, 0.5, 0.95])
## `y_t` is indexed on the reported series, which runs 1:n; `obs_days` puts
## those reports back on the infection process's day axis for the plot.
bands = reduce(hcat, map(quantiles, 1:n))

fig4 = Figure(; size = (760, 300))
ax5 = Axis(
    fig4[1, 1]; xlabel = "Day", ylabel = "Reported cases",
    yscale = log10
)
band!(
    ax5, obs_days, max.(bands[1, :], 1), max.(bands[3, :], 1);
    color = (:teal, 0.25)
)
lines!(
    ax5, obs_days, max.(bands[2, :], 1); color = :teal, linewidth = 2,
    label = "posterior predictive"
)
scatter!(
    ax5, obs_days, max.(y_obs, 1); color = :black,
    markersize = 7, label = "simulated"
)
axislegend(ax5; position = :lt)
fig4

md"""
## Fitting the stochastic renewal

Nothing about the fit changes when the infections are stochastic.
The model is passed to `sample` the same way, and the noise specification is the one used to simulate.
The difference is in what comes back.
Infections are no longer a deterministic function of the other parameters, so they carry a posterior of their own and are recovered by name from the chain.
"""

y_stochastic = sim_stochastic.generated_y_t
chain_stochastic = sample(
    as_turing_model(stochastic, y_stochastic, n),
    NUTS(0.95; adtype = AutoMooncake(; config = nothing)),
    MCMCThreads(), 250, 2; progress = false
)
nothing # hide

infection_draws(i) = Float64.(vec(chain_stochastic[@varname(I_t[i])]))
infection_bands = reduce(
    hcat, [quantile(infection_draws(i), [0.05, 0.5, 0.95]) for i in inf_days]
)

fig5 = Figure(; size = (760, 300))
ax6 = Axis(
    fig5[1, 1]; xlabel = "Day", ylabel = "Infections Iₜ", yscale = log10
)
band!(
    ax6, inf_days, infection_bands[1, :], infection_bands[3, :];
    color = (:teal, 0.25)
)
lines!(
    ax6, inf_days, infection_bands[2, :]; color = :teal, linewidth = 2,
    label = "posterior"
)
lines!(
    ax6, inf_days, sim_stochastic.I_t; color = :black, linewidth = 2,
    label = "simulated"
)
axislegend(ax6; position = :lt)
fig5

md"""
The simulated path sits inside the 90% interval on most days, which is what a calibrated fit to one series looks like.
"""

covered = count(
    i -> infection_bands[1, i] <= sim_stochastic.I_t[i] <= infection_bands[3, i],
    eachindex(inf_days)
)
(; covered, days = length(inf_days))

md"""
## Summary

Each extension is one positional argument on [`Renewal`](@ref), and none of them changes the observation model, the latent process, or the fitting code.
Each draws whatever it does not know through the same seam every other component uses, so the importation rate can equally be a fixed constant, an unknown constant, or a time-varying process such as a [`RandomWalk`](@ref).
In every case the prior is on the log scale, like every other unknown positive quantity here.

Stochastic infections are the one extension with a second form.
As a modifier they are non-centred, and [`StochasticRenewal`](@ref) draws the same specification centred, which is the parameterisation to prefer when the data inform the infection path.
"""
