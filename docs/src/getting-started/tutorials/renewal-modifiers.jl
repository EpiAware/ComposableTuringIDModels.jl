md"""
# [Renewal modifiers: susceptible depletion and importation](@id renewal-modifiers)

A [`Renewal`](@ref) process can be extended by *modifiers* that change how the
next incidence is formed.
This page takes one delayed renewal process and adds two of them —
[`SusceptibleDepletion`](@ref), which bounds the epidemic by a finite
population, and [`ImportedCases`](@ref), which seeds infections from outside it
— so the contribution of each is visible against the same baseline.

Everything else is held fixed: the same generation interval, the same
reproduction number, the same reporting delay, and the same initial incidence.
Only the modifier list changes.
"""

using ComposableTuringIDModels, Distributions, Random, Turing, CairoMakie
using DynamicPPL: fix
using Statistics: median, quantile

Random.seed!(20260727)

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
    Renewal(gen_int, modifiers...;
        rt = FixedIntercept(log(1.3)), initialisation = Normal())
end

md"""
Reported cases are the infections convolved with a reporting delay and observed
with negative-binomial noise — a [`LatentDelay`](@ref) wrapped around a
[`NegativeBinomialError`](@ref).
"""

delay = [0.1, 0.4, 0.3, 0.2]
obs = LatentDelay(NegativeBinomialError(cluster_factor = HalfNormal(0.1)), delay)

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
        renewal(SusceptibleDepletion(pop_size),
            ImportedCases(FixedIntercept(log(2.0)))),
        obs))

md"""
Simulating from each model with `missing` observations and the same fixed
initial incidence returns its infections and reported cases.
"""

simulate(model) = fix(as_turing_model(model, missing, n),
    (init_incidence = log(1.0),))()

sims = map(simulate, models)

md"""
## What each modifier does
"""

fig = Figure(; size = (760, 560))
ax1 = Axis(fig[1, 1]; ylabel = "Infections Iₜ", yscale = log10)
ax2 = Axis(fig[2, 1]; xlabel = "Day", ylabel = "Reported cases", yscale = log10)
colours = (plain = :grey30, depleting = :purple, seeded = :teal)
labels = (plain = "renewal", depleting = "+ susceptible depletion",
    seeded = "+ importation")
for k in keys(sims)
    lines!(ax1, 1:n, max.(sims[k].I_t, 1e-1); color = colours[k],
        linewidth = 2, label = labels[k])
    scatter!(ax2, 1:n, max.(sims[k].generated_y_t, 1e-1); color = colours[k],
        markersize = 6, label = labels[k])
end
axislegend(ax1; position = :lt)
fig

md"""
The plain renewal grows exponentially without bound.
Susceptible depletion turns that growth over: as susceptibles are used up the
effective reproduction number falls below one, so incidence peaks and declines.
Importation keeps seeding the epidemic instead of leaving it to grow from the
initial incidence alone, which brings the peak forward by more than two weeks
and roughly doubles it.
Once the susceptible pool is spent, incidence settles at the importation rate
rather than decaying away.

That last point is what importation is mainly for.
With ``R_t = 0.8`` and no susceptible depletion the plain process dies out,
while a seeded one levels off where importation and decay balance.
"""

function subcritical(modifiers...)
    fix(
        as_turing_model(
            IDModel(
                Renewal(gen_int, modifiers...; rt = FixedIntercept(log(0.8)),
                    initialisation = Normal()),
                obs),
            missing, n),
        (init_incidence = log(50.0),))()
end

sub_plain = subcritical()
sub_seeded = subcritical(ImportedCases(FixedIntercept(log(2.0))))

fig2 = Figure(; size = (760, 300))
ax3 = Axis(fig2[1, 1]; xlabel = "Day", ylabel = "Infections Iₜ",
    yscale = log10)
lines!(ax3, 1:n, max.(sub_plain.I_t, 1e-3); color = :grey30, linewidth = 2,
    label = "Rₜ = 0.8")
lines!(ax3, 1:n, max.(sub_seeded.I_t, 1e-3); color = :teal, linewidth = 2,
    label = "Rₜ = 0.8 + importation")
axislegend(ax3; position = :rt)
fig2

md"""
## Fitting a model with modifiers

Modifiers are part of the model, so a model carrying them is fitted exactly like
any other.
An [`ImportedCases`](@ref) rate is a prior slot, so it is *estimated* rather than
assumed: here it is a single unknown constant on the log scale, drawn once
before the renewal recursion runs.
We fit that model to the reported cases simulated above.
"""

model = IDModel(
    renewal(SusceptibleDepletion(pop_size),
        ImportedCases(Normal(0.0, 1.0))),
    obs)
y_obs = sims.seeded.generated_y_t
chain = sample(as_turing_model(model, y_obs, n), NUTS(), 250; progress = false)

md"""
The rate is namespaced by the modifier's position on the renewal step, so
several modifiers carrying priors can never collide.
Exponentiating the draws puts it back on the scale of infections per day, where
the true value was two.
"""

import_draws = exp.(vec(chain[@varname(modifier_2.import_rates)]))
(posterior = round.(quantile(import_draws, [0.05, 0.5, 0.95]), digits = 2),
    truth = 2.0)

md"""
The posterior predictive then tracks the simulated series.
"""

pred = predict(as_turing_model(model, fill(missing, n), n), chain)
## The reporting delay leaves the first days unscored, so band only the days
## the predictive chain carries.
y_draws(i) =
    try
        Float64.(vec(pred[@varname(y_t[i])]))
    catch
        Float64[]
    end
scored = filter(i -> !isempty(y_draws(i)), 1:n)
bands = reduce(hcat, map(i -> quantile(y_draws(i), [0.05, 0.5, 0.95]), scored))

fig3 = Figure(; size = (760, 300))
ax4 = Axis(fig3[1, 1]; xlabel = "Day", ylabel = "Reported cases",
    yscale = log10)
band!(ax4, scored, max.(bands[1, :], 1e-1), max.(bands[3, :], 1e-1);
    color = (:teal, 0.25))
lines!(ax4, scored, max.(bands[2, :], 1e-1); color = :teal, linewidth = 2,
    label = "posterior predictive")
scatter!(ax4, 1:n, max.(y_obs, 1e-1); color = :black, markersize = 6,
    label = "simulated")
axislegend(ax4; position = :lt)
fig3

md"""
## Summary

Both extensions are one positional argument on [`Renewal`](@ref), and neither
changes the observation model, the latent process, or the fitting code.
Susceptible depletion is a deterministic modifier; importation carries a prior
and draws it through the same seam every other component uses, so its rate can
equally be a fixed constant, an unknown constant, or a time-varying process such
as a [`RandomWalk`](@ref).
"""
