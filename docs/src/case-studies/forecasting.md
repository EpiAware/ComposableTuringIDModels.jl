# [Out-of-sample forecasting](@id case-study-forecasting)

Fitting a model to a series of length ``T`` estimates the past; a forecast asks
what comes next.
This case study fits a renewal model to the first part of an epidemic and then
predicts the reported cases over a future horizon ``t = T+1, \dots, T+h``,
propagating both parameter and latent uncertainty into the forecast.

Forecasting is a single call, [`forecast`](@ref), because the package's latent
processes are non-centred.
A [`RandomWalk`](@ref) or [`AR`](@ref) process accumulates an i.i.d. sequence of
standard innovations, so extending the reproduction-number path past the data is
just a matter of drawing more innovations from the prior.
[`forecast`](@ref) carries each posterior draw forward — holding the fitted
parameters and the in-sample latent path fixed, and continuing the latent
process over the horizon with fresh prior innovations — then draws the future
observations from the observation model.

## The model

The model is the renewal process with a second-order autoregressive ``\log R_t``
and negative-binomial reporting from the [renewal case study](@ref
case-study-renewal), so here it is assembled without further comment.

```@example forecasting
using ComposableTuringIDModels, Distributions, Random, Turing, Mooncake
using ADTypes: AutoMooncake
Random.seed!(1234)

latent = AR(
    damp = [truncated(Normal(0.8, 0.05), 0, 1),
        truncated(Normal(0.1, 0.05), 0, 1)],
    init = [Normal(0.0, 0.2), Normal(0.0, 0.2)],
    ϵ_t = HierarchicalNormal(std = HalfNormal(0.1)))
renewal = Renewal(gen_distribution = Gamma(6.5, 0.62); rt = latent,
    initialisation = Normal(log(1.0), 0.1))
model = IDModel(renewal, NegativeBinomialError(cluster_factor = HalfNormal(0.1)))
nothing # hide
```

## Split the series into a training window and a held-out horizon

We use the same South Korean COVID-19 series as the renewal case study, but hold
out the last `h` days so the forecast can be checked against observations the
model never saw.

```@example forecasting
using CSV, DataFrames
datapath = joinpath(pkgdir(ComposableTuringIDModels),
    "docs", "src", "case-studies", "data", "south_korea_data.csv")
south_korea = CSV.read(datapath, DataFrame)

tspan = (45, 85)
cases = south_korea.cases_new[first(tspan):last(tspan)]
h = 10
T = length(cases) - h
y_train = cases[1:T]
y_holdout = cases[(T + 1):end]
(T = T, h = h, train_total = sum(y_train))
```

## Fit to the training window

Conditioning on the training cases and sampling with NUTS recovers the posterior
over the model parameters and the in-sample ``\log R_t`` path.
We differentiate with [Mooncake](https://chalk-lab.github.io/Mooncake.jl/), the
recommended backend (see [Automatic differentiation backend](@ref ad-backend)).

```@example forecasting
posterior = as_turing_model(model, y_train, T)
chain = sample(
    posterior, NUTS(0.9; adtype = AutoMooncake(; config = nothing)),
    MCMCThreads(), 150, 2; progress = false)
nothing # hide
```

## Forecast the horizon

[`forecast`](@ref) takes the fitted model, the training series, the posterior
`chain`, and the horizon `h`.
It returns a chain of the same shape carrying the predicted observations
`y_t[T+1] … y_t[T+h]`; the in-sample points stay conditioned on the data.

```@example forecasting
fc = forecast(model, y_train, chain, h)
size(fc)
```

Because the latent process was continued rather than re-drawn, re-running the
horizon-length model over the forecast chain with `returned` recovers a
reproduction-number path that runs unbroken from the fitted window into the
forecast.

```@example forecasting
fc_model = as_turing_model(model, vcat(y_train, fill(missing, h)), T + h)
gens = vec(returned(fc_model, fc))
nothing # hide
```

## Bands and plot

A couple of small helpers reduce the per-draw trajectories to credible bands and
draw a median line with 50% and 95% ribbons.

```@setup forecasting
using Statistics

const CI_QS = [0.025, 0.25, 0.5, 0.75, 0.975]

# time × 5 credible bands from a time × draws matrix
function credible_bands(mat; qs = CI_QS)
    reduce(hcat, (map(eachrow(mat)) do row
        vals = collect(skipmissing(row))
        isempty(vals) ? missing : quantile(vals, q)
    end for q in qs))
end

# median line with 50% and 95% ribbons
function ci_ribbon!(ax, ts, bands; color, label)
    keep = findall(!ismissing, view(bands, :, 3))
    x, b = ts[keep], Float64.(bands[keep, :])
    band!(ax, x, b[:, 1], b[:, 5]; color = (color, 0.15))
    band!(ax, x, b[:, 2], b[:, 4]; color = (color, 0.3))
    lines!(ax, x, b[:, 3]; color = color, linewidth = 2, label = label)
end

# posterior-predictive y_t bands over the forecast horizon T+1 … T+h
function forecast_bands(fc, T, h)
    rows = map((T + 1):(T + h)) do i
        permutedims(vec(fc[@varname(y_t[i])]))
    end
    credible_bands(reduce(vcat, rows))
end
```

The reproduction number ``R_t = \exp(Z_t)`` runs across the whole span; the
forecast cases occupy the held-out window and are compared against the
observations the model never saw.

```@example forecasting
using CairoMakie

Rt = credible_bands(reduce(hcat, (exp.(g.Z_t) for g in gens)))
yt = forecast_bands(fc, T, h)

fig = Figure(; size = (760, 620))
ax1 = Axis(fig[1, 1]; ylabel = "Reproduction number Rₜ")
ci_ribbon!(ax1, 1:size(Rt, 1), Rt; color = :purple, label = "posterior median")
vlines!(ax1, [T + 0.5]; color = :grey, linestyle = :dash)
hlines!(ax1, [1.0]; color = :grey, linestyle = :dot)
axislegend(ax1; position = :rt)

ax2 = Axis(fig[2, 1]; xlabel = "Day", ylabel = "Reported cases")
ci_ribbon!(ax2, (T + 1):(T + h), yt; color = :teal, label = "forecast")
scatter!(ax2, 1:T, y_train; color = :black, markersize = 7, label = "training")
scatter!(ax2, (T + 1):(T + h), y_holdout; color = :red, markersize = 9,
    marker = :diamond, label = "held out")
vlines!(ax2, [T + 0.5]; color = :grey, linestyle = :dash)
axislegend(ax2; position = :lt)
fig
```

The dashed line marks the last day of training data.
Left of it the ``R_t`` path is informed by the cases; right of it the forecast
fans out as the autoregressive process reverts towards its mean and the credible
interval widens with the horizon.
The forecast band for reported cases brackets the held-out observations, so the
model extrapolates the tail of the wave without having seen it.

## Forecasting from an `IDProblem`

The same call works on an [`IDProblem`](@ref), taking the infection and
observation models from the problem and extending its `tspan`:

```@example forecasting
problem = IDProblem(infection = renewal,
    observation_model = NegativeBinomialError(cluster_factor = HalfNormal(0.1)),
    tspan = (1, T))
prob_chain = sample(as_turing_model(problem, (; y_t = y_train)),
    NUTS(0.9; adtype = AutoMooncake(; config = nothing)), 150; progress = false)
size(forecast(problem, y_train, prob_chain, h))
```
