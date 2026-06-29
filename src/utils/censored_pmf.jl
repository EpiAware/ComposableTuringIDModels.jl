# Right-truncated, double-interval-censored discrete PMF built with
# CensoredDistributions.jl (replaces upstream's bespoke `censored_pmf`).

# Right-truncated, double-interval-censored discrete PMF of a continuous
# distribution, built with CensoredDistributions.jl rather than a bespoke
# quadrature (the original EpiAware shipped its own `censored_pmf`/`censored_cdf`;
# here that is delegated to the org's CensoredDistributions package).
#
# `double_interval_censored(dist; upper = D, interval = Δd)` applies primary
# (uniform-window) censoring, right-truncation at `D`, then secondary interval
# censoring of width `Δd`. Evaluating its `pdf` on the bin left-edges
# `0, Δd, …, D-Δd` and normalising gives the discrete PMF the models consume.
# When `D` is `nothing` it defaults to the `upper`th quantile rounded to a
# multiple of `Δd`, matching the original behaviour.
function _discretised_pmf(dist::Distribution; Δd = 1.0, D = nothing, upper = 0.99)
    @assert minimum(dist)>=0.0 "Distribution must be non-negative."
    @assert Δd>0.0 "Δd must be positive."
    if isnothing(D)
        D = round(Int64, invlogcdf(dist, log(upper)) / Δd) * Δd
    end
    @assert D>=Δd "D can't be shorter than Δd."
    censored = double_interval_censored(dist; upper = D, interval = Δd)
    ts = 0.0:Δd:(D - Δd)
    probs = [pdf(censored, t) for t in ts]
    return probs ./ sum(probs)
end
