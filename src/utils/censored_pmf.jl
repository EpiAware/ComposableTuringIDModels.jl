# Right-truncated, double-interval-censored discrete PMF of a continuous
# distribution, built with CensoredDistributions.jl.
#
# `double_interval_censored(dist; upper = D, interval = Δd)` applies primary
# (uniform-window) censoring, right-truncation at `D`, then secondary interval
# censoring of width `Δd`. Evaluating its `pdf` on the bin left-edges
# `0, Δd, …, D-Δd` and normalising gives the discrete PMF the models consume.
# When `D` is `nothing`, a distribution with finite support (e.g. one the
# caller has already `truncated`) sets the horizon from that bound, which
# already says where the mass ends; a quantile would give a shorter PMF than
# the caller's own truncation implies. The bound rounds *up* to a multiple of
# `Δd` because every consumer treats the PMF as covering `[0, D)` in whole
# `Δd`-wide bins: rounding down would drop the support between the rounded
# point and the true bound, while rounding up only widens the last bin into a
# region with no mass. An unbounded distribution uses the `upper`th quantile.
function _discretised_pmf(dist::Distribution; Δd = 1.0, D = nothing, upper = 0.99)
    @assert minimum(dist) >= 0.0 "Distribution must be non-negative."
    @assert Δd > 0.0 "Δd must be positive."
    if isnothing(D)
        D = isfinite(maximum(dist)) ? ceil(maximum(dist) / Δd) * Δd :
            round(Int64, invlogcdf(dist, log(upper)) / Δd) * Δd
    end
    @assert D >= Δd "D can't be shorter than Δd."
    censored = double_interval_censored(dist; upper = D, interval = Δd)
    ts = 0.0:Δd:(D - Δd)
    probs = [pdf(censored, t) for t in ts]
    return probs ./ sum(probs)
end
