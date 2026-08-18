# Right-truncated, double-interval-censored discrete PMF of a continuous
# distribution, built with CensoredDistributions.jl.
#
# `double_interval_censored(dist; upper = D, interval = Δd)` applies primary
# (uniform-window) censoring, right-truncation at `D`, then secondary interval
# censoring of width `Δd`. Evaluating its `pdf` on the bin left-edges
# `0, Δd, …, D-Δd` and normalising gives the discrete PMF the models consume.
# When `D` is `nothing`, a distribution with finite support (e.g. one the
# caller has already `truncated`) uses that bound directly as the horizon —
# it already says where the mass ends, so re-deriving a horizon from a
# quantile would silently produce a shorter PMF than the caller's own
# truncation implies. That bound is rounded *up* to the nearest multiple of
# `Δd`: the downstream bins `0, Δd, …, D-Δd` must land exactly on `D`, and
# every consumer (`LatentDelay`, `RightTruncate`, `ReportTriangle`,
# `Renewal`) treats the PMF as covering `[0, D)` in whole `Δd`-wide steps.
# Rounding down would silently discard the slice of the caller's own support
# between the rounded point and the true bound; rounding up only widens the
# final bin into a region where the distribution's density is already zero
# (nothing past `maximum(dist)` has mass), so no mass is fabricated or lost.
# An unbounded distribution falls back to the `upper`th quantile rounded to a
# multiple of `Δd`, unchanged from before.
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
