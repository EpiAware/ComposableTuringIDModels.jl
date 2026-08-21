Changes are documented in Github releases.

## Unreleased

  - A gap in an observation series is now marginalised out rather than imputed
    as a latent. A missing-at-random observation contributes no likelihood
    term, so there is nothing to infer at that point: it is skipped. Missing
    entries cost no parameters and no gradient, and a count stream with gaps
    can now be sampled by NUTS (previously it failed, because a sampled count
    has no bijector to link through). Predictive values at unobserved points
    are generated after fitting, by replaying the posterior through the model
    and drawing from the error model at the expected values that replay
    produces.

    `forecast` now builds its horizon that way rather than by reading sampled
    blanks out of the chain. It returns the same per-time-point predictions,
    and additionally carries the in-sample posterior predictive, since
    generating the horizon generates every other point at the same time.

    **Breaking**: a fitted model no longer carries `y_t[i]` parameters, so they
    are absent from the chain and from posterior summaries. Read predictions
    from `predict` on a fully missing series, or from `forecast`.
