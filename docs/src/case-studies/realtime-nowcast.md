# [Real-time nowcasting: correcting right-truncation](@id case-study-nowcast)

In real time, the most recent days of a surveillance series are **incomplete**:
a case with reference day ``t`` is only reported after a reporting delay, so on
day ``\text{now}`` we have seen just a fraction of the cases that will eventually
be attributed to recent days. Fitting a renewal model to this *right-truncated*
tail without correction produces the familiar artefact of real-time estimation —
an apparent late down-turn in ``R_t`` that is really an artefact of not-yet-reported
counts [abbott2020estimating](@citep).

This case study simulates an outbreak, truncates it at ``\text{now}`` to mimic a
real-time snapshot, and contrasts two fits: a **naive** one that treats the
truncated counts as complete, and one that wraps the observation model in
[`RightTruncate`](@ref) to scale each reference day's expected count by the
fraction of its eventual total reported so far. The correction is the
EpiNow2-style CDF-scaling nowcast, expressed here as a one-line observation
modifier.

!!! note "Prototype"
    Short sampler runs and a single simulated path keep this page fast to build.
    For real analyses use more iterations, check convergence, and supply observed
    data.

## The idea

The infection pipeline produces ``Y_t = \mu_t``, the expected *eventual* total
for reference day ``t``. At time ``\text{now}`` a reference day of age
``a = \text{now} - t`` has only had a fraction ``F[a+1]`` of its eventual total
reported, where ``F`` is the reporting-delay CDF. The expected observed-so-far is
therefore ``\mu_t \cdot F[a+1]``. [`RightTruncate`](@ref) conditions the
observation error on that down-weighted mean, so the model's ``Y_t`` stays the
eventual total and the nowcast is just ``Y_t`` read back out.

```math
\begin{aligned}
R_t &= \exp(Z_t), & Z_t &= \rho\, Z_{t-1} + \epsilon_t, \\
I_t &= R_t \sum_{s \ge 1} g_s\, I_{t-s}, \\
y_t &\sim \mathrm{NegBinomial}\!\big(I_t \cdot F[(\text{now}-t)+1],\ \phi\big).
\end{aligned}
```

A naive fit drops the ``F[\cdot]`` factor (equivalently assumes ``F \equiv 1``).

## The full-data model

We build the same composed renewal model as the
[renewal case study](@ref case-study-renewal): an autoregressive ``\log R_t``
process folded into a [`Renewal`](@ref) infection process, observed with a
[`NegativeBinomialError`](@ref).

```@example nowcast
using EpiAwarePrototype, Distributions, Random, Turing
Random.seed!(20240625)

latent = AR(
    damp_priors = [truncated(Normal(0.8, 0.05), 0, 1)],
    init_priors = [Normal(0.0, 0.25)],
    ϵ_t = HierarchicalNormal(std_prior = HalfNormal(0.1)))
data = EpiData(gen_distribution = Gamma(6.5, 0.62))
renewal = Renewal(data; rt = latent, initialisation_prior = Normal(log(1.4), 0.1))
error = NegativeBinomialError(cluster_factor_prior = HalfNormal(0.1))
nothing # hide
```

## Simulate an outbreak and truncate it

We simulate a full series of eventual totals, then impose a reporting delay and
truncate at ``\text{now}`` to produce the partially-reported tail a real-time
analyst actually sees. The reporting delay is a continuous distribution
discretised to a CDF with the same released
[CensoredDistributions.jl](https://github.com/EpiAware/CensoredDistributions.jl)
path the rest of the package uses; [`ReportingCDF`](@ref) builds the completeness
curve ``F``.

```@example nowcast
n = 50
reporting_delay = LogNormal(1.6, 0.5)        # mean ≈ 5.6 days
cdf_curve = ReportingCDF(reporting_delay)
nothing # hide
```

Simulate the eventual totals from the model (no truncation), and treat one path
as the ground truth.

```@example nowcast
truth_model = EpiAwareModel(renewal, error)
truth = as_turing_model(truth_model, fill(missing, n), n)()
eventual = truth.generated_y_t              # the eventual totals
R_true = exp.(truth.Z_t)
extrema(R_true)
```

Now form the right-truncated snapshot: thin each reference day's eventual total
to the fraction reported by ``\text{now} = n``. The completeness by reference day
is ``F`` reversed onto the time axis (the most recent day is least complete).

```@example nowcast
completeness = as_turing_model(cdf_curve, n)()    # F by age, F[1] = age 0
scale = reverse(completeness)                     # by reference day t = 1..n
observed_so_far = @. rand(Binomial(eventual, scale))
(complete_tail = eventual[(end - 4):end], truncated_tail = observed_so_far[(end - 4):end])
```

The last few days are visibly thinned: the most recent day shows only a fraction
of its eventual count.

## Fit 1 — naive (no truncation correction)

Condition the plain renewal model on the truncated counts as though they were
complete.

```@example nowcast
naive_model = EpiAwareModel(renewal, error)
naive_post = as_turing_model(naive_model, observed_so_far, n)
naive_chain = sample(naive_post, NUTS(0.9), 100; progress = false)
nothing # hide
```

## Fit 2 — corrected with [`RightTruncate`](@ref)

Wrap the *same* error model in [`RightTruncate`](@ref) with the *same* reporting
delay. Nothing else changes — the infection process and its latent ``R_t`` are
untouched; only how the expected counts are compared to the truncated data.

```@example nowcast
corrected_obs = RightTruncate(error, reporting_delay)
corrected_model = EpiAwareModel(renewal, corrected_obs)
corrected_post = as_turing_model(corrected_model, observed_so_far, n)
corrected_chain = sample(corrected_post, NUTS(0.9), 100; progress = false)
nothing # hide
```

## The right-truncation bias and its correction

The reproduction number ``R_t = \exp(Z_t)`` is a generated quantity of the model,
not a sampled parameter. Re-running the fitted model over the chain with
[`returned`](https://turinglang.org/Turing.jl/) (re-exported by `Turing`) recovers
the latent ``Z_t`` trajectory per draw, from which we take the posterior-mean
``R_t`` and average it over the most recent window. Right-truncation biases the
naive fit's recent ``R_t`` *downward* (the not-yet-reported counts look like a
decline); the [`RightTruncate`](@ref) fit removes that bias.

```@example nowcast
using Statistics

function recent_Rt(posterior, chain; window = 7)
    # Per-draw R_t = exp(Z_t); average over draws, then over the recent window.
    gens = vec(returned(posterior, chain))
    Rt_mean = mean(exp.(g.Z_t) for g in gens)
    mean(Rt_mean[(end - window + 1):end])
end

R_true_recent = mean(R_true[(end - 6):end])
R_naive_recent = recent_Rt(naive_post, naive_chain)
R_corrected_recent = recent_Rt(corrected_post, corrected_chain)

(truth = round(R_true_recent, digits = 2),
    naive = round(R_naive_recent, digits = 2),
    corrected = round(R_corrected_recent, digits = 2))
```

The naive recent-``R_t`` sits **below** the truth — the artefactual late
down-turn produced by treating the not-yet-reported tail as complete — whereas the
[`RightTruncate`](@ref) fit, which knows the recent days are incomplete, pulls the
recent ``R_t`` back **up** off that spurious decline. (These are short,
illustrative runs on a single simulated path, so the exact corrected value is
noisy; the robust, repeatable signal is the *direction* — the naive fit
under-estimates recent ``R_t``, and the correction removes that downward pull.)
The nowcast of the eventual totals is the corrected model's ``Y_t``, recovered the
same way from the per-draw generated quantities; this page stops at the ``R_t``
summary, consistent with the other case studies (the docs build runs no plotting
stack).

## Reading the shared parameters

`sample` returns a [FlexiChains](https://github.com/penelopeysm/FlexiChains.jl)
chain, which we index by variable name directly — no conversion step. Wrapping the
error model in [`RightTruncate`](@ref) does not touch the renewal process, so the
corrected fit recovers the *same* shared parameters: the autoregressive damping
``\rho`` (`damp_AR[1]`), the innovation scale ``\sigma`` (`std`), the observation
overdispersion (`cluster_factor`), and the initial infections (`init_incidence`).
They keep their flat names (prefixing is disabled throughout the package).

```@example nowcast
using Turing: @varname

posterior_draws(chain, vn) = vec(chain[vn])
posterior_summary(chain, vn) = (mean = mean(posterior_draws(chain, vn)),
    std = std(posterior_draws(chain, vn)))

(damp = posterior_summary(corrected_chain, @varname(damp_AR[1])),
    sigma = posterior_summary(corrected_chain, @varname(std)),
    cluster_factor = posterior_summary(corrected_chain, @varname(cluster_factor)),
    init_incidence = posterior_summary(corrected_chain, @varname(init_incidence)))
```

## From the marginal to the joint

[`RightTruncate`](@ref) corrects right-truncation by conditioning on each
reference day's observed-so-far *total* — the **marginal** of the full
reference-day × reporting-delay structure. When the delay structure itself is of
interest (e.g. reporting that drifts over the outbreak), the
[`ReportTriangle`](@ref) observation model keeps the full 2D reporting triangle
and scores it cell by cell; its observed row-sums reconcile with this CDF-scaling
to machine precision. The marginal correction here is the cheaper, released-code
first step.

## References

```@bibliography
Pages = ["realtime-nowcast.md"]
Canonical = false
```
