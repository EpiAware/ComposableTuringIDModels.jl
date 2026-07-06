# [Overview](@id overview)

`EpiAwarePrototype` builds an epidemiological model by assembling small parts
rather than writing one bespoke model.
Each part plays one of three roles — a **latent** process, an **infection**
process, or an **observation** process — and each part becomes a
[Turing](https://turinglang.org) / `DynamicPPL` model through a single generic
constructor, [`as_turing_model`](@ref).
Because every part speaks that one interface, parts nest inside one another and
a whole model is *composed* from the pieces.

<figure style="margin:1.5rem 0">
<svg viewBox="0 0 820 445" role="img" aria-labelledby="ovw-t ovw-d" style="width:100%;height:auto;max-width:820px;font-family:system-ui,Segoe UI,Helvetica,Arial,sans-serif">
<title id="ovw-t">Composable design of EpiAwarePrototype</title>
<desc id="ovw-d">A latent process nested inside an infection model feeds an observation model through the single as_turing_model interface; any of the three roles can be swapped.</desc>
<rect x="2" y="2" width="816" height="441" rx="16" fill="#fbfafc" stroke="#e6e3ec"/>
<defs>
<marker id="ovw-arrow" viewBox="0 0 10 10" refX="8" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse"><path d="M0 0 L10 5 L0 10 z" fill="#6b6b72"/></marker>
</defs>
<rect x="28" y="52" width="312" height="150" rx="14" fill="#e7eef8" stroke="#2f6fb0" stroke-width="2"/>
<text x="44" y="78" font-size="15" font-weight="700" fill="#204a75">Infection model</text>
<text x="44" y="96" font-size="11" fill="#3f6fa0">owns and drives its own latent process</text>
<rect x="44" y="110" width="284" height="82" rx="10" fill="#efe7f6" stroke="#8a4faf" stroke-width="1.6"/>
<text x="58" y="134" font-size="13" font-weight="700" fill="#6a3d8f">Latent process  Z<tspan baseline-shift="sub" font-size="9">t</tspan></text>
<text x="58" y="154" font-size="10.5" fill="#86659e">random walk · AR · differenced AR …</text>
<text x="58" y="174" font-size="10.5" font-style="italic" fill="#86659e">the infection model draws it, then maps it to I<tspan baseline-shift="sub" font-size="8">t</tspan></text>
<line x1="340" y1="127" x2="396" y2="127" stroke="#6b6b72" stroke-width="1.6" marker-end="url(#ovw-arrow)"/>
<text x="356" y="118" font-size="12" font-style="italic" fill="#4a4a52">I<tspan baseline-shift="sub" font-size="8">t</tspan></text>
<rect x="400" y="52" width="300" height="150" rx="14" fill="#e7f3df" stroke="#3f9a2c" stroke-width="2"/>
<text x="416" y="78" font-size="15" font-weight="700" fill="#2c6a1e">Observation model</text>
<text x="416" y="96" font-size="11" fill="#47823a">maps infections to the observed data</text>
<rect x="416" y="110" width="268" height="82" rx="10" fill="#f2f9ec" stroke="#79bd5d" stroke-width="1.4"/>
<text x="430" y="134" font-size="10.5" fill="#3d7a2b">reporting delay · ascertainment</text>
<text x="430" y="154" font-size="10.5" fill="#3d7a2b">day-of-week · right-truncation</text>
<text x="430" y="174" font-size="10.5" font-style="italic" fill="#3d7a2b">counts: Poisson / negative-binomial</text>
<line x1="700" y1="127" x2="748" y2="127" stroke="#6b6b72" stroke-width="1.6" marker-end="url(#ovw-arrow)"/>
<text x="710" y="118" font-size="12" font-style="italic" fill="#4a4a52">y<tspan baseline-shift="sub" font-size="8">t</tspan></text>
<rect x="750" y="104" width="60" height="46" rx="8" fill="#edeaf1" stroke="#9a93a8"/>
<text x="780" y="124" font-size="10.5" text-anchor="middle" fill="#4a4553">data</text>
<text x="780" y="140" font-size="12" font-weight="700" text-anchor="middle" fill="#4a4553">y<tspan baseline-shift="sub" font-size="8">t</tspan></text>
<line x1="184" y1="202" x2="184" y2="246" stroke="#6d5b8a" stroke-width="1.4" stroke-dasharray="3 4"/>
<line x1="550" y1="202" x2="550" y2="246" stroke="#6d5b8a" stroke-width="1.4" stroke-dasharray="3 4"/>
<rect x="28" y="248" width="672" height="52" rx="12" fill="#ece7f3" stroke="#6d5b8a" stroke-width="1.8" stroke-dasharray="6 4"/>
<text x="364" y="272" font-size="14" font-weight="700" text-anchor="middle" fill="#4c3d6b" font-family="ui-monospace,SFMono-Regular,Menlo,monospace">as_turing_model</text>
<text x="364" y="290" font-size="11" text-anchor="middle" fill="#6b5c86">the one interface every part implements — parts compose as submodels</text>
<text x="419" y="336" font-size="12.5" font-weight="700" text-anchor="middle" fill="#4a4553">Swap any part — change one assumption without touching the rest</text>
<text x="42" y="356" font-size="12" font-weight="700" fill="#6a3d8f">Latent</text>
<rect x="36" y="362" width="226" height="20" rx="6" fill="#f4eef9" stroke="#8a4faf"/><text x="46" y="376" font-size="11" fill="#6a3d8f">RandomWalk</text>
<rect x="36" y="386" width="226" height="20" rx="6" fill="#f4eef9" stroke="#8a4faf"/><text x="46" y="400" font-size="11" fill="#6a3d8f">AR · MA · ARIMA</text>
<rect x="36" y="410" width="226" height="20" rx="6" fill="#f4eef9" stroke="#8a4faf"/><text x="46" y="424" font-size="11" fill="#6a3d8f">DiffLatentModel</text>
<text x="304" y="356" font-size="12" font-weight="700" fill="#204a75">Infection</text>
<rect x="298" y="362" width="244" height="20" rx="6" fill="#eef3fa" stroke="#2f6fb0"/><text x="308" y="376" font-size="11" fill="#204a75">DirectInfections</text>
<rect x="298" y="386" width="244" height="20" rx="6" fill="#eef3fa" stroke="#2f6fb0"/><text x="308" y="400" font-size="11" fill="#204a75">Renewal  (drives Rt)</text>
<rect x="298" y="410" width="244" height="20" rx="6" fill="#eef3fa" stroke="#2f6fb0"/><text x="308" y="424" font-size="11" fill="#204a75">ExpGrowthRate</text>
<text x="584" y="356" font-size="12" font-weight="700" fill="#2c6a1e">Observation</text>
<rect x="578" y="362" width="232" height="20" rx="6" fill="#eef6e8" stroke="#3f9a2c"/><text x="588" y="376" font-size="11" fill="#2c6a1e">PoissonError</text>
<rect x="578" y="386" width="232" height="20" rx="6" fill="#eef6e8" stroke="#3f9a2c"/><text x="588" y="400" font-size="11" fill="#2c6a1e">NegativeBinomialError</text>
<rect x="578" y="410" width="232" height="20" rx="6" fill="#eef6e8" stroke="#3f9a2c"/><text x="588" y="424" font-size="11" fill="#2c6a1e">LatentDelay wrapper</text>
</svg>
<figcaption style="font-size:0.85rem;color:#6b6b72;text-align:center;margin-top:0.4rem">The three roles plug into one <code>as_turing_model</code> interface. The infection model owns its latent process, and any part can be swapped for a compatible one.</figcaption>
</figure>

## Three roles, one interface

A model is put together from parts filling three roles.

  - A **latent** process describes an unobserved series ``Z_t`` over time, e.g.
    a log reproduction number or a growth rate.
  - An **infection** process turns that latent series into unobserved infections
    ``I_t`` (directly, through exponential growth, or through the renewal
    equation).
  - An **observation** process turns infections into the observed data ``y_t``,
    adding reporting delays, ascertainment, day-of-week effects,
    right-truncation, and count noise.

Each part is a plain struct with a single method of
[`as_turing_model`](@ref), which returns a `DynamicPPL.Model`.
There is no deep type hierarchy: a part is identified by the method it
implements, not by its place in a tree.
A part that contains another part builds the inner model and samples it as a
submodel, so the parts nest through the same interface they expose.

## The infection model owns its latent

The latent process is not a separate top-level component.
The reason is that the latent (e.g. ``\log R_t``) is not always the quantity you
care about, so it is handed to the infection model that consumes it rather than
threaded through the composer.
An infection model takes a latent slot — `Z` for [`DirectInfections`](@ref),
`rt` for [`ExpGrowthRate`](@ref) and [`Renewal`](@ref) — draws that process
internally, and maps it to infections.
Only [`Renewal`](@ref) needs a generation interval, so it alone carries an
[`EpiData`](@ref); the others take a transformation directly.

## Swap a part to change an assumption

Because the parts share one interface, you compare modelling assumptions by
swapping one struct for another and leaving the rest untouched.

```@example overview
using EpiAwarePrototype, Distributions

# One latent process: an ARIMA-style differenced AR.
latent = DiffLatentModel(; model = AR(), init_priors = [Normal(), Normal()])

# Fold it into a direct-infections process, then swap only the observation
# model. Everything else stays the same.
poisson_model = EpiAwareModel(
    DirectInfections(; Z = latent, initialisation_prior = Normal()),
    PoissonError())

negbin_model = EpiAwareModel(
    DirectInfections(; Z = latent, initialisation_prior = Normal()),
    NegativeBinomialError())

# Each assembly is turned into one Turing model. `missing` data simulates from
# the prior; the composed model exposes its generated quantities.
turing_model = as_turing_model(poisson_model, missing, 20)
(; generated_y_t, I_t, Z_t) = turing_model()
length(generated_y_t), length(I_t), length(Z_t)
```

## Where to go next

  - [Composable design](@ref) explains the `as_turing_model` protocol and how
    parts nest as submodels in more detail.
  - The [case studies](@ref case-studies-overview) build complete models and fit
    them to real surveillance data, from a renewal model to a compartmental SIR.
  - The [Public API](@ref public-api) lists every component you can compose.
