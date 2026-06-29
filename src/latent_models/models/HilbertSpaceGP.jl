# Hilbert-space approximate Gaussian-process latent model.

@doc raw"
A **Hilbert-space approximate Gaussian process** (HSGP) latent process.

A Gaussian process places a prior over functions and is a natural latent process
for a smoothly varying quantity such as ``\log R_t``. An exact GP is impractical
inside a sampler: it needs an ``n \times n`` covariance factorisation that costs
``O(n^3)`` per leapfrog step. This model uses the **Hilbert-space basis-function
approximation** of [riutortmayol2023practical](@citep), which writes the GP as a
short weighted sum of fixed basis functions

```math
f(x) \approx \sum_{j=1}^{m} \phi_j(x)\, \sqrt{S(\sqrt{\lambda_j})}\; \beta_j,
\qquad \beta_j \sim \mathrm{Normal}(0, 1),
```

where the eigenfunctions ``\phi_j`` and eigenvalues ``\lambda_j`` of the Laplacian
on the interval ``[-L, L]`` are

```math
\phi_j(x) = \sqrt{\tfrac{1}{L}}\, \sin\!\Big(\sqrt{\lambda_j}\,(x + L)\Big),
\qquad \sqrt{\lambda_j} = \frac{\pi j}{2 L},
```

and ``S`` is the spectral density of the chosen covariance kernel. For the
squared-exponential kernel with marginal standard deviation ``\sigma`` and length
scale ``\ell``,

```math
S(\omega) = \sigma^2 \sqrt{2\pi}\, \ell \, \exp\!\Big(-\tfrac{1}{2}\ell^2\omega^2\Big).
```

Only ``\ell``, ``\sigma`` and the ``m`` weights ``\beta`` are sampled; the basis
``\phi_j`` and eigenvalues ``\lambda_j`` depend only on `n`, `m` and the boundary
factor `c`, **not** on any sampled parameter, so they fall entirely outside the
differentiated path. The latent path is therefore a cheap matrix–vector product
of a fixed basis against a small set of standard-normal weights — a non-centred
parameterisation that is fast and samples well under NUTS, including with
[Mooncake](https://chalk-lab.github.io/Mooncake.jl/) reverse-mode AD.

The accuracy/speed trade-off is controlled by two numbers
[riutortmayol2023practical](@citep): the number of basis functions `m` (more
basis functions resolve shorter length scales, at linear cost) and the boundary
factor `c` (the domain is extended to ``L = c\,S`` beyond the half-range ``S`` of
the rescaled inputs, so that boundary effects do not distort the fit). The
defaults (`m = 20`, `c = 1.5`) are a reasonable starting point for a smooth
latent process of moderate length; short length scales relative to the series may
need a larger `m`.

## Fields

  - `length_scale_prior`: prior for the length scale ``\ell``.
  - `marginal_std_prior`: prior for the marginal standard deviation ``\sigma``.
  - `m`: number of basis functions.
  - `c`: boundary factor; the GP is approximated on ``[-L, L]`` with ``L = c S``.

# Examples
```@example HilbertSpaceGP
using EpiAwarePrototype, Distributions
gp = HilbertSpaceGP()
mdl = as_turing_model(gp, 30)
rand(mdl)
```
"
@kwdef struct HilbertSpaceGP{L <: Sampleable, S <: Sampleable, M <: Int, C <: Real} <:
              AbstractLatentModel
    "Prior distribution for the length scale ``\\ell``."
    length_scale_prior::L = truncated(Normal(0.0, 0.4), 0, Inf)
    "Prior distribution for the marginal standard deviation ``\\sigma``."
    marginal_std_prior::S = truncated(Normal(0.0, 1.0), 0, Inf)
    "Number of basis functions."
    m::M = 20
    "Boundary factor: the GP is approximated on ``[-L, L]`` with ``L = c S``."
    c::C = 1.5

    function HilbertSpaceGP(length_scale_prior::Sampleable,
            marginal_std_prior::Sampleable, m::Int, c::Real)
        @assert m>0 "m (the number of basis functions) must be greater than 0"
        @assert c>1 "c (the boundary factor) must be greater than 1"
        new{typeof(length_scale_prior), typeof(marginal_std_prior), typeof(m),
            typeof(c)}(length_scale_prior, marginal_std_prior, m, c)
    end
end

@doc raw"
Build the Hilbert-space GP basis for `n` evenly spaced inputs.

Returns `(Φ, sqrt_λ)` where `Φ` is the ``n \times m`` matrix of eigenfunctions
``\phi_j`` evaluated at the rescaled inputs and `sqrt_λ` is the length-`m` vector
of ``\sqrt{\lambda_j}``. The inputs ``t = 1, \ldots, n`` are centred and rescaled
to unit spacing about their midpoint, so the half-range is ``S = (n-1)/2`` and the
GP is approximated on ``[-L, L]`` with ``L = c S``. Both outputs depend only on
`n`, `m` and `c` — none of the sampled parameters — so this computation lies
entirely outside the differentiated path of the [`HilbertSpaceGP`](@ref) model.
Requires `n > 1` so the half-range ``S`` is positive.
"
function hsgp_basis(n::Int, m::Int, c::Real)
    @assert n>1 "n must be greater than 1 for a well-defined basis (S = (n-1)/2 > 0)"
    x = collect(1:n) .- (n + 1) / 2          # centre the integer index about 0
    S = (n - 1) / 2                           # half-range of the rescaled inputs
    L = c * S
    j = collect(1:m)'
    sqrt_λ = (π .* j) ./ (2L)                  # √eigenvalues, 1×m
    Φ = sqrt(1 / L) .* sin.(sqrt_λ .* (x .+ L))  # n×m eigenfunctions
    return Φ, vec(sqrt_λ)
end

@doc raw"
Spectral density of the squared-exponential kernel at frequency ``\omega``.

```math
S(\omega) = \sigma^2 \sqrt{2\pi}\, \ell \, \exp\!\Big(-\tfrac{1}{2}\ell^2\omega^2\Big)
```

with marginal standard deviation ``\sigma`` and length scale ``\ell``. Used by
[`HilbertSpaceGP`](@ref) to weight each basis function.
"
function se_spectral_density(ω, σ, ℓ)
    return σ^2 * sqrt(2π) * ℓ .* exp.(-(ℓ^2 / 2) .* ω .^ 2)
end

@model function as_turing_model(model::HilbertSpaceGP, n)
    @assert n>1 "n must be greater than 1"
    Φ, sqrt_λ = hsgp_basis(n, model.m, model.c)
    ℓ ~ model.length_scale_prior
    σ ~ model.marginal_std_prior
    β ~ filldist(Normal(), model.m)
    spectral_weights = sqrt.(se_spectral_density(sqrt_λ, σ, ℓ))
    gp = Φ * (spectral_weights .* β)
    return gp
end
