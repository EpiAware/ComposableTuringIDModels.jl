# Hilbert-space approximate Gaussian-process latent model.

@doc raw"
Supertype for the covariance kernels available to [`HilbertSpaceGP`](@ref).

A kernel enters the Hilbert-space approximation only through its **spectral
density** ``S(\omega)`` — the Fourier transform of the stationary covariance
function. The eigenfunction basis is shared across all kernels; switching kernel
just reweights the basis functions by ``\sqrt{S(\sqrt{\lambda_j})}``. Concrete
kernels implement

```julia
spectral_density(kernel, ω, σ, ℓ)  # ⇒ S(ω) for marginal sd σ, length scale ℓ
```

Members: [`SquaredExponentialKernel`](@ref), [`Matern32Kernel`](@ref),
[`Matern52Kernel`](@ref).
"
abstract type AbstractGPKernel end

@doc raw"
Squared-exponential (radial basis function) kernel for [`HilbertSpaceGP`](@ref).

Its one-dimensional spectral density is

```math
S(\omega) = \sigma^2 \sqrt{2\pi}\, \ell \, \exp\!\Big(-\tfrac{1}{2}\ell^2\omega^2\Big),
```

giving infinitely differentiable (very smooth) sample paths.
"
struct SquaredExponentialKernel <: AbstractGPKernel end

@doc raw"
Matérn-3/2 kernel for [`HilbertSpaceGP`](@ref).

Its one-dimensional spectral density is

```math
S(\omega) = \sigma^2 \, \frac{4 \,\nu_\ell^3}{(\nu_\ell^2 + \omega^2)^2},
\qquad \nu_\ell = \frac{\sqrt 3}{\ell},
```

giving once-differentiable, rougher sample paths than the squared-exponential
kernel — useful when the latent process is less smooth.
"
struct Matern32Kernel <: AbstractGPKernel end

@doc raw"
Matérn-5/2 kernel for [`HilbertSpaceGP`](@ref).

Its one-dimensional spectral density is

```math
S(\omega) = \sigma^2 \, \frac{16}{3}\, \frac{\nu_\ell^5}{(\nu_\ell^2 + \omega^2)^3},
\qquad \nu_\ell = \frac{\sqrt 5}{\ell},
```

giving twice-differentiable sample paths — an intermediate smoothness between
Matérn-3/2 and the squared-exponential kernel.
"
struct Matern52Kernel <: AbstractGPKernel end

@doc raw"
Spectral density ``S(\omega)`` of a [`HilbertSpaceGP`](@ref) covariance kernel at
frequency `ω`, for marginal standard deviation `σ` and length scale `ℓ`.

`ω` may be a scalar or a vector (the call broadcasts). The squared-exponential,
Matérn-3/2 and Matérn-5/2 one-dimensional spectral densities are
[riutortmayol2023practical](@citep)

```math
S_{\mathrm{SE}}(\omega) = \sigma^2 \sqrt{2\pi}\, \ell \,
    \exp\!\Big(-\tfrac{1}{2}\ell^2\omega^2\Big), \qquad
S_{3/2}(\omega) = \sigma^2 \frac{4\nu_\ell^3}{(\nu_\ell^2 + \omega^2)^2}, \qquad
S_{5/2}(\omega) = \sigma^2 \frac{16}{3}\frac{\nu_\ell^5}{(\nu_\ell^2 + \omega^2)^3},
```

with ``\nu_\ell = \sqrt{2p+1}/\ell`` for Matérn order ``p``. [`HilbertSpaceGP`](@ref)
weights each basis function by ``\sqrt{S(\sqrt{\lambda_j})}``.
"
function spectral_density(::SquaredExponentialKernel, ω, σ, ℓ)
    return σ^2 * sqrt(2π) * ℓ .* exp.(-(ℓ^2 / 2) .* ω .^ 2)
end

function spectral_density(::Matern32Kernel, ω, σ, ℓ)
    ν = sqrt(3) / ℓ
    return σ^2 * 4 * ν^3 ./ (ν^2 .+ ω .^ 2) .^ 2
end

function spectral_density(::Matern52Kernel, ω, σ, ℓ)
    ν = sqrt(5) / ℓ
    return σ^2 * (16 / 3) * ν^5 ./ (ν^2 .+ ω .^ 2) .^ 3
end

@doc raw"
Spectral density of the squared-exponential kernel at frequency ``\omega``.

```math
S(\omega) = \sigma^2 \sqrt{2\pi}\, \ell \, \exp\!\Big(-\tfrac{1}{2}\ell^2\omega^2\Big)
```

with marginal standard deviation ``\sigma`` and length scale ``\ell``. A thin
convenience wrapper over `spectral_density(SquaredExponentialKernel(), ω, σ, ℓ)`,
kept for the worked examples that check the approximation against the
squared-exponential kernel directly.
"
function se_spectral_density(ω, σ, ℓ)
    return spectral_density(SquaredExponentialKernel(), ω, σ, ℓ)
end

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

and ``S`` is the spectral density of the chosen covariance `kernel`. The kernel
controls the smoothness of the prior: a [`SquaredExponentialKernel`](@ref) (the
default) gives very smooth paths, while [`Matern32Kernel`](@ref) /
[`Matern52Kernel`](@ref) give progressively rougher ones. Only the spectral
density changes between kernels; the basis is shared.

Only ``\ell``, ``\sigma`` and the ``m`` weights ``\beta`` are sampled; the basis
``\phi_j`` and eigenvalues ``\lambda_j`` depend only on `n`, `m` and the boundary
factor `c`, **not** on any sampled parameter. They are therefore built **once**
when [`as_turing_model`](@ref) is called (outside the `@model`, so outside the
per-gradient-evaluation path) and captured by the model; each log-density
evaluation only reweights and combines them. The latent path is a cheap
matrix–vector product of a fixed basis against a small set of standard-normal
weights — a non-centred parameterisation that is fast and samples well under
NUTS, including with [Mooncake](https://chalk-lab.github.io/Mooncake.jl/)
reverse-mode AD.

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
  - `kernel`: the covariance kernel (an [`AbstractGPKernel`](@ref)).

# Examples
```@example HilbertSpaceGP
using EpiAwarePrototype, Distributions
gp = HilbertSpaceGP()
mdl = as_turing_model(gp, 30)
rand(mdl)
```

A rougher prior with a Matérn-3/2 kernel:
```@example HilbertSpaceGP
gp_matern = HilbertSpaceGP(kernel = Matern32Kernel())
length(as_turing_model(gp_matern, 30)())
```
"
struct HilbertSpaceGP{L <: Sampleable, S <: Sampleable, K <: AbstractGPKernel} <:
       AbstractLatentModel
    "Prior distribution for the length scale ``\\ell``."
    length_scale_prior::L
    "Prior distribution for the marginal standard deviation ``\\sigma``."
    marginal_std_prior::S
    "Number of basis functions."
    m::Int
    "Boundary factor: the GP is approximated on ``[-L, L]`` with ``L = c S``."
    c::Float64
    "Covariance kernel (an [`AbstractGPKernel`](@ref))."
    kernel::K

    function HilbertSpaceGP(length_scale_prior::Sampleable,
            marginal_std_prior::Sampleable, m::Int, c::Real,
            kernel::AbstractGPKernel)
        @assert m>0 "m (the number of basis functions) must be greater than 0"
        @assert c>1 "c (the boundary factor) must be greater than 1"
        new{typeof(length_scale_prior), typeof(marginal_std_prior), typeof(kernel)}(
            length_scale_prior, marginal_std_prior, m, Float64(c), kernel)
    end
end

function HilbertSpaceGP(;
        length_scale_prior::Sampleable = truncated(
            Normal(0.0, 0.4), 0, Inf),
        marginal_std_prior::Sampleable = truncated(Normal(0.0, 1.0), 0, Inf),
        m::Int = 20, c::Real = 1.5,
        kernel::AbstractGPKernel = SquaredExponentialKernel())
    return HilbertSpaceGP(length_scale_prior, marginal_std_prior, m, c, kernel)
end

@doc raw"
Build the Hilbert-space GP basis for `n` evenly spaced inputs.

Returns `(Φ, sqrt_λ)` where `Φ` is the ``n \times m`` matrix of eigenfunctions
``\phi_j`` evaluated at the rescaled inputs and `sqrt_λ` is the length-`m` vector
of ``\sqrt{\lambda_j}``. The inputs ``t = 1, \ldots, n`` are centred and rescaled
to unit spacing about their midpoint, so the half-range is ``S = (n-1)/2`` and the
GP is approximated on ``[-L, L]`` with ``L = c S``. Both outputs depend only on
`n`, `m` and `c` — none of the sampled parameters — so [`HilbertSpaceGP`](@ref)
calls this once per model construction, outside the differentiated per-evaluation
path. Requires `n > 1` so the half-range ``S`` is positive.
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

# Inner Turing model over a PREBUILT basis. Keeping the basis out of the `@model`
# body means it is computed once (in `as_turing_model` below) rather than on every
# log-density / gradient evaluation: only `ℓ`, `σ`, `β` and the matrix–vector
# product remain inside the differentiated path.
@model function _hsgp_model(kernel::AbstractGPKernel, Φ, sqrt_λ, m,
        length_scale_prior, marginal_std_prior)
    ℓ ~ length_scale_prior
    σ ~ marginal_std_prior
    β ~ filldist(Normal(), m)
    spectral_weights = sqrt.(spectral_density(kernel, sqrt_λ, σ, ℓ))
    gp = Φ * (spectral_weights .* β)
    return gp
end

function as_turing_model(model::HilbertSpaceGP, n)
    @assert n>1 "n must be greater than 1"
    Φ, sqrt_λ = hsgp_basis(n, model.m, model.c)
    return _hsgp_model(model.kernel, Φ, sqrt_λ, model.m,
        model.length_scale_prior, model.marginal_std_prior)
end
