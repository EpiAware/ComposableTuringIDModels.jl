# Hilbert-space approximate Gaussian-process latent model.
#
# The covariance kernels are the ecosystem-standard types from
# [KernelFunctions.jl](https://juliagaussianprocesses.github.io/KernelFunctions.jl/):
# `SqExponentialKernel`, `Matern32Kernel` and `Matern52Kernel`. KernelFunctions
# defines the kernels (and their Gram matrices, which the reconstruction tests use
# as ground truth) but not the *spectral densities* the Hilbert-space
# approximation needs, so this file adds a `spectral_density` method for each. New
# kernels plug in by adding a `spectral_density(::MyKernel, ω, σ, ℓ)` method — no
# other change to `HilbertSpaceGP` is required.

@doc raw"
Spectral density ``S(\omega)`` of a [`HilbertSpaceGP`](@ref) covariance `kernel` at
frequency `ω`, for marginal standard deviation `σ` and length scale `ℓ`.

The kernels are [KernelFunctions.jl](https://juliagaussianprocesses.github.io/KernelFunctions.jl/)
types. A kernel enters the Hilbert-space approximation only through this spectral
density — the Fourier transform of the stationary covariance — so switching kernel
just reweights the shared basis by ``\sqrt{S(\sqrt{\lambda_j})}``. Adding a new
kernel means adding a `spectral_density` method; nothing else changes.

`ω` may be a scalar or a vector (the call broadcasts). The one-dimensional
squared-exponential, Matérn-3/2 and Matérn-5/2 spectral densities are
[riutortmayol2023practical](@citep)

```math
S_{\mathrm{SE}}(\omega) = \sigma^2 \sqrt{2\pi}\, \ell \,
    \exp\!\Big(-\tfrac{1}{2}\ell^2\omega^2\Big), \qquad
S_{3/2}(\omega) = \sigma^2 \frac{4\nu_\ell^3}{(\nu_\ell^2 + \omega^2)^2}, \qquad
S_{5/2}(\omega) = \sigma^2 \frac{16}{3}\frac{\nu_\ell^5}{(\nu_\ell^2 + \omega^2)^3},
```

with ``\nu_\ell = \sqrt{2p+1}/\ell`` for Matérn order ``p``. The
squared-exponential (`SqExponentialKernel`) gives infinitely differentiable, very
smooth paths; `Matern32Kernel` (once-differentiable) and `Matern52Kernel`
(twice-differentiable) give progressively rougher ones. These three cover the
kernels offered by the EpiNow2 Gaussian-process implementation.
[`HilbertSpaceGP`](@ref) weights each basis function by ``\sqrt{S(\sqrt{\lambda_j})}``.

# Examples
```@example spectral_density
using ComposableTuringIDModels
spectral_density(SqExponentialKernel(), [0.0, 1.0, 2.0], 1.0, 0.8)
```
"
function spectral_density(::SqExponentialKernel, ω, σ, ℓ)
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
convenience wrapper over `spectral_density(SqExponentialKernel(), ω, σ, ℓ)`, kept
for the worked examples that check the approximation against the
squared-exponential kernel directly.
"
function se_spectral_density(ω, σ, ℓ)
    return spectral_density(SqExponentialKernel(), ω, σ, ℓ)
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

and ``S`` is the spectral density of the chosen covariance `kernel`. Kernels are
[KernelFunctions.jl](https://juliagaussianprocesses.github.io/KernelFunctions.jl/)
types, so the model reuses the ecosystem-standard kernels rather than defining its
own: `SqExponentialKernel` (the default) gives very smooth paths, while
`Matern32Kernel` / `Matern52Kernel` give progressively rougher ones. Only the
[`spectral_density`](@ref) changes between kernels; the basis is shared. See the
Gaussian-process case study for how this relates to `AbstractGPs.jl` /
`TemporalGPs.jl` as ecosystem alternatives.

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
the standardised inputs, so that boundary effects do not distort the fit). Because
the inputs are standardised to unit standard deviation (see [`hsgp_basis`](@ref)),
``\ell`` is scale-free — measured in standard deviations of the inputs, not raw
time steps — so a fixed `m` stays adequate as the series length changes. The
defaults (`m = 20`, `c = 1.5`) are a reasonable starting point for a smooth
latent process; short length scales relative to the standardised range may still
need a larger `m`.

## Fields

  - `length_scale_prior`: prior for the length scale ``\ell``.
  - `marginal_std_prior`: prior for the marginal standard deviation ``\sigma``.
  - `m`: number of basis functions.
  - `c`: boundary factor; the GP is approximated on ``[-L, L]`` with ``L = c S``.
  - `kernel`: the covariance kernel, a KernelFunctions.jl `Kernel` (default
    `SqExponentialKernel()`).

# Examples
```@example HilbertSpaceGP
using ComposableTuringIDModels, Distributions
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
struct HilbertSpaceGP{L <: Sampleable, S <: Sampleable, K <: Kernel} <:
       AbstractLatentModel
    "Prior distribution for the length scale ``\\ell``."
    length_scale_prior::L
    "Prior distribution for the marginal standard deviation ``\\sigma``."
    marginal_std_prior::S
    "Number of basis functions."
    m::Int
    "Boundary factor: the GP is approximated on ``[-L, L]`` with ``L = c S``."
    c::Float64
    "Covariance kernel, a KernelFunctions.jl `Kernel`."
    kernel::K

    function HilbertSpaceGP(length_scale_prior::Sampleable,
            marginal_std_prior::Sampleable, m::Int, c::Real,
            kernel::Kernel)
        @assert m>0 "m (the number of basis functions) must be greater than 0"
        @assert c>1 "c (the boundary factor) must be greater than 1"
        new{typeof(length_scale_prior), typeof(marginal_std_prior), typeof(kernel)}(
            length_scale_prior, marginal_std_prior, m, Float64(c), kernel)
    end
end

# Small positive floor on the default length scale. The prior is truncated at
# `ℓ_floor` rather than 0 because the Matérn spectral density stiffens as ℓ→0
# (ν = √(2p+1)/ℓ → ∞): a hard floor keeps ν finite and the sampler well-behaved.
# In standardised input units 0.05 is short (well below the ~√3 half-range) yet
# safely above the singular limit.
const _DEFAULT_LENGTH_SCALE_FLOOR = 0.05

function HilbertSpaceGP(;
        length_scale_prior::Sampleable = truncated(
            Normal(0.0, 0.4), _DEFAULT_LENGTH_SCALE_FLOOR, Inf),
        marginal_std_prior::Sampleable = truncated(Normal(0.0, 1.0), 0, Inf),
        m::Int = 20, c::Real = 1.5,
        kernel::Kernel = SqExponentialKernel())
    return HilbertSpaceGP(length_scale_prior, marginal_std_prior, m, c, kernel)
end

# Standardise the integer index 1:n to zero mean and unit standard deviation so
# the length scale ℓ is scale-free: the half-range then approaches √3 as n grows
# rather than scaling like (n-1)/2, keeping a short ℓ representable by a fixed
# number of basis functions m regardless of series length. Internal, but shared
# with the reconstruction tests so they build the exact kernel on the same
# coordinates the basis uses.
function _hsgp_standardised_index(n::Int)
    (collect(1:n) .- Statistics.mean(1:n)) ./ Statistics.std(1:n)
end

@doc raw"
Build the Hilbert-space GP basis for `n` evenly spaced inputs.

Returns `(Φ, sqrt_λ)` where `Φ` is the ``n \times m`` matrix of eigenfunctions
``\phi_j`` evaluated at the standardised inputs and `sqrt_λ` is the length-`m`
vector of ``\sqrt{\lambda_j}``. The integer indices ``t = 1, \ldots, n`` are
**standardised** to zero mean and unit standard deviation, so the length scale
``\ell`` is scale-free: the half-range ``S = \max_i |x_i|`` approaches
``\sqrt 3`` as ``n`` grows rather than scaling like ``(n-1)/2``, and the GP is
approximated on ``[-L, L]`` with ``L = c S``. Standardising keeps a fixed
``\ell`` prior (and a fixed `m`) meaningful across series lengths: ``\ell`` is
measured in standard deviations of the inputs, not raw time steps. Both outputs
depend only on `n`, `m` and `c` — none of the sampled parameters — so
[`HilbertSpaceGP`](@ref) calls this once per model construction, outside the
differentiated per-evaluation path. Requires `n > 1` so the standard deviation
(and hence ``S``) is positive.
"
function hsgp_basis(n::Int, m::Int, c::Real)
    @assert n>1 "n must be greater than 1 for a well-defined basis (S > 0)"
    x = _hsgp_standardised_index(n)
    S = maximum(abs, x)          # half-range of the standardised inputs
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
@model function _hsgp_model(kernel::Kernel, Φ, sqrt_λ, m,
        length_scale_prior, marginal_std_prior)
    ℓ ~ length_scale_prior
    σ ~ marginal_std_prior
    β ~ filldist(Normal(), m)
    spectral_weights = sqrt.(spectral_density(kernel, sqrt_λ, σ, ℓ))
    gp = Φ * (spectral_weights .* β)
    return gp
end

# Architecture note: CLAUDE.md's directive is "one `@model function
# as_turing_model(m::MyModel, ...)` per struct". Here `as_turing_model` is
# deliberately a *plain* function that builds the fixed basis once and then
# delegates to the inner `@model _hsgp_model`. This keeps the basis construction
# out of the differentiated per-evaluation path while preserving the single
# `as_turing_model(model, n)` entry point; the `@model` is an implementation
# detail of that one method, not a second public model per struct.
function as_turing_model(model::HilbertSpaceGP, n)
    @assert n>1 "n must be greater than 1"
    Φ, sqrt_λ = hsgp_basis(n, model.m, model.c)
    return _hsgp_model(model.kernel, Φ, sqrt_λ, model.m,
        model.length_scale_prior, model.marginal_std_prior)
end
