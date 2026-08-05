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
(twice-differentiable) give progressively rougher ones — the three choices most
often used for a smooth epidemiological latent process.
[`HilbertSpaceGP`](@ref) weights each basis function by ``\sqrt{S(\sqrt{\lambda_j})}``.

Requires ``\ell > 0`` and ``\sigma \ge 0``. Both are asserted rather than
allowed through: at ``\ell = 0`` the Matérn densities are ``\infty/\infty`` and
return `NaN`, which would silently propagate into the basis weights. A new
method should assert the same.

# Examples
```@example spectral_density
using ComposableTuringIDModels
spectral_density(SqExponentialKernel(), [0.0, 1.0, 2.0], 1.0, 0.8)
```
"
function spectral_density(::SqExponentialKernel, ω, σ, ℓ)
    _check_spectral_args(σ, ℓ)
    return σ^2 * sqrt(2π) * ℓ .* exp.(-(ℓ^2 / 2) .* ω .^ 2)
end

function spectral_density(::Matern32Kernel, ω, σ, ℓ)
    _check_spectral_args(σ, ℓ)
    ν = sqrt(3) / ℓ
    return σ^2 * 4 * ν^3 ./ (ν^2 .+ ω .^ 2) .^ 2
end

function spectral_density(::Matern52Kernel, ω, σ, ℓ)
    _check_spectral_args(σ, ℓ)
    ν = sqrt(5) / ℓ
    return σ^2 * (16 / 3) * ν^5 ./ (ν^2 .+ ω .^ 2) .^ 3
end

# Shared guard for the `spectral_density` methods. A zero length scale makes the
# Matérn densities NaN (ν = √(2p+1)/ℓ = Inf, then Inf^p/Inf^p), so reject it
# with a diagnostic rather than let a NaN basis reach the sampler. The two GP
# models reject a hyperprior that could reach here at construction (see
# `_check_hyperprior_support`), so for them this is a backstop for a direct call
# rather than something a chain can hit.
function _check_spectral_args(σ, ℓ)
    @assert ℓ>0 "the length scale ℓ must be greater than 0"
    @assert σ>=0 "the marginal standard deviation σ must not be negative"
    return nothing
end

# Shared construction-time guard for the two GP latent models. Both need ℓ > 0
# (no spectral density, and a singular covariance, at ℓ ≤ 0) and σ ≥ 0. Checking
# the *support* at construction rejects an unusable hyperprior — say
# `length_scale = Normal()` — up front, rather than letting the first negative
# proposal abort a chain a long way from its cause. The bound is ≥ 0 rather than
# > 0 so that a prior with an open lower limit at zero (`Gamma`, `LogNormal`) is
# accepted: it puts no mass on ℓ = 0, and `_check_spectral_args` is the backstop
# for the measure-zero case.
function _check_hyperprior_support(length_scale, marginal_std)
    ℓ_msg = "the length scale prior must not put mass on ℓ < 0"
    σ_msg = "the marginal standard deviation prior must not put mass on σ < 0"
    @assert minimum(length_scale)>=0 ℓ_msg
    @assert minimum(marginal_std)>=0 σ_msg
    return nothing
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
[`spectral_density`](@ref) changes between kernels; the basis is shared. The
Gaussian-process case study checks the basis against the Gram matrix
[AbstractGPs.jl](https://juliagaussianprocesses.github.io/AbstractGPs.jl/) builds
from the same kernel, and compares this model against [`ExactGP`](@ref).

Only ``\ell``, ``\sigma`` and the ``m`` weights ``\beta`` are sampled; the basis
``\phi_j`` and eigenvalues ``\lambda_j`` depend only on `n`, `m` and the boundary
factor `c`, **not** on any sampled parameter. [`as_turing_model`](@ref) therefore
builds them outside the `@model` body and captures them, so nothing in the basis
is differentiated and each log-density evaluation only reweights and combines a
fixed matrix. (Composed inside another model — a [`Renewal`](@ref) whose `rt`
slot is this GP — the enclosing `@model` reconstructs its submodels on every
evaluation, so the basis is rebuilt *inside* the traced call. It still depends
on no sampled parameter, but it is executed and taped like the rest of the body,
so it is not free; the Gaussian-process case study measures what that costs.)
The latent path is a cheap matrix–vector product of a fixed basis against a
small set of standard-normal weights — a non-centred parameterisation that is
fast and samples well under NUTS, including with
[Mooncake](https://chalk-lab.github.io/Mooncake.jl/) reverse-mode AD.

The accuracy/speed trade-off is controlled by two numbers
[riutortmayol2023practical](@citep): the number of basis functions `m` (more
basis functions resolve shorter length scales, at linear cost) and the boundary
factor `c` (the domain is extended to ``L = c\,S`` beyond the half-range ``S`` of
the standardised inputs, so that boundary effects do not distort the fit).
Because the inputs are standardised to unit standard deviation (see
[`standardised_index`](@ref)), ``\ell`` is scale-free — measured in standard
deviations of the inputs, not raw time steps — so a fixed `m` stays adequate as
the series length changes.

The two act at opposite ends of the length-scale range, and each sets the floor
on accuracy where the other cannot help. A **short** ``\ell`` needs a larger `m`:
the basis has to resolve wiggles finer than ``2L/m``. A **long** ``\ell`` —
comparable to the standardised half-range ``S \approx \sqrt 3`` — needs a larger
`c`, because the approximation is periodic on ``[-L, L]`` and a slowly varying
path feels that boundary; adding basis functions does nothing for it. Below
`c = 1.2` that boundary error cannot be cleared by any `m`, so the constructor
rejects it.

!!! warning \"The defaults are tuned for the squared-exponential kernel\"
    With `m = 20` and `c = 1.5` the squared-exponential covariance is
    reconstructed to a few parts in ten thousand for ``\ell`` between roughly
    0.3 and 0.5, degrading to 1.6% at ``\ell = 0.2`` and 26% at ``\ell = 0.1``
    (where `m` is what fixes it: `m = 60` brings ``\ell = 0.1`` back to three
    parts in ten thousand) and to 1.5% at ``\ell = 0.8`` (where `c` sets the
    floor instead). **The Matérn kernels are markedly worse at the same `m`** —
    their spectral density has an algebraic rather than Gaussian tail, so the
    truncated basis discards more of it. At ``\ell = 0.3`` the `Matern32Kernel`
    error is about 2.5%, roughly 110 times the squared-exponential figure, and
    at ``\ell = 0.5`` it is still 0.6%; `Matern52Kernel` sits between the two.

    A truncated basis loses variance rather than adding it, so where `m` is too
    small the process is systematically **under-dispersed**: at `m = 20` the
    implied marginal standard deviation against a nominal ``\sigma = 1`` is
    0.99/0.96/0.97 (squared-exponential / Matérn-3/2 / Matérn-5/2) at
    ``\ell = 0.2`` and 0.89/0.84/0.86 at ``\ell = 0.1``, which biases the
    inferred ``\sigma`` upward without being visible in a plot of the fit. The
    default ``\ell`` prior puts roughly a fifth of its mass below 0.15, so raise
    `m` for a fit that settles on a short length scale, and raise it further for
    a Matérn kernel. The prior's floor of 0.05 is a numerical guard on the
    spectral density, not a claim that a 20-function basis resolves that scale.

## Fields

  - `length_scale`: prior for the length scale ``\ell``; it must put no mass
    below zero, since ``\ell \le 0`` has no spectral density. Checked at
    construction.
  - `marginal_std`: prior for the marginal standard deviation ``\sigma``; it
    must put no mass below zero. Checked at construction.
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
struct HilbertSpaceGP{L <: UnivariateDistribution, S <: UnivariateDistribution,
    K <: Kernel} <: AbstractPriorModel
    "Prior for the length scale ``\\ell``; puts no mass below zero."
    length_scale::L
    "Prior for ``\\sigma``; puts no mass below zero."
    marginal_std::S
    "Number of basis functions."
    m::Int
    "Boundary factor: the GP is approximated on ``[-L, L]`` with ``L = c S``."
    c::Float64
    "Covariance kernel, a KernelFunctions.jl `Kernel`."
    kernel::K

    function HilbertSpaceGP(length_scale::UnivariateDistribution,
            marginal_std::UnivariateDistribution, m::Int, c::Real,
            kernel::Kernel)
        @assert m>0 "m (the number of basis functions) must be greater than 0"
        # Below c = 1.2 the boundary effect leaves an error floor no number of
        # basis functions can clear (24% at c = 1.05, 6% at c = 1.2, measured as
        # relative Frobenius error against the Gram matrix), so the constructor
        # rejects it rather than accepting an unusable configuration.
        @assert c>=1.2 "c (the boundary factor) must be at least 1.2"
        _check_hyperprior_support(length_scale, marginal_std)
        new{typeof(length_scale), typeof(marginal_std), typeof(kernel)}(
            length_scale, marginal_std, m, Float64(c), kernel)
    end
end

# Small positive floor on the default length scale. The prior is truncated at
# `ℓ_floor` rather than 0 because the Matérn spectral density stiffens as ℓ→0
# (ν = √(2p+1)/ℓ → ∞): a hard floor keeps ν finite and the sampler well-behaved.
# In standardised input units 0.05 is short (well below the ~√3 half-range) yet
# safely above the singular limit.
const _DEFAULT_LENGTH_SCALE_FLOOR = 0.05

function HilbertSpaceGP(;
        length_scale::UnivariateDistribution = truncated(
            Normal(0.0, 0.4), _DEFAULT_LENGTH_SCALE_FLOOR, Inf),
        marginal_std::UnivariateDistribution = truncated(Normal(0.0, 1.0), 0, Inf),
        m::Int = 20, c::Real = 1.5,
        kernel::Kernel = SqExponentialKernel())
    return HilbertSpaceGP(length_scale, marginal_std, m, c, kernel)
end

@doc raw"
Standardise the integer index ``1:n`` to zero mean and unit standard deviation.

This is the input grid both [`HilbertSpaceGP`](@ref) and [`ExactGP`](@ref) hand
to their covariance kernel, so a given length scale ``\ell`` means the same
thing for both. It also makes ``\ell`` scale-free: the half-range approaches
``\sqrt 3`` as ``n`` grows rather than scaling like ``(n-1)/2``, so a short
``\ell`` stays representable by a fixed number of basis functions ``m``
regardless of series length.

It is public so that a comparison against another Gaussian-process
implementation — an [AbstractGPs.jl](https://juliagaussianprocesses.github.io/AbstractGPs.jl/)
`GP`, say — can be built on the same coordinates rather than on a
reimplementation of this formula. Requires `n > 1`, since the standard deviation
of a single point is zero.

# Arguments

  - `n`: the series length.

# Examples
```@example standardised_index
using ComposableTuringIDModels
standardised_index(5)
```
"
function standardised_index(n::Int)
    @assert n>1 "n must be greater than 1 to standardise the index"
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
[`HilbertSpaceGP`](@ref) calls this once per model construction rather than from
inside the `@model` body; see [`HilbertSpaceGP`](@ref) for what that leaves when
the GP is composed inside another model. Requires `n > 1` so the standard
deviation (and hence ``S``) is positive; `n = 1` would give ``S = L = 0`` and a
basis of `NaN`.
"
function hsgp_basis(n::Int, m::Int, c::Real)
    @assert n>1 "n must be greater than 1 for a well-defined basis (S > 0)"
    x = standardised_index(n)
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
        length_scale, marginal_std)
    ℓ ~ length_scale
    σ ~ marginal_std
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
function as_turing_model(model::HilbertSpaceGP, n::Int)
    @assert n>1 "n must be greater than 1"
    Φ, sqrt_λ = hsgp_basis(n, model.m, model.c)
    return _hsgp_model(model.kernel, Φ, sqrt_λ, model.m,
        model.length_scale, model.marginal_std)
end
