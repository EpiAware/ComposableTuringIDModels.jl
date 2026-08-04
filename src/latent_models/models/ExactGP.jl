# Exact Gaussian-process latent model.
#
# The exact counterpart of [`HilbertSpaceGP`](@ref): rather than the fixed
# basis-function approximation, this forms the full covariance matrix from the
# same ecosystem-standard [KernelFunctions.jl](https://juliagaussianprocesses.github.io/KernelFunctions.jl/)
# kernel and factorises it. It is the reference an approximation is judged
# against — accurate but ``O(n^3)`` per evaluation. The two share the standardised
# input grid (`standardised_index`), so a given length scale means the same
# thing for both.

@doc raw"
An **exact Gaussian-process** latent process.

The exact counterpart of [`HilbertSpaceGP`](@ref). Where the Hilbert-space model
approximates the GP by a short weighted sum of fixed basis functions, this model
forms the full ``n \times n`` covariance matrix ``K`` from the covariance
`kernel` and draws the path from it directly, so it is the *exact* GP the
Hilbert-space model approximates:

```math
K_{ij} = \sigma^2\, k(x_i, x_j; \ell), \qquad
f = L z, \quad L L^\top = K + \tau\sigma^2 I,
\quad z_i \sim \mathrm{Normal}(0, 1).
```

The nugget ``\tau`` is *relative*: it scales with ``\sigma^2``, the diagonal of
``K``. A fixed absolute nugget is swamped once the sampler visits a large
``\sigma``, and the Cholesky factorisation then throws on a matrix that is only
numerically indefinite, ending the chain; a nugget with an absolute floor
instead dominates the covariance at small ``\sigma``. Scaling it leaves the
relative variance inflation at ``\tau`` throughout. A floor of
``\tau\,\varepsilon`` is added so the factorisation stays defined at
``\sigma = 0`` exactly, where ``K`` vanishes.

The path is drawn non-centred: standard-normal weights ``z`` are pushed through
the Cholesky factor ``L`` of the covariance. As with [`HilbertSpaceGP`](@ref)
this keeps only ``\ell``, ``\sigma`` and the length-`n` weights ``z`` sampled, a
parameterisation NUTS handles well. Unlike the Hilbert-space model, the
covariance and its Cholesky factorisation depend on the sampled ``\ell`` and
``\sigma``, so they are rebuilt on **every** log-density evaluation at ``O(n^3)``
cost. That is the price of exactness, and the reason the Hilbert-space
approximation exists; this model is the accuracy reference to compare it against,
best suited to short series.

Kernels are [KernelFunctions.jl](https://juliagaussianprocesses.github.io/KernelFunctions.jl/)
types (`SqExponentialKernel`, `Matern32Kernel`, `Matern52Kernel`, ...) — the same
kernels [`HilbertSpaceGP`](@ref) uses, and the ones
[AbstractGPs.jl](https://juliagaussianprocesses.github.io/AbstractGPs.jl/) builds
exact GPs from. The kernel sees the standardised index
[`standardised_index`](@ref), exactly as in [`HilbertSpaceGP`](@ref), so the
length scale ``\ell`` is scale-free and means the same thing for both models.
That grid depends only on `n`, so — again as in [`HilbertSpaceGP`](@ref) —
[`as_turing_model`](@ref) builds it once and captures it rather than rebuilding
it inside the `@model` body; only the covariance and its factorisation are
per-evaluation work.

## Fields

  - `length_scale`: prior for the length scale ``\ell``; it must put no mass
    below zero, since the covariance is not positive definite at
    ``\ell \le 0``. Checked at construction.
  - `marginal_std`: prior for the marginal standard deviation ``\sigma``; it
    must put no mass below zero. Checked at construction.
  - `kernel`: the covariance kernel, a KernelFunctions.jl `Kernel` (default
    `SqExponentialKernel()`).
  - `jitter`: relative diagonal nugget ``\tau`` for a stable Cholesky factor
    (default `1e-6`); the amount added is ``\tau\sigma^2``.

# Examples
```@example ExactGP
using ComposableTuringIDModels, Distributions
gp = ExactGP()
mdl = as_turing_model(gp, 30)
rand(mdl)
```

A rougher prior with a Matérn-3/2 kernel:
```@example ExactGP
gp_matern = ExactGP(kernel = Matern32Kernel())
length(as_turing_model(gp_matern, 30)())
```
"
struct ExactGP{L <: UnivariateDistribution, S <: UnivariateDistribution,
    K <: Kernel} <: AbstractPriorModel
    "Prior for the length scale ``\\ell``; puts no mass below zero."
    length_scale::L
    "Prior for ``\\sigma``; puts no mass below zero."
    marginal_std::S
    "Covariance kernel, a KernelFunctions.jl `Kernel`."
    kernel::K
    "Relative Cholesky nugget: the amount added is ``\\tau\\sigma^2``."
    jitter::Float64

    function ExactGP(length_scale::UnivariateDistribution,
            marginal_std::UnivariateDistribution, kernel::Kernel, jitter::Real)
        @assert jitter>0 "jitter must be greater than 0"
        _check_hyperprior_support(length_scale, marginal_std)
        new{typeof(length_scale), typeof(marginal_std), typeof(kernel)}(
            length_scale, marginal_std, kernel, Float64(jitter))
    end
end

function ExactGP(;
        length_scale::UnivariateDistribution = truncated(
            Normal(0.0, 0.4), _DEFAULT_LENGTH_SCALE_FLOOR, Inf),
        marginal_std::UnivariateDistribution = truncated(Normal(0.0, 1.0), 0, Inf),
        kernel::Kernel = SqExponentialKernel(), jitter::Real = 1e-6)
    return ExactGP(length_scale, marginal_std, kernel, jitter)
end

# Inner Turing model over a PREBUILT input grid. The grid depends only on `n`,
# so it is computed once in `as_turing_model` below rather than on every
# log-density / gradient evaluation, exactly as `HilbertSpaceGP` does with its
# basis. What is left inside the differentiated path is the part that genuinely
# depends on the sampled parameters: the covariance and its Cholesky factor.
@model function _exact_gp_model(kernel::Kernel, x, jitter,
        length_scale, marginal_std)
    ℓ ~ length_scale
    σ ~ marginal_std
    z ~ filldist(Normal(), length(x))
    # Scale the Gram matrix, not the kernel: `σ^2 * kernel` builds a
    # KernelFunctions `ScaledKernel`, which asserts σ² > 0 at construction.
    # NUTS can sample σ = 0.0 exactly (the prior floor, or an underflowed
    # extreme excursion), so scaling the matrix keeps that a valid, if
    # degenerate, draw instead of throwing mid-chain.
    K = σ^2 .* kernelmatrix(with_lengthscale(kernel, ℓ), x)
    # The nugget tracks the diagonal of K, which is σ². A fixed absolute nugget
    # is swamped once σ is large and the factorisation then fails on a matrix
    # that is only numerically indefinite, which ends the chain; a nugget with a
    # floor of `jitter` instead swamps the covariance when σ is small (at
    # σ = 1e-3 a `jitter * (σ² + 1)` nugget inflates the variance by 100%). The
    # absolute floor here is `jitter * eps()`, far below σ² for any σ the
    # sampler can represent, and is only there to keep the factorisation
    # defined at σ = 0 exactly, where K itself vanishes.
    L = cholesky(Symmetric(K + (jitter * (σ^2 + eps())) * I)).L
    gp = L * z
    return gp
end

# See the architecture note in HilbertSpaceGP.jl: `as_turing_model` is a plain
# function that hoists the parameter-independent work and delegates to a single
# inner `@model`, keeping the one `as_turing_model(model, n)` entry point.
function as_turing_model(model::ExactGP, n)
    @assert n>1 "n must be greater than 1"
    x = standardised_index(n)
    return _exact_gp_model(model.kernel, x, model.jitter,
        model.length_scale, model.marginal_std)
end
