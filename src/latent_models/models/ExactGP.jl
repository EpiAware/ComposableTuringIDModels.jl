# Exact Gaussian-process latent model.
#
# The exact counterpart of [`HilbertSpaceGP`](@ref): rather than the fixed
# basis-function approximation, this forms the full covariance matrix from the
# same ecosystem-standard [KernelFunctions.jl](https://juliagaussianprocesses.github.io/KernelFunctions.jl/)
# kernel and factorises it. It is the reference an approximation is judged
# against — accurate but ``O(n^3)`` per evaluation. The two share the standardised
# input grid ([`_hsgp_standardised_index`](@ref)), so a given length scale means
# the same thing for both.

@doc raw"
An **exact Gaussian-process** latent process.

The exact counterpart of [`HilbertSpaceGP`](@ref). Where the Hilbert-space model
approximates the GP by a short weighted sum of fixed basis functions, this model
forms the full ``n \times n`` covariance matrix ``K`` from the covariance
`kernel` and draws the path from it directly, so it is the *exact* GP the
Hilbert-space model approximates:

```math
K_{ij} = \sigma^2\, k(x_i, x_j; \ell), \qquad
f = L z, \quad L L^\top = K + \tau I, \quad z_i \sim \mathrm{Normal}(0, 1).
```

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
exact GPs from. The inputs are standardised with the shared
[`_hsgp_standardised_index`](@ref), so the length scale ``\ell`` is scale-free and
means the same thing here and in the Hilbert-space model.

## Fields

  - `length_scale`: prior for the length scale ``\ell``.
  - `marginal_std`: prior for the marginal standard deviation ``\sigma``.
  - `kernel`: the covariance kernel, a KernelFunctions.jl `Kernel` (default
    `SqExponentialKernel()`).
  - `jitter`: diagonal ``\tau`` added to ``K`` for a stable Cholesky factor
    (default `1e-6`).

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
struct ExactGP{L <: Sampleable, S <: Sampleable, K <: Kernel} <:
       AbstractLatentModel
    "Prior distribution for the length scale ``\\ell``."
    length_scale::L
    "Prior distribution for the marginal standard deviation ``\\sigma``."
    marginal_std::S
    "Covariance kernel, a KernelFunctions.jl `Kernel`."
    kernel::K
    "Diagonal jitter added to the covariance for a stable Cholesky factor."
    jitter::Float64

    function ExactGP(length_scale::Sampleable, marginal_std::Sampleable,
            kernel::Kernel, jitter::Real)
        @assert jitter>0 "jitter must be greater than 0"
        new{typeof(length_scale), typeof(marginal_std), typeof(kernel)}(
            length_scale, marginal_std, kernel, Float64(jitter))
    end
end

function ExactGP(;
        length_scale::Sampleable = truncated(
            Normal(0.0, 0.4), _DEFAULT_LENGTH_SCALE_FLOOR, Inf),
        marginal_std::Sampleable = truncated(Normal(0.0, 1.0), 0, Inf),
        kernel::Kernel = SqExponentialKernel(), jitter::Real = 1e-6)
    return ExactGP(length_scale, marginal_std, kernel, jitter)
end

@model function as_turing_model(model::ExactGP, n)
    @assert n>1 "n must be greater than 1"
    ℓ ~ model.length_scale
    σ ~ model.marginal_std
    z ~ filldist(Normal(), n)
    x = _hsgp_standardised_index(n)
    K = kernelmatrix(σ^2 * with_lengthscale(model.kernel, ℓ), x)
    L = cholesky(Symmetric(K + model.jitter * I)).L
    gp = L * z
    return gp
end
