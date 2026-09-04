# Safe count distributions wrap `Distributions.Poisson` and `NegativeBinomial` to
# avoid `InexactError` at very large means.
# They declare a `SafeIntValued` support so `eltype` stays integer-typed inside a
# Turing model.

const SafeInt = Union{Int, BigInt}

@doc raw"
A value-support tag for real-valued count distributions whose `eltype` must stay
integer-typed inside a `Turing` model even when `rand` is called.
"
struct SafeIntValued <: Distributions.ValueSupport end
function Base.eltype(::Type{<:Distributions.Sampleable{F, SafeIntValued}}) where {F}
    return SafeInt
end

const SafeDiscreteUnivariateDistribution = Distributions.Distribution{
    Distributions.Univariate, SafeIntValued,
}

@doc raw"
Guard shared by every entry point that floors a `SafePoisson` /
`SafeNegativeBinomial` parameter, quantile or summary statistic to an
integer. Raises a clear `DomainError` naming the offending value for a
non-finite input, instead of letting the eventual `floor`/`Int` conversion
fail with an opaque `InexactError` deep in `Distributions.jl` or `Base`.
"
function _require_finite_for_int(x::Real, caller::AbstractString = "_safe_int_floor")
    isfinite(x) || throw(
        DomainError(
            x,
            "$caller: cannot floor a non-finite value to an integer. " *
                "This usually means an upstream rate or mean was Inf or NaN " *
                "(e.g. an explosive prior draw feeding a SafePoisson or " *
                "SafeNegativeBinomial rate)."
        )
    )
    return nothing
end

function _safe_int_floor(x::Real)
    _require_finite_for_int(x)
    Tf = typeof(x)
    if (Tf(typemin(Int)) - one(Tf)) < x < (Tf(typemax(Int)) + one(Tf))
        return floor(Int, x)
    else
        return floor(BigInt, x)
    end
end
