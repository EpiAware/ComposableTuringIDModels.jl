# Safe count distributions wrap `Distributions.Poisson` / `NegativeBinomial` to
# avoid `InexactError` at very large means and declare a `SafeIntValued` support
# so `eltype` stays integer-typed inside a Turing model.

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

function _safe_int_floor(x::Real)
    isfinite(x) || throw(
        DomainError(
            x,
            "_safe_int_floor: cannot floor a non-finite value to an integer. " *
                "This usually means an upstream rate or mean was Inf or NaN " *
                "(e.g. an explosive prior draw feeding a SafePoisson or " *
                "SafeNegativeBinomial rate)."
        )
    )
    Tf = typeof(x)
    if (Tf(typemin(Int)) - one(Tf)) < x < (Tf(typemax(Int)) + one(Tf))
        return floor(Int, x)
    else
        return floor(BigInt, x)
    end
end
