# Latent-model utility: expand a vector of distributions into a product
# distribution.

@doc raw"
Expand a vector of distributions into a single product distribution.

If every element of `dist` is equal, a `filldist` is returned for efficiency;
otherwise an `arraydist` over the heterogeneous vector is returned.
"
function _expand_dist(dist::Vector{D} where {D <: Distribution})
    d = length(dist)
    product_dist = all(first(dist) .== dist) ?
                   filldist(first(dist), d) : arraydist(dist)
    return product_dist
end
