# Reshape a `Chains` object into a `(draws x chains)` array of `NamedTuple`s.

@doc raw"
Reshape an `MCMCChains.Chains` object into a `(draws × chains)` array of
per-sample `NamedTuple`s.

# Arguments

  - `chn`: the `Chains` object.

# Examples
```@example get_param_array
using ComposableTuringIDModels
nothing
```
"
function get_param_array(chn::Chains)
    return rowtable(chn) |> x -> reshape(x, size(chn, 1), size(chn, 3))
end
