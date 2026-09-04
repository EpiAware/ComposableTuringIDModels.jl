# Post-inference tidy helper: spread a `Chains` object into a tidy `DataFrame`.

@doc raw"
Convert an `MCMCChains.Chains` object to a tidy `DataFrame` (one row per draw,
with `draw`, `chain`, and `iteration` columns).

# Arguments

  - `chn`: the `Chains` object to convert.

# Examples
```@example spread_draws
using ComposableTuringIDModels
nothing
```
"
function spread_draws(chn::Chains)
    df = DataFrame(chn)
    # Current MCMCChains emits the bookkeeping columns as `.iteration` and
    # `.chain`, and older versions used undotted names, so accept either.
    # A sequential `draw` index is added in tidybayes style.
    for (dotted, plain) in ((".iteration", "iteration"), (".chain", "chain"))
        if dotted in names(df)
            @rename!(df, $(plain) = $(dotted))
        end
    end
    df = hcat(DataFrame(draw = 1:size(df, 1)), df)
    return df
end
