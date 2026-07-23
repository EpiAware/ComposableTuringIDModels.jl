struct ImportedCases{I} <: AbstractRenewalModifier
    importation_rate::I
    function ImportedCases(importation_rate)
        return new{typeof(importation_rate)}(importation_rate)
    end
end
modifier_init_state(::ImportedCases) = 0.0
apply_modifier(::ImportedCases, inc, _) = (inc, 0.0)
