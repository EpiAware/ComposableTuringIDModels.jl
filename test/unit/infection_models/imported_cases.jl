# Tests for the ImportedCases renewal modifier (issue #189).

@testitem "ImportedCases struct carries the importation prior" begin
    using ComposableTuringIDModels, Distributions
    # Constant unknown rate (Distribution)
    ic = ImportedCases(Normal(0.0, 0.1))
    @test ic.importation_rate isa Normal
    @test ic isa ComposableTuringIDModels.AbstractRenewalModifier

    # Time-varying rate (latent process)
    ic_tv = ImportedCases(RandomWalk())
    @test ic_tv.importation_rate isa RandomWalk
end

@testitem "ImportedCases composes into a Renewal step" begin
    using ComposableTuringIDModels, Distributions
    ic = ImportedCases(Normal(0.0, 0.1))
    gen_int = [0.2, 0.3, 0.5]
    # Via the positional modifier constructor
    r = Renewal(gen_int, ic; rt = RandomWalk(), initialisation = Normal())
    @test r.recurrent_step isa ComposableTuringIDModels.RenewalStep
    @test only(r.recurrent_step.modifiers) isa ImportedCases
    @test r.recurrent_step.modifiers[1].importation_rate isa Normal

    # Composed with another modifier
    ic2 = ImportedCases(Normal(0.0, 0.1))
    r2 = Renewal(gen_int, ic2,
        ComposableTuringIDModels.SusceptibleDepletion(1000.0);
        rt = RandomWalk(), initialisation = Normal())
    @test length(r2.recurrent_step.modifiers) == 2
    @test r2.recurrent_step.modifiers[1] isa ImportedCases
    @test r2.recurrent_step.modifiers[2] isa
          ComposableTuringIDModels.SusceptibleDepletion
end

@testitem "ImportedCases with constant rate adds to incidence" begin
    using ComposableTuringIDModels, Distributions, Random
    using DynamicPPL: fix
    Random.seed!(1891)
    gen_int = [0.2, 0.3, 0.5]
    logR = log(1.0)  # Rt = 1, no growth

    # No importation: flat incidence.
    plain = Renewal(gen_int; rt = FixedIntercept(logR),
        initialisation = Normal())
    fixinit = (init_incidence = log(1.0),)
    I_plain = fix(as_turing_model(plain, 30), fixinit)().I_t

    # With importation: incidence grows.
    ic = ImportedCases(Normal(0.5, 0.01))
    imported = Renewal(gen_int, ic; rt = FixedIntercept(logR),
        initialisation = Normal())
    I_imported = fix(as_turing_model(imported, 30), fixinit)().I_t

    @test all(isfinite, I_imported)
    @test all(>=(0), I_imported)
    @test I_imported[end] > I_plain[end]  # Importation adds to incidence
end

@testitem "ImportedCases with time-varying rate samples via as_turing_submodel" begin
    using ComposableTuringIDModels, Distributions, Random
    Random.seed!(1892)
    gen_int = [0.2, 0.3, 0.5]
    # Time-varying importation via RandomWalk
    ic_tv = ImportedCases(RandomWalk())
    r = Renewal(gen_int, ic_tv; rt = RandomWalk(),
        initialisation = Normal())
    out = as_turing_model(r, 20)()
    @test length(out.I_t) == 20
    @test all(isfinite, out.I_t)
    @test all(>=(0), out.I_t)

    # The importation rate is a sampled parameter.
    draw = rand(as_turing_model(r, 20))
    @test any(k -> startswith(string(k), "import_rates"), keys(draw))
end

@testitem "ImportedCases with fixed rate samples under NUTS" tags=[:sample] begin
    using ComposableTuringIDModels, Distributions, Turing, Random
    Random.seed!(1893)
    gen_int = [0.2, 0.3, 0.5]
    ic = ImportedCases(Normal(0.0, 0.1))
    model = IDModel(
        Renewal(gen_int, ic; rt = RandomWalk(),
            initialisation = Normal()),
        PoissonError())
    y = as_turing_model(model, missing, 20)().generated_y_t
    chn = sample(as_turing_model(model, y, 20), NUTS(), 30;
        progress = false)
    @test chn !== nothing
end