# Tests for the component-tree `show`. A composed model must render as a
# recursive, indented tree of role + concrete component name, never the nested
# parametric type signature.

@testitem "show renders a composed model as a component tree" begin
    using ComposableTuringIDModels, Distributions

    data = IDData([0.2, 0.3, 0.5], exp)
    model = IDModel(
        Renewal(data; rt = RandomWalk(), initialisation_prior = Normal()),
        NegativeBinomialError())
    out = sprint(show, MIME"text/plain"(), model)

    # Root and each nested component appear by concrete name.
    for name in ("IDModel", "Renewal", "RandomWalk", "NegativeBinomialError")
        @test occursin(name, out)
    end
    # Roles are labelled and the two top-level branches use tree connectors.
    @test occursin("├─ infection: Renewal", out)
    @test occursin("└─ observation: NegativeBinomialError", out)
    # The latent slot is recursed into, indented under its infection branch.
    @test occursin("rt: RandomWalk", out)
    @test occursin("│  ", out)   # continuation indentation

    # The raw parametric type signature must NOT leak into the display.
    @test !occursin("{", out)
    @test !occursin("BroadcastPrior", out)
    @test !occursin("ConstantRenewalStep", out)
end

@testitem "show recurses through nested observation modifiers" begin
    using ComposableTuringIDModels, Distributions

    obs = Ascertainment(PoissonError(), FixedIntercept(0.1))
    out = sprint(show, MIME"text/plain"(), obs)

    @test startswith(out, "Ascertainment")
    @test occursin("├─ model: PoissonError", out)
    # The latent modifier is wrapped and recursed into.
    @test occursin("latent: PrefixLatentModel", out)
    @test occursin("FixedIntercept", out)
    @test !occursin("{", out)
end

@testitem "compact 2-arg show prints only the component name" begin
    using ComposableTuringIDModels, Distributions

    model = IDModel(
        DirectInfections(; Z = RandomWalk(), initialisation_prior = Normal()),
        PoissonError())
    # The compact form (used inside arrays / `repr`) is a single clean line with
    # no nested parametric type dump.
    @test sprint(show, model) == "IDModel"
    @test sprint(show, RandomWalk()) == "RandomWalk"
    @test !occursin("{", repr(model))
end
