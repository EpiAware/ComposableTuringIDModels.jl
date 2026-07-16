# Tests for the component-tree `show`. A composed model must render as a
# recursive, indented tree of role + concrete component name, never the nested
# parametric type signature.

@testitem "show renders a composed model as a component tree" begin
    using ComposableTuringIDModels, Distributions

    gen_int = [0.2, 0.3, 0.5]
    model = IDModel(
        Renewal(gen_int; rt = RandomWalk(), initialisation = Normal()),
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
    @test !occursin("ConstantRenewalStep", out)
end

@testitem "show recurses through a vector of component children" begin
    using ComposableTuringIDModels, Distributions

    # A manipulator holding a vector of latent models exercises the vector branch
    # of the child collector: each element is indexed and recursed into.
    model = ConcatLatentModels([Intercept(Normal(2, 0.2)), RandomWalk()])
    out = sprint(show, MIME"text/plain"(), model)
    @test occursin("ConcatLatentModels", out)
    # Each vector element is indexed and recursed into (the models are prefixed).
    @test occursin("models[1]: PrefixLatentModel", out)
    @test occursin("models[2]: PrefixLatentModel", out)
    @test occursin("model: Intercept", out)
    @test occursin("model: RandomWalk", out)
    # A raw-distribution prior slot stays a leaf even inside a vector child.
    @test !occursin("Normal{", out)
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

@testitem "show renders each Split stream on its own line" begin
    using ComposableTuringIDModels, Distributions

    split = Split((
        cases = LatentDelay(NegativeBinomialError(), [0.4, 0.3, 0.2, 0.1]),
        deaths = LatentDelay(PoissonError(), [0.1, 0.2, 0.3, 0.4])))
    out = sprint(show, MIME"text/plain"(), split)

    # Each stream is a labelled branch on its own line, keyed by its name, and the
    # inner model of each stream is recursed into beneath it.
    @test startswith(out, "Split")
    @test occursin("├─ cases: LatentDelay", out)
    @test occursin("└─ deaths: LatentDelay", out)
    @test occursin("model: NegativeBinomialError", out)
    @test occursin("model: PoissonError", out)
    # Two streams means two separate lines, not one collapsed leaf.
    @test count("LatentDelay", out) == 2
    @test !occursin("{", out)
end

@testitem "Split show composes with nesting in both directions" begin
    using ComposableTuringIDModels, Distributions

    # A Split nested inside a modifier, with a further submodel nested inside one
    # of its streams: indentation must compose through both layers.
    model = LatentDelay(
        Split((
            cases = PoissonError(),
            deaths = LatentDelay(
                Ascertainment(PoissonError(), FixedIntercept(log(0.1))),
                [0.2, 0.3, 0.5]))),
        [0.5, 0.3, 0.2])
    out = sprint(show, MIME"text/plain"(), model)

    # The Split sits under the outer modifier, and its streams are indented one
    # level deeper again (three leading spaces of continuation prefix).
    @test occursin("└─ model: Split", out)
    @test occursin("   ├─ cases: PoissonError", out)
    @test occursin("   └─ deaths: LatentDelay", out)
    # The downstream stream keeps recursing, deeper still.
    @test occursin("Ascertainment", out)
    @test occursin("FixedIntercept", out)
    @test !occursin("{", out)
end

@testitem "data-driven strata Split shows its template" begin
    using ComposableTuringIDModels, Distributions

    split = Split(PoissonError())
    out = sprint(show, MIME"text/plain"(), split)

    @test startswith(out, "Split")
    @test occursin("└─ template: PoissonError", out)
    @test !occursin("{", out)
end

@testitem "compact 2-arg show prints only the component name" begin
    using ComposableTuringIDModels, Distributions

    model = IDModel(
        DirectInfections(; Z = RandomWalk(), initialisation = Normal()),
        PoissonError())
    # The compact form (used inside arrays / `repr`) is a single clean line with
    # no nested parametric type dump.
    @test sprint(show, model) == "IDModel"
    @test sprint(show, RandomWalk()) == "RandomWalk"
    @test !occursin("{", repr(model))
end
