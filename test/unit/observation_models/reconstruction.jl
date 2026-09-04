@testitem "an identity rebuild returns the same component" begin
    using ComposableTuringIDModels, Distributions
    using Accessors: Accessors

    # A wrapper whose constructor derives a field is rebuilt by re-deriving it
    # from the stored fields, so a no-op has to return the same model. Nothing
    # throws when it does not, so the check is on the concrete type: a second
    # wrapping shows up there.
    err = PoissonError()
    wrappers = [
        LatentDelay(err, [0.5, 0.3, 0.2]),
        LatentDelay(err, Gamma(2.0, 1.0); D = 5.0),
        Ascertainment(err, FixedIntercept(0.1)),
        Ascertainment(err, Normal(0.0, 0.1)),
        Ascertainment(err, FixedIntercept(0.1); latent_prefix = ""),
        Aggregate(err, [0, 0, 3]),
        RecordExpectedObs(err),
        TransformObservationModel(err),
        PrefixObservationModel(err, "a"),
        RightTruncate(err, [0.5, 1.0]),
        ReportTriangle(err, [0.6, 0.4]),
        Split((cases = err, deaths = NegativeBinomialError())),
        Split(err, [1.0 1.0]),
    ]
    for w in wrappers
        rebuilt = Accessors.modify(identity, w, Accessors.Properties())
        @test typeof(rebuilt) === typeof(w)
        for f in fieldnames(typeof(w))
            @test typeof(getfield(rebuilt, f)) === typeof(getfield(w, f))
        end
    end
end

@testitem "swapping a component leaves the rest of a wrapper alone" begin
    using ComposableTuringIDModels, Distributions
    using Accessors: Accessors, @set

    # `Ascertainment` stores its prior already wrapped in a `PrefixLatentModel`,
    # so a rebuild that re-applied the constructor would wrap it a second time.
    asc = Ascertainment(NormalError(), FixedIntercept(0.1))
    swapped = @set asc.model = PoissonError()
    @test swapped.model isa PoissonError
    @test typeof(getfield(swapped, :latent_model)) ===
        typeof(getfield(asc, :latent_model))
    @test swapped.latent_prefix == asc.latent_prefix

    # `LatentDelay` holds its PMF as it was given, so a rebuild leaves it alone.
    ld = LatentDelay(PoissonError(), [0.5, 0.3, 0.2])
    ld2 = @set ld.model = NegativeBinomialError()
    @test ld2.model isa NegativeBinomialError
    @test ld2.delay == ld.delay

    # `Aggregate` derives its presence mask from its window lengths, so a
    # rebuild must keep the two consistent.
    ag = Aggregate(PoissonError(), [0, 0, 3])
    ag2 = @set ag.model = NormalError()
    @test ag2.aggregation == ag.aggregation
    @test ag2.present == ag.present
end

@testitem "a rebuilt chain generates the same expected series" begin
    using ComposableTuringIDModels, Distributions, Random
    using Accessors: Accessors

    # The consequence a type check stands in for: a corrupted rebuild draws a
    # different prior for the doubly-prefixed effect, so the expected series it
    # generates differs from the original's.
    obs = LatentDelay(
        Ascertainment(PoissonError(), FixedIntercept(log(0.1))), [0.5, 0.3, 0.2]
    )
    model = IDModel(
        DirectInfections(; Z = RandomWalk(), initialisation = Normal(1.0, 0.1)),
        obs
    )
    rebuilt = IDModel(
        model.infection_model,
        Accessors.modify(identity, obs, Accessors.Properties())
    )
    Random.seed!(42)
    a = as_turing_model(model, missing, 20)()
    Random.seed!(42)
    b = as_turing_model(rebuilt, missing, 20)()
    @test a.expected_y_t == b.expected_y_t
end

@testitem "a chain derived with @set samples the same variables" begin
    using ComposableTuringIDModels, Distributions
    using Accessors: @set
    using DynamicPPL: VarInfo

    # The ascertainment's prior is namespaced by its prefix wrapper. Losing the
    # wrapper renames the sampled variable and throws nothing, so the variable
    # names are what is pinned here.
    chain(prior, delay) = LatentDelay(
        Ascertainment(
            NegativeBinomialError(; cluster_factor = HalfNormal(0.1)), prior
        ), delay
    )
    infection = DirectInfections(;
        Z = RandomWalk(), initialisation = Normal(1.0, 0.1)
    )
    varnames(obs) = keys(VarInfo(as_turing_model(IDModel(infection, obs), missing, 30)))

    cases = chain(FixedIntercept(log(0.6)), LogNormal(1.6, 0.5))
    built = chain(Intercept(Normal(log(0.015), 0.25)), LogNormal(2.8, 0.4))

    derived = @set cases.model.latent_model = Intercept(Normal(log(0.015), 0.25))
    derived = @set derived.delay = LogNormal(2.8, 0.4)

    @test varnames(derived) == varnames(built)
    @test Symbol("Ascertainment.intercept") in Symbol.(varnames(built))
    @test derived.model.latent_model isa PrefixLatentModel
    @test derived.delay == built.delay
end

@testitem "a prior slot widens whatever the prefix leaves behind" begin
    using ComposableTuringIDModels, Distributions

    # A `PrefixLatentModel` holds a `PriorLike`, so re-prefixing one to an empty
    # prefix unwraps it and can uncover a bare `Distribution`. The slot still has
    # to widen that into a constant path, which its declared field type would
    # otherwise reject at construction.
    wrapped = PrefixLatentModel(Normal(0.0, 0.1), "own")
    asc = Ascertainment(PoissonError(), wrapped; latent_prefix = "")
    @test asc.latent_model isa Intercept

    combined = CombineLatentModels([wrapped, Normal()], ["", "b"])
    @test combined.models[1] isa Intercept
end
