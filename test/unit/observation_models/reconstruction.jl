@testitem "an identity rebuild returns the same component" begin
    using ComposableTuringIDModels, Distributions
    using Accessors: Accessors

    # A wrapper whose constructor transforms an argument before storing it is
    # rebuilt by re-applying the transform to the already-transformed value, so
    # a no-op returns a different — and silently wrong — model. Nothing throws,
    # so the check is on the concrete type: a second wrapping shows up there.
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

    # `LatentDelay` stores its PMF reversed for the convolution, so a rebuild
    # that re-applied the constructor would reverse it a second time — and a
    # reversed PMF is still a valid one, so nothing would complain.
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
