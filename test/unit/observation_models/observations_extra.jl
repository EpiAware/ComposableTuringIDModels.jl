@testitem "Ascertainment scales expected observations by a latent model" begin
    using EpiAwarePrototype, Distributions, Random
    Random.seed!(51)
    Y = fill(100.0, 14)
    asc = Ascertainment(model = NegativeBinomialError(), latent_model = FixedIntercept(0.1))
    sim = as_turing_model(asc, missing, Y)()
    @test length(sim) == length(Y)
    @test all(>=(0), sim)

    adw = ascertainment_dayofweek(PoissonError())
    @test length(as_turing_model(adw, missing, Y)()) == length(Y)
end

@testitem "Aggregate sums expected observations over windows" begin
    using EpiAwarePrototype, Random
    Random.seed!(52)
    agg = Aggregate(PoissonError(), [0, 0, 0, 0, 7, 0, 0])
    out = as_turing_model(agg, missing, fill(1.0, 28))()
    @test length(out) == 28
    # Only the present (weekly) positions are non-zero.
    @test count(!=(0), out) == 4
end

@testitem "CDFScaledObs construction and interface" begin
    using EpiAwarePrototype, Distributions

    # Role + interface conformance (the nowcasting modifier is an observation model).
    o = CDFScaledObs(PoissonError(), [0.2, 0.6, 1.0])
    @test o isa AbstractObservationModel
    @test implements_observation_interface(o)

    # The precomputed-CDF constructor stores the CDF as-is.
    @test o.delay_cdf == [0.2, 0.6, 1.0]

    # An invalid CDF (decreasing / out of [0, 1]) is rejected at construction.
    @test_throws Exception CDFScaledObs(PoissonError(), [0.6, 0.2, 1.0])
    @test_throws Exception CDFScaledObs(PoissonError(), [0.2, 0.6, 1.5])
    @test_throws Exception CDFScaledObs(PoissonError(), [-0.1, 0.6, 1.0])

    # The distribution constructor builds the CDF from the released-CD
    # double-interval-censored PMF (the LatentDelay / EpiData path): cumulative,
    # non-decreasing, ending at 1.
    od = CDFScaledObs(NegativeBinomialError(), truncated(Normal(5.0, 2.0), 0.0, Inf))
    @test issorted(od.delay_cdf)
    @test isapprox(od.delay_cdf[end], 1.0)
    @test all(0 .<= od.delay_cdf .<= 1 + 1e-8)
end

@testitem "CDFScaledObs scales expected eventual totals by the delay CDF" begin
    using EpiAwarePrototype, Distributions, Random
    Random.seed!(60)

    # Huge expected means make the Poisson draw ≈ its mean, so the realised
    # observed-so-far reveals the applied scaling. With a 3-bin CDF and a length-6
    # series, reference days older than the delay support are fully reported.
    n = 6
    Y = fill(1.0e6, n)
    F = [0.2, 0.5, 1.0]                       # completeness by age 0, 1, 2
    obs = CDFScaledObs(PoissonError(), F)
    sim = as_turing_model(obs, missing, Y)()
    # Completeness by reference day t = reverse of [F; ones(n - length(F))]:
    # the most recent day (t = n, age 0) is scaled by F[1], the oldest by 1.
    expected_scale = reverse(vcat(F, ones(n - length(F))))
    @test all(abs.(sim ./ Y .- expected_scale) .< 0.01)

    # The oldest day (t = 1) is fully reported; the most recent (t = n) least so.
    @test isapprox(sim[1] / Y[1], 1.0; atol = 0.01)
    @test isapprox(sim[n] / Y[n], F[1]; atol = 0.01)
end

@testitem "CDFScaledObs with a complete delay reduces to the inner model" begin
    using EpiAwarePrototype, Random
    Random.seed!(61)
    n = 8
    Y = fill(50.0, n)

    # A fully-reported delay (CDF all ones, or a single bin F = [1.0] padded to
    # ones) means every reference day is complete, so the modifier is a no-op:
    # it reduces exactly to the wrapped error model.
    for F in (ones(n), [1.0])
        Random.seed!(123)
        scaled = as_turing_model(CDFScaledObs(PoissonError(), F), missing, Y)()
        Random.seed!(123)
        inner = as_turing_model(PoissonError(), missing, Y)()
        @test scaled == inner
    end
end

@testitem "CDFScaledObs simulate-then-condition and length handling" begin
    using EpiAwarePrototype, Random
    Random.seed!(62)
    n = 10
    Y = fill(100.0, n)
    F = collect(range(0.1, 1.0; length = n))
    obs = CDFScaledObs(PoissonError(), F)

    sim = as_turing_model(obs, missing, Y)()
    @test length(sim) == n
    @test all(>=(0), sim)
    # Conditioning on the simulated data returns it.
    @test as_turing_model(obs, sim, Y)() == sim

    # A CDF shorter than the series is fine: older days are taken complete.
    short = CDFScaledObs(PoissonError(), [0.3, 0.7, 1.0])
    @test length(as_turing_model(short, missing, fill(10.0, 20))()) == 20
    # A CDF longer than the series is also fine (only its head is used).
    long = CDFScaledObs(PoissonError(), collect(range(0.05, 1.0; length = 30)))
    @test length(as_turing_model(long, missing, fill(10.0, 5))()) == 5
end

@testitem "CDFScaledObs composes with a renewal model end-to-end" begin
    using EpiAwarePrototype, Distributions, Random
    Random.seed!(63)
    data = EpiData([0.2, 0.3, 0.5], exp)
    n = 25
    model = EpiAwareModel(
        Renewal(data; rt = RandomWalk(), initialisation_prior = Normal()),
        CDFScaledObs(NegativeBinomialError(),
            truncated(Normal(4.0, 1.5), 0.0, Inf)))

    sim = as_turing_model(model, missing, n)()
    y = sim.generated_y_t
    @test length(y) == n
    # The composed posterior builds and evaluates when conditioned on the data,
    # reproducing the conditioned observations as a generated quantity.
    @test as_turing_model(model, y, n)().generated_y_t == y

    # It also stacks under StackObservationModels.
    stk = StackObservationModels((reported = CDFScaledObs(PoissonError(),
        truncated(Normal(4.0, 1.5), 0.0, Inf)),))
    out = as_turing_model(stk, (reported = missing,), fill(100.0, n))()
    @test length(out) == 1
end

@testitem "PrefixObservationModel prefixes observation parameters" begin
    using EpiAwarePrototype, Random
    Random.seed!(53)
    pom = PrefixObservationModel(model = NegativeBinomialError(), prefix = "Test")
    names = string.(collect(keys(rand(as_turing_model(pom, missing, fill(10.0, 5))))))
    @test any(startswith("Test."), names)
end

@testitem "RecordExpectedObs and TransformObservationModel wrap an error model" begin
    using EpiAwarePrototype, Random
    Random.seed!(54)
    Y = fill(10.0, 30)
    reo = RecordExpectedObs(NegativeBinomialError())
    @test length(as_turing_model(reo, missing, Y)()) == length(Y)

    tom = TransformObservationModel(NegativeBinomialError())
    @test length(as_turing_model(tom, missing, Y)()) == length(Y)
end

@testitem "StackObservationModels prefixes and stacks several models" begin
    using EpiAwarePrototype, Distributions, Random
    Random.seed!(55)
    stk = StackObservationModels((cases = PoissonError(),
        deaths = NegativeBinomialError()))
    yt = (cases = missing, deaths = missing)
    sm = as_turing_model(stk, yt, fill(10.0, 10))
    names = string.(collect(keys(rand(sm))))
    @test any(startswith("cases."), names)
    @test any(startswith("deaths."), names)
    out = sm()
    @test length(out) == 2
    @test length(out[1]) == 10
end
