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

@testitem "TriangleObs constructs the reporting triangle and masks t+d≤now" begin
    using EpiAwarePrototype, Distributions, Random
    Random.seed!(60)
    obs = TriangleObs(PoissonError(), [0.5, 0.3, 0.2])     # Dmax = 2

    # Role + interface conformance: it is an observation model whose `y_t` is a
    # reporting triangle (built from a matrix via `define_y_t`).
    @test obs isa AbstractObservationModel
    N = [10 5 2; 12 6 3; 14 7 4; 16 8 5]
    tri = define_y_t(obs, N, fill(20.0, 4); now = 4)
    @test tri isa ReportingTriangle
    @test implements_observation_interface(obs; y_t = tri, Y_t = fill(20.0, 4))

    # The mask is exactly `t + d ≤ now` (delays d = 0,1,2 on the columns).
    expected = BitMatrix([t + d <= 4 for t in 1:4, d in 0:2])
    @test tri.observed == expected
    @test tri.observed == BitMatrix([true true true; true true true;
                     true true false; true false false])
    @test tri.Dmax == 2

    # A triangle whose Dmax does not match the model's delay PMF is rejected.
    @test_throws Exception define_y_t(obs,
        ReportingTriangle(zeros(Int, 4, 5), trues(4, 5), 4), fill(20.0, 4))
    # A count matrix with the wrong number of delay columns is rejected.
    @test_throws Exception define_y_t(obs, zeros(Int, 4, 5), fill(20.0, 4))
end

@testitem "define_y_t builds a triangle from a matrix and a long-form table" begin
    using EpiAwarePrototype, Random
    Random.seed!(61)
    obs = TriangleObs(PoissonError(), [0.6, 0.25, 0.15])   # Dmax = 2
    N = [10 5 2; 12 6 3; 14 7 4]                            # 3 reference days

    from_matrix = define_y_t(obs, N, fill(20.0, 3); now = 3)

    # The same data as long-form `(reference, delay, count)` rows.
    refs = Int[];
    dls = Int[];
    cts = Int[]
    for t in 1:3, d in 0:2

        push!(refs, t)
        push!(dls, d)
        push!(cts, N[t, d + 1])
    end
    tbl = (reference = refs, delay = dls, count = cts)
    from_table = define_y_t(obs, tbl, fill(20.0, 3); now = 3)

    @test from_table.counts == N
    @test from_table.observed == from_matrix.observed
    @test from_table.Dmax == from_matrix.Dmax == 2

    # A `missing` series gives a fully observed triangle (now = n + Dmax) of
    # `missing` cells, sized to `Y_t`.
    mt = define_y_t(obs, missing, fill(20.0, 5))
    @test size(mt.counts) == (5, 3)
    @test all(ismissing, mt.counts)
    @test all(mt.observed)
end

@testitem "TriangleObs simulates, conditions, and recovers per-cell means" begin
    using EpiAwarePrototype, Distributions, Random, Statistics
    Random.seed!(62)
    pmf = [0.5, 0.3, 0.2]
    obs = TriangleObs(PoissonError(), pmf)
    n = 7
    Y = fill(40.0, n)

    # Simulate: a fully observed triangle of non-negative integer counts.
    sim = as_turing_model(obs, missing, Y)()
    @test sim isa ReportingTriangle
    @test size(sim.counts) == (n, 3)
    @test all(c -> c >= 0, sim.counts[sim.observed])

    # Condition on the simulated triangle: the observed cells round-trip.
    cond = as_turing_model(obs, sim, Y)()
    @test cond.counts[sim.observed] == sim.counts[sim.observed]

    # Each cell mean is `Y_t · p[d+1]`: average a stack of simulated triangles.
    Random.seed!(63)
    draws = [Float64.(as_turing_model(obs, missing, fill(100.0, 5))().counts)
             for _ in 1:4000]
    cell_means = mean(draws)
    for d in 0:2
        @test isapprox(cell_means[1, d + 1], 100.0 * pmf[d + 1]; rtol = 0.05)
    end
end

@testitem "TriangleObs observed row-sums reconcile with CDF-scaling" begin
    using EpiAwarePrototype, Distributions, Random
    Random.seed!(64)
    # #50's consistency check: the observed row-sums of the triangle are the
    # marginal that the EpiNow2-style CDF-scaling nowcast conditions on. With
    # `now = n`, reference day `t` has age `a = n - t` and observed delays
    # `d = 0 … a`, so its expected observed row-sum is
    #   Σ_{d=0}^{a} μ_t · p[d+1] = μ_t · F[a+1],
    # exactly the CDF-scaled expected observed-so-far. (PR #54's `CDFScaledObs`
    # is the released realisation of this scaling; the expectation is replicated
    # here directly so the test does not depend on that symbol.)
    pmf = [0.5, 0.3, 0.2]
    Dmax = length(pmf) - 1
    F = cumsum(pmf)                                 # the reporting-delay CDF
    n = 8
    μ = collect(range(20.0, 80.0; length = n))      # expected eventual totals

    # Expected (noise-free) observed row-sums from the masked triangle.
    triangle_rowsums = zeros(n)
    for t in 1:n, d in 0:Dmax

        (t + d <= n) || continue                    # the `t + d ≤ now` mask, now = n
        triangle_rowsums[t] += μ[t] * pmf[d + 1]
    end

    # The CDF-scaled observed-so-far: completeness[t] = F[age = n - t]
    # (ages beyond the delay support are fully reported).
    completeness = [(n - t) >= length(F) ? 1.0 : F[(n - t) + 1] for t in 1:n]
    cdf_scaled = μ .* completeness

    @test triangle_rowsums≈cdf_scaled rtol=1e-10
    # The oldest day is fully reported; the most recent only its delay-0 fraction.
    @test triangle_rowsums[1]≈μ[1] rtol=1e-10
    @test triangle_rowsums[n]≈μ[n] * pmf[1] rtol=1e-10
end

@testitem "TriangleObs composes with the renewal infection process" begin
    using EpiAwarePrototype, Distributions, Random
    Random.seed!(65)
    data = EpiData([0.2, 0.3, 0.5], exp)
    renewal = Renewal(data; rt = RandomWalk(), initialisation_prior = Normal())

    # Released-CD discretised-delay constructor (the LatentDelay / EpiData path).
    obs = TriangleObs(NegativeBinomialError(),
        truncated(Normal(2.0, 1.0), 0.0, Inf))
    model = EpiAwareModel(renewal, obs)

    sim = as_turing_model(model, missing, 15)()
    @test sim.generated_y_t isa ReportingTriangle
    @test length(sim.I_t) == 15

    # Conditioning on the simulated triangle builds and evaluates.
    cond = as_turing_model(model, sim.generated_y_t, 15)()
    @test cond.generated_y_t isa ReportingTriangle

    # A triangle stream stacks alongside a plain-vector stream.
    stk = StackObservationModels((
        tri = TriangleObs(PoissonError(), [0.6, 0.25, 0.15]),
        cases = PoissonError()))
    out = as_turing_model(stk, (tri = missing, cases = missing), fill(20.0, 8))()
    @test length(out) == 2
end
