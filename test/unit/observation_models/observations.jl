@testitem "observation error models simulate and condition" begin
    using ComposableTuringIDModels, Distributions, Random
    Random.seed!(11)
    Y_t = fill(10.0, 15)
    for obs in (PoissonError(), NegativeBinomialError())
        # Simulate from the prior (missing observations). Every observation model
        # returns the uniform `(; y_t, expected)` contract.
        res = as_turing_model(obs, missing, Y_t)()
        @test keys(res) == (:y_t, :expected)
        @test res.expected == Y_t
        sim = res.y_t
        @test length(sim) == length(Y_t)
        @test all(>=(0), sim)
        # Condition on simulated data: the model still builds and evaluates.
        cond = as_turing_model(obs, sim, Y_t)
        @test cond().y_t == sim
    end
end

@testitem "NormalError is a continuous observation-error model" begin
    using ComposableTuringIDModels, Distributions, Random
    Random.seed!(111)

    ne = NormalError()
    # Role + interface conformance (the post-role-hierarchy contract).
    @test ne isa AbstractObservationErrorModel
    @test ne isa AbstractObservationModel
    @test implements_observation_interface(ne)

    # Default std prior is positive (a HalfNormal).
    @test minimum(ne.std) >= 0

    Y_t = fill(10.0, 20)
    # Simulate from the prior: continuous (real) observations. Simulating from a
    # `missing` series returns a `Union{Missing,Float64}` vector, so check the
    # sampled values are real rather than the container eltype.
    sim = as_turing_model(ne, missing, Y_t)().y_t
    @test length(sim) == length(Y_t)
    @test all(x -> x isa Real, sim)
    # Condition on the simulated data: the model still builds and evaluates,
    # returning the same observations.
    cond = as_turing_model(ne, sim, Y_t)
    @test cond().y_t == sim
    # The standard deviation is an inferred parameter.
    draw = rand(as_turing_model(ne, sim, Y_t))
    @test any(k -> occursin("σ", string(k)), keys(draw))

    # A custom std prior is honoured.
    ne2 = NormalError(; std = truncated(Normal(0, 2), 0, Inf))
    @test length(as_turing_model(ne2, missing, fill(5.0, 6))().y_t) == 6

    # Expected-mean alignment: conditioning on data with a matching expected
    # series evaluates to that data (the Gaussian likelihood is centred on Y_t).
    μ = 50.0
    obs = μ .+ 0.5 .* randn(200)
    m = as_turing_model(
        NormalError(; std = truncated(Normal(0, 1), 0, Inf)),
        obs, fill(μ, 200)
    )
    @test m().y_t == obs
end

@testitem "BinomialError reads trials from NamedTuple data" begin
    using ComposableTuringIDModels, Distributions, Random
    Random.seed!(222)

    be = BinomialError()
    # The struct carries no data.
    @test isempty(fieldnames(BinomialError))
    # Role + interface conformance. The expected series is a probability, and the
    # number of trials is supplied via the NamedTuple data `y_t`.
    @test be isa AbstractObservationErrorModel
    @test be isa AbstractObservationModel
    @test implements_observation_interface(
        be; y_t = (y = missing, N = 20),
        Y_t = fill(0.3, 10)
    )

    # Scalar N (in the data) is broadcast across the series; successes lie in 0..N.
    p = fill(0.3, 12)
    sim = as_turing_model(be, (y = missing, N = 20), p)().y_t
    @test length(sim) == length(p)
    @test all(x -> 0 <= x <= 20, sim)
    # Conditioning on the data returns it.
    @test as_turing_model(be, (y = sim, N = 20), p)().y_t == sim

    # A per-time-point trials vector in the data is honoured.
    Nvec = collect(5:14)            # length 10
    simv = as_turing_model(be, (y = missing, N = Nvec), fill(0.8, 10))().y_t
    @test all(i -> 0 <= simv[i] <= Nvec[i], eachindex(simv))
    @test as_turing_model(be, (y = simv, N = Nvec), fill(0.8, 10))().y_t == simv

    # `y_t` must be a NamedTuple carrying `N`: a plain vector or a NamedTuple
    # without `N` is rejected.
    @test_throws Exception as_turing_model(be, fill(3, 10), fill(0.3, 10))()
    @test_throws Exception as_turing_model(be, (y = missing,), fill(0.3, 10))()

    # A trials vector whose length does not match the series is rejected.
    @test_throws Exception as_turing_model(
        be, (y = missing, N = [5, 5, 5]),
        fill(0.2, 10)
    )()

    # The success probability is clamped away from 0/1 (no degenerate likelihood).
    edge = as_turing_model(be, (y = missing, N = 8), [0.0, 1.0, 0.5, 0.5])().y_t
    @test all(x -> 0 <= x <= 8, edge)
end

@testitem "NamedTuple data keeps blanks latent and leaves the caller's data alone" begin
    using ComposableTuringIDModels, Distributions, Turing, Random
    using DynamicPPL: @varname
    Random.seed!(64)

    # A `missing` entry reached through a NamedTuple field must behave exactly
    # as one reached through a plain vector: sampled as a latent, present in
    # the chain, redrawn on every evaluation, and never written back into the
    # array the caller handed over.
    trials = fill(1000, 3)
    y = Vector{Union{Missing, Int}}([10, 25, missing])
    y_before = copy(y)
    mdl = as_turing_model(BinomialError(), (y = y, N = trials), fill(0.02, 3))
    chn = sample(mdl, Prior(), 20; progress = false)

    @test isequal(y, y_before)
    draws = vec(chn[@varname(y_t[3])])
    @test length(draws) == 20
    @test length(unique(draws)) > 1
    # The observed entries stay data: they are not sampled.
    @test !haskey(chn, @varname(y_t[1]))

    # The same holds for the shared error-model loop (Poisson, NegBin, Normal),
    # which also accepts its counts in a NamedTuple `y` field.
    yp = Vector{Union{Missing, Int}}([5, missing, 7])
    yp_before = copy(yp)
    chp = sample(
        as_turing_model(PoissonError(), (y = yp,), fill(6.0, 3)), Prior(), 20;
        progress = false
    )
    @test isequal(yp, yp_before)
    @test length(unique(vec(chp[@varname(y_t[2])]))) > 1

    # A fully missing series in a NamedTuple field is left alone too.
    ym = Vector{Union{Missing, Int}}(missing, 3)
    chm = sample(
        as_turing_model(BinomialError(), (y = ym, N = trials), fill(0.02, 3)),
        Prior(), 5; progress = false
    )
    @test all(ismissing, ym)
    @test length(unique(vec(chm[@varname(y_t[1])]))) > 1

    # Wrapping the error model in a modifier changes none of this.
    yw = Vector{Union{Missing, Int}}([10, 25, missing])
    yw_before = copy(yw)
    wrapped = TransformObservationModel(BinomialError(), x -> x ./ 100)
    chw = sample(
        as_turing_model(wrapped, (y = yw, N = trials), fill(2.0, 3)), Prior(),
        20; progress = false
    )
    @test isequal(yw, yw_before)
    @test length(unique(vec(chw[@varname(y_t[3])]))) > 1

    # Reached through a composed `IDModel` the data takes the other route into
    # the same helper: `concrete_observations` narrows the NamedTuple
    # field-wise when the model is built, so the series arrives already a
    # carrier. The guarantee has to hold on that path too.
    yc = Vector{Union{Missing, Int}}([10, 25, missing])
    yc_before = copy(yc)
    composed = IDModel(
        DirectInfections(; Z = RandomWalk(), initialisation = Normal()),
        BinomialError()
    )
    chc = sample(
        as_turing_model(composed, (y = yc, N = trials), 3), Prior(), 20;
        progress = false
    )
    @test isequal(yc, yc_before)
    @test length(unique(vec(chc[@varname(y_t[3])]))) > 1
    @test !haskey(chc, @varname(y_t[1]))
end

@testitem "define_y_t unpacks counts for vector or NamedTuple data" begin
    using ComposableTuringIDModels
    Y_t = fill(10.0, 5)
    # Plain vector passes through.
    @test define_y_t(PoissonError(), [1, 2, 3, 4, 5], Y_t) == [1, 2, 3, 4, 5]
    # NamedTuple: the `y` field is unpacked.
    @test define_y_t(PoissonError(), (y = [1, 2, 3, 4, 5],), Y_t) == [1, 2, 3, 4, 5]
    # `missing` (plain or in the `y` field) becomes a length-Y_t missing vector.
    @test all(ismissing, define_y_t(PoissonError(), missing, Y_t))
    @test length(define_y_t(PoissonError(), missing, Y_t)) == 5
    @test all(ismissing, define_y_t(PoissonError(), (y = missing,), Y_t))
    # BinomialError shares the default unpacking for its `y` field.
    @test define_y_t(BinomialError(), (y = [3, 4], N = 10), fill(0.5, 2)) == [3, 4]

    # A `MissingObservations` carrier (what `concrete_observations` returns for
    # a partially-missing vector) is rebuilt into the ragged vector, at the top
    # level and unpacked from a `y` field alike.
    using ComposableTuringIDModels: MissingObservations
    carrier = MissingObservations([1, 0, 3, 4, 5], Bool[1, 0, 1, 1, 1])
    rebuilt = define_y_t(PoissonError(), carrier, Y_t)
    @test isequal(rebuilt, [1, missing, 3, 4, 5])
    @test isequal(
        define_y_t(PoissonError(), (y = carrier,), Y_t),
        [1, missing, 3, 4, 5]
    )
end

@testitem "LatentDelay shortens expectations and wraps an error model" begin
    using ComposableTuringIDModels, Distributions, Random
    Random.seed!(12)
    obs = LatentDelay(PoissonError(), truncated(Normal(5.0, 2.0), 0.0, Inf))
    Y_t = fill(10.0, 30)
    sim = as_turing_model(obs, missing, Y_t)().y_t
    @test length(sim) == length(Y_t)
    # LatentDelay deliberately leaves the head of the series unobserved
    # (partially observed data), so only the non-missing tail is filled.
    observed = filter(!ismissing, sim)
    @test !isempty(observed)
    @test all(>=(0), observed)
end

@testitem "safe count distributions tolerate very large means" begin
    using ComposableTuringIDModels, Distributions, Random
    Random.seed!(13)
    # exp(48) overflows Int; the safe samplers must not throw.
    bigλ = exp(48.0)
    @test rand(SafePoisson(bigλ)) >= 0
    σ² = bigλ + 0.05 * bigλ^2
    p = bigλ / σ²
    r = bigλ * p / (1 - p)
    @test rand(SafeNegativeBinomial(r, p)) >= 0
end
