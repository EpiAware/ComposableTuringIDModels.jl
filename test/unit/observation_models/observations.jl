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

    # A trials vector too short to cover every scored step is rejected.
    @test_throws Exception as_turing_model(
        be, (y = missing, N = [5, 5, 5]),
        fill(0.2, 10)
    )()

    # The success probability is clamped away from 0/1 (no degenerate likelihood).
    edge = as_turing_model(be, (y = missing, N = 8), [0.0, 1.0, 0.5, 0.5])().y_t
    @test all(x -> 0 <= x <= 8, edge)
end

@testitem "BinomialError right-aligns its trials with the data" begin
    using ComposableTuringIDModels, Distributions, Random
    using ComposableTuringIDModels: MissingObservations
    using DynamicPPL: VarInfo, logjoint
    Random.seed!(306)

    be = BinomialError()
    # An expected series longer than the data: what a `Split` branch with a
    # shorter lead-in than its neighbours produces. The trials are data, so
    # they are right-aligned exactly like the successes and may be given at
    # the length of either series.
    p = collect(range(0.2, 0.6; length = 10))
    y = [2, 3, 4, 5]
    N_expected = collect(11:20)             # padded to the expected series
    N_data = N_expected[(end - 3):end]      # only the days the caller observed

    ref = sum(logpdf.(Binomial.(N_data, p[(end - 3):end]), y))
    padded = as_turing_model(be, (y = y, N = N_expected), p)
    unpadded = as_turing_model(be, (y = y, N = N_data), p)
    @test logjoint(padded, VarInfo(padded)) ≈ ref
    @test logjoint(unpadded, VarInfo(unpadded)) ≈ ref

    # A scalar is broadcast across the scored steps either way.
    scalar = as_turing_model(be, (y = y, N = 20), p)
    @test logjoint(scalar, VarInfo(scalar)) ≈
        sum(logpdf.(Binomial.(20, p[(end - 3):end]), y))

    # The other direction — the data runs longer than the expected series (the
    # chain's lead-in) — takes the same alignment.
    y_long = collect(2:13)
    N_long = collect(21:32)
    ref_long = sum(logpdf.(Binomial.(N_long[3:12], p), y_long[3:12]))
    lead_in = as_turing_model(be, (y = y_long, N = N_long), p)
    trimmed = as_turing_model(be, (y = y_long, N = N_long[3:12]), p)
    @test logjoint(lead_in, VarInfo(lead_in)) ≈ ref_long
    @test logjoint(trimmed, VarInfo(trimmed)) ≈ ref_long

    # A partially observed series takes the same alignment: the gap is drawn
    # against the trials of the day it falls on.
    carrier = MissingObservations([2, 0, 4, 5], Bool[1, 0, 1, 1])
    out = as_turing_model(be, (y = carrier, N = N_data), p)()
    @test length(out.y_t) == 4
    @test out.y_t[[1, 3, 4]] == [2, 4, 5]
    @test 0 <= out.y_t[2] <= N_data[2]

    # Trials that do not cover every scored step are rejected.
    @test_throws Exception as_turing_model(be, (y = y, N = [5, 5, 5]), p)()
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

@testitem "an expected series longer than the data is right-aligned" begin
    using ComposableTuringIDModels, Distributions
    using ComposableTuringIDModels: MissingObservations
    using DynamicPPL: VarInfo, logjoint

    # The model runs longer than the data. The extra leading expected values
    # are unobserved run-in and every observation is still scored against the
    # day it belongs to.
    Y_t = collect(10.0:1.0:19.0)
    y = [12, 13, 14, 15]

    mdl = as_turing_model(PoissonError(), y, Y_t)
    ref = sum(logpdf.(SafePoisson.(Y_t[(end - 3):end] .+ 1.0e-6), y))
    @test logjoint(mdl, VarInfo(mdl)) ≈ ref

    # A `MissingObservations` carrier takes the same alignment: the observed
    # entries come back as given and the gap is drawn.
    carrier = MissingObservations([12, 0, 14, 15], Bool[1, 0, 1, 1])
    out = as_turing_model(PoissonError(), carrier, Y_t)()
    @test length(out.y_t) == 4
    @test out.y_t[[1, 3, 4]] == [12, 14, 15]
end
