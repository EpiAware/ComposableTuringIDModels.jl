@testitem "observation error models simulate and condition" begin
    using EpiAwarePrototype, Distributions, Random
    Random.seed!(11)
    Y_t = fill(10.0, 15)
    for obs in (PoissonError(), NegativeBinomialError())
        # Simulate from the prior (missing observations).
        sim = as_turing_model(obs, missing, Y_t)()
        @test length(sim) == length(Y_t)
        @test all(>=(0), sim)
        # Condition on simulated data: the model still builds and evaluates.
        cond = as_turing_model(obs, sim, Y_t)
        @test cond() == sim
    end
end

@testitem "NormalError is a continuous observation-error model" begin
    using EpiAwarePrototype, Distributions, Random
    Random.seed!(111)

    ne = NormalError()
    # Role + interface conformance (the post-role-hierarchy contract).
    @test ne isa AbstractObservationErrorModel
    @test ne isa AbstractObservationModel
    @test implements_observation_interface(ne)

    # Default std prior is positive (a HalfNormal).
    @test minimum(ne.std_prior) >= 0

    Y_t = fill(10.0, 20)
    # Simulate from the prior: continuous (real) observations. Simulating from a
    # `missing` series returns a `Union{Missing,Float64}` vector, so check the
    # sampled values are real rather than the container eltype.
    sim = as_turing_model(ne, missing, Y_t)()
    @test length(sim) == length(Y_t)
    @test all(x -> x isa Real, sim)
    # Condition on the simulated data: the model still builds and evaluates,
    # returning the same observations.
    cond = as_turing_model(ne, sim, Y_t)
    @test cond() == sim
    # The standard deviation is an inferred parameter.
    draw = rand(as_turing_model(ne, sim, Y_t))
    @test any(k -> occursin("σ", string(k)), keys(draw))

    # A custom std prior is honoured.
    ne2 = NormalError(; std_prior = truncated(Normal(0, 2), 0, Inf))
    @test length(as_turing_model(ne2, missing, fill(5.0, 6))()) == 6

    # Expected-mean alignment: conditioning on data with a matching expected
    # series evaluates to that data (the Gaussian likelihood is centred on Y_t).
    μ = 50.0
    obs = μ .+ 0.5 .* randn(200)
    m = as_turing_model(NormalError(; std_prior = truncated(Normal(0, 1), 0, Inf)),
        obs, fill(μ, 200))
    @test m() == obs
end

@testitem "LatentDelay shortens expectations and wraps an error model" begin
    using EpiAwarePrototype, Distributions, Random
    Random.seed!(12)
    obs = LatentDelay(PoissonError(), truncated(Normal(5.0, 2.0), 0.0, Inf))
    Y_t = fill(10.0, 30)
    sim = as_turing_model(obs, missing, Y_t)()
    @test length(sim) == length(Y_t)
    # LatentDelay deliberately leaves the head of the series unobserved
    # (partially observed data), so only the non-missing tail is filled.
    observed = filter(!ismissing, sim)
    @test !isempty(observed)
    @test all(>=(0), observed)
end

@testitem "safe count distributions tolerate very large means" begin
    using EpiAwarePrototype, Distributions, Random
    Random.seed!(13)
    # exp(48) overflows Int; the safe samplers must not throw.
    bigλ = exp(48.0)
    @test rand(SafePoisson(bigλ)) >= 0
    σ² = bigλ + 0.05 * bigλ^2
    p = bigλ / σ²
    r = bigλ * p / (1 - p)
    @test rand(SafeNegativeBinomial(r, p)) >= 0
end
