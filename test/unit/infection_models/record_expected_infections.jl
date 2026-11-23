# Recording an infection model's noise-free expectation.
#
# `RecordExpectedInfections` wraps any infection model and tracks `exp_I_t`,
# the incidence entering the model's noise layer. A model with no noise layer
# commits its own expectation, so wrapping one records the committed series.

@testitem "RecordExpectedInfections records a model that draws nothing" begin
    using ComposableTuringIDModels, Distributions, Random
    gen_int = [0.2, 0.3, 0.5]
    models = (
        DirectInfections(; Z = RandomWalk(), initialisation = Normal()),
        Renewal(;
            generation_time = gen_int, rt = RandomWalk(),
            initialisation = Normal()
        ),
        Renewal(
            gen_int, SusceptibleDepletion(500.0);
            rt = RandomWalk(), initialisation = Normal()
        ),
    )
    for model in models
        out = as_turing_model(RecordExpectedInfections(model), 20)(Xoshiro(11))
        @test out.exp_I_t == out.I_t
        @test length(out.exp_I_t) == 20
    end
end

@testitem "RecordExpectedInfections records the pre-noise incidence" begin
    using ComposableTuringIDModels, Distributions, Random
    gen_int = [0.2, 0.3, 0.5]
    noisy = (
        Renewal(
            gen_int, SusceptibleDepletion(5.0e3), InfectionNoise();
            rt = FixedIntercept(0.0),
            initialisation = Normal(log(100.0), 0.0)
        ),
        StochasticRenewal(
            gen_int; rt = FixedIntercept(0.0),
            initialisation = Normal(log(100.0), 0.0)
        ),
    )
    for model in noisy
        out = as_turing_model(RecordExpectedInfections(model), 20)(Xoshiro(12))
        @test out.exp_I_t != out.I_t
        @test all(isfinite, out.exp_I_t)
        @test all(>(0), out.exp_I_t)
    end
end

@testitem "wrapping renames no sampled variable" begin
    using ComposableTuringIDModels, Distributions
    using DynamicPPL: VarInfo
    gen_int = [0.2, 0.3, 0.5]
    models = (
        DirectInfections(; Z = RandomWalk(), initialisation = Normal()),
        Renewal(
            gen_int, SusceptibleDepletion(500.0), InfectionNoise();
            rt = RandomWalk(), initialisation = Normal()
        ),
        StochasticRenewal(
            gen_int; rt = RandomWalk(), initialisation = Normal()
        ),
    )
    for model in models
        bare = keys(VarInfo(as_turing_model(model, 15)))
        wrapped = keys(VarInfo(as_turing_model(RecordExpectedInfections(model), 15)))
        @test bare == wrapped
    end
end

@testitem "a stratified renewal reports one expectation per stratum" begin
    using ComposableTuringIDModels, Distributions, Random
    model = Renewal(
        [0.2, 0.3, 0.5], SusceptibleDepletion(500.0);
        rt = Stratify(RandomWalk(), FixedIntercept(0.0)),
        initialisation = Normal(log(50.0), 0.0)
    )
    out = as_turing_model(RecordExpectedInfections(model), (2, 20))(Xoshiro(13))
    @test size(out.exp_I_t) == (2, 20)
    @test out.exp_I_t == out.I_t
end

@testitem "exp_I_t comes back from a chain" tags = [:sample] begin
    using ComposableTuringIDModels, Distributions, Turing, Random
    Random.seed!(482)
    model = IDModel(
        RecordExpectedInfections(
            Renewal(
                [0.2, 0.3, 0.5], SusceptibleDepletion(1000.0);
                rt = RandomWalk(), initialisation = Normal()
            )
        ),
        PoissonError()
    )
    y = as_turing_model(model, missing, 20)().generated_y_t
    chn = sample(as_turing_model(model, y, 20), NUTS(), 15; progress = false)
    # `exp_I_t` is tracked with `:=`, so it comes out of the chain, one series
    # per draw.
    draws = vec(chn[:exp_I_t])
    @test length(draws) == 15
    @test all(d -> length(d) == 20, draws)
    @test all(isfinite, reduce(vcat, draws))
end
