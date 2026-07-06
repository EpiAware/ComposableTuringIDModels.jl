# Tests for the numeric `Hierarchy` partial-pooling construct: it returns a
# length-`n_groups` numeric per-group level path (shared mean + group deviations),
# with the cross-group relationship supplied through the prior interface
# (iid-Normal ⇒ classic exchangeable pooling, RandomWalk ⇒ correlated neighbours).

@testitem "Hierarchy constructs and is an AbstractLatentModel" begin
    using EpiAwarePrototype, Distributions
    h = Hierarchy(Normal(), IID(Normal(0.0, 1.0)))
    @test h isa Hierarchy
    @test h isa AbstractLatentModel
    # Both slots are coerced to the prior interface.
    @test h.mean isa AbstractPriorModel
    @test h.across isa IID
    # Keyword constructor; across defaults to an iid Normal.
    hk = Hierarchy()
    @test hk isa Hierarchy
    @test hk.across isa IID
    # A bare Distribution in `mean` is coerced via as_prior.
    @test Hierarchy(; mean = Normal(0, 2)).mean isa BroadcastPrior
end

@testitem "Hierarchy returns a numeric length-n_groups path" begin
    using EpiAwarePrototype, Distributions, Random
    Random.seed!(381)
    h = Hierarchy(Normal(), IID(Normal(0.0, 1.0)))
    # n_groups is NOT stored on the struct; it is read at build time.
    @test !(:n_groups in fieldnames(typeof(h)))
    for G in (2, 4, 7)
        vals = as_turing_model(h, G)()
        @test length(vals) == G
        @test all(x -> x isa Real, vals)      # numeric, not model variants
    end
end

@testitem "iid-Normal across recovers classic (exchangeable) partial pooling" begin
    using EpiAwarePrototype, Distributions, Random, Statistics
    Random.seed!(383)
    # Zero shared mean so the group values are the (exchangeable) deviations.
    h = Hierarchy(Dirac(0.0), IID(Normal(0.0, 1.0)))
    G = 6
    reps = 3000
    M = reduce(hcat, [as_turing_model(h, G)() for _ in 1:reps])   # G × reps
    # Exchangeable draws ⇒ near-zero correlation between neighbouring groups.
    neighbour_cors = [cor(M[g, :], M[g + 1, :]) for g in 1:(G - 1)]
    @test mean(abs.(neighbour_cors)) < 0.1
    @test abs(mean(M)) < 0.1
    @test isapprox(std(vec(M)), 1.0; atol = 0.1)
end

@testitem "RandomWalk across gives correlated neighbouring group effects" begin
    using EpiAwarePrototype, Distributions, Random, Statistics
    Random.seed!(384)
    h_rw = Hierarchy(Dirac(0.0), RandomWalk())
    h_iid = Hierarchy(Dirac(0.0), IID(Normal(0.0, 1.0)))
    G = 6
    reps = 3000
    rw = reduce(hcat, [as_turing_model(h_rw, G)() for _ in 1:reps])
    iid = reduce(hcat, [as_turing_model(h_iid, G)() for _ in 1:reps])
    rw_neighbour = mean([cor(rw[g, :], rw[g + 1, :]) for g in 1:(G - 1)])
    iid_neighbour = mean([cor(iid[g, :], iid[g + 1, :]) for g in 1:(G - 1)])
    # The random walk relates neighbours: strong positive lag-1 correlation,
    # well above the (near-zero) iid case.
    @test rw_neighbour > 0.5
    @test rw_neighbour > iid_neighbour + 0.4
end

@testitem "Hierarchy conforms to the latent interface and composes numerically" begin
    using EpiAwarePrototype, Distributions, Random
    Random.seed!(385)
    h = Hierarchy(Normal(), IID(Normal(0.0, 1.0)))
    @test implements_latent_interface(h; n = 4)
    @test !implements_observation_interface(h)
    @test !implements_infection_interface(h)
    # Because it returns numeric values, it threads straight into an infection
    # model's latent slot (no group-axis contract change).
    inf = DirectInfections(; Z = h, initialisation = Normal())
    out = as_turing_model(inf, 5)()
    @test length(out.I_t) == 5
    @test all(x -> x isa Real, out.Z_t)
end

@testitem "a partially-pooled model samples under NUTS (ForwardDiff)" tags=[:sample] begin
    using EpiAwarePrototype, Distributions, Turing, Random
    using DynamicPPL: to_submodel
    Random.seed!(386)
    # Each group has its own partially-pooled level mapped to Poisson counts.
    # n_groups is read from the data (the length of the observed group vector).
    h = Hierarchy(Normal(), IID(Normal(0.0, 1.0)))

    @model function pooled_counts(h, y, n_groups)
        levels ~ to_submodel(as_turing_model(h, n_groups), false)
        for g in 1:n_groups
            y[g] ~ Poisson(exp(levels[g]))
        end
        return y
    end

    G = 5
    ydata = Int.(pooled_counts(h, fill(missing, G), G)())
    chn = sample(pooled_counts(h, ydata, G),
        NUTS(0.8; adtype = Turing.AutoForwardDiff()), 60; progress = false)
    @test chn !== nothing
    @test size(chn, 1) == 60
end
