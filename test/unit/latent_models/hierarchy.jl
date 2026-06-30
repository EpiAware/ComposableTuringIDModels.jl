# Tests for the `Hierarchy` partial-pooling construct: it varies one field of a
# base model across a data-derived number of groups, with the cross-group
# relationship given by a latent model (iid-Normal ⇒ classic partial pooling,
# RandomWalk ⇒ correlated neighbouring groups). The field is written with an
# Accessors.jl lens via `set`.

@testitem "Hierarchy constructs and is an AbstractLatentModel" begin
    using EpiAwarePrototype, Distributions, Accessors
    base = FixedIntercept(0.0)
    lens = @optic _.intercept
    h = Hierarchy(base, lens, IID(Normal(0.0, 1.0)))
    @test h isa Hierarchy
    @test h isa AbstractLatentModel
    @test h.base === base
    @test h.latent isa IID
    # Keyword constructor; latent defaults to an iid Normal.
    hk = Hierarchy(; base = base, lens = lens)
    @test hk isa Hierarchy
    @test hk.latent isa IID
end

@testitem "Hierarchy rejects a lens that misses the base at construction" begin
    using EpiAwarePrototype, Distributions, Accessors
    # An optic onto a non-existent field fails fast at construction, not at
    # sample time.
    @test_throws Exception Hierarchy(
        FixedIntercept(0.0), (@optic _.not_a_field), IID(Normal()))
end

@testitem "Hierarchy derives n_groups from the build call, not a field" begin
    using EpiAwarePrototype, Distributions, Accessors, Random
    Random.seed!(381)
    h = Hierarchy(FixedIntercept(0.0), (@optic _.intercept), IID(Normal(0.0, 1.0)))
    # n_groups is NOT stored on the struct.
    @test !(:n_groups in fieldnames(typeof(h)))
    # The same struct produces however many group variants the build call asks
    # for — the number of groups is read at build time.
    for G in (2, 4, 7)
        variants = as_turing_model(h, G)()
        @test length(variants) == G
        @test all(v -> v isa FixedIntercept, variants)
    end
end

@testitem "Hierarchy writes each group value onto its base variant via the lens" begin
    using EpiAwarePrototype, Distributions, Accessors, Random
    Random.seed!(382)
    base = FixedIntercept(0.0)
    lens = @optic _.intercept
    h = Hierarchy(base, lens, IID(Normal(0.0, 1.0)))
    variants = as_turing_model(h, 5)()
    # Each variant is a copy of the base with only the lens field overwritten;
    # the base itself is untouched (Accessors.set returns a new struct).
    @test base.intercept == 0.0
    @test all(v -> v isa FixedIntercept, variants)
    # The drawn values differ between groups (iid draws), and each variant
    # broadcasts its own value as a latent path.
    vals = [v.intercept for v in variants]
    @test length(unique(vals)) == 5
    path = as_turing_model(variants[1], 4)()
    @test path == fill(variants[1].intercept, 4)
end

@testitem "set_per_group builds one variant per value, including nested lenses" begin
    using EpiAwarePrototype, Distributions, Accessors
    # Flat field.
    base = FixedIntercept(0.0)
    vs = EpiAwarePrototype.set_per_group(base, (@optic _.intercept), [0.1, 0.2, 0.3])
    @test length(vs) == 3
    @test [v.intercept for v in vs] == [0.1, 0.2, 0.3]
    # Deeply-nested field reached by a composed lens (the key Accessors upside):
    # the cluster-factor prior buried inside a LatentDelay-wrapped error model.
    delayed = LatentDelay(NegativeBinomialError(), [0.5, 0.5])
    nlens = @optic _.model.cluster_factor_prior
    nvs = EpiAwarePrototype.set_per_group(delayed, nlens, [Dirac(0.1), Dirac(0.2)])
    @test length(nvs) == 2
    @test all(v -> v isa LatentDelay, nvs)
    @test [only(Accessors.getall(v, nlens)).value for v in nvs] == [0.1, 0.2]
end

@testitem "iid-Normal latent recovers classic (exchangeable) partial pooling" begin
    using EpiAwarePrototype, Distributions, Accessors, Random, Statistics
    Random.seed!(383)
    base = FixedIntercept(0.0)
    lens = @optic _.intercept
    h = Hierarchy(base, lens, IID(Normal(0.0, 1.0)))
    G = 6
    reps = 3000
    # Matrix of group values: G rows × reps columns.
    M = reduce(hcat, [[v.intercept for v in as_turing_model(h, G)()] for _ in 1:reps])
    # Exchangeable draws ⇒ near-zero correlation between neighbouring groups.
    neighbour_cors = [cor(M[g, :], M[g + 1, :]) for g in 1:(G - 1)]
    @test mean(abs.(neighbour_cors)) < 0.1
    # Each group is marginally standard-normal-ish (shared hyperparameters).
    @test abs(mean(M)) < 0.1
    @test isapprox(std(vec(M)), 1.0; atol = 0.1)
end

@testitem "RandomWalk latent gives correlated neighbouring group effects" begin
    using EpiAwarePrototype, Distributions, Accessors, Random, Statistics
    Random.seed!(384)
    base = FixedIntercept(0.0)
    lens = @optic _.intercept
    h_rw = Hierarchy(base, lens, RandomWalk())
    h_iid = Hierarchy(base, lens, IID(Normal(0.0, 1.0)))
    G = 6
    reps = 3000
    rw = reduce(hcat, [[v.intercept for v in as_turing_model(h_rw, G)()] for _ in 1:reps])
    iid = reduce(hcat, [[v.intercept for v in as_turing_model(h_iid, G)()] for _ in 1:reps])
    rw_neighbour = mean([cor(rw[g, :], rw[g + 1, :]) for g in 1:(G - 1)])
    iid_neighbour = mean([cor(iid[g, :], iid[g + 1, :]) for g in 1:(G - 1)])
    # The random walk relates neighbours: strong positive lag-1 correlation,
    # well above the (near-zero) iid case.
    @test rw_neighbour > 0.5
    @test rw_neighbour > iid_neighbour + 0.4
end

@testitem "Hierarchy conforms to the latent interface and composes" begin
    using EpiAwarePrototype, Distributions, Accessors
    h = Hierarchy(FixedIntercept(0.0), (@optic _.intercept), IID(Normal(0.0, 1.0)))
    # The reusable latent-interface checker passes (Hierarchy is a latent process
    # over the grouping dimension): as_turing_model(h, n_groups) is a model.
    @test implements_latent_interface(h; n = 4)
    @test !implements_observation_interface(h)
    @test !implements_infection_interface(h)
end

@testitem "a partially-pooled model samples under NUTS (ForwardDiff)" tags=[:sample] begin
    using EpiAwarePrototype, Distributions, Accessors, Turing, Random
    using DynamicPPL: to_submodel
    Random.seed!(385)
    # Each group has its own intercept (level), partially pooled across groups by
    # the iid-Normal hierarchy; each level maps to Poisson counts. n_groups is
    # read from the data (the length of the observed group vector).
    h = Hierarchy(FixedIntercept(0.0), (@optic _.intercept), IID(Normal(0.0, 1.0)))

    @model function pooled_counts(h, y, n_groups)
        levels ~ to_submodel(as_turing_model(h, n_groups), false)
        for g in 1:n_groups
            y[g] ~ Poisson(exp(levels[g].intercept))
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
