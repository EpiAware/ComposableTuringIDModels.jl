@testitem "composed model: prior simulation and generated quantities" begin
    using ComposableTuringIDModels, Distributions, Random
    Random.seed!(21)
    model = IDModel(
        DirectInfections(; Z = RandomWalk(), initialisation = Normal()),
        PoissonError()
    )
    n = 20
    mdl = as_turing_model(model, missing, n)

    draw = rand(mdl)
    @test !isempty(keys(draw))

    gen = mdl()
    @test length(gen.generated_y_t) == n
    @test length(gen.I_t) == n
    @test length(gen.Z_t) == n
    @test all(>=(0), gen.I_t)
end

@testitem "concrete_observations narrows only when nothing is missing" begin
    using ComposableTuringIDModels: concrete_observations, MissingObservations

    # Fully concrete data (no `missing` left) is narrowed to a concrete eltype.
    y = Vector{Union{Missing, Int}}([1, 2, 3])
    narrowed = concrete_observations(y)
    @test narrowed == [1, 2, 3]
    @test eltype(narrowed) == Int

    # A genuine gap in a vector is split into a `MissingObservations` carrier:
    # a concrete value vector plus a presence mask, so no `Missing` survives in
    # either field's type (narrowing the gap away would drop it entirely).
    y_gap = Vector{Union{Missing, Int}}([1, missing, 3])
    split = concrete_observations(y_gap)
    @test split isa MissingObservations
    @test !(eltype(split.value) >: Missing)
    @test split.present == [true, false, true]
    @test split.value[[1, 3]] == [1, 3]

    # A fully missing vector survives unchanged: there is nothing concrete to
    # split out.
    y_allgap = Vector{Union{Missing, Int}}(missing, 3)
    @test concrete_observations(y_allgap) === y_allgap

    # A `NamedTuple` of data streams (e.g. `BinomialError`'s `(y, N)`, or a
    # `Split` observation's per-stream series) is narrowed field-wise.
    nt = (y = Vector{Union{Missing, Int}}([4, 5]), N = [10, 10])
    narrowed_nt = concrete_observations(nt)
    @test eltype(narrowed_nt.y) == Int

    # A scalar `missing` (the "simulate from the prior" sentinel) passes
    # through unchanged.
    @test concrete_observations(missing) === missing

    # An already-concrete array is returned as-is, not copied, so a caller
    # holding the array keeps its identity.
    y_concrete = [1, 2, 3]
    @test concrete_observations(y_concrete) === y_concrete

    # A partially-missing MATRIX (the data-driven `Split` strata shape) is left
    # unchanged: only vectors are split into a `MissingObservations` carrier.
    y_mat_gap = Matrix{Union{Missing, Int}}([1 2; missing 4])
    @test concrete_observations(y_mat_gap) === y_mat_gap
end

@testitem "composed model: as_turing_model narrows y_t automatically" begin
    using ComposableTuringIDModels, Distributions, Random
    using DynamicPPL: DynamicPPL, logjoint
    Random.seed!(24)
    model = IDModel(
        DirectInfections(; Z = RandomWalk(), initialisation = Normal()),
        PoissonError()
    )
    n = 10
    # Data simulated from the prior carries DynamicPPL's `Union{Missing, T}`
    # predictive eltype even though every entry is concrete.
    y = as_turing_model(model, missing, n)().generated_y_t
    @test eltype(y) >: Missing

    # `as_turing_model` narrows it before storing it as a model argument.
    posterior = as_turing_model(model, y, n)
    @test !(eltype(posterior.args.y_t) >: Missing)
    @test posterior isa DynamicPPL.Model
    @test isfinite(logjoint(posterior, rand(posterior)))
end

@testitem "composed model: fix and condition" begin
    using ComposableTuringIDModels, Distributions, Random
    using DynamicPPL: fix, condition
    Random.seed!(22)
    model = IDModel(
        DirectInfections(; Z = RandomWalk(), initialisation = Normal()),
        PoissonError()
    )
    n = 20
    mdl = as_turing_model(model, missing, n)

    # fix removes a parameter from the sampled set.
    fixed = fix(mdl, (init_incidence = 0.0,))
    fixed_names = string.(collect(keys(rand(fixed))))
    @test !("init_incidence" in fixed_names)

    # condition and the | syntax produce equivalent conditioned models.
    c1 = condition(mdl, (std = 0.1,))
    c2 = mdl | (std = 0.1,)
    @test typeof(c1) == typeof(c2)
end

@testitem "composed model: short NUTS sample runs" tags = [:sample] begin
    using ComposableTuringIDModels, Distributions, Turing, Random
    Random.seed!(23)
    model = IDModel(
        DirectInfections(; Z = RandomWalk(), initialisation = Normal()),
        NegativeBinomialError()
    )
    n = 20
    y = as_turing_model(model, missing, n)()
    cond_model = as_turing_model(model, y.generated_y_t, n)
    chn = sample(cond_model, NUTS(), 50; progress = false)
    @test chn !== nothing
end

# --- the Stratify + Split panel replacement -------------------------------
#
# A grouped/panel model is a shared `R_t` with partially pooled per-group
# deviations (`Stratify` + `Hierarchy`) observed as a data-driven strata
# `Split`, resolved from the data alone via the two-argument
# `as_turing_model(model, Y)`.

@testitem "Stratify + Split simulates a panel, grouping dim from the data" begin
    using ComposableTuringIDModels, Distributions, Random
    Random.seed!(451)
    model = IDModel(
        DirectInfections(;
            Z = Stratify(
                RandomWalk(), Hierarchy(; across = IID(Normal(0, 0.5)))
            ),
            initialisation = Normal(log(50.0), 0.2)
        ),
        Split(PoissonError())
    )
    for (n_time, n_groups) in ((24, 4), (12, 6))
        Ymiss = Matrix{Union{Missing, Float64}}(missing, n_groups, n_time)
        sim = as_turing_model(model, Ymiss)()
        @test size(sim.I_t) == (n_groups, n_time)
        @test length(sim.generated_y_t) == n_groups
        @test all(length(sim.generated_y_t[g]) == n_time for g in 1:n_groups)
        @test all(>=(0), reduce(vcat, collect(sim.generated_y_t)))
    end
end

@testitem "Stratify + Split: variable names are namespaced per group" begin
    using ComposableTuringIDModels, Distributions
    model = IDModel(
        DirectInfections(;
            Z = Stratify(
                RandomWalk(), Hierarchy(; across = IID(Normal(0, 0.5)))
            ),
            initialisation = Normal()
        ),
        Split(PoissonError())
    )
    Ymiss = Matrix{Union{Missing, Float64}}(missing, 3, 8)
    mdl = as_turing_model(model, Ymiss)
    names = string.(collect(keys(rand(mdl))))
    @test any(startswith(n, "group1.") for n in names)
    @test any(startswith(n, "group2.") for n in names)
end

@testitem "Stratify + Split panel samples under NUTS" tags = [:sample] begin
    using ComposableTuringIDModels, Distributions, Turing, Random
    Random.seed!(388)
    model = IDModel(
        DirectInfections(;
            Z = Stratify(
                RandomWalk(), Hierarchy(;
                    mean = Normal(0.0, 0.5),
                    across = IID(Normal(0.0, 0.5))
                )
            ),
            initialisation = Normal(log(50.0), 0.2)
        ),
        Split(PoissonError())
    )

    n_time, n_groups = 16, 5
    Ymiss = Matrix{Union{Missing, Float64}}(missing, n_groups, n_time)
    sim = as_turing_model(model, Ymiss)()
    Ydata = Float64.(
        reduce(
            vcat, [
                permutedims(sim.generated_y_t[g])
                    for g in 1:n_groups
            ]
        )
    )

    posterior = as_turing_model(model, Ydata)
    chain = sample(
        posterior, NUTS(0.8; adtype = Turing.AutoForwardDiff()), 30;
        progress = false
    )
    @test size(chain, 1) == 30
end
