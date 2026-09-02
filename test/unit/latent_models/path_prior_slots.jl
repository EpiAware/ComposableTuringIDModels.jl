# Every construction path must widen a PATH slot identically. A bare
# `Distribution` in a length-`n` PATH slot is wrapped in an `Intercept`, so the
# positional and keyword forms have to agree on the stored field type and on
# the generated path length. Widening belongs in an inner constructor: one
# reached only through an outer keyword constructor leaves the positional form
# storing an unwidened slot, which draws a scalar and sizes the path wrongly.

@testitem "PATH slots widen identically on every construction path" begin
    using ComposableTuringIDModels, Distributions
    using ComposableTuringIDModels: path_prior

    d = Normal(0, 0.05)
    damp = truncated(Normal(0.0, 0.05), 0, 1)
    θ = [truncated(Normal(0.0, 0.05), -1, 1)]

    # (name, positional form, keyword form, PATH slots, shape, path length)
    cases = [
        (
            "RandomWalk", RandomWalk(Normal(), d),
            RandomWalk(; init = Normal(), ϵ_t = d), (:ϵ_t,), 20, 20,
        ),
        ("MA", MA(θ, 1, d), MA(; θ = θ, ϵ_t = d), (:ϵ_t,), 20, 20),
        (
            "AR", AR(damp, Normal(), 1, d, identity),
            AR(; damp = damp, init = Normal(), ϵ_t = d), (:ϵ_t,), 20, 20,
        ),
        (
            "Stratify", Stratify(d, d, +), Stratify(d, d),
            (:shared, :across), (3, 20), 60,
        ),
    ]

    # The infection models take the same treatment on their latent slot, but
    # return `(; I_t, Z_t)` rather than a bare path, so they are measured on
    # `I_t`.
    infection_cases = [
        (
            "DirectInfections", DirectInfections(d, exp, Normal()),
            DirectInfections(; Z = d), :Z, 20,
        ),
        (
            "ExpGrowthRate",
            ExpGrowthRate(d, ComposableTuringIDModels._oneexpy, Normal()),
            ExpGrowthRate(; rt = d), :rt, 20,
        ),
    ]

    for (name, pos, kw, slots, shape, len) in cases
        for slot in slots
            # The positional form must widen, not store the bare distribution.
            @test getfield(pos, slot) isa Intercept
            @test typeof(getfield(pos, slot)) == typeof(getfield(kw, slot))
        end
        # Both forms generate the same-shaped path. An unwidened slot draws a
        # single scalar, which sizes the path by the recursion rather than `n`.
        @test length(as_turing_model(pos, shape)()) == len
        @test length(as_turing_model(kw, shape)()) == len
    end

    for (name, pos, kw, slot, n) in infection_cases
        @test getfield(pos, slot) isa Intercept
        @test typeof(getfield(pos, slot)) == typeof(getfield(kw, slot))
        # An unwidened latent slot sizes `I_t` by the scalar draw rather than
        # by `n`, so the length is the assertion that catches it.
        @test length(as_turing_model(pos, n)().I_t) == n
        @test length(as_turing_model(kw, n)().I_t) == n
        @test path_prior(getfield(pos, slot)) === getfield(pos, slot)
    end

    # `path_prior` is idempotent, so rebuilding a component from its own stored
    # fields is a fixed point. `Accessors.@set` and any other field-wise
    # reconstruction depend on this.
    @test path_prior(path_prior(d)) === path_prior(d)
    for (_, pos, _, slots, _, _) in cases
        for slot in slots
            @test path_prior(getfield(pos, slot)) === getfield(pos, slot)
        end
    end
end

@testitem "rebuilding a latent model from its stored fields is a fixed point" begin
    using ComposableTuringIDModels, Distributions
    using Accessors: @set

    d = Normal(0, 0.05)
    rw = RandomWalk(; init = Normal(), ϵ_t = d)
    # Reconstructing from the stored fields must not re-wrap the already
    # widened slot.
    rebuilt = RandomWalk(rw.init, rw.ϵ_t)
    @test typeof(rebuilt) == typeof(rw)
    @test rebuilt.ϵ_t === rw.ϵ_t
    # `Accessors` rebuilds through the positional constructor, so a set on one
    # field must leave the other exactly as stored.
    moved = @set rw.init = Normal(1, 1)
    @test moved.ϵ_t === rw.ϵ_t
    @test length(as_turing_model(moved, 20)()) == 20

    st = Stratify(d, d)
    @test typeof(Stratify(st.shared, st.across, st.combine)) == typeof(st)
    st2 = @set st.combine = *
    @test st2.shared === st.shared
    @test st2.across === st.across

    di = DirectInfections(; Z = d)
    @test typeof(
        DirectInfections(di.Z, di.transformation, di.initialisation)
    ) == typeof(di)
    di2 = @set di.initialisation = Normal(1, 1)
    @test di2.Z === di.Z
    @test length(as_turing_model(di2, 20)().I_t) == 20
end

@testitem "the role guard survives moving widening into the constructor" begin
    using ComposableTuringIDModels, Distributions
    # The widening happens before `new`, so the slot's `PriorLike` bound still
    # rejects a wrong-role component at construction rather than at sampling.
    @test_throws TypeError RandomWalk(PoissonError(), Normal())
    @test_throws TypeError MA(PoissonError(), 1, Normal())
    @test_throws TypeError AR(PoissonError(), Normal(), 1, Normal(), identity)
    @test_throws TypeError Stratify(PoissonError(), Normal(), +)
    @test_throws TypeError DirectInfections(PoissonError(), exp, Normal())
    @test_throws TypeError ExpGrowthRate(PoissonError(), exp, Normal())
end
