@testitem "AbstractLatentModel is a deprecated alias of AbstractPriorModel" begin
    using ComposableTuringIDModels
    @test AbstractLatentModel === AbstractPriorModel
    @test RandomWalk() isa AbstractLatentModel
    @test RandomWalk() isa AbstractPriorModel
    @test implements_latent_interface(RandomWalk())
end

@testitem "widened slots accept Distribution / vector / process" begin
    using ComposableTuringIDModels, Distributions
    # innovation slots
    @test AR(; ϵ_t = Normal()) isa AbstractPriorModel
    @test RandomWalk(; ϵ_t = Normal()) isa AbstractPriorModel
    @test MA(; ϵ_t = RandomWalk()) isa AbstractPriorModel
    # infection rt / Z slots
    data = IDData([0.2, 0.3, 0.5], exp)
    @test Renewal(data; rt = Normal()) isa AbstractInfectionModel
    @test DirectInfections(; Z = Normal()) isa AbstractInfectionModel
    @test ExpGrowthRate(; rt = RandomWalk()) isa AbstractInfectionModel
    # a fixed vector of per-element distributions still sets the AR order
    @test AR(;
        damp = [truncated(Normal(0.5, 0.1), 0, 1),
            truncated(Normal(0.2, 0.1), 0, 1)],
        init = [Normal(), Normal()]).p == 2
end

@testitem "bare-distribution dynamic slots build and sample" begin
    using ComposableTuringIDModels, Distributions, Turing, Random
    Random.seed!(11)
    # white-noise (iid) innovation: AR(ϵ_t = Normal())
    @test rand(as_turing_model(AR(; ϵ_t = Normal()), 8)) !== nothing
    # iid Rt: Renewal(rt = Normal()) inside a composed model
    data = IDData([0.2, 0.3, 0.5], exp)
    idmodel = IDModel(Renewal(data; rt = Normal(), initialisation = Normal()),
        PoissonError())
    y = as_turing_model(idmodel, missing, 8)().generated_y_t
    @test length(y) == 8
    chain = sample(as_turing_model(idmodel, y, 8),
        NUTS(0.8; adtype = Turing.AutoForwardDiff()), 30; progress = false)
    @test size(chain, 1) == 30
end

@testitem "a mismatched vector prior on a dynamic slot errors clearly" begin
    using ComposableTuringIDModels, Distributions
    # a length-2 vector innovation cannot produce the n-1 = 7 innovations of a
    # length-8 AR (a single Distribution or a process would adapt to length n-1).
    ar = AR(; ϵ_t = [Normal(), Normal()])
    @test_throws Exception as_turing_model(ar, 8)()
end

@testitem "manipulators accept a bare-Distribution member" begin
    using ComposableTuringIDModels, Distributions, Random
    Random.seed!(12)
    # a Distribution member composes as iid draws alongside a process
    comb = CombineLatentModels([Normal(2, 0.2), RandomWalk()])
    @test length(as_turing_model(comb, 10)()) == 10
    conc = ConcatLatentModels([Normal(2, 0.2), RandomWalk()])
    @test length(as_turing_model(conc, 10)()) == 10
    dm = DiffLatentModel(; model = Normal(), init = [Normal(), Normal()])
    @test length(as_turing_model(dm, 10)()) == 10
end
