@testitem "AbstractLatentModel is a deprecated alias of AbstractPriorModel" begin
    using ComposableTuringIDModels
    @test AbstractLatentModel === AbstractPriorModel
    @test RandomWalk() isa AbstractLatentModel
    @test RandomWalk() isa AbstractPriorModel
    @test implements_latent_interface(RandomWalk())
end

@testitem "widened slots accept Distribution / IID / vector / process" begin
    using ComposableTuringIDModels, Distributions
    # innovation slots take a process (or IID for white noise)
    @test AR(; ϵ_t = IID(Normal())) isa AbstractPriorModel
    @test RandomWalk(; ϵ_t = HierarchicalNormal()) isa AbstractPriorModel
    @test MA(; ϵ_t = RandomWalk()) isa AbstractPriorModel
    # infection rt / Z process slots take a process (IID for iid)
    gen_int = [0.2, 0.3, 0.5]
    @test Renewal(; generation_time = gen_int, rt = IID(Normal())) isa
          AbstractInfectionModel
    @test DirectInfections(; Z = RandomWalk()) isa AbstractInfectionModel
    @test ExpGrowthRate(; rt = RandomWalk()) isa AbstractInfectionModel
    # a fixed vector of per-element distributions still sets the AR order
    @test AR(;
        damp = [truncated(Normal(0.5, 0.1), 0, 1),
            truncated(Normal(0.2, 0.1), 0, 1)],
        init = [Normal(), Normal()]).p == 2
    # a bare Distribution in a per-step PARAMETER slot is a scalar constant
    @test AR(; damp = Normal()).p == 1
end

@testitem "IID / process dynamic slots build and sample" begin
    using ComposableTuringIDModels, Distributions, Turing, Random
    Random.seed!(11)
    # white-noise (iid) innovation via IID()
    @test rand(as_turing_model(AR(; ϵ_t = IID(Normal())), 8)) !== nothing
    # iid Rt via IID() inside a composed model
    gen_int = [0.2, 0.3, 0.5]
    idmodel = IDModel(
        Renewal(; generation_time = gen_int, rt = IID(Normal()),
            initialisation = Normal()),
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
    # length-8 AR (an IID() or a process would adapt to length n-1).
    ar = AR(; ϵ_t = [Normal(), Normal()])
    @test_throws Exception as_turing_model(ar, 8)()
end

@testitem "manipulators accept an IID / process member" begin
    using ComposableTuringIDModels, Distributions, Random
    Random.seed!(12)
    # an IID() member composes as iid draws alongside a process
    comb = CombineLatentModels([IID(Normal(2, 0.2)), RandomWalk()])
    @test length(as_turing_model(comb, 10)()) == 10
    conc = ConcatLatentModels([IID(Normal(2, 0.2)), RandomWalk()])
    @test length(as_turing_model(conc, 10)()) == 10
    dm = DiffLatentModel(; model = IID(Normal()), init = [Normal(), Normal()])
    @test length(as_turing_model(dm, 10)()) == 10
end
