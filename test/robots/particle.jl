@testset "Robots: particle" begin 
    # TODO: add tests
    q0 = nominal_configuration(particle)
    # visualizer 
    vis = RoboDojo.Visualizer();
    @test visualize!(vis, particle, [q0], Δt=0.1);
end