@testset "Robots: box" begin 
    # TODO: add tests
    q0 = nominal_configuration(box)
    # visualizer 
    vis = RoboDojo.Visualizer();
    @test visualize!(vis, box, [q0], Δt=0.1);
end