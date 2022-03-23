@testset "Robots: centroidal quadruped" begin 
    # TODO: add tests
    q0 = nominal_configuration(centroidal_quadruped)
    # visualizer 
    vis = RoboDojo.Visualizer();
    @test visualize!(vis, box, [centroidal_quadruped], Δt=0.1);
end