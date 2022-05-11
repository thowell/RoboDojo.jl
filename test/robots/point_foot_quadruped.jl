@testset "Robots: point foot quadruped" begin
    # TODO: add tests
    q0 = nominal_configuration(point_foot_quadruped)
    # visualizer
    vis = RoboDojo.Visualizer();
    @test visualize!(vis, box, [point_foot_quadruped], Δt=0.1);
end
