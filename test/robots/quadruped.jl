@testset "Robots: quadruped" begin 
    q0 = [0.0; 0.0; 0.0; 0.5 * π * ones(8)] 
    v0 = zeros(11)
    @test norm(lagrangian(quadruped, q0, v0)) < 1.0e-8

    # visualize
    vis = RobotDojo.Visualizer();
    @test visualize!(vis, quadruped, [q0], Δt=0.1);
end
