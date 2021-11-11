@testset "Robots: biped" begin 
    q0 = [0.0; 0.0; 0.5 * π * ones(7)]
    v0 = zeros(9)
    @test norm(lagrangian(biped, q0, v0)) < 1.0e-8

    # visualize
    vis = RobotDojo.Visualizer();
    @test visualize!(vis, biped, [q0], Δt=0.1);
end