@testset "Simulator: simulate" begin 
    q1 = nominal_configuration(biped)
    v1 = zeros(biped.nq)

    h = 0.01
    T = 100

    s = Simulator(biped, T, h=h, diff_sol=true, sim_opts=RobotDojo.SimulatorOptions(record=true))
    status = simulate!(s, q1, v1, reset_traj=true)
    @test status

    RobotDojo.process!(s.stats, 1)
    @test all(s.stats.policy_time .> 0.0)
    @test s.stats.policy_mean[1] > 0.0

    @test all(s.stats.sim_time .> 0.0)
    @test s.stats.sim_mean[1] > 0.0
end
