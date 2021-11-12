@testset "Simulator: simulate" begin 
    ## hopper
    q1 = nominal_configuration(hopper)
    v1 = zeros(hopper.nq)

    h = 0.01
    T = 100

    # simulator 
    s = Simulator(hopper, T, h=h, diff_sol=true, sim_opts=RobotDojo.SimulatorOptions(record=true))

    # step
    q2 = step!(s, q1, v1, zeros(hopper.nu), 1)
    @test sum(q2) != 0.0

    # simulate
    status = simulate!(s, q1, v1, reset_traj=true)
    @test status

    ## quadruped
    q1 = nominal_configuration(quadruped)
    v1 = zeros(quadruped.nq)

    h = 0.01
    T = 100

    # simulator 
    s = Simulator(quadruped, T, h=h, diff_sol=true, sim_opts=RobotDojo.SimulatorOptions(record=true))

    # step
    q2 = step!(s, q1, v1, zeros(quadruped.nu), 1)
    @test sum(q2) != 0.0

    # simulate
    status = simulate!(s, q1, v1, reset_traj=true)
    @test status

    ## biped 
    q1 = nominal_configuration(biped)
    v1 = zeros(biped.nq)

    h = 0.01
    T = 100

    # simulator 
    s = Simulator(biped, T, h=h, diff_sol=true, sim_opts=RobotDojo.SimulatorOptions(record=true))

    # step
    q2 = step!(s, q1, v1, zeros(biped.nu), 1)
    @test sum(q2) != 0.0

    # simulate
    status = simulate!(s, q1, v1, reset_traj=true)
    @test status

    # process timing results
    RobotDojo.process!(s.stats, 1)
    @test sum(s.stats.policy_time) > 0.0
    @test s.stats.policy_mean[1] > 0.0

    @test sum(s.stats.sim_time) > 0.0
    @test s.stats.sim_mean[1] > 0.0
end
