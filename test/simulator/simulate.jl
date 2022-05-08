# @testset "Simulator: simulate" begin
#     ## hopper
#     q1 = nominal_configuration(hopper)
#     v1 = zeros(hopper.nq)
#
#     h = 0.01
#     T = 100
#
#     # simulator
#     s = Simulator(hopper, T, h=h, diff_sol=true, sim_opts=RoboDojo.SimulatorOptions(record=true))
#
#     # step
#     q2 = step!(s, q1, v1, zeros(hopper.nu), 1)
#     @test sum(q2) != 0.0
#
#     # simulate
#     status = simulate!(s, q1, v1, reset_traj=true)
#     @test status
#
#
#     ## quadruped
#     q1 = nominal_configuration(quadruped)
#     v1 = zeros(quadruped.nq)
#
#     h = 0.01
#     T = 100
#
#     # simulator
#     s = Simulator(quadruped, T, h=h, diff_sol=true, sim_opts=RoboDojo.SimulatorOptions(record=true))
#
#     # step
#     q2 = step!(s, q1, v1, zeros(quadruped.nu), 1)
#     @test sum(q2) != 0.0
#
#     # simulate
#     status = simulate!(s, q1, v1, reset_traj=true)
#     @test status
#
#
#     ## halfquadruped
#     q1 = nominal_configuration(halfquadruped)
#     v1 = zeros(halfquadruped.nq)
#
#     h = 0.01
#     T = 100
#
#     # simulator
#     s = Simulator(halfquadruped, T, h=h, diff_sol=true, sim_opts=RoboDojo.SimulatorOptions(record=true))
#
#     # step
#     q2 = step!(s, q1, v1, zeros(halfquadruped.nu), 1)
#     @test sum(q2) != 0.0
#
#     # simulate
#     status = simulate!(s, q1, v1, reset_traj=true)
#     @test status
#
#
#     ## halfcheetah
#     q1 = nominal_configuration(halfcheetah)
#     v1 = zeros(halfcheetah.nq)
#
#     h = 0.01
#     T = 100
#
#     # simulator
#     s = Simulator(halfcheetah, T, h=h, diff_sol=true, sim_opts=RoboDojo.SimulatorOptions(record=true))
#
#     # step
#     q2 = step!(s, q1, v1, zeros(halfcheetah.nu), 1)
#     @test sum(q2) != 0.0
#
#     # simulate
#     status = simulate!(s, q1, v1, reset_traj=true)
#     @test status
#
#
#     ## biped
#     q1 = nominal_configuration(biped)
#     v1 = zeros(biped.nq)
#
#     h = 0.01
#     T = 100
#
#     # simulator
#     s = Simulator(biped, T, h=h, diff_sol=true, sim_opts=RoboDojo.SimulatorOptions(record=true))
#
#     # step
#     q2 = step!(s, q1, v1, zeros(biped.nu), 1)
#     @test sum(q2) != 0.0
#
#     # simulate
#     status = simulate!(s, q1, v1, reset_traj=true)
#     @test status
#
#     # process timing results
#     RoboDojo.process!(s.stats, 1)
#     @test sum(s.stats.policy_time) > 0.0
#     @test s.stats.policy_mean[1] > 0.0
#
#     @test sum(s.stats.sim_time) > 0.0
#     @test s.stats.sim_mean[1] > 0.0
#
#
#     ## box
#     q1 = nominal_configuration(box)
#     v1 = zeros(box.nq)
#
#     h = 0.01
#     T = 100
#
#     # simulator
#     s = Simulator(box, T, h=h, diff_sol=true, sim_opts=RoboDojo.SimulatorOptions(record=true))
#
#     # step
#     q2 = step!(s, q1, v1, zeros(box.nu), 1)
#     @test sum(q2) != 0.0
#
#     # simulate
#     status = simulate!(s, q1, v1, reset_traj=true)
#     @test status
#
#     # process timing results
#     RoboDojo.process!(s.stats, 1)
#     @test sum(s.stats.policy_time) > 0.0
#     @test s.stats.policy_mean[1] > 0.0
#
#     @test sum(s.stats.sim_time) > 0.0
#     @test s.stats.sim_mean[1] > 0.0
#
#
#     ## particle
#     q1 = nominal_configuration(particle)
#     v1 = zeros(particle.nq)
#
#     h = 0.01
#     T = 100
#
#     # simulator
#     s = Simulator(particle, T, h=h, diff_sol=true, sim_opts=RoboDojo.SimulatorOptions(record=true))
#
#     # step
#     q2 = step!(s, q1, v1, zeros(particle.nu), 1)
#     @test sum(q2) != 0.0
#
#     # simulate
#     status = simulate!(s, q1, v1, reset_traj=true)
#     @test status
#
#     # process timing results
#     RoboDojo.process!(s.stats, 1)
#     @test sum(s.stats.policy_time) > 0.0
#     @test s.stats.policy_mean[1] > 0.0
#
#     @test sum(s.stats.sim_time) > 0.0
#     @test s.stats.sim_mean[1] > 0.0
#
#
#     ## particle (time budget)
#     q1 = nominal_configuration(particle)
#     v1 = zeros(particle.nq)
#
#     h = 0.01
#     T = 100
#
#     # simulator
#     s = Simulator(particle, T, h=h, diff_sol=true,
#         sim_opts=RoboDojo.SimulatorOptions(record=true),
#         solver_opts=solver_opts=RoboDojo.InteriorPointOptions(
#             max_time=0.49,
#             undercut=Inf,
#             γ_reg=0.1,
#             r_tol=1e-8,
#             κ_tol=1e-8,
#             max_ls=25,
#             ϵ_min=0.25,
#             diff_sol=false,
#             verbose=false),)
#
#     # step
#     q2 = step!(s, q1, v1, zeros(particle.nu), 1)
#     @test sum(q2) != 0.0
#
#     # simulate
#     status = simulate!(s, q1, v1, reset_traj=true)
#     @test status
#
#      # simulator
#      s = Simulator(particle, T, h=h, diff_sol=true,
#      sim_opts=RoboDojo.SimulatorOptions(record=true),
#      solver_opts=solver_opts=RoboDojo.InteriorPointOptions(
#          max_time=1.0e-6,
#          undercut=Inf,
#          γ_reg=0.1,
#          r_tol=1e-8,
#          κ_tol=1e-8,
#          max_ls=25,
#          ϵ_min=0.25,
#          diff_sol=false,
#          verbose=false),)
#
#     # simulate
#     status = simulate!(s, q1, v1, reset_traj=true)
#     @test !status
#
#
#     ## centroidal quadruped
#     q1 = nominal_configuration(centroidal_quadruped)
#     v1 = zeros(centroidal_quadruped.nq)
#
#     h = 0.01
#     T = 100
#
#     # simulator
#     s = Simulator(centroidal_quadruped, T, h=h, diff_sol=true, sim_opts=RoboDojo.SimulatorOptions(record=true))
#
#     # step
#     q2 = step!(s, q1, v1, zeros(centroidal_quadruped.nu), 1)
#     @test sum(q2) != 0.0
#
#     # simulate
#     status = simulate!(s, q1, v1, reset_traj=true)
#     @test status
#
#     # process timing results
#     RoboDojo.process!(s.stats, 1)
#     @test sum(s.stats.policy_time) > 0.0
#     @test s.stats.policy_mean[1] > 0.0
#
#     @test sum(s.stats.sim_time) > 0.0
#     @test s.stats.sim_mean[1] > 0.0
# end



## hopper
q1 = nominal_configuration(hopper)
v1 = zeros(hopper.nq)

h = 0.01
T = 100

# simulator
s = Simulator(hopper, T, h=h, diff_sol=true, sim_opts=RoboDojo.SimulatorOptions(record=true))

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
s = Simulator(quadruped, T, h=h, diff_sol=true, sim_opts=RoboDojo.SimulatorOptions(record=true))

# step
q2 = step!(s, q1, v1, zeros(quadruped.nu), 1)
@test sum(q2) != 0.0

# simulate
status = simulate!(s, q1, v1, reset_traj=true)
@test status


## halfquadruped
q1 = nominal_configuration(halfquadruped)
v1 = zeros(halfquadruped.nq)

h = 0.01
T = 100

# simulator
s = Simulator(halfquadruped, T, h=h, diff_sol=true, sim_opts=RoboDojo.SimulatorOptions(record=true))

# step
q2 = step!(s, q1, v1, zeros(halfquadruped.nu), 1)
@test sum(q2) != 0.0

# simulate
status = simulate!(s, q1, v1, reset_traj=true)
@test status


## halfcheetah
q1 = nominal_configuration(halfcheetah)
v1 = zeros(halfcheetah.nq)

h = 0.01
T = 100

# simulator
s = Simulator(halfcheetah, T, h=h, diff_sol=true, sim_opts=RoboDojo.SimulatorOptions(record=true))

# step
q2 = step!(s, q1, v1, zeros(halfcheetah.nu), 1)
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
s = Simulator(biped, T, h=h, diff_sol=true, sim_opts=RoboDojo.SimulatorOptions(record=true))

# step
q2 = step!(s, q1, v1, zeros(biped.nu), 1)
@test sum(q2) != 0.0

# simulate
status = simulate!(s, q1, v1, reset_traj=true)
@test status

# process timing results
RoboDojo.process!(s.stats, 1)
@test sum(s.stats.policy_time) > 0.0
@test s.stats.policy_mean[1] > 0.0

@test sum(s.stats.sim_time) > 0.0
@test s.stats.sim_mean[1] > 0.0


## box
q1 = nominal_configuration(box)
v1 = zeros(box.nq)

h = 0.01
T = 100

# simulator
s = Simulator(box, T, h=h, diff_sol=true, sim_opts=RoboDojo.SimulatorOptions(record=true))

# step
q2 = step!(s, q1, v1, zeros(box.nu), 1)
@test sum(q2) != 0.0

# simulate
status = simulate!(s, q1, v1, reset_traj=true)
@test status

# process timing results
RoboDojo.process!(s.stats, 1)
@test sum(s.stats.policy_time) > 0.0
@test s.stats.policy_mean[1] > 0.0

@test sum(s.stats.sim_time) > 0.0
@test s.stats.sim_mean[1] > 0.0


## particle
q1 = nominal_configuration(particle)
v1 = zeros(particle.nq)

h = 0.01
T = 100

# simulator
s = Simulator(particle, T, h=h, diff_sol=true, sim_opts=RoboDojo.SimulatorOptions(record=true))

# step
q2 = step!(s, q1, v1, zeros(particle.nu), 1)
@test sum(q2) != 0.0

# simulate
status = simulate!(s, q1, v1, reset_traj=true)
@test status

# process timing results
RoboDojo.process!(s.stats, 1)
@test sum(s.stats.policy_time) > 0.0
@test s.stats.policy_mean[1] > 0.0

@test sum(s.stats.sim_time) > 0.0
@test s.stats.sim_mean[1] > 0.0


## particle (time budget)
q1 = nominal_configuration(particle)
v1 = zeros(particle.nq)

h = 0.01
T = 100

# simulator
s = Simulator(particle, T, h=h, diff_sol=true,
    sim_opts=RoboDojo.SimulatorOptions(record=true),
    solver_opts=solver_opts=RoboDojo.InteriorPointOptions(
        max_time=0.49,
        undercut=Inf,
        γ_reg=0.1,
        r_tol=1e-8,
        κ_tol=1e-8,
        max_ls=25,
        ϵ_min=0.25,
        diff_sol=false,
        verbose=false),)

# step
q2 = step!(s, q1, v1, zeros(particle.nu), 1)
@test sum(q2) != 0.0

# simulate
status = simulate!(s, q1, v1, reset_traj=true)
@test status

 # simulator
 s = Simulator(particle, T, h=h, diff_sol=true,
 sim_opts=RoboDojo.SimulatorOptions(record=true),
 solver_opts=solver_opts=RoboDojo.InteriorPointOptions(
     max_time=1.0e-6,
     undercut=Inf,
     γ_reg=0.1,
     r_tol=1e-8,
     κ_tol=1e-8,
     max_ls=25,
     ϵ_min=0.25,
     diff_sol=false,
     verbose=false),)

# simulate
status = simulate!(s, q1, v1, reset_traj=true)
@test !status


## centroidal quadruped
q1 = nominal_configuration(centroidal_quadruped)
v1 = zeros(centroidal_quadruped.nq)

h = 0.01
T = 100

# simulator
s = Simulator(centroidal_quadruped, T, h=h, diff_sol=true, sim_opts=RoboDojo.SimulatorOptions(record=true))

# step
q2 = step!(s, q1, v1, zeros(centroidal_quadruped.nu), 1)
@test sum(q2) != 0.0

# simulate
status = simulate!(s, q1, v1, reset_traj=true)
@test status

# process timing results
RoboDojo.process!(s.stats, 1)
@test sum(s.stats.policy_time) > 0.0
@test s.stats.policy_mean[1] > 0.0

@test sum(s.stats.sim_time) > 0.0
@test s.stats.sim_mean[1] > 0.0
