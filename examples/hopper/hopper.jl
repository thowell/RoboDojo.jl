
q1 = nominal_configuration(hopper)
v1 = zeros(hopper.nq)

T = 100
h = 0.01

# p_raibert = raibert_policy(hopper, h=h)
s = Simulator(hopper, T, h=h)#, policy=p_raibert)
s.ip.opts.diff_sol = true
simulate!(s, q1, v1, reset_traj=false)
using BenchmarkTools
@benchmark simulate!($s, $q1, $v1)

# fill!(s.ip.rz, 0.0)
# s.ip.rz
# RobotDojo.rz!(s.ip, s.ip.rz, s.ip.z, s.ip.θ)
# s.ip.methods.rθ!(s.ip.rθ, s.ip.z, s.ip.θ)
# # linear_solve!(s.ip.solver, s.ip.δzs, s.ip.rz, s.ip.rθ)
# s.ip.solver.A
# fill!(s.ip.solver.A, 0.0)
# s.ip.rz
# s.ip.solver.A .= s.ip.rz
# fill!(s.ip.solver.ipiv, 0)
# RobotDojo.getrf!(s.ip.solver.A, s.ip.solver.ipiv, s.ip.solver.lda, s.ip.solver.info)
# RobotDojo.factorize!(s.ip.solver, copy(s.ip.rz))

# RobotDojo.differentiate_solution!(s.ip)
# s.ip.solver.A

vis = Visualizer()
open(vis)
visualize!(vis, s)