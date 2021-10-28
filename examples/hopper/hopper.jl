q1 = nominal_configuration(hopper)
v1 = zeros(hopper.nq)

T = 100
h = 0.01

# p_raibert = raibert_policy(hopper, h=h)
s = Simulator(hopper, T, h=h)#, policy=p_raibert)

num_var(hopper)
num_data(hopper)

s.ip.opts.diff_sol = false
simulate!(s, q1, v1, reset_traj=false)
@code_warntype simulate!(s, q1, v1, reset_traj=false)

@benchmark simulate!($s, $q1, $v1)

rank(s.ip.rz)
s.ip.solver.A
rz!(s.ip, s.ip.rz, s.ip.z, s.ip.θ, reg = 0.0)
rθ!(s.ip, s.ip.rθ, s.ip.z, s.ip.θ)

differentiate_solution!(s.ip) 
s.ip.solver.A #.= 0.0
s.ip.solver.A
s.ip.solver.ipiv
s.ip.solver.info
s.ip.δz

gradient!(s.grad, s.ip.δz, s.idx_z, s.idx_θ, h, 1)
@benchmark gradient!($s.grad, $s.ip.δz, $s.idx_z, $s.idx_θ, $h, 1)


using InteractiveUtils
t = 1
gradient!(s.grad, s.ip.δz, s.idx_z, s.idx_θ, h, t)
@code_warntype gradient!(s.grad, s.ip.δz, s.idx_z, s.idx_θ, h, t)

vis = Visualizer()
open(vis)
visualize!(vis, s)