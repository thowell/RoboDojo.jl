
q1 = nominal_configuration(hopper)
v1 = zeros(hopper.nq)

T = 100
h = 0.01

# p_raibert = raibert_policy(hopper, h=h)
s = Simulator(hopper, T, h=h)#, policy=p_raibert)
s.ip.opts.diff_sol = true
simulate!(s, q1, v1, reset_traj=false)
@benchmark simulate!($s, $q1, $v1)

vis = Visualizer()
open(vis)
visualize!(vis, s)