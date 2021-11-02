q1 = nominal_configuration(quadruped)
v1 = zeros(quadruped.nq)

h = 0.01
T = 100

s = Simulator(quadruped, T, h=h)
s.ip.opts.diff_sol = true
simulate!(s, q1, v1)
@benchmark simulate!($s, $q1, $v1)

vis = Visualizer()
open(vis)
visualize!(vis, s)