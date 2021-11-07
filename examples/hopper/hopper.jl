q1 = nominal_configuration(hopper)
v1 = zeros(hopper.nq)

T = 100
h = 0.01

s = Simulator(hopper, T, h=h)
simulate!(s, q1, v1)

vis = Visualizer()
open(vis)
visualize!(vis, s)