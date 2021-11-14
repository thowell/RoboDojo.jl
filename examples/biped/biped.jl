
q1 = nominal_configuration(biped)
v1 = zeros(biped.nq)

h = 0.01
T = 100

s = Simulator(biped, T, h=h)
status = simulate!(s, q1, v1)

vis = Visualizer()
open(vis)
visualize!(vis, biped, [q1])

q1 = [0.0; 0.8; -0.01; -0.15 * π; 0.1 * π; 0.0; -0.15 * π; 0.1 * π; 0.0]
