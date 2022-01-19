# PREAMBLE

# PKG_SETUP

# ## Setup

using RoboDojo 

# ## Initial conditions
q1 = nominal_configuration(box)
# q1[2] += 0.5
# q1[3] += 0.2 * π

v1 = zeros(box.nq)
v1[1] = 5.0

# ## Time 
T = 100
h = 0.01

# ## Simulator 
s = Simulator(box, T, h=h)

# ## Simulate
simulate!(s, q1, v1)

# ## Visualizer
vis = Visualizer()
render(vis)

# ## Visualize
visualize!(vis, s)