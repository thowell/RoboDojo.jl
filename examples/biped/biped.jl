# PREAMBLE

# PKG_SETUP

# ## Setup

using RoboDojo 

# ## Initial conditions
q1 = nominal_configuration(biped)
v1 = zeros(biped.nq)

# ## Time
h = 0.01
T = 100

# ## Simulator
s = Simulator(biped, T, h=h)

# ## Simulate
status = simulate!(s, q1, v1)

# ## Visualizer
vis = Visualizer()
render(vis)

# ## Visualize
visualize!(vis, s)