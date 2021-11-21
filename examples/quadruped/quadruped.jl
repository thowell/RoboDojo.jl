# PREAMBLE

# PKG_SETUP

# ## Setup

using RoboDojo 

# ## Initial conditions
q1 = nominal_configuration(quadruped)
v1 = zeros(quadruped.nq)

# ## Time 
h = 0.01
T = 100

# ## Simulator
s = Simulator(quadruped, T, h=h)

# ## Simulate
simulate!(s, q1, v1)

# ## Visualizer
vis = Visualizer()
render(vis)

# ## Visualize
visualize!(vis, s)
