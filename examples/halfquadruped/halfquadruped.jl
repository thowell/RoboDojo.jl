# PREAMBLE

# PKG_SETUP

# ## Setup

using RoboDojo

# ## Initial conditions
q1 = nominal_configuration(halfquadruped)
v1 = zeros(halfquadruped.nq)

# ## Time
h = 0.01
T = 100

# ## Simulator
s = Simulator(halfquadruped, T, h=h)

# ## Simulate
simulate!(s, q1, v1)

# ## Visualizer
vis = Visualizer()
render(vis)

# ## Visualize
visualize!(vis, s)
