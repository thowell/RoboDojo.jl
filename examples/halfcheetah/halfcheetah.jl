# PREAMBLE

# PKG_SETUP

# ## Setup

using RoboDojo

# ## Initial conditions
q1 = nominal_configuration(halfcheetah)
v1 = zeros(halfcheetah.nq)

# ## Time
h = 0.01
T = 100

# ## Simulator
s = Simulator(halfcheetah, T, h=h)

# ## Simulate
simulate!(s, q1, v1)

# ## Visualizer
vis = Visualizer()
render(vis)

# ## Visualize
visualize!(vis, s)
