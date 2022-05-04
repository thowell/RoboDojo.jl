# PREAMBLE

# PKG_SETUP

# ## Setup

using RoboDojo

# ## Initial conditions
q1 = nominal_configuration(halfcheetah4)
v1 = zeros(halfcheetah4.nq)
q1[2] += 0.25
q1[3] += 0.0 * π
# ## Time
h = 0.01
T = 1000

# ## Simulator
s = Simulator(halfcheetah4, T, h=h)

# ## Simulate
simulate!(s, q1, v1)

# ## Visualizer
vis = Visualizer()
render(vis)

# ## Visualize
visualize!(vis, s)
