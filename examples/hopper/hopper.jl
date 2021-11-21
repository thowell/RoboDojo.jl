# PREAMBLE

# PKG_SETUP

# ## Setup

using RoboDojo 

# ## Initial conditions
q1 = nominal_configuration(hopper)
v1 = zeros(hopper.nq)

# ## Time 
T = 100
h = 0.01

# ## Simulator 
s = Simulator(hopper, T, h=h)

# ## Simulate
simulate!(s, q1, v1)

# ## Visualizer
vis = Visualizer()
render(vis)

# ## Visualize
visualize!(vis, s)