# PREAMBLE

# PKG_SETUP

# ## Setup

using RoboDojo 

# ## Initial conditions
q1 = nominal_configuration(hopper1)
v1 = zeros(hopper1.nq)

# ## Time 
T = 200
h = 0.01

# ## Simulator 
s = Simulator(hopper1, T, h=h)

# ## Simulate
simulate!(s, q1, v1)

# ## Visualizer
vis = Visualizer()
render(vis)

# ## Visualize
visualize!(vis, s)