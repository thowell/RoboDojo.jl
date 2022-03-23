# PREAMBLE

# PKG_SETUP

# ## Setup

using RoboDojo 

# ## Initial conditions
q1 = nominal_configuration(particle)
q1[3] = 0.25
v1 = zeros(particle.nq)
v1[1] = 10.0

# ## Time 
T = 100
h = 0.01

# ## Simulator 
s = Simulator(particle, T, h=h)

# ## Simulate
simulate!(s, q1, v1)

# ## Visualizer
vis = Visualizer()
render(vis)

# ## Visualize
visualize!(vis, s, 
    fixed_camera=false)