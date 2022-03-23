# PREAMBLE

# PKG_SETUP

# ## Setup

using RoboDojo 

# ## Initial conditions
q1 = nominal_configuration(centroidal_quadruped)
v1 = zeros(centroidal_quadruped.nq)

# ## Time 
h = 0.1
T = 100

# ## Simulator
s = Simulator(centroidal_quadruped, T, 
    h=h)

# ## Simulate
simulate!(s, q1, v1, 
    verbose=true)

# ## Visualizer
vis = Visualizer()
render(vis)

# ## Visualize
visualize!(vis, s, 
    fixed_camera=false)
