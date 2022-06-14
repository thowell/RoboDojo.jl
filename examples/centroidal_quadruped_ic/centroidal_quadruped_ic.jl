# PREAMBLE

# PKG_SETUP

# ## Setup

using RoboDojo 

# ## Initial conditions
q1 = nominal_configuration(centroidal_quadruped_ic) 
height = 0.5 
q1[3] += height 
q1[6 + 3] += height 
q1[9 + 3] += height 
q1[12 + 3] += height 
q1[15 + 3] += height 

v1 = zeros(centroidal_quadruped_ic.nq)
# v1 = randn(centroidal_quadruped_ic.nq)

# ## Time 
h = 0.01
T = 100

# ## Simulator
s = Simulator(centroidal_quadruped_ic, T, 
    h=h)

# ## Simulate
simulate!(s, q1, v1)

using BenchmarkTools 
@benchmark simulate!(s, q1, v1) setup=(s=s, q1=q1, v1=v1)

# ## Visualizer
vis = Visualizer()
render(vis)

# ## Visualize
visualize!(vis, s, 
    fixed_camera=false)
