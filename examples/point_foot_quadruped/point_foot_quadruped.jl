# PREAMBLE

# PKG_SETUP

# ## Setup

using RoboDojo

# ## Initial conditions
q1 = nominal_configuration(point_foot_quadruped)
height = 0.5
q1[3] += height
q1[6 + 3] += height
q1[9 + 3] += height
q1[12 + 3] += height
q1[15 + 3] += height

v1 = zeros(point_foot_quadruped.nq)
# v1= randn(point_foot_quadruped.nq)
# ## Time
h = 0.01
T = 100

# ## Simulator
s = Simulator(point_foot_quadruped, T,
    h=h)

# ## Simulate
simulate!(s, q1, v1)

using BenchmarkTools
@benchmark simulate!(s, q1, v1) setup=(s=s, q1=q1, v1=v1)

# ## Visualizer
vis = Visualizer()
open(vis)

# ## Visualize
visualize!(vis, s,
    fixed_camera=false)
