include(joinpath("../src/utils.jl"))
include(joinpath("../src/simulator/model.jl"))
include(joinpath("../src/robots/biped/model.jl"))
include(joinpath("../src/robots/biped/visuals.jl"))
include(joinpath("../src/robots/visual_utils.jl"))
include(joinpath("../src/simulator/dimensions.jl"))
include(joinpath("../src/simulator/indices.jl"))
include(joinpath("../src/simulator/disturbances.jl"))
include(joinpath("../src/simulator/trajectory.jl"))
include(joinpath("../src/simulator/policy.jl"))
include(joinpath("../src/simulator/simulator.jl"))

model = biped 

q0 = [0.0; 1.0; 0.01 * π; -0.5 * π; 0.0 * π; -0.01 * π; 0.5 * π; 0.0 * π; -0.0 * π]

f = [biped.friction_body_world; biped.friction_foot_world]
h = 0.01

T = 100

s = Simulator(model, T, h=h, f=f, residual=r_biped!, jacobian_z=rz_biped!, jacobian_θ=rθ_biped!)
simulate!(s, q0, q0)
@benchmark simulate!($s, $q0, $q0)

vis = Visualizer()
open(vis)
visualize!(vis, model, 
     s.traj.q[1:1:end],
     Δt = 1 * h)

