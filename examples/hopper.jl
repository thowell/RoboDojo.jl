include(joinpath("../src/utils.jl"))
include(joinpath("../src/simulator/model.jl"))
include(joinpath("../src/robots/hopper/model.jl"))
include(joinpath("../src/robots/hopper/simulator.jl"))
include(joinpath("../src/robots/hopper/visuals.jl"))
include(joinpath("../src/robots/visual_utils.jl"))
include(joinpath("../src/simulator/dimensions.jl"))
include(joinpath("../src/simulator/indices.jl"))
include(joinpath("../src/simulator/disturbances.jl"))
include(joinpath("../src/simulator/trajectory.jl"))
include(joinpath("../src/simulator/policy.jl"))
include(joinpath("../src/simulator/simulator.jl"))

model = hopper 

q0 = configuration(model)
v0 = zeros(model.nq)
h = 0.01

T = 100

f = [hopper.friction_body_world; hopper.friction_foot_world]
h = 0.01

T = 100

s = Simulator(model, T, h=h, f=f, residual=r_hopper!, jacobian_z=rz_hopper!, jacobian_θ=rθ_hopper!)
simulate!(s, q0, v0)
@benchmark simulate!($s, $q0, $v0)

vis = Visualizer()
open(vis)
visualize!(vis, s)
