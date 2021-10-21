include(joinpath("../src/utils.jl"))
include(joinpath("../src/simulator/model.jl"))
include(joinpath("../src/robots/quadruped/model.jl"))
include(joinpath("../src/robots/quadruped/visuals.jl"))
include(joinpath("../src/robots/visual_utils.jl"))
include(joinpath("../src/robots/visual_utils.jl"))
include(joinpath("../src/simulator/dimensions.jl"))
include(joinpath("../src/simulator/indices.jl"))
include(joinpath("../src/simulator/disturbances.jl"))
include(joinpath("../src/simulator/trajectory.jl"))
include(joinpath("../src/simulator/policy.jl"))
include(joinpath("../src/simulator/simulator.jl"))

model = quadruped 
q0 = configuration(model)
v0 = zeros(model.nq)
h = 0.01
f = [model.friction_body_world; model.friction_foot_world]
T = 100
s = Simulator(model, T, h=h, f=f, residual=r_quadruped!, jacobian_z=rz_quadruped!, jacobian_θ=rθ_quadruped!)
simulate!(s, q0, v0)
@benchmark simulate!($s, $q0, $v0)

vis = Visualizer()
open(vis)
visualize!(vis, model, 
     s.traj.q[1:1:end],
     Δt=1 * h)

settransform!(vis["/Cameras/default"],
        compose(Translation(0.0, -50.0, -1.0),LinearMap(RotZ(-pi / 2.0))))

setprop!(vis["/Cameras/default/rotated/<object>"], "zoom", 25)
