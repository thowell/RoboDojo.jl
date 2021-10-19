include(joinpath("../src/utils.jl"))
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

# ip_opts = InteriorPointOptions(
#     undercut = Inf,
#     γ_reg = 0.1,
#     r_tol = 1e-8,
#     κ_tol = 1e-8,  
#     max_ls = 25,
#     ϵ_min = 0.25,
#     diff_sol = false,
#     verbose = false)


q1 = copy(q0)
u = zeros(biped.nu) 
w = zeros(biped.nw)
f = [biped.friction_body_world; biped.friction_foot_world]

idx_x = indices_z(biped)
# idx_θ = indices_θ(biped, nf=length(f))
# idx_opt = indices_optimization(biped)

# nz = num_var(biped) 
# nθ = num_data(biped, nf=length(f))
# z0 = zeros(nz) 
# θ0 = zeros(nθ) 
# initialize_z!(z0, idx_x, q1)
# initialize_θ!(θ0, idx_θ, q0, q1, u, w, f, h)

# ip = interior_point(
#          z0,
#          θ0,
#          idx = idx_opt,
#          r! = r_biped!,
#          rz! = rz_biped!,
#          rθ! = rθ_biped!,
#          rz = zeros(nz, nz),
#          rθ = zeros(nz, nθ),
#          opts = ip_opts)

# interior_point_solve!(ip) 
T = 100
h = 0.01 

v0 = SVector{model.nq}(zeros(model.nq))
q0 = SVector{model.nq}([0.0; 1.0; 0.01 * π; -0.5 * π; 0.0 * π; -0.01 * π; 0.5 * π; 0.0 * π; -0.0 * π])

s = Simulator(model, T, h=h, residual=r_biped!, jacobian_z=rz_biped!, jacobian_θ=rθ_biped!)
@benchmark simulate!($s, $q0, $q0)

vis = Visualizer()
open(vis)
visualize!(vis, model, 
     traj.q[1:1:end],
     Δt = 1 * h)

settransform!(vis["/Cameras/default"],
        compose(Translation(0.0, -50.0, -1.0),LinearMap(RotZ(-pi / 2.0))))

setprop!(vis["/Cameras/default/rotated/<object>"], "zoom", 25)
