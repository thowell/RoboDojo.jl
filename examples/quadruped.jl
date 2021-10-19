include(joinpath("../src/utils.jl"))
include(joinpath("../src/robots/quadruped/model.jl"))
include(joinpath("../src/robots/quadruped/visuals.jl"))
include(joinpath("../src/robots/visual_utils.jl"))
include(joinpath("../src/simulator/dimensions.jl"))
include(joinpath("../src/simulator/indices.jl"))

model = quadruped 

ip_opts = InteriorPointOptions(
    undercut = Inf,
    γ_reg = 0.1,
    r_tol = 1e-8,
    κ_tol = 1e-8,  
    max_ls = 25,
    ϵ_min = 0.25,
    diff_sol = true,
    verbose = false)

q0 = [0.0; 0.5; 0.0 * π; 0.25 * π; 0.5 * π; 0.1 * π; 0.3 * π; -0.25 * π; 0.1 * π; -0.5 * π; -0.1 * π]
q1 = copy(q0)
u = zeros(nu) 
w = zeros(nw)
h = 0.01
f = [model.friction_body_world; model.friction_foot_world]

idx_z = indices_z(quadruped)
idx_θ = indices_θ(quadruped, nf=length(f))
idx_opt = indices_optimization(quadruped)

nz = num_var(quadruped) 
nθ = num_data(quadruped, nf=length(f))
z0 = zeros(nz) 
θ0 = zeros(nθ) 
initialize_z!(z0, idx_z, q1)
initialize_θ!(θ0, idx_θ, q0, q1, u, w, f, h)

ip = interior_point(
         z0,
         θ0,
         idx = idx_opt,
         r! = r_quadruped!,
         rz! = rz_quadruped!,
         rθ! = rθ_quadruped!,
         rz = zeros(nz, nz),
         rθ = zeros(nz, nθ),
         opts = ip_opts)

interior_point_solve!(ip) 

T = 100
traj = Trajectory(model, T)
traj.q[1] = q0 
traj.q[2] = q1
traj.v[1] = (q1 - q0) ./ h

p = empty_policy(model)
w = empty_disturbances(model) 


simulate!(traj, ip, p, w, f, h, idx_x, idx_θ, T)
@benchmark simulate!($traj, $ip, $p, $w, $f, $h, $idx_x, $idx_θ, $T)

vis = Visualizer()
open(vis)
visualize!(vis, model, 
     traj.q[1:1:end],
     Δt = 1 * h)

settransform!(vis["/Cameras/default"],
        compose(Translation(0.0, -50.0, -1.0),LinearMap(RotZ(-pi / 2.0))))

setprop!(vis["/Cameras/default/rotated/<object>"], "zoom", 25)
