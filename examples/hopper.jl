include(joinpath("../src/utils.jl"))
include(joinpath("../src/robots/hopper/model.jl"))
include(joinpath("../src/robots/hopper/visuals.jl"))
include(joinpath("../src/robots/visual_utils.jl"))

model = hopper 

space = Euclidean(nz)

ip_opts = InteriorPointOptions(
    undercut = Inf,
    γ_reg = 0.1,
    r_tol = 1e-8,
    κ_tol = 1e-8,  
    max_ls = 25,
    ϵ_min = 0.25,
    diff_sol = true,
    verbose = false
    )

q0 = [0.0, 1.0, 0.0 * π, 0.5]
q1 = copy(q0)
u = zeros(2) 
w = zeros(nw)
h = 0.01

idx_x = indices_z(hopper)
idx_θ = indices_θ(hopper)
idx_opt = indices_optimization(hopper)

z0 = zeros(num_var(hopper)) 
θ0 = zeros(num_data(hopper)) 
initialize_z!(z0, idx_x, q1)
initialize_θ!(θ0, idx_θ, q0, q1, u, w, f1, f2, h)

ip = interior_point(
         z0,
         θ0,
         s = space,
         idx = idx_opt,
         r! = r_hopper!,
         rz! = rz_hopper!,
         rθ! = rθ_hopper!,
         rz = zeros(nz, nz),
         rθ = zeros(nz, nθ),
         opts = ip_opts)

T = 100
q_hist = [zeros(nq) for t = 1:(T+2)]
q_hist[1] = q0
q_hist[2] = q1

function simulate!(q_hist, ip, u, w, f1, f2, h, idx_x, idx_θ, T)
    for t = 1:T
        # initialize
        initialize_z!(ip.z, idx_x, q_hist[t+1])
        initialize_θ!(ip.θ, idx_θ, q_hist[t], q_hist[t+1], u, w, f1, f2, h)

        # solve
        status = interior_point_solve!(ip)

        # status check
        if !status 
            @show norm(ip.r) 
            @warn "residual failure"
            break
        end

        # cache solution
        q = @views ip.z[idx_x.q]
        q_hist[t+2] .= q
    end
end

simulate!(q_hist, ip, u, w, f1, f2, h, idx_x, idx_θ, T)

@benchmark simulate!($q_hist, $ip, $u, $w, $f1, $f2, $h, $idx_x, $idx_θ, $T)

vis = Visualizer()
open(vis)
visualize!(vis, hopper, q_hist[1:1:end], Δt = 1*h)

settransform!(vis["/Cameras/default"],
        compose(Translation(0.0, -50.0, -1.0),LinearMap(RotZ(-pi / 2.0))))

setprop!(vis["/Cameras/default/rotated/<object>"], "zoom", 25)
