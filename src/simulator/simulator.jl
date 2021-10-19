struct Simulator{T}
    model 
    policy::Policy 
    dist::Disturbances 
    traj::Trajectory 
    ip::InteriorPoint 
    idx_z::IndicesZ
    idx_θ::Indicesθ
    f::Vector{T}
    h::T
end

function Simulator(model, T; 
        h=0.01, policy=empty_policy(model), dist=empty_disturbances(model), f=Float64[],
        residual=r!, jacobian_z=rz!, jacobian_θ=rθ!,
        opts=InteriorPointOptions(
            undercut=Inf,
            γ_reg=0.1,
            r_tol=1e-8,
            κ_tol=1e-8,  
            max_ls=25,
            ϵ_min=0.25,
            diff_sol=false,
            verbose=false))

    idx_z = indices_z(model)
    idx_θ = indices_θ(model, nf=length(f))
    idx_opt = indices_optimization(model)
    
    nz = num_var(model) 
    nθ = num_data(model, nf=length(f))
    z0 = zeros(nz) 
    θ0 = zeros(nθ) 
    q0 = configuration(model)
    initialize_z!(z0, idx_z, q0)
    initialize_θ!(θ0, idx_θ, q0, q0, zeros(model.nu), zeros(model.nw), f, h)

    ip = interior_point(
         z0,
         θ0,
         idx=idx_opt,
         r! = residual,
         rz! = jacobian_z,
         rθ! = jacobian_θ,
         rz=zeros(nz, nz),
         rθ=zeros(nz, nθ),
         opts=opts)

    traj = Trajectory(model, T)

    Simulator(model, policy, dist, traj, ip, idx_z, idx_θ, f, h)
end

function step!(traj::Trajectory, ip::InteriorPoint, p::Policy, w::Disturbances, f, h, idx_z, idx_θ, t)
    # policy 
    traj.u[t] = policy(p, traj, t)

    # disturbances 
    traj.w[t] = disturbances(w, traj.q[t+1], t)

    # initialize
    initialize_z!(ip.z, idx_z, traj.q[t+1])
    initialize_θ!(ip.θ, idx_θ, traj.q[t], traj.q[t+1], traj.u[t], traj.w[t], f, h)

    # solve
    status = interior_point_solve!(ip)

    # status check
    if !status 
        @show norm(ip.r) 
        @warn "residual failure"
        return false
    end

    # cache solution
    q = @views ip.z[idx_x.q]
    γ = @views ip.z[idx_x.γ] 
    b = @views ip.z[idx_x.b]
    traj.q[t+2] = q
    traj.v[t+1] = (traj.q[t+2] - traj.q[t+1]) ./ h
    traj.γ[t] = γ 
    traj.b[t] = b

    return true
end

function simulate!(traj::Trajectory, ip::InteriorPoint, p::Policy, w::Disturbances, f, h, idx_z, idx_θ, T)
    for t = 1:T
        status = step!(traj, ip, p, w, f, h, idx_z, idx_θ, t)
        !status && break
    end
end

function simulate!(s::Simulator, q0, q1)
    s.traj.q[1] = q0
    # s.traj.q[2] = q1
    # s.traj.v[1] = v
    # simulate!(s.traj, s.ip, s.p, s.w, s.f, s.h, s.idx_z, s.idx_θ, s.T)
end
