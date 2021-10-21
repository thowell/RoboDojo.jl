
struct Simulator{T,R,RZ,Rθ,M<:Model{T},P<:Policy{T},D<:Disturbances{T}}
    model::M
    policy::P
    dist::D
    traj::Trajectory{T}
    ip::InteriorPoint{T,R,RZ,Rθ}
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

function step!(traj::Trajectory{T}, ip::InteriorPoint{T}, p::Policy{T}, w::Disturbances{T}, f::Vector{T}, h::T, idx_z::IndicesZ, idx_θ::Indicesθ, t::Int) where T
    # policy 
    traj.u[t] .= policy(p, traj, t)

    # disturbances 
    traj.w[t] .= disturbances(w, traj.q[t+1], t)

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
    q = @views ip.z[idx_z.q]
    γ = @views ip.z[idx_z.γ] 
    b = @views ip.z[idx_z.b]
    traj.q[t+2] .= q
    # traj.v[t+1] .= (traj.q[t+2] - traj.q[t+1]) ./ h #TODO: make allocation free
    traj.v[t+1] .= traj.q[t+2]
    traj.v[t+1] .-= traj.q[t+1] 
    traj.v[t+1] ./= h
    traj.γ[t] .= γ 
    traj.b[t] .= b

    return true
end

function simulate!(traj::Trajectory{T}, ip::InteriorPoint{T}, 
        p::Policy{T}, w::Disturbances{T}, f::Vector{T}, h::T, idx_z::IndicesZ, idx_θ::Indicesθ, N::Int) where T
    for t = 1:N
        status = step!(traj, ip, p, w, f, h, idx_z, idx_θ, t)
        !status && break
    end
end

function simulate!(s::Simulator{T}, q::Vector{T}, v::Vector{T}) where T
    s.traj.q[2] .= q
    s.traj.v[1] .= v
    # s.traj.q[1] .= s.traj.q[2] - v * h
    s.traj.q[1] .= v 
    s.traj.q[1] .*= s.h
    s.traj.q[1] .*= -1.0
    s.traj.q[1] .+= q
    N = length(s.traj.u)
    simulate!(s.traj, s.ip, s.policy, s.dist, s.f, s.h, s.idx_z, s.idx_θ, N)
end

function visualize!(vis, s::Simulator; skip=1, fixed_camera=true)
    visualize!(vis, s.model, s.traj.q[1:skip:end], Δt=skip * s.h, fixed_camera=fixed_camera)
end
