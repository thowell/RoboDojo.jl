
@with_kw struct SimulatorOptions{T}
    warmstart::Bool = true
    z_warmstart::T = 0.001
    κ_warmstart::T = 0.001
	failure_abort::Int = 50
end

struct SimulatorStatistics{T}
    time::Vector{T}
    mean::Vector{T}
    std::Vector{T}
end

function SimulatorStatistics()
    time = zeros(0)
    mean = zeros(0)
    std = zeros(0)
    return SimulatorStatistics(time, mean, std)
end

function process!(stats::SimulatorStatistics, N_sample::Int)
    H_sim = length(stats.time)
    H = Int(H_sim / N_sample)
    time = sum(reshape(stats.time, (H, N_sample)), dims=2)
    stats.mean[1] = mean(time)
    stats.std[1] = sqrt(mean( (time .- stats.mean[1]).^2 ))
    return nothing
end

struct Simulator{T,R,RZ,Rθ,M<:Model{T},P<:Policy{T},D<:Disturbances{T}}
    model::M
    policy::P
    dist::D
    traj::Trajectory{T}
    grad::GradientTrajectory{T}
    ip::InteriorPoint{T,R,RZ,Rθ}
    idx_z::IndicesZ
    idx_θ::Indicesθ
    f::Vector{T}
    h::T
    stats::SimulatorStatistics{T}
    opts::SimulatorOptions{T}
end

function Simulator(model, T; 
        h=0.01, 
        policy=empty_policy(model), 
        dist=empty_disturbances(model), 
        f=friction_coefficients(model),
        residual=eval(residual_expr(model)), 
        jacobian_z=eval(jacobian_var_expr(model)), 
        jacobian_θ=eval(jacobian_data_expr(model)),
        diff_sol=false,
        solver_opts=InteriorPointOptions(
            undercut=Inf,
            γ_reg=0.1,
            r_tol=1e-8,
            κ_tol=1e-8,  
            max_ls=25,
            ϵ_min=0.25,
            diff_sol=diff_sol,
            verbose=false),
        stats=SimulatorStatistics(),
        sim_opts=SimulatorOptions()
        )

    idx_z = indices_z(model)
    idx_θ = indices_θ(model, nf=length(f))
    idx_opt = indices_optimization(model)
    
    nz = num_var(model) 
    nθ = num_data(model, nf=length(f))
    z0 = zeros(nz) 
    θ0 = zeros(nθ) 
    q0 = nominal_configuration(model)
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
         opts=solver_opts)

    traj = Trajectory(model, T)
    grad = GradientTrajectory(model, T)

    Simulator(model, policy, dist, traj, grad, ip, idx_z, idx_θ, f, h, stats, sim_opts)
end

function step!(traj::Trajectory{T}, grad::GradientTrajectory{T}, ip::InteriorPoint{T}, p::Policy{T}, w::Disturbances{T}, f::Vector{T}, h::T, idx_z::IndicesZ, idx_θ::Indicesθ, t::Int) where T
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
    solution!(traj, ip.z, idx_z, h, t)

    # cache gradients
    ip.opts.diff_sol && gradient!(grad, ip.δz, idx_z, idx_θ, t)

    return true
end

function solution!(traj::Trajectory{T}, z::Vector{T}, idx_z::IndicesZ, h::T, t::Int) where T 
    q = @views z[idx_z.q]
    γ = @views z[idx_z.γ] 
    b = @views z[idx_z.b]

    # current configuration
    traj.q[t+2] .= q

    # current velocity
    traj.v[t+1] .= traj.q[t+2]
    traj.v[t+1] .-= traj.q[t+1] 
    traj.v[t+1] ./= h

    # contact impulses
    traj.γ[t] .= γ 
    traj.b[t] .= b
    return nothing
end

function gradient!(grad::GradientTrajectory{T}, δz::Matrix{T}, idx_z::IndicesZ, idx_θ::Indicesθ, t::Int) where T
    ∂q3∂q1 = @views δz[idx_z.q, idx_θ.q1]
    ∂q3∂q2 = @views δz[idx_z.q, idx_θ.q2]
    ∂q3∂u1 = @views δz[idx_z.q, idx_θ.u]
    ∂γ1∂q1 = @views δz[idx_z.γ, idx_θ.q1]
    ∂γ1∂q2 = @views δz[idx_z.γ, idx_θ.q2]
    ∂γ1∂u1 = @views δz[idx_z.γ, idx_θ.u]
    ∂b1∂q1 = @views δz[idx_z.b, idx_θ.q1]
    ∂b1∂q2 = @views δz[idx_z.b, idx_θ.q2]
    ∂b1∂u1 = @views δz[idx_z.b, idx_θ.u]
    
    grad.∂q3∂q1[t] .= ∂q3∂q1
    grad.∂q3∂q2[t] .= ∂q3∂q2
    grad.∂q3∂u1[t] .= ∂q3∂u1
    grad.∂γ1∂q1[t] .= ∂γ1∂q1
    grad.∂γ1∂q2[t] .= ∂γ1∂q2
    grad.∂γ1∂u1[t] .= ∂γ1∂u1
    grad.∂b1∂q1[t] .= ∂b1∂q1
    grad.∂b1∂q2[t] .= ∂b1∂q2
    grad.∂b1∂u1[t] .= ∂b1∂u1

    return nothing
end

function simulate!(traj::Trajectory{T}, grad::GradientTrajectory{T}, ip::InteriorPoint{T}, 
        p::Policy{T}, w::Disturbances{T}, f::Vector{T}, h::T, idx_z::IndicesZ, idx_θ::Indicesθ, N::Int) where T
    for t = 1:N
        status = step!(traj, grad, ip, p, w, f, h, idx_z, idx_θ, t)
        !status && break
    end
end

function simulate!(s::Simulator{T}, q::Vector{T}, v::Vector{T}; reset_traj=false) where T
    # reset trajectory
    reset_traj && reset!(s.traj) 
    reset_traj && reset!(s.grad)

    # reset solver
    

    # initial configuration and velocity
    s.traj.q[2] .= q # q2
    s.traj.v[1] .= v # v1

    #   q1
    s.traj.q[1] .= v 
    s.traj.q[1] .*= s.h
    s.traj.q[1] .*= -1.0
    s.traj.q[1] .+= q

    # simulation length
    N = length(s.traj.u)

    # simulate
    simulate!(s.traj, s.grad, s.ip, s.policy, s.dist, s.f, s.h, s.idx_z, s.idx_θ, N)
end

function visualize!(vis, s::Simulator; skip=1, fixed_camera=true)
    visualize!(vis, s.model, s.traj.q[1:skip:end], Δt=skip * s.h, fixed_camera=fixed_camera)
end
