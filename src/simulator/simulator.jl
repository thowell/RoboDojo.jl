@with_kw struct SimulatorOptions{T}
    warmstart::Bool=false
    record::Bool=true
    z_warmstart::T=0.001
    κ_warmstart::T=0.001
	failure_abort::Int=50
end

struct SimulatorStatistics{T}
    policy_time::Vector{T}
    policy_mean::Vector{T}
    policy_std::Vector{T}
    sim_time::Vector{T}
    sim_mean::Vector{T}
    sim_std::Vector{T}
end

function SimulatorStatistics(H)
    policy_time = zeros(H)
    policy_mean = zeros(1)
    policy_std = zeros(1)
    sim_time = zeros(H)
    sim_mean = zeros(1)
    sim_std = zeros(1)
    return SimulatorStatistics(policy_time, policy_mean, policy_std, 
        sim_time, sim_mean, sim_std)
end

function process!(stats::SimulatorStatistics, N_sample::Int)
    # policy
    H_sim = length(stats.policy_time)
    H = Int(H_sim / N_sample)
    policy_time = sum(reshape(stats.policy_time, (H, N_sample)), dims=2)
    stats.policy_mean[1] = mean(policy_time)
    stats.policy_std[1] = sqrt(mean((policy_time .- stats.policy_mean[1]).^2.0))
    
    # simulation
    H_sim = length(stats.sim_time)
    H = Int(H_sim / N_sample)
    sim_time = sum(reshape(stats.sim_time, (H, N_sample)), dims=2)
    stats.sim_mean[1] = mean(sim_time)
    stats.sim_std[1] = sqrt(mean((sim_time .- stats.sim_mean[1]).^2.0))
    
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
        stats=SimulatorStatistics(T),
        sim_opts=SimulatorOptions())
    
    idx_z = indices_z(model)
    idx_θ = indices_θ(model, nf=length(f))
    idx_opt = indices_optimization(model)
    
    nz = num_var(model) 
    nθ = num_data(model, nf=length(f))
    z0 = zeros(nz) 
    θ0 = zeros(nθ) 
    q0 = nominal_configuration(model)
    initialize_z!(z0, model, idx_z, q0)
    initialize_θ!(θ0, model, idx_θ, q0, q0, zeros(model.nu), zeros(model.nw), f, h)

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

function step!(s::Simulator{T}, t::Int) where T
    model = s.model
    traj = s.traj
    grad = s.grad
    ip = s.ip
    idx_z = s.idx_z
    idx_θ = s.idx_θ
    f = s.f
    h = s.h

    # initialize
    initialize_z!(ip.z, model, idx_z, traj.q[t+1])
    initialize_θ!(ip.θ, model, idx_θ, traj.q[t], traj.q[t+1], traj.u[t], traj.w[t], f, h)

    # solve
    sim_time = @elapsed status = interior_point_solve!(ip)
    s.opts.record && (s.stats.sim_time[t] = sim_time)

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

    return status
end

function step!(s::Simulator{T}, q::Vector{T}, v::Vector{T}, u::Vector{T}, t::Int) where T
    # set state
    set_state!(s, q, v, t)

    # set control 
    s.traj.u[t] .= u 

    # step 
    step!(s, t) 

    return s.traj.q[t+2] #TODO: return gradients
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

    #TODO: velocity gradients

    return nothing
end

function set_state!(s::Simulator{T}, q::Vector{T}, v::Vector{T}, t::Int) where T 
    # initial configuration and velocity
    s.traj.q[t+1] .= q # q2
    s.traj.v[t] .= v   # v1

    # q1
    s.traj.q[t] .= v 
    s.traj.q[t] .*= s.h
    s.traj.q[t] .*= -1.0
    s.traj.q[t] .+= q
    
    return nothing 
end

function simulate!(s::Simulator{T}, q::Vector{T}, v::Vector{T}; reset_traj=false) where T
    # reset trajectory
    reset_traj && reset!(s.traj) 
    reset_traj && reset!(s.grad)

    # reset solver
    
    # set initial state
    set_state!(s, q, v, 1)

    # simulate
    eval_simulate!(s)
end

function eval_simulate!(s::Simulator{T}) where T
    status = false

    N = length(s.traj.u)
    p = s.policy
    w = s.dist
    traj = s.traj

    for t = 1:N
        # policy 
        policy_time = @elapsed traj.u[t] .= policy(p, traj, t)
        s.opts.record && (s.stats.policy_time[t] = policy_time)

        # disturbances 
        traj.w[t] .= disturbances(w, traj.q[t+1], t)

        # step
        status = step!(s, t)
        !status && break
    end
    return status
end

function visualize!(vis, s::Simulator; skip=1, fixed_camera=true)
    visualize!(vis, s.model, s.traj.q[1:skip:end], Δt=skip * s.h, fixed_camera=fixed_camera)
end
