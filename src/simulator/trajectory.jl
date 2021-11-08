struct Trajectory{T} 
    q::Vector{Vector{T}} # generalized coordinates
    v::Vector{Vector{T}} # generalized velocities
    u::Vector{Vector{T}} # control inputs
    γ::Vector{Vector{T}} # impact impulses 
    b::Vector{Vector{T}} # friction impulses
    w::Vector{Vector{T}} # disturbances
end

function Trajectory(model, T; nv=model.nq, nc=model.nc, nb=model.nc) 
    nq = model.nq 
    nu = model.nu
    nw = model.nw 

    q = [zeros(nq) for t = 1:T+2]
    v = [zeros(nv) for t = 1:T+1] # midpoint velocities; TODO: knotpoint velocities
    u = [zeros(nu) for t = 1:T]
    γ = [zeros(nc) for t = 1:T] 
    b = [zeros(nb) for t = 1:T] 
    w = [zeros(nw) for t = 1:T] 
    Trajectory(q, v, u, γ, b, w)
end

function reset!(traj::Trajectory) 
    T = length(traj.u) 
    for t = 1:T
        fill!(traj.q[t], 0.0) 
        fill!(traj.v[t], 0.0) 
        fill!(traj.u[t], 0.0) 
        fill!(traj.γ[t], 0.0) 
        fill!(traj.b[t], 0.0) 
        fill!(traj.w[t], 0.0) 
    end
    fill!(traj.q[T+1], 0.0) 
    fill!(traj.q[T+2], 0.0) 
    fill!(traj.v[T+1], 0.0) 
    return nothing
end

struct GradientTrajectory{T} 
    ∂q3∂q1::Vector{Matrix{T}}
    ∂q3∂q2::Vector{Matrix{T}}
    ∂q3∂u1::Vector{Matrix{T}}
    ∂γ1∂q1::Vector{Matrix{T}} 
    ∂γ1∂q2::Vector{Matrix{T}}
    ∂γ1∂u1::Vector{Matrix{T}}
    ∂b1∂q1::Vector{Matrix{T}}
    ∂b1∂q2::Vector{Matrix{T}}
    ∂b1∂u1::Vector{Matrix{T}}
end

function GradientTrajectory(model, T; nc=model.nc, nb=model.nc) 
    nq = model.nq 
    nu = model.nu

    ∂q3∂q1 = [zeros(nq, nq) for t = 1:T]
    ∂q3∂q2 = [zeros(nq, nq) for t = 1:T]
    ∂q3∂u1 = [zeros(nq, nu) for t = 1:T]
    ∂γ1∂q1 = [zeros(nc, nq) for t = 1:T]
    ∂γ1∂q2 = [zeros(nc, nq) for t = 1:T]
    ∂γ1∂u1 = [zeros(nc, nu) for t = 1:T]
    ∂b1∂q1 = [zeros(nb, nq) for t = 1:T]
    ∂b1∂q2 = [zeros(nb, nq) for t = 1:T]
    ∂b1∂u1 = [zeros(nb, nu) for t = 1:T]

    GradientTrajectory(∂q3∂q1, ∂q3∂q2, ∂q3∂u1, 
                       ∂γ1∂q1, ∂γ1∂q2, ∂γ1∂u1,
                       ∂b1∂q1, ∂b1∂q2, ∂b1∂u1)
end

function reset!(traj::GradientTrajectory) 
    T = length(traj.∂q3∂q1) 
    for t = 1:T
        fill!(traj.∂q3∂q1[t], 0.0)
        fill!(traj.∂q3∂q2[t], 0.0)
        fill!(traj.∂q3∂u1[t], 0.0)
        fill!(traj.∂γ1∂q1[t], 0.0)
        fill!(traj.∂γ1∂q2[t], 0.0)
        fill!(traj.∂γ1∂u1[t], 0.0)
        fill!(traj.∂b1∂q1[t], 0.0)
        fill!(traj.∂b1∂q2[t], 0.0)
        fill!(traj.∂b1∂u1[t], 0.0)
    end
    return nothing
end