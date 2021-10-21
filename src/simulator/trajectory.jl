struct Trajectory{T} 
    q::Vector{Vector{T}} # generalized coordinates
    v::Vector{Vector{T}} # generalized velocities
    u::Vector{Vector{T}} # control inputs
    γ::Vector{Vector{T}} # impact impulses 
    b::Vector{Vector{T}} # friction impulses
    w::Vector{Vector{T}} # disturbances
end

function Trajectory(model, T; nv=model.nq, nb=model.nc) 
    nq = model.nq 
    nu = model.nu
    nw = model.nw 
    nc = model.nc
    q = [zeros(nq) for t = 1:T+2]
    v = [zeros(nv) for t = 1:T+1] 
    u = [zeros(nu) for t = 1:T]
    γ = [zeros(nc) for t = 1:T] 
    b = [zeros(nb) for t = 1:T] 
    w = [zeros(nw) for t = 1:T] 
    Trajectory(q, v, u, γ, b, w)
end