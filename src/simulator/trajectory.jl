struct Trajectory{nq,nv,nu,nγ,nb,T} 
    q::Vector{SVector{nq,T}} # generalized coordinates
    v::Vector{SVector{nv,T}} # generalized velocities
    u::Vector{SVector{nu,T}} # control inputs
    γ::Vector{SVector{nγ,T}} # impact impulses 
    b::Vector{SVector{nb,T}} # friction impulses
    w::Vector{SVector{nw,T}} # disturbances
end

function Trajectory(model, T; nv=model.nq, nb=model.nc) 
    nq = model.nq 
    nu = model.nu
    nw = model.nw 
    nc = model.nc
    q = [SVector{nq}(zeros(nq)) for t = 1:T+2]
    v = [SVector{nv}(zeros(nv)) for t = 1:T+1] 
    u = [SVector{nu}(zeros(nu)) for t = 1:T]
    γ = [SVector{nc}(zeros(nc)) for t = 1:T] 
    b = [SVector{nb}(zeros(nb)) for t = 1:T] 
    w = [SVector{nw}(zeros(nw)) for t = 1:T] 
    Trajectory(q, v, u, γ, b, w)
end