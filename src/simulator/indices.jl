struct IndicesZ 
    q::Vector{Int}
    γ::Vector{Int}  
    sγ::Vector{Int}
    ψ::Vector{Int}
    b::Vector{Int}
    sψ::Vector{Int}
    sb::Vector{Int}  
end

function indices_z(model) 
    nq = model.nq 
    nc = model.nc
    q = collect(1:nq) 
    γ = collect(nq .+ (1:nc)) 
    sγ = collect(nq + nc .+ (1:nc))
    ψ = collect(nq + nc + nc .+ (1:nc)) 
    b = collect(nq + nc + nc + nc .+ (1:nc)) 
    sψ = collect(nq + nc + nc + nc + nc .+ (1:nc)) 
    sb = collect(nq + nc + nc + nc + nc + nc .+ (1:nc)) 
    IndicesZ(q, γ, sγ, ψ, b, sψ, sb)
end

function initialize_z!(z, idx::IndicesZ, q)
    z[idx.q] .= q
    z[idx.γ] .= 1.0
    z[idx.sγ] .= 1.0
    z[idx.ψ] .= 1.0
    z[idx.b] .= 0.1
    z[idx.sψ] .= 1.0
    z[idx.sb] .= 0.1
end

# @benchmark initialize_z!($z0, $idx_x, $q1)

struct Indicesθ 
    q1::Vector{Int}
    q2::Vector{Int} 
    u::Vector{Int} 
    w::Vector{Int} 
    f::Vector{Int} 
    h::Vector{Int} 
end

function indices_θ(model; nf=1) 
    nq = model.nq 
    nu = model.nu 
    nw = model.nw 

    q1 = collect(1:nq)
    q2 = collect(nq .+ (1:nq))
    u = collect(2nq .+ (1:nu))
    w = collect(2nq + nu .+ (1:nw))
    f = collect(2nq + nu + nw .+ (1:nf))
    h = collect(2nq + nu + nw + nf .+ (1:1))

    Indicesθ(q1, q2, u, w, f, h) 
end

function initialize_θ!(θ, idx, q1, q2, u, w, f, h) 
    θ[idx.q1] .= q1
    θ[idx.q2] .= q2 
    θ[idx.u] .= u
    θ[idx.w] .= w
    θ[idx.f] .= f
    θ[idx.h] .= h
end

# @benchmark initialize_θ!($θ0, $idx_θ, $q1, $q2, $u0, $w0, $model.friction_body_world, $model.friction_foot_world, $h) 

function indices_optimization(model) 
    nq = model.nq
    nz = num_var(model) 
    ny = nz - nq 
    nc = model.nc
    IndicesOptimization(
        nz, 
        nz, 
        ny,
        collect(1:nq),
        collect(1:nq), 
        [collect(nq .+ (1:nc)), collect(nq + nc .+ (1:nc))],
        [collect(nq .+ (1:nc)), collect(nq + nc .+ (1:nc))],
        [
        [collect([nq + nc + nc + i, nq + nc + nc + nc + i]) for i = 1:nc], 
        [collect([nq + nc + nc + nc + nc + i, nq + nc + nc + nc + nc + nc + i]) for i = 1:nc]
        ],
        [
        [collect([nq + nc + nc + i, nq + nc + nc + nc + i]) for i = 1:nc], 
        [collect([nq + nc + nc + nc + nc + i, nq + nc + nc + nc + nc + nc + i]) for i = 1:nc]
        ],
        collect(1:(nq + nc + nc + nc)),
        collect(nq + nc + nc + nc .+ (1:nc)),
        collect(nq + nc + nc + nc + nc .+ (1:(2 * nc))),
        [collect(nq + nc + nc + nc + nc + (i - 1) * 2 .+ (1:2)) for i = 1:nc],
        collect(1:nq),
        collect(nq .+ (1:(nc + nc + nc))),
        collect(nq + nc + nc + nc .+ (1:(3 * nc))))
end