function num_var(model::Hopper) 
    nq = model.nq
    nq + 4 + 4 + 2 + 2 + 2 + 2 
end 

function indices_z(model::Hopper) 
    nq = model.nq 
    q = collect(1:nq) 
    γ = collect(nq .+ (1:4)) 
    sγ = collect(nq + 4 .+ (1:4))
    ψ = collect(nq + 4 + 4 .+ (1:2)) 
    b = collect(nq + 4 + 4 + 2 .+ (1:2)) 
    sψ = collect(nq + 4 + 4 + 2 + 2 .+ (1:2)) 
    sb = collect(nq + 4 + 4 + 2 + 2 + 2 .+ (1:2))   
    IndicesZ(q, γ, sγ, ψ, b, sψ, sb)
end

function indices_optimization(model::Hopper) 
    nz = num_var(model) 
    ny = nz - model.nq 

    IndicesOptimization(
        nz, 
        nz, 
        ny,
        collect(1:nq),
        collect(1:nq), 
        [collect(nq .+ (1:4)), collect(nq + 4 .+ (1:4))],
        [collect(nq .+ (1:4)), collect(nq + 4 .+ (1:4))],
        [[collect([nq + 4 + 4 + 1, nq + 4 + 4 + 2 + 1]), collect([nq + 4 + 4 + 2, nq + 4 + 4 + 2 + 2])], [collect([nq + 4 + 4 + 2 + 2 + 1, nq + 4 + 4 + 2 + 2 + 2 + 1]), collect([nq + 4 + 4 + 2 + 2 + 2, nq + 4 + 4 + 2 + 2 + 2 + 2])]],
        [[collect([nq + 4 + 4 + 1, nq + 4 + 4 + 2 + 1]), collect([nq + 4 + 4 + 2, nq + 4 + 4 + 2 + 2])], [collect([nq + 4 + 4 + 2 + 2 + 1, nq + 4 + 4 + 2 + 2 + 2 + 1]), collect([nq + 4 + 4 + 2 + 2 + 2, nq + 4 + 4 + 2 + 2 + 2 + 2])]],
        collect(1:(nq + 4 + 1 + 1 + 1 + 1)),
        collect(nq + 4 + 1 + 1 + 1 + 1 .+ (1:4)),
        collect(nq + 4 + 1 + 1 + 1 + 1 + 4 .+ (1:4)),
        [collect(nq + 4 + 1 + 1 + 1 + 1 + 4 + (i - 1) * 2 .+ (1:2)) for i = 1:2],
        collect(1:nq),
        collect(nq .+ (1:(4 + 1 + 1 + 1 + 1))),
        collect(nq + 4 + 1 + 1 + 1 + 1 .+ (1:8)))
end

function Trajectory(model::Hopper, T; nv=model.nq, nb=2) 
    nq = model.nq 
    nu = model.nu
    nw = model.nw 
    nc = 4
    q = [zeros(nq) for t = 1:T+2]
    v = [zeros(nv) for t = 1:T+1] 
    u = [zeros(nu) for t = 1:T]
    γ = [zeros(nc) for t = 1:T] 
    b = [zeros(nb) for t = 1:T] 
    w = [zeros(nw) for t = 1:T] 
    Trajectory(q, v, u, γ, b, w)
end