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
    nq = model.nq
    nz = num_var(model) 

    IndicesOptimization(
        nz, 
        nz,  
        [collect(nq .+ (1:4)), collect(nq + 4 .+ (1:4))],
        [collect(nq .+ (1:4)), collect(nq + 4 .+ (1:4))],
        [[collect([nq + 4 + 4 + 1, nq + 4 + 4 + 2 + 1]), collect([nq + 4 + 4 + 2 + 2 + 1, nq + 4 + 4 + 2 + 2 + 2 + 1])],
         [collect([nq + 4 + 4 + 2, nq + 4 + 4 + 2 + 2]), collect([nq + 4 + 4 + 2 + 2 + 2, nq + 4 + 4 + 2 + 2 + 2 + 2])]],
        [[collect([nq + 4 + 4 + 1, nq + 4 + 4 + 2 + 1]), collect([nq + 4 + 4 + 2 + 2 + 1, nq + 4 + 4 + 2 + 2 + 2 + 1])],
         [collect([nq + 4 + 4 + 2, nq + 4 + 4 + 2 + 2]), collect([nq + 4 + 4 + 2 + 2 + 2, nq + 4 + 4 + 2 + 2 + 2 + 2])]],
        collect(1:(nq + 4 + 1 + 1 + 1 + 1)),
        collect(nq + 4 + 1 + 1 + 1 + 1 .+ (1:4)),
        collect(nq + 4 + 1 + 1 + 1 + 1 + 4 .+ (1:4)),
        [collect(nq + 4 + 1 + 1 + 1 + 1 + 4 + (i - 1) * 2 .+ (1:2)) for i = 1:2],
        collect(nq + 4 + 1 + 1 + 1 + 1 .+ (1:8)))
end

Trajectory(model::Hopper, T) = Trajectory(model, T, nc=4, nb=2)
GradientTrajectory(model::Hopper, T) = GradientTrajectory(model, T, nc=4, nb=2) 
 


