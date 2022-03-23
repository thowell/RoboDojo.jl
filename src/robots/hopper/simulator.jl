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
 
function residual(model::Hopper, mass_matrix, dynamics_bias, kinematics, kinematics_jacobians, z, θ, μ)
    
    # dimensions 
    nq = model.nq
    nu = model.nu
    nw = model.nw

    # unpack data
    q0 = θ[1:nq] 
    q1 = θ[nq .+ (1:nq)] 
    u1 = θ[2nq .+ (1:nu)] 
    w1 = θ[2nq + nu .+ (1:nw)]
    friction_body_world = θ[2nq + nu + nw .+ (1:1)] 
    friction_foot_world = θ[2nq + nu + nw + 1 .+ (1:1)] 
    h = θ[2nq + nu + nw + 2 .+ (1:1)] 

    # unpack variables
    q2 = z[1:nq] 
    γ1 = z[nq .+ (1:4)] 
    sγ1 = z[nq + 4 .+ (1:4)]
    ψ1 = z[nq + 4 + 4 .+ (1:2)] 
    b1 = z[nq + 4 + 4 + 2 .+ (1:2)] 
    sψ1 = z[nq + 4 + 4 + 2 + 2 .+ (1:2)] 
    sb1 = z[nq + 4 + 4 + 2 + 2 + 2 .+ (1:2)]   

    rotation_body = vector_rotation_matrix([0.0; 1.0]) 
    rotation_foot = vector_rotation_matrix([0.0; 1.0])

    # tangential velocity
    v1 = (q2 - q1) / h[1] 
    vT_body = (rotation_body * v1[1:2])[1] + model.body_radius * v1[3]
    vT_foot = (rotation_foot * kinematics_foot_jacobian(model, q2) * v1)[1]

    # contact forces
    J = contact_jacobian(model, q1)
    λ1 = transpose(J) * [transpose(rotation_body) * [b1[1]; γ1[1]];
                         transpose(rotation_foot) * [b1[2]; γ1[2]];
                        γ1[3:4]]
    λ1[3] += model.body_radius * b1[1] # friction on body creates a moment

    [
    dynamics(model, mass_matrix, dynamics_bias, h, q0, q1, u1, w1, λ1, q2);
    signed_distance(model, q2) - sγ1;
    vT_body - sb1[1];
    vT_foot - sb1[2];
    friction_body_world[1] * γ1[1] - ψ1[1];
    friction_foot_world[1] * γ1[2] - ψ1[2];
    γ1 .* sγ1 .- μ[1];
    cone_product([ψ1[1]; b1[1]], [sψ1[1]; sb1[1]]) - [μ[1]; 0.0]; 
    cone_product([ψ1[2]; b1[2]], [sψ1[2]; sb1[2]]) - [μ[1]; 0.0]; 
    ]
end

