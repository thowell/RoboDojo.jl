function num_var(model::Particle) 
    nq = model.nq
    nc = model.nc 
    nb = 2
    nq + nc + nc + nc + nc * nb + nc + nc * nb
end 

function indices_z(model::Particle) 
    nq = model.nq 
    nc = model.nc 
    nb = 2

    q = collect(1:nq) 
    γ = collect(nq .+ (1:nc)) 
    sγ = collect(nq + nc .+ (1:nc))
    ψ = collect(nq + nc + nc .+ (1:nc)) 
    b = collect(nq + nc + nc + nc .+ (1:(nc * nb))) 
    sψ = collect(nq + nc + nc + nc + nc * nb .+ (1:nc)) 
    sb = collect(nq + nc + nc + nc + nc * nb + nc .+ (1:(nc * nb)))   
    IndicesZ(q, γ, sγ, ψ, b, sψ, sb)
end

function indices_optimization(model::Particle) 
    nq = model.nq
    nc = model.nc 
    nb = 2
    nz = num_var(model) 

    IndicesOptimization(
        nz, 
        nz,  
        [collect(nq .+ (1:nc)), collect(nq + nc .+ (1:nc))],
        [collect(nq .+ (1:nc)), collect(nq + nc .+ (1:nc))],
        [
            [collect([nq + nc + nc + 1, (nq + nc + nc + nc .+ (1:nb))...]), collect([nq + nc + nc + nc + nb + 1, (nq + nc + nc + nc + nb + nc .+ (1:nb))...])],
        ],
        [
            [collect([nq + nc + nc + 1, (nq + nc + nc + nc .+ (1:nb))...]), collect([nq + nc + nc + nc + nb + 1, (nq + nc + nc + nc + nb + nc .+ (1:nb))...])],
        ],
        collect(1:(nq + nc + nc * nb + nc)),
        collect(nq + nc + nc * nb + nc .+ (1:nc)),
        collect(nq + nc + nc * nb + nc + nc .+ (1:3)),
        [collect(nq + nc + nc * nb + nc + nc .+ (1:3))],
        collect(nq + nc + nc * nb + nc .+ (1:(nc + 3))))
end

Trajectory(model::Particle, T) = Trajectory(model, T, nc=model.nc, nb=2)
GradientTrajectory(model::Particle, T) = GradientTrajectory(model, T, nc=model.nc, nb=2) 
 
function residual(model::Particle, mass_matrix, dynamics_bias, kinematics, kinematics_jacobians, z, θ, μ)
    
    # dimensions 
    nq = model.nq
    nu = model.nu
    nw = model.nw
    nc = model.nc 
    nf = length(friction_coefficients(model))
    nb = 2

    # unpack data
    q0 = θ[1:nq] 
    q1 = θ[nq .+ (1:nq)] 
    u1 = θ[2nq .+ (1:nu)] 
    w1 = θ[2nq + nu .+ (1:nw)]
    friction_body_world = θ[2nq + nu + nw .+ (1:nf)] 
    h = θ[2nq + nu + nw + nf .+ (1:1)] 

    # unpack variables
    q2 = z[1:nq] 
    γ1 = z[nq .+ (1:nc)] 
    sγ1 = z[nq + nc .+ (1:nc)]
    ψ1 = z[nq + nc + nc .+ (1:nc)] 
    b1 = z[nq + nc + nc + nc .+ (1:(nc * nb))] 
    sψ1 = z[nq + nc + nc + nc + nc * nb .+ (1:nc)] 
    sb1 = z[nq + nc + nc + nc + nc * nb + nc .+ (1:(nc * nb))]   

    # tangential velocity
    v1 = (q2 - q1) / h[1] 
    vT = v1[1:2]

    # contact forces
    J = contact_jacobian(model, q1)
    λ1 = transpose(J) * [b1; γ1]

    [
    dynamics(model, mass_matrix, dynamics_bias, h, q0, q1, u1, w1, λ1, q2);
    signed_distance(model, q2) - sγ1;
    vT - sb1;
    friction_body_world[1] * γ1[1] - ψ1[1];
    γ1 .* sγ1 .- μ[1];
    cone_product([ψ1; b1], [sψ1; sb1]) - [μ[1]; 0.0; 0.0]; 
    ]
end

