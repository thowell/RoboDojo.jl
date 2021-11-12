function residual(model::Model, mass_matrix, dynamics_bias, kinematics, kinematics_jacobians, z, θ, μ)
    # dimensions 
    nq = model.nq
    nu = model.nu
    nw = model.nw
    nc = model.nc
    nf = length(friction_coefficients(model))

    # unpack data
    q0 = θ[1:nq] 
    q1 = θ[nq .+ (1:nq)] 
    u1 = θ[2nq .+ (1:nu)] 
    w1 = θ[2nq + nu .+ (1:nw)]
    friction = θ[2nq + nu + nw .+ (1:nf)] 
    h = θ[2nq + nu + nw + nf .+ (1:1)] 

    # unpack variables
    q2 = z[1:nq] 
    γ1 = z[nq .+ (1:nc)] 
    sγ1 = z[nq + nc .+ (1:nc)]
    ψ1 = z[nq + nc + nc .+ (1:nc)] 
    b1 = z[nq + nc + nc + nc .+ (1:nc)] 
    sψ1 = z[nq + nc + nc + nc + nc .+ (1:nc)] 
    sb1 = z[nq + nc + nc + nc + nc + nc .+ (1:nc)]   

    # rotation matrix at surface
    rotation = [vector_rotation_matrix([0.0; 1.0]) for i = 1:nc] #TODO non-flat surfaces

    # tangential velocity
    v1 = (q2 - q1) / h[1] 
    vT = [(rotation[i] * kinematics_jacobians[i](q2) * v1)[1] for i = 1:nc]

    # contact forces
    J = contact_jacobian(model, q1)
    λ1 = transpose(J) * vcat([transpose(rotation[i]) * [b1[i]; γ1[i]] for i = 1:nc]...)

    # residual
    [
    dynamics(model, mass_matrix, dynamics_bias, h, q0, q1, u1, w1, λ1, q2);
    signed_distance(model, q2) - sγ1;
    vcat([vT[i] - sb1[i] for i = 1:nc]...);
    vcat([friction[i] * γ1[i] - ψ1[i] for i = 1:nc]...);
    γ1 .* sγ1 .- μ[1];
    vcat([cone_product([ψ1[i]; b1[i]], [sψ1[i]; sb1[i]]) - [μ[1]; 0.0] for i = 1:nc]...); 
    ]
end