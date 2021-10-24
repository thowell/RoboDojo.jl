path = joinpath(@__DIR__, "expr/expr.jld2")

if BIPED_CODEGEN == :load 
    @load path r_biped rz_biped rθ_biped 
else
    nq = biped.nq
    nu = biped.nu
    nw = biped.nw
    # variables
    @variables q[1:nq] q̇[1:nq]

    # Lagrangian
    L = lagrangian(biped, q, q̇)

    ddL = Symbolics.hessian(L, [q; q̇])#, simplify=true)
    dLq = Symbolics.gradient(L, q)#, simplify=true)
    ddLq̇q = ddL[nq .+ (1:nq), 1:nq]

    M = ddL[nq .+ (1:nq), nq .+ (1:nq)]
    C = ddLq̇q * q̇ - dLq

    mass_matrix_sym = eval(Symbolics.build_function(M, q)[1])
    dynamics_bias_sym = eval(Symbolics.build_function(C, q, q̇)[1])

    function lagrangian_derivatives(model::Biped, q, v)
        D1L = -1.0 * dynamics_bias_sym(q, v)
        D2L = mass_matrix_sym(q) * v
        return D1L, D2L
    end

    function dynamics(model::Biped, h, q0, q1, u1, w1, λ1, q2)
        # evalutate at midpoint
        qm1 = 0.5 * (q0 + q1)
        vm1 = (q1 - q0) / h[1]
        qm2 = 0.5 * (q1 + q2)
        vm2 = (q2 - q1) / h[1]

        D1L1, D2L1 = lagrangian_derivatives(model, qm1, vm1)
        D1L2, D2L2 = lagrangian_derivatives(model, qm2, vm2)

        return (0.5 * h[1] * D1L1 + D2L1 + 0.5 * h[1] * D1L2 - D2L2
            + transpose(input_jacobian(model)) * u1 
            - [0.0; 0.0; 0.0; model.friction_joint .* vm2[4:9]] # joint friction
            + λ1)
    end

    # @variables h[1:1] q0[1:nq] q1[1:nq] u1[1:nu] w1[1:nw] λ1[1:nq] q2[1:nq]
    # d = dynamics(biped, h, q0, q1, u1, w1, λ1, q2)

    function residual(model::Biped, z, θ, μ)
        # dimensions 
        nq = model.nq
        nu = model.nu
        nw = model.nw
        nc = model.nc

        # unpack data
        q0 = θ[1:nq] 
        q1 = θ[nq .+ (1:nq)] 
        u1 = θ[2nq .+ (1:nu)] 
        w1 = θ[2nq + nu .+ (1:nw)]
        friction_body_world = θ[2nq + nu + nw .+ (1:4)] 
        friction_foot_world = θ[2nq + nu + nw + 4 .+ (1:4)] 
        h = θ[2nq + nu + nw + 8 .+ (1:1)] 

        # unpack variables
        q2 = z[1:nq] 
        γ1 = z[nq .+ (1:nc)] 
        sγ1 = z[nq + nc .+ (1:nc)]
        ψ1 = z[nq + nc + nc .+ (1:nc)] 
        b1 = z[nq + nc + nc + nc .+ (1:nc)] 
        sψ1 = z[nq + nc + nc + nc + nc .+ (1:nc)] 
        sb1 = z[nq + nc + nc + nc + nc + nc .+ (1:nc)]   

        rotation_toe1 = vector_rotation_matrix([0.0; 1.0])
        rotation_heel1 = vector_rotation_matrix([0.0; 1.0])
        rotation_toe2 = vector_rotation_matrix([0.0; 1.0])
        rotation_heel2 = vector_rotation_matrix([0.0; 1.0])

        rotation_knee1 = vector_rotation_matrix([0.0; 1.0])
        rotation_knee2 = vector_rotation_matrix([0.0; 1.0])

        rotation_body1 = vector_rotation_matrix([0.0; 1.0])
        rotation_body2 = vector_rotation_matrix([0.0; 1.0])

        # tangential velocity
        v1 = (q2 - q1) / h[1] 

        vT_toe1 = (rotation_toe1 * kinematics_jacobian_foot(model, q2, leg=:leg1, mode=:toe) * v1)[1]
        vT_heel1 = (rotation_heel1 * kinematics_jacobian_foot(model, q2, leg=:leg1, mode=:heel) * v1)[1]
        vT_toe2 = (rotation_toe2 * kinematics_jacobian_foot(model, q2, leg=:leg2, mode=:toe) * v1)[1]
        vT_heel2 = (rotation_heel2 * kinematics_jacobian_foot(model, q2, leg=:leg2, mode=:heel) * v1)[1]

        vT_knee1 = (rotation_knee1 * kinematics_jacobian_thigh(model, q2, leg=:leg1, mode=:ee) * v1)[1]
        vT_knee2 = (rotation_knee2 * kinematics_jacobian_thigh(model, q2, leg=:leg2, mode=:ee) * v1)[1]
        
        vT_body1 = (rotation_body1 * kinematics_jacobian_body(model, q2, mode=:hip) * v1)[1]
        vT_body2 = (rotation_body2 * kinematics_jacobian_body(model, q2, mode=:ee) * v1)[1]

        # contact forces
        J = contact_jacobian(model, q2)
        λ1 = transpose(J) * [
                            transpose(rotation_toe1)  * [b1[1]; γ1[1]];
                            transpose(rotation_heel1) * [b1[2]; γ1[2]];
                            transpose(rotation_toe2)  * [b1[3]; γ1[3]];
                            transpose(rotation_heel2) * [b1[4]; γ1[4]];
                            transpose(rotation_knee1) * [b1[5]; γ1[5]];
                            transpose(rotation_knee2) * [b1[6]; γ1[6]];
                            transpose(rotation_body1) * [b1[7]; γ1[7]];
                            transpose(rotation_body2) * [b1[8]; γ1[8]];
                            ]

        [
        dynamics(model, h, q0, q1, u1, w1, λ1, q2);
        signed_distance(model, q2) - sγ1;
        vT_toe1 - sb1[1];
        vT_heel1 - sb1[2];
        vT_toe2 - sb1[3];
        vT_heel2 - sb1[4];
        vT_knee1 - sb1[5];
        vT_knee2 - sb1[6];
        vT_body1 - sb1[7];
        vT_body2 - sb1[8];
        friction_foot_world[1] * γ1[1] - ψ1[1];
        friction_foot_world[2] * γ1[2] - ψ1[2];
        friction_foot_world[3] * γ1[3] - ψ1[3];
        friction_foot_world[4] * γ1[4] - ψ1[4];
        friction_body_world[1] * γ1[5] - ψ1[5];
        friction_body_world[2] * γ1[6] - ψ1[6];
        friction_body_world[3] * γ1[7] - ψ1[7];
        friction_body_world[4] * γ1[8] - ψ1[8];
        γ1 .* sγ1 .- μ[1];
        cone_product([ψ1[1]; b1[1]], [sψ1[1]; sb1[1]]) - [μ[1]; 0.0]; 
        cone_product([ψ1[2]; b1[2]], [sψ1[2]; sb1[2]]) - [μ[1]; 0.0]; 
        cone_product([ψ1[3]; b1[3]], [sψ1[3]; sb1[3]]) - [μ[1]; 0.0]; 
        cone_product([ψ1[4]; b1[4]], [sψ1[4]; sb1[4]]) - [μ[1]; 0.0]; 
        cone_product([ψ1[5]; b1[5]], [sψ1[5]; sb1[5]]) - [μ[1]; 0.0]; 
        cone_product([ψ1[6]; b1[6]], [sψ1[6]; sb1[6]]) - [μ[1]; 0.0]; 
        cone_product([ψ1[7]; b1[7]], [sψ1[7]; sb1[7]]) - [μ[1]; 0.0]; 
        cone_product([ψ1[8]; b1[8]], [sψ1[8]; sb1[8]]) - [μ[1]; 0.0]; 
        ]
    end

    nz = num_var(biped)
    nθ = num_data(biped, nf=length(friction_coefficients(biped)))

    @variables z[1:nz] θ[1:nθ] μ[1:1]
    r = residual(biped, z, θ, μ)
    rz = Symbolics.jacobian(r, z)
    rθ = Symbolics.jacobian(r, θ)

    r_biped = Symbolics.build_function(r, z, θ, μ)[2]
    rz_biped = Symbolics.build_function(rz, z, θ)[2]
    rθ_biped = Symbolics.build_function(rθ, z, θ)[2]

    @save path r_biped rz_biped rθ_biped
end

r_biped! = eval(r_biped)
rz_biped! = eval(rz_biped)
rθ_biped! = eval(rθ_biped)

residual_name(::Biped) = :r_biped!
jacobian_var_name(::Biped) = :rz_biped! 
jacobian_data_name(::Biped) = :rθ_biped! 

# using BenchmarkTools
# @benchmark r_biped!($r0, $z0, $θ0, $μ0)
# @benchmark rz_biped!($rz0, $z0, $θ0)
# @benchmark rθ_biped!($rθ0, $z0, $θ0)
