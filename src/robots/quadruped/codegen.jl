path = joinpath(@__DIR__, "expr/expr.jld2")

if false#QUADRUPED_CODEGEN == :load 
    @load path r_quadruped rz_quadruped rθ_quadruped
else
    nq = quadruped.nq
    nu = quadruped.nu
    nw = quadruped.nw

    # variables
    @variables q[1:nq] q̇[1:nq]

    # Lagrangian
    L = lagrangian(quadruped, q, q̇);

    ddL = Symbolics.hessian(L, [q; q̇])#, simplify=true)
    dLq = Symbolics.gradient(L, q)#, simplify=true)
    ddLq̇q = ddL[nq .+ (1:nq), 1:nq]

    M = ddL[nq .+ (1:nq), nq .+ (1:nq)]
    C = ddLq̇q * q̇ - dLq

    mass_matrix_sym = eval(Symbolics.build_function(M, q)[1])
    dynamics_bias_sym = eval(Symbolics.build_function(C, q, q̇)[1])

    # ϕ = signed_distance(quadruped, q)
    # dϕ = Symbolics.jacobian(ϕ, q)
    # signed_distance_jacobian = eval(Symbolics.build_function(dϕ, q)[1])

    function lagrangian_derivatives(model::Quadruped, q, v)
        D1L = -1.0 * dynamics_bias_sym(q, v)
        D2L = mass_matrix_sym(q) * v
        return D1L, D2L
    end

    function dynamics(model::Quadruped, h, q0, q1, u1, w1, λ1, q2)
        # evalutate at midpoint
        qm1 = 0.5 * (q0 + q1)
        vm1 = (q1 - q0) / h[1]
        qm2 = 0.5 * (q1 + q2)
        vm2 = (q2 - q1) / h[1]

        D1L1, D2L1 = lagrangian_derivatives(model, qm1, vm1)
        D1L2, D2L2 = lagrangian_derivatives(model, qm2, vm2)

        return (0.5 * h[1] * D1L1 + D2L1 + 0.5 * h[1] * D1L2 - D2L2
            + transpose(input_jacobian(model)) * u1 
            - [0.0; 0.0; 0.0; model.friction_joint .* vm2[4:11]] # joint friction
            + λ1)
    end

    # @variables h[1:1] q0[1:nq] q1[1:nq] u1[1:nu] w1[1:nw] λ1[1:nq] q2[1:nq]
    # d = dynamics(quadruped, h, q0, q1, u1, w1, λ1, q2)

    function residual(model::Quadruped, z, θ, μ)
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
        friction_body_world = θ[2nq + nu + nw .+ (1:6)] 
        friction_foot_world = θ[2nq + nu + nw + 6 .+ (1:4)] 
        h = θ[2nq + nu + nw + 10 .+ (1:1)] 

        # unpack variables
        q2 = z[1:nq] 
        γ1 = z[nq .+ (1:nc)] 
        sγ1 = z[nq + nc .+ (1:nc)]
        ψ1 = z[nq + nc + nc .+ (1:nc)] 
        b1 = z[nq + nc + nc + nc .+ (1:nc)] 
        sψ1 = z[nq + nc + nc + nc + nc .+ (1:nc)] 
        sb1 = z[nq + nc + nc + nc + nc + nc .+ (1:nc)]   

        rotation_foot1 = vector_rotation_matrix([0.0; 1.0])
        rotation_foot2 = vector_rotation_matrix([0.0; 1.0])
        rotation_foot3 = vector_rotation_matrix([0.0; 1.0])
        rotation_foot4 = vector_rotation_matrix([0.0; 1.0])

        rotation_knee1 = vector_rotation_matrix([0.0; 1.0])
        rotation_knee2 = vector_rotation_matrix([0.0; 1.0])
        rotation_knee3 = vector_rotation_matrix([0.0; 1.0])
        rotation_knee4 = vector_rotation_matrix([0.0; 1.0])

        rotation_hip1 = vector_rotation_matrix([0.0; 1.0])
        rotation_hip2 = vector_rotation_matrix([0.0; 1.0])

        # tangential velocity
        v1 = (q2 - q1) / h[1] 

        vT_foot1 = (rotation_foot1 * kinematics_jacobian_calf(model, q2, leg=:leg1, mode=:ee) * v1)[1]
        vT_foot2 = (rotation_foot2 * kinematics_jacobian_calf(model, q2, leg=:leg2, mode=:ee) * v1)[1]
        vT_foot3 = (rotation_foot3 * kinematics_jacobian_calf(model, q2, leg=:leg3, mode=:ee) * v1)[1]
        vT_foot4 = (rotation_foot4 * kinematics_jacobian_calf(model, q2, leg=:leg4, mode=:ee) * v1)[1]

        vT_knee1 = (rotation_knee1 * kinematics_jacobian_thigh(model, q2, leg=:leg1, mode=:ee) * v1)[1]
        vT_knee2 = (rotation_knee2 * kinematics_jacobian_thigh(model, q2, leg=:leg2, mode=:ee) * v1)[1]
        vT_knee3 = (rotation_knee3 * kinematics_jacobian_thigh(model, q2, leg=:leg3, mode=:ee) * v1)[1]
        vT_knee4 = (rotation_knee4 * kinematics_jacobian_thigh(model, q2, leg=:leg4, mode=:ee) * v1)[1]

        vT_hip1 = (rotation_hip1 * kinematics_jacobian_hip(model, q2, hip=:hip1) * v1)[1]
        vT_hip2 = (rotation_hip2 * kinematics_jacobian_hip(model, q2, hip=:hip2) * v1)[1]

        # contact forces
        J = contact_jacobian(model, q2)
        λ1 = transpose(J) * [
                            transpose(rotation_foot1) * [b1[1]; γ1[1]];
                            transpose(rotation_foot2) * [b1[2]; γ1[2]];
                            transpose(rotation_foot3) * [b1[3]; γ1[3]];
                            transpose(rotation_foot4) * [b1[4]; γ1[4]];
                            transpose(rotation_knee1) * [b1[5]; γ1[5]];
                            transpose(rotation_knee2) * [b1[6]; γ1[6]];
                            transpose(rotation_knee3) * [b1[7]; γ1[7]];
                            transpose(rotation_knee4) * [b1[8]; γ1[8]];
                            transpose(rotation_hip1)  * [b1[9]; γ1[9]];
                            transpose(rotation_hip2)  * [b1[10]; γ1[10]];
                            ]

        [
        dynamics(model, h, q0, q1, u1, w1, λ1, q2);
        signed_distance(model, q2) - sγ1;
        vT_foot1 - sb1[1];
        vT_foot2 - sb1[2];
        vT_foot3 - sb1[3];
        vT_foot4 - sb1[4];
        vT_knee1 - sb1[5];
        vT_knee2 - sb1[6];
        vT_knee3 - sb1[7];
        vT_knee4 - sb1[8];
        vT_hip1  - sb1[9]; 
        vT_hip2  - sb1[10];
        friction_foot_world[1] * γ1[1] - ψ1[1];
        friction_foot_world[2] * γ1[2] - ψ1[2];
        friction_foot_world[3] * γ1[3] - ψ1[3];
        friction_foot_world[4] * γ1[4] - ψ1[4];
        friction_body_world[1] * γ1[5] - ψ1[5];
        friction_body_world[2] * γ1[6] - ψ1[6];
        friction_body_world[3] * γ1[7] - ψ1[7];
        friction_body_world[4] * γ1[8] - ψ1[8];
        friction_body_world[5] * γ1[9] - ψ1[9];
        friction_body_world[6] * γ1[10] - ψ1[10];
        γ1 .* sγ1 .- μ[1];
        cone_product([ψ1[1]; b1[1]], [sψ1[1]; sb1[1]]) - [μ[1]; 0.0]; 
        cone_product([ψ1[2]; b1[2]], [sψ1[2]; sb1[2]]) - [μ[1]; 0.0]; 
        cone_product([ψ1[3]; b1[3]], [sψ1[3]; sb1[3]]) - [μ[1]; 0.0]; 
        cone_product([ψ1[4]; b1[4]], [sψ1[4]; sb1[4]]) - [μ[1]; 0.0]; 
        cone_product([ψ1[5]; b1[5]], [sψ1[5]; sb1[5]]) - [μ[1]; 0.0]; 
        cone_product([ψ1[6]; b1[6]], [sψ1[6]; sb1[6]]) - [μ[1]; 0.0]; 
        cone_product([ψ1[7]; b1[7]], [sψ1[7]; sb1[7]]) - [μ[1]; 0.0]; 
        cone_product([ψ1[8]; b1[8]], [sψ1[8]; sb1[8]]) - [μ[1]; 0.0]; 
        cone_product([ψ1[9]; b1[9]], [sψ1[9]; sb1[9]]) - [μ[1]; 0.0]; 
        cone_product([ψ1[10]; b1[10]], [sψ1[10]; sb1[10]]) - [μ[1]; 0.0]; 
        ]
    end

    nz = num_var(quadruped)
    nθ = num_data(quadruped, nf=length(friction_coefficients(quadruped)))

    @variables z[1:nz] θ[1:nθ] μ[1:1]
    r = residual(quadruped, z, θ, μ)
    rz = Symbolics.jacobian(r, z)
    rθ = Symbolics.jacobian(r, θ)
  
    r_quadruped = Symbolics.build_function(r, z, θ, μ)[2]
    rz_quadruped = Symbolics.build_function(rz, z, θ)[2]
    rθ_quadruped = Symbolics.build_function(rθ, z, θ)[2]

    @save path r_quadruped rz_quadruped rθ_quadruped
end

r_quadruped! = eval(r_quadruped)
rz_quadruped! = eval(rz_quadruped)
rθ_quadruped! = eval(rθ_quadruped)

residual_name(::Quadruped) = :r_quadruped!
jacobian_var_name(::Quadruped) = :rz_quadruped! 
jacobian_data_name(::Quadruped) = :rθ_quadruped! 

# using BenchmarkTools
# @benchmark r_quadruped!($r0, $z0, $θ0, $μ0)
# @benchmark rz_quadruped!($rz0, $z0, $θ0)
# @benchmark rθ_quadruped!($rθ0, $z0, $θ0)