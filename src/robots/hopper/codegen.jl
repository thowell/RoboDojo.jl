path = joinpath(@__DIR__, "expr/expr.jld2")

if HOPPER_CODEGEN == :load 
    JLD2.@load path expr
elseif HOPPER_CODEGEN == :generate
    nq = hopper.nq
    nu = hopper.nu
    nw = hopper.nw

    # variables
    @variables q[1:nq] q̇[1:nq]

    # Lagrangian
    L = lagrangian(hopper, q, q̇)

    ddL = Symbolics.hessian(L, [q; q̇])#, simplify=true)
    dLq = Symbolics.gradient(L, q)#, simplify=true)
    ddLq̇q = ddL[nq .+ (1:nq), 1:nq]

    M = ddL[nq .+ (1:nq), nq .+ (1:nq)]
    C = ddLq̇q * q̇ - dLq

    mass_matrix_sym = eval(Symbolics.build_function(M, q)[1])
    dynamics_bias_sym = eval(Symbolics.build_function(C, q, q̇)[1])

    # ϕ = signed_distance(model, q)
    # dϕ = Symbolics.jacobian(ϕ, q)
    # signed_distance_jacobian = eval(Symbolics.build_function(dϕ, q)[1])

    function lagrangian_derivatives(model::Hopper, q, v)
        D1L = -1.0 * dynamics_bias_sym(q, v)
        D2L = mass_matrix_sym(q) * v
        return D1L, D2L
    end

    function dynamics(model::Hopper, h, q0, q1, u1, w1, λ1, q2)
        # evalutate at midpoint
        qm1 = 0.5 * (q0 + q1)
        vm1 = (q1 - q0) / h[1]
        qm2 = 0.5 * (q1 + q2)
        vm2 = (q2 - q1) / h[1]

        D1L1, D2L1 = lagrangian_derivatives(model, qm1, vm1)
        D1L2, D2L2 = lagrangian_derivatives(model, qm2, vm2)

        return (0.5 * h[1] * D1L1 + D2L1 + 0.5 * h[1] * D1L2 - D2L2
            + input_jacobian(model, qm2) * u1 
            - model.friction_joint * [0.0; 0.0; 0.0; vm2[4]] # joint friction
            + λ1)
    end

    # @variables h[1:1] q0[1:nq] q1[1:nq] u1[1:nu] w1[1:nw] λ1[1:4] q2[1:nq]
    # d = dynamics(hopper, h, q0, q1, u1, w1, λ1, q2)

    function residual(model::Hopper, z, θ, μ)
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
        J = contact_jacobian(model, q2)
        λ1 = transpose(J) * [transpose(rotation_body) * [b1[1]; γ1[1]];
                            transpose(rotation_body) * [b1[2]; γ1[2]];
                            γ1[3:4]]
        λ1[3] += model.body_radius * b1[1] # friction on body creates a moment

        [
        dynamics(model, h, q0, q1, u1, w1, λ1, q2);
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

    nz = num_var(hopper)
    nθ = num_var(hopper)

    @variables z[1:nz] θ[1:nθ] μ[1:1]
    r = residual(hopper, z, θ, μ)
    rz = Symbolics.jacobian(r, z)
    rθ = Symbolics.jacobian(r, θ)

    r_hopper = Symbolics.build_function(r, z, θ, μ)[2]
    rz_hopper = Symbolics.build_function(rz, z, θ)[2]
    rθ_hopper = Symbolics.build_function(rθ, z, θ)[2]

    expr = Dict{Symbol, Expr}()
    expr[:r] = r_hopper 
    expr[:rz] = rz_hopper 
    expr[:rθ] = rθ_hopper

    JLD2.@save path expr
end

r_hopper! = eval(expr[:r])
rz_hopper! = eval(expr[:rz])
rθ_hopper! = eval(expr[:rθ])

residual_name(::Hopper) = :r_hopper!
jacobian_var_name(::Hopper) = :rz_hopper! 
jacobian_data_name(::Hopper) = :rθ_hopper! 

# using BenchmarkTools
# @benchmark r_hopper!($r0, $z0, $θ0, $μ0)
# @benchmark rz_hopper!($rz0, $z0, $θ0)
# @benchmark rθ_hopper!($rθ0, $z0, $θ0)