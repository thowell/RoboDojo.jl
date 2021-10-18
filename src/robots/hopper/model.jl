"""
    hopper
    	model inspired by Raibert hopper 
            from "Dynamically Stable Legged Locomotion"
		q = (x, z, t, r)
			x - lateral position
			z - vertical position
			t - body orientation
			r - leg length
"""
struct Hopper{T}
    # dimensions
    nq::Int # generalized coordinates
    nu::Int # controls
    nw::Int # parameters

    # parameters
    mass_body::T
    mass_foot::T 
    inertia_body::T 
    body_radius::T 
    foot_radius::T 
    leg_len_max::T 
    leg_len_min::T
    friction_joint::T

    # environment
    friction_body_world::T
    friction_foot_world::T
	gravity::T
end

# kinematics
function kinematics_body(::Hopper, q)
	[q[1],
	 q[2]]
end

function kinematics_body_jacobian(::Hopper, q) 
    [1.0 0.0 0.0 0.0;
	 0.0 1.0 0.0 0.0] 
end

function kinematics_foot(::Hopper, q)
	[q[1] + q[4] * sin(q[3]),
	 q[2] - q[4] * cos(q[3])]
end

function kinematics_foot_jacobian(::Hopper, q) 
    [1.0 0.0 (q[4] * cos(q[3])) sin(q[3]);
	 0.0 1.0 (q[4] * sin(q[3])) (-1.0 * cos(q[3]))] 
end

# Lagrangian
function lagrangian(model::Hopper, q, q̇) 
    L = 0.0 

    # body kinetic energy
    v_body = q̇[1:2] 
    ω_body = q̇[3]
    L += 0.5 * model.mass_body * dot(v_body, v_body) 
    L += 0.5 * model.inertia_body * ω_body^2.0 

    # body potential energy
    z_body = q[2]
    L -= model.mass_body * model.gravity * z_body 

    # foot kinetic energy 
    v_foot = kinematics_foot_jacobian(model, q) * q̇
    L += 0.5 * model.mass_foot * dot(v_foot, v_foot)

    # foot potential energy
    z_foot =  kinematics_foot(model, q)[2]
    L -= model.mass_foot * model.gravity * z_foot

    return L
end

# signed distance
function signed_distance(model::Hopper, q)
    [
     q[2] - model.body_radius,                    # body-world distance
     kinematics_foot(model, q)[2] - model.foot_radius, # foot-world distance 
     q[4] - model.leg_len_min,                    # leg length min.
     model.leg_len_max - q[4]                     # leg length max.
    ]
end

function signed_distance(model::Hopper, q, pt1, pt2)
    # error between robot kinematics and point on nearest surface
    e1 = kinematics_body(model, q) - pt1
    e2 = kinematics_foot(model, q) - pt2
    @show pt1 
    @show pt2
    [
     sqrt(dot(e1, e1)) - model.body_radius,  # body-world distance
     sqrt(dot(e2, e2)) - model.foot_radius,  # foot-world distance 
     q[4] - model.leg_len_min,             # leg length min.
     model.leg_len_max - q[4]              # leg length max.
    ]
end

# input matrix
function input_jacobian(model, q)
	[
     0.0 -sin(q[3]); 
     0.0 cos(q[3]); 
     1.0 0.0; 
     0.0 1.0
    ]
end

# dimensions
nq = 4
nu = 2
nw = 0

# parameters
mass_body = 1.0 
mass_foot = 0.1 
inertia_body = 0.1 
body_radius = 0.1 
foot_radius = 0.05
leg_len_max = 1.0 
leg_len_min = 0.25
friction_joint = 0.1

# environment
friction_body_world = 0.5
friction_foot_world = 0.5
gravity = 9.81 

# model
hopper = Hopper(nq, nu, nw,
            mass_body, mass_foot, inertia_body, 
            body_radius, foot_radius, 
            leg_len_max, leg_len_min,
            friction_joint,
            friction_body_world, friction_foot_world, 
            gravity)
			   
# code-gen
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

function contact_jacobian(model, q)
    [1.0 0.0 0.0 0.0; # note that body torque due to friction is added in residual
     0.0 1.0 0.0 0.0;
     1.0 0.0 (q[4] * cos(q[3])) sin(q[3]);
     0.0 1.0 (q[4] * sin(q[3])) (-1.0 * cos(q[3]));
     0.0 0.0 0.0 1.0; 
     0.0 0.0 0.0 -1.0] 
end

function lagrangian_derivatives(model, q, v)
	D1L = -1.0 * dynamics_bias_sym(q, v)
    D2L = mass_matrix_sym(q) * v
	return D1L, D2L
end

function dynamics(model, h, q0, q1, u1, w1, λ1, q2)
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

function residual(model, z, θ, μ)
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

nz = nq + 4 + 4 + 2 + 2 + 2 + 2 
ny = 4 + 4 + 2 + 2 + 2 + 2
nθ = 2nq + nu + nw + 2 + 1 

@variables z[1:nz] θ[1:nθ] μ[1:1]
r = residual(hopper, z, θ, μ)
rz = Symbolics.jacobian(r, z)
rθ = Symbolics.jacobian(r, θ)
rz_sp = Symbolics.sparsejacobian(r, z)
rθ_sp = Symbolics.sparsejacobian(r, θ)

r_hopper! = eval(Symbolics.build_function(r, z, θ, μ)[2])
rz_hopper! = eval(Symbolics.build_function(rz, z, θ)[2])
rθ_hopper! = eval(Symbolics.build_function(rθ, z, θ)[2])
# rz_sp_hopper! = eval(Symbolics.build_function(rz_sp.nzval, z, θ)[2])
# rθ_sp_hopper! = eval(Symbolics.build_function(rθ_sp.nzval, z, θ)[2])
# idx_rz_sp = [findnz(rz_sp)[1:2]...]
# idx_rθ_sp = [findnz(rθ_sp)[1:2]...]

# r0 = zeros(length(r))
# rz0 = zeros(size(rz))
# rθ0 = zeros(size(rθ))
# rz0_sp_vec = similar(rz_sp.nzval, Float64)
# rθ0_sp_vec = similar(rθ_sp.nzval, Float64)
# rz0_sp_mat = similar(rz_sp, Float64)
# rθ0_sp_mat = similar(rθ_sp, Float64)

# rz_sp

# # rz0_sp = zeros(size(rz_sp.nzval))
# # rθ0_sp = zeros(size(rθ_sp.nzval))

# z0 = rand(nz) 
# θ0 = rand(nθ) 
# μ0 = [1.0]

# rz_sp_hopper!(rz0_sp_vec, z0, θ0)
# rθ_sp_hopper!(rθ0_sp_vec, z0, θ0)
# rz0_sp_mat


# rz0_sp_mat[CartesianIndices(idx_rz_sp[1], idx_rz_sp[2])]
# rθ0_sp_mat[idx_rθ_sp[1], idx_rθ_sp[2]] .= rθ0_sp


# @benchmark rz_sp_hopper!($rz0_sp, $z0, $θ0)
# @benchmark rθ_sp_hopper!($rθ0_sp, $z0, $θ0)


# using BenchmarkTools
# @benchmark r_hopper!($r0, $z0, $θ0, $μ0)
# @benchmark rz_hopper!($rz0, $z0, $θ0)
# @benchmark rθ_hopper!($rθ0, $z0, $θ0)


# A = rand(50, 50)
# A = A' * A
# b = rand(50)

# x = A \ b

# xsp = zero(x)
# Asp = sparse(A)
# fact = ilu0(Asp)

# function ilu_solve!(fact, A, b, x)
#     ilu0!(fact, A)
#     ldiv!(x, fact, b)
# end
# ilu_solve!(fact, Asp, b, xsp)
# @benchmark ilu_solve!($fact, $Asp, $b, $xsp)

# norm(x - xsp, Inf) #< 1.0e-8

# @benchmark $x = $A \ $b

# solv = lu_solver(copy(A))
# xlu = zero(x) 
# linear_solve!(solv, xlu, A, b)
# @benchmark linear_solve!($solv, $xlu, $A, $b)

# norm(x - xlu, Inf)

# fact_sp = lu(Asp)
