"""
    Egg(Bot)
		q = (x, z, t)
			x - lateral position
			z - vertical position
			t - body orientation
"""
struct Egg{T}
    # dimensions
    nq::Int # generalized coordinates
    nu::Int # controls
    nw::Int # parameters

    # parameters
    mass_body::T
    inertia_body::T 
    radius_a::T 
    radius_b::T 

    # environment
    friction_body_world::T
	gravity::T
end

# Lagrangian
function lagrangian(model::Egg, q, q̇) 
    L = 0.0 

    # body kinetic energy
    v_body = q̇[1:2] 
    ω_body = q̇[3]
    L += 0.5 * model.mass_body * dot(v_body, v_body) 
    L += 0.5 * model.inertia_body * ω_body^2.0 

    # body potential energy
    z_body = q[2]
    L -= model.mass_body * model.gravity * z_body 

    return L
end

# signed distance
function signed_distance(model::Egg, q)
    [
     q[2] - max(model.radius_a, model.radius_b)  # body-world distance
    ]
end

function signed_distance(model::Egg, q, pt1, pt2)
    # error between egg and point on nearest surface
    e1 = pt1 - pt2
    [
     sqrt(dot(e1, e1)) # body-world distance
    ]
end

# input matrix
function input_jacobian(model, q)
	[
     0.0; 
     0.0; 
     1.0
    ]
end

# dimensions
nq = 3
nu = 1
nw = 0

# parameters
mass_body = 1.0 
inertia_body = 0.1 
radius_a = 0.1 
radius_b = 0.05

# environment
friction_body_world = 0.5 
gravity = 9.81 

# model
model = Egg(nq, nu, nw,
            mass_body, inertia_body, 
            radius_a, radius_b,
            friction_body_world, 
            gravity)
			   
# code-gen

nq = model.nq
nu = model.nu
nw = model.nw

# variables
@variables q[1:nq] q̇[1:nq]

# Lagrangian
L = lagrangian(model, q, q̇)

ddL = Symbolics.hessian(L, [q; q̇])#, simplify=true)
dLq = Symbolics.gradient(L, q)#, simplify=true)
ddLq̇q = ddL[nq .+ (1:nq), 1:nq]

M = ddL[nq .+ (1:nq), nq .+ (1:nq)]
C = ddLq̇q * q̇ - dLq

mass_matrix_sym = eval(Symbolics.build_function(M, q)[1])
dynamics_bias_sym = eval(Symbolics.build_function(C, q, q̇)[1])

ϕ = signed_distance(model, q)
dϕ = Symbolics.jacobian(ϕ, q)
signed_distance_jacobian = eval(Symbolics.build_function(dϕ, q)[1])

function contact_jacobian(model, q)
    [
     1.0 0.0 0.0; # note that body torque due to friction is added in residual
     0.0 1.0 0.0
    ] 
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
		+ input_jacobian(model, qm2) * u1[1]
        + λ1)
end

@variables h[1:1] q0[1:nq] q1[1:nq] u1[1:nu] w1[1:nw] λ1[1:3] q2[1:nq]
dyn = dynamics(model, h, q0, q1, u1, w1, λ1, q2)

function residual(model, z, θ, η, μ)
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
    h = θ[2nq + nu + nw + 1 .+ (1:1)] 

    # unpack contact points and directions
    pt_contact_body = η[1:2] 
    pt_contact_surface = η[3:4]
    dir_contact = η[5:6] 

    # unpack variables
    q2 = z[1:nq] 
    γ1 = z[nq .+ (1:1)] 
    sγ1 = z[nq + 1 .+ (1:1)]
    ψ1 = z[nq + 1 + 1 .+ (1:1)] 
    b1 = z[nq + 1 + 1 + 1 .+ (1:1)] 
    sψ1 = z[nq + 1 + 1 + 1 + 1 .+ (1:1)] 
    sb1 = z[nq + 1 + 1 + 1 + 1 + 1 .+ (1:1)]   

    rotation_body = rotation(dir_contact) 
    r_contact = [pt_contact_body - q2[1:2]; 0.0]

    # tangential velocity
    v1 = (q2 - q1) / h[1] 
    vT_body = (rotation_body * (v1[1:2] + cross([0.0; 0.0; v1[3]], r_contact)[1:2]))[1]

    # contact forces
    J = contact_jacobian(model, q2)
    λ1 = transpose(J) * transpose(rotation_body) * [b1[1]; γ1[1]]
                        
    λ1[3] += cross(r_contact, [transpose(rotation_body) * [b1[1]; 0.0]; 0.0])[3] # friction on body creates a moment

	[
     dynamics(model, h, q0, q1, u1, w1, λ1, q2);
	 signed_distance(model, q2, pt_contact_body, pt_contact_surface) - sγ1;
     vT_body - sb1[1];
     friction_body_world[1] * γ1[1] - ψ1[1];
     γ1 .* sγ1 .- μ[1];
	 cone_product([ψ1[1]; b1[1]], [sψ1[1]; sb1[1]]) - [μ[1]; 0.0]; 
    ]
end

nz = nq + 1 + 1 + 1 + 1 + 1 + 1
ny = 1 + 1 + 1 + 1 + 1 + 1
nθ = 2nq + nu + nw + 1 + 1 
nη = 2 + 2 + 2

@variables z[1:nz] θ[1:nθ] η[1:nη] μ[1:1]
r = residual(model, z, θ, η, μ)
rz = Symbolics.jacobian(r, z)
rη = Symbolics.jacobian(r, η)
rθ = Symbolics.jacobian(r, θ)
rz_sp = Symbolics.sparsejacobian(r, z)
rη_sp = Symbolics.sparsejacobian(r, η)
rθ_sp = Symbolics.sparsejacobian(r, θ)

r_Egg! = eval(Symbolics.build_function(r, z, θ, η, μ)[2])
rz_Egg! = eval(Symbolics.build_function(rz, z, θ, η)[2])
rη_Egg! = eval(Symbolics.build_function(rη, z, θ, η)[2])
rθ_Egg! = eval(Symbolics.build_function(rθ, z, θ, η)[2])

r_ = eval(Symbolics.build_function(r, z, θ, η, μ)[1])
rz_ = eval(Symbolics.build_function(rz, z, θ, η)[1])
rη_ = eval(Symbolics.build_function(rη, z, θ, η)[1])
rθ_ = eval(Symbolics.build_function(rθ, z, θ, η)[1])
rz_sp! = eval(Symbolics.build_function(rz_sp.nzval, z, θ, η)[2])
rη_sp! = eval(Symbolics.build_function(rη_sp.nzval, z, θ, η)[2])
rθ_sp! = eval(Symbolics.build_function(rθ_sp.nzval, z, θ, η)[2])

r0 = zeros(length(r))
rz0 = zeros(size(rz))
rθ0 = zeros(size(rθ))
rη0 = zeros(size(rη))

rz0_sp = zeros(size(rz_sp.nzval))
rθ0_sp = zeros(size(rθ_sp.nzval))
rη0_sp = zeros(size(rη_sp.nzval))

z0 = rand(nz) 
θ0 = rand(nθ) 
η0 = rand(nη) 

μ0 = [1.0]

# using BenchmarkTools
# @benchmark r!($r0, $z0, $θ0, $η0, $μ0)
# @benchmark rz!($rz0, $z0, $θ0, $η0)
# @benchmark rη!($rη0, $z0, $θ0, $η0)
# @benchmark rθ!($rθ0, $z0, $θ0, $η0)

# @benchmark $r0 .= r_($z0, $θ0, $η0, $μ0)
# @benchmark $rz0 .= rz_($z0, $θ0, $η0)
# @benchmark $rη0 .= rη_($z0, $θ0, $η0)
# @benchmark $rθ0 .= rθ_($z0, $θ0, $η0)

# @benchmark rz_sp!($rz0_sp, $z0, $θ0, $η0)
# @benchmark rη_sp!($rη0_sp, $z0, $θ0, $η0)
# @benchmark rθ_sp!($rθ0_sp, $z0, $θ0, $η0)

## 
q0 = [0.0, 1.0, 0.0 * π]
z0[1:nq] = copy(q0)
z0[nq .+ (1:1)] .= 1.0
z0[nq + 1 .+ (1:1)] .= 1.0
z0[nq + 1 + 1 .+ (1:1)] .= 1.0
z0[nq + 1 + 1 + 1 .+ (1:1)] .= 0.1
z0[nq + 1 + 1 + 1 + 1 .+ (1:1)] .= 1.0
z0[nq + 1 + 1 + 1 + 1 + 1 .+ (1:1)] .= 0.1

k_body = kinematics_body(model, q0) 
k_foot = kinematics_foot(model, q0)

p_body = contact_point(surf, k_body, r=2.0 * model.body_radius)
p_foot = contact_point(surf, k_foot, r=2.0 * model.body_radius)

η0 = [p_body[1]; p_foot[1]; p_body[3]; p_foot[3]]
rz0 = rz_(z0, θ0, η0)

r1 = rη_(z0, θ0, η0)[:, 1:2] * p_body[6] * kinematics_body_jacobian(model, q0)
r2 = rη_(z0, θ0, η0)[:, 2 .+ (1:2)] * p_foot[6] * kinematics_foot_jacobian(model, q0)
r3 = rη_(z0, θ0, η0)[:, 4 .+ (1:2)] * p_body[7] * p_body[6] * kinematics_body_jacobian(model, q0)

rη_update = r1 + r2 + r3

rz0[:, 1:nq] += rη_update

function res!(r, z, θ, μ) 
    q2 = z[1:nq]
    k_body = kinematics_body(model, q2) 
    k_foot = kinematics_foot(model, q2)

    p_body = contact_point(surf, k_body, r=2.0 * model.body_radius)
    p_foot = contact_point(surf, k_foot, r=2.0 * model.foot_radius)

    η = [p_body[1]; p_foot[1]; p_body[3]; p_foot[3]]

    r_Egg!(r, z, θ, η, μ)
end

function res_z!(rz, z, θ) 
    q2 = z[1:nq]
    k_body = kinematics_body(model, q2) 
    k_foot = kinematics_foot(model, q2)

    p_body = contact_point(surf, k_body, r=2.0 * model.body_radius)
    p_foot = contact_point(surf, k_foot, r=2.0 * model.foot_radius)

    # @show p_body[1] 
    # @show p_foot[1]

    η = [p_body[1]; p_foot[1]; p_body[3]; p_foot[3]]

    rz_Egg!(rz, z, θ, η)
    rη = rη_(z, θ, η)

    r1 = rη[:, 1:2] * p_body[6] * kinematics_body_jacobian(model, q2)
    r2 = rη[:, 2 .+ (1:2)] * p_foot[6] * kinematics_foot_jacobian(model, q2)
    r3 = rη[:, 4 .+ (1:2)] * p_body[7] * p_body[6] * kinematics_body_jacobian(model, q2)
    r4 = rη[:, 6 .+ (1:2)] * p_foot[7] * p_foot[6] * kinematics_foot_jacobian(model, q2)

    rη_update = r1 + r2 + r3 + r4

    rz[:, 1:nq] += rη_update
end

# res(z0, θ0, μ0)
# res_z(z0, θ0, μ0)

# res!(r0, z0, θ0, μ0)
# res_z!(rz0, z0, θ0)

# r0
