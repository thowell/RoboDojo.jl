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
struct Hopper1{T} <: Model{T}
    # dimensions
    nq::Int # generalized coordinates
    nu::Int # controls
    nc::Int # contacts
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
function kinematics_body(::Hopper1, q)
	[q[1],
	 q[2]]
end

function kinematics_body_jacobian(::Hopper1, q) 
    [1.0 0.0 0.0 0.0;
	 0.0 1.0 0.0 0.0] 
end

function kinematics_foot(::Hopper1, q)
	[q[1] + q[4] * sin(q[3]),
	 q[2] - q[4] * cos(q[3])]
end

function kinematics_foot_jacobian(::Hopper1, q) 
    [1.0 0.0 (q[4] * cos(q[3])) sin(q[3]);
	 0.0 1.0 (q[4] * sin(q[3])) (-1.0 * cos(q[3]))] 
end

# Lagrangian
function lagrangian(model::Hopper1, q, q̇) 
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
function signed_distance(model::Hopper1, q)
    [
    #  q[2] - model.body_radius,                    # body-world distance
     kinematics_foot(model, q)[2] - model.foot_radius, # foot-world distance 
    #  q[4] - model.leg_len_min,                    # leg length min.
    #  model.leg_len_max - q[4]                     # leg length max.
    ]
end

function signed_distance(model::Hopper1, q, pt1, pt2)
    # error between robot kinematics and point on nearest surface
    e1 = kinematics_body(model, q) - pt1
    e2 = kinematics_foot(model, q) - pt2
    @show pt1 
    @show pt2
    [
    #  sqrt(dot(e1, e1)) - model.body_radius,  # body-world distance
     dot(e2, e2) - model.foot_radius^2,  # foot-world distance 
    #  q[4] - model.leg_len_min,             # leg length min.
    #  model.leg_len_max - q[4]              # leg length max.
    ]
end

# input matrix
function input_jacobian(::Hopper1, q)
	transpose([
     0.0 -sin(q[3]); 
     0.0 cos(q[3]); 
     1.0 0.0; 
     0.0 1.0
    ])
end

# contact Jacobian
function contact_jacobian(::Hopper1, q)
    [
    #  1.0 0.0 0.0 0.0; # note that body torque due to friction is added in residual
    #  0.0 1.0 0.0 0.0;
     1.0 0.0 (q[4] * cos(q[3])) sin(q[3]);
     0.0 1.0 (q[4] * sin(q[3])) (-1.0 * cos(q[3]));
    #  0.0 0.0 0.0 1.0; 
    #  0.0 0.0 0.0 -1.0;
     ] 
end

# nominal configuratoin
function nominal_configuration(model::Hopper1) 
    return [0.0, 1.0, 0.0 * π, 0.5]
end

function friction_coefficients(model::Hopper1) 
    return [model.friction_body_world; model.friction_foot_world]
end

# dimensions
nq = 4
nu = 2
nc = 1
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
hopper1 = Hopper1(nq, nu, nc, nw,
            mass_body, mass_foot, inertia_body, 
            body_radius, foot_radius, 
            leg_len_max, leg_len_min,
            friction_joint,
            friction_body_world, friction_foot_world, 
            gravity)

# contact kinematics
hopper1_contact_kinematics = [
    # q -> kinematics_body(hopper1, q),
    q -> kinematics_foot(hopper, q),
    ] 

hopper1_contact_kinematics_jacobians = [
    # q -> kinematics_body_jacobian(hopper1, q), 
    q -> kinematics_foot_jacobian(hopper, q),
    ]

name(::Hopper1) = :hopper1