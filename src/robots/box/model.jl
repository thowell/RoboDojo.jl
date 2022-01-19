"""
    box
        -fully actuated
"""
struct Box{T} <: Model{T}
	# dimensions
	nq::Int # generalized coordinates
	nu::Int # controls
	nw::Int # parameters
	nc::Int # contact points

    mass::T
    inertia::T
    gravity::T

    friction_body_world::Vector{T}
    contact_corner_offset::Vector{Vector{T}}
    friction_joint::Vector{T} 
end

# Methods
function lagrangian(model::Box, q, q̇) 
    L = 0.0 
    L += 0.5 * model.mass * dot(q̇[1:2], q̇[1:2]) 
    L -= model.mass * model.gravity * q[2] 

    L += 0.5 * model.inertia * q̇[3]^2

    return L 
end

function kinematics(model::Box, q)
    pos = q[1:2]
    θ = q[3]
    R = rotation_matrix(θ)
    [
        pos + R * model.contact_corner_offset[1]; 
        pos + R * model.contact_corner_offset[2];
        pos + R * model.contact_corner_offset[3];
        pos + R * model.contact_corner_offset[4];
    ]
end

function signed_distance(model::Box, q)
    kinematics(model, q)[collect([2, 4, 6, 8])]
end

function input_jacobian(model::Box, q)
	[
	 1.0 0.0 0.0;
	 0.0 1.0 0.0;
     0.0 0.0 1.0;
    ]
end

function contact_jacobian(model::Box, q)
	j = kinematics(model, q)
	Symbolics.jacobian(j, q)
end

# nominal configuration 
function nominal_configuration(model::Box)
	[0.0; 0.1; 0.0]
end

# friction coefficients 
friction_coefficients(model::Box) = model.friction_body_world

# Dimensions
nq = 3 # configuration dimension
nu = 3 # control dimension
nw = 0 # parameters
nc = 4 # number of contact points

# Parameters
r_dim = 0.5
gravity = 10.0#9.81
mass = 1.0 # mass
inertia = 1.0 / 12.0 * mass * ((2.0 * r_dim)^2 + (2.0 * r_dim)^2)

# Contact corner
cc1 = [r_dim, r_dim]
cc2 = [-r_dim, r_dim]
cc3 = [r_dim, -r_dim]
cc4 = [-r_dim, -r_dim]

contact_corner_offset = [cc1, cc2, cc3, cc4]
friction_body_world = [1.0, 1.0, 1.0, 1.0]  # coefficient of friction

# Model
box = Box(nq, nu, nw, nc,
		  mass, inertia, gravity,
		  friction_body_world, contact_corner_offset, zeros(0))

box_contact_kinematics = [q -> kinematics(box, q)[1:2], 
                          q -> kinematics(box, q)[3:4],
                          q -> kinematics(box, q)[5:6], 
                          q -> kinematics(box, q)[7:8]]
    
box_contact_kinematics_jacobians = [q -> contact_jacobian(box, q)[1:2, :], 
                                    q -> contact_jacobian(box, q)[3:4, :],
                                    q -> contact_jacobian(box, q)[5:6, :], 
                                    q -> contact_jacobian(box, q)[7:8, :]]

name(::Box) = :box