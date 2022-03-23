"""
    particle 3D
        -fully actuated
"""
struct Particle{T} <: Model{T}
	# dimensions
	nq::Int # generalized coordinates
	nu::Int # controls
	nw::Int # parameters
	nc::Int # contact points

    mass::T
    gravity::T

    friction_body_world::Vector{T}
    friction_joint::Vector{T} 
end

# Methods
function lagrangian(model::Particle, q, q̇) 
    L = 0.0 
    L += 0.5 * model.mass * dot(q̇[1:3], q̇[1:3]) 
    L -= model.mass * model.gravity * q[3] 

    return L 
end

function kinematics(model::Particle, q)
    q[1:3]
end

function signed_distance(model::Particle, q)
    q[3:3]
end

function input_jacobian(model::Particle, q)
	[
	 1.0 0.0 0.0;
	 0.0 1.0 0.0;
     0.0 0.0 1.0;
    ]
end

function contact_jacobian(model::Particle, q)
	[
	 1.0 0.0 0.0;
	 0.0 1.0 0.0;
     0.0 0.0 1.0;
    ]
end

# nominal configuration 
function nominal_configuration(model::Particle)
	[0.0; 0.0; 1.0]
end

# friction coefficients 
friction_coefficients(model::Particle) = model.friction_body_world

# Dimensions
nq = 3 # configuration dimension
nu = 3 # control dimension
nw = 0 # parameters
nc = 1 # number of contact points

# Parameters
gravity = 10.0#9.81
mass = 1.0 # mass

# Contact corner

friction_body_world = [1.0]  # coefficient of friction

# Model
particle = Particle(nq, nu, nw, nc,
		  mass, gravity,
		  friction_body_world, zeros(0))

particle_contact_kinematics = [q -> kinematics(particle, q)]
    
particle_contact_kinematics_jacobians = [q -> contact_jacobian(particle, q)]

name(::Particle) = :particle
floating_base_dim(::Particle) = 3
