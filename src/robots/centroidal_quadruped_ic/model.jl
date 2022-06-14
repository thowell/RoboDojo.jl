"""
    centroidal quadruped 
    q = (p, r, f1, f2, f3, f4) 
        p - body position
        r - body orientation (modified Rodriques parameters)
        f1 - foot 1 position 
        f2 - foot 2 position 
        f3 - foot 3 position 
        f4 - foot 4 position 
"""
mutable struct CentroidalQuadrupedIC{T} <: Model{T}
    # dimensions
	nq::Int # generalized coordinates 
    nu::Int # controls 
    nw::Int # parameters
    nc::Int # contact points

    # parameters
    mass_body::T 
    inertia_body::Matrix{T} 
    mass_foot::T 

    # environment 
    friction_joint::Vector{T}
    friction_body_world::Vector{T} 
    friction_foot_world::Vector{T} 
    gravity::Vector{T}
end

# function skew(x)
#     return [0.0  -x[3]  x[2];
#             x[3]   0.0 -x[1];
#            -x[2]  x[1]   0.0]
# end

# function L_mult(x)
#     [x[1] -transpose(x[2:4]); 
#      x[2:4] x[1] * I(3) + skew(x[2:4])]
# end

# # right quaternion multiply as matrix
# function R_mult(x)
#     [x[1] -transpose(x[2:4]); x[2:4] x[1] * I(3) - skew(x[2:4])]
# end

# # rotation matrix
# function quaternion_rotation_matrix(q) 
#     H = [zeros(1, 3); I(3)]
#     transpose(H) * L_mult(q) * transpose(R_mult(q)) * H
# end

# function quaternion_from_mrp(p)
#     """Quaternion (scalar first) from MRP"""
#     return (1.0 / (1.0 + dot(p, p))) * [(1 - dot(p, p)); 2.0 * p]
# end

# function mrp_rotation_matrix(x) 
#     quaternion_rotation_matrix(quaternion_from_mrp(x))
# end

# #https://en.wikipedia.org/wiki/Rotation_matrix#In_three_dimensions

# function euler_rotation_matrix(θ)
# 	a = θ[1]
# 	b = θ[2]
# 	c = θ[3]

# 	[cos(a) * cos(b) (cos(a) * sin(b) * sin(c) - sin(a) * cos(c)) (cos(a) * sin(b) * cos(c) + sin(a) * sin(c));
# 	 sin(a) * cos(b) (sin(a) * sin(b) * sin(c) + cos(a) * cos(c)) (sin(a) * sin(b) * cos(c) - cos(a) * sin(c));
# 	 -sin(b) cos(b) * sin(c) cos(b) * cos(c)]
# end

function mass_matrix(model::CentroidalQuadrupedIC, q)
    cat(
        model.mass_body * Diagonal(ones(3)),     # body position
        model.inertia_body,                      # body orienation
        model.mass_foot * Diagonal(ones(3 * 4)), # feet position
        dims=(1, 2)
        )
end

function dynamics_bias(model::CentroidalQuadrupedIC, q, q̇)
    [
        model.mass_body * model.gravity;            # body position
        skew(q̇[4:6]) * model.inertia_body * q̇[4:6]; # body orienation 
        model.mass_foot * model.gravity;
        model.mass_foot * model.gravity;
        model.mass_foot * model.gravity;
        model.mass_foot * model.gravity;
    ]
end

function signed_distance(model::CentroidalQuadrupedIC, q)

    position_body = q[1:3] 
    orientation_body = q[3 .+ (1:3)]

    position_foot1 = q[6 .+ (1:3)]
    position_foot2 = q[9 .+ (1:3)]
    position_foot3 = q[12 .+ (1:3)]
	position_foot4 = q[15 .+ (1:3)]

    return [position_foot1[3]; position_foot2[3]; position_foot3[3]; position_foot4[3]; position_body[3]]
end

function input_jacobian(model::CentroidalQuadrupedIC, q)
    position_body = q[1:3]
    orientation_body = q[3 .+ (1:3)]
    # R = mrp_rotation_matrix(orientation_body)
    R = euler_rotation_matrix(orientation_body)
    
    # kinematics in world frame
	r1 = q[6 .+ (1:3)] - position_body
	r2 = q[9 .+ (1:3)] - position_body
	r3 = q[12 .+ (1:3)] - position_body
	r4 = q[15 .+ (1:3)] - position_body

	z3 = zeros(3, 3)

	transpose([
        I(3) z3 I(3) I(3) I(3) I(3);
        z3 I(3) transpose(R) * skew(r1) transpose(R) * skew(r2) transpose(R) * skew(r3) transpose(R) * skew(r4);
        z3 z3 -I(3)    z3    z3   z3;
        z3 z3 z3    -I(3)    z3   z3;
        z3 z3 z3       z3 -I(3)   z3;
        z3 z3 z3       z3    z3 -I(3)
    ])
end

function contact_jacobian(model::CentroidalQuadrupedIC, q) 
    z3 = zeros(3, 3)

    [
        z3   z3 I(3)   z3   z3   z3;
        z3   z3   z3 I(3)   z3   z3;
        z3   z3   z3   z3 I(3)   z3;
        z3   z3   z3   z3   z3 I(3);
        I(3) z3   z3   z3   z3   z3;
    ]
end

# nominal configuration 
function nominal_configuration(model::CentroidalQuadrupedIC) 
    [
        0.0; 0.0; 0.15; 
        0.0; 0.0; 0.0;
        0.1; 0.1; 0.0;
        0.1;-0.1; 0.0;
       -0.1; 0.1; 0.0;
       -0.1;-0.1; 0.0;
    ]
end

# friction coefficients 
function friction_coefficients(model::CentroidalQuadrupedIC) 
    return [model.friction_foot_world; model.friction_body_world]
end

# dimensions
nq = 3 + 3 + 3 * 4       # generalized coordinates
nu = 3 * 4 + 6           # controls
nw = 0                   # parameters
nc = 5                   # contact points 

# parameters
gravity = [0.0; 0.0; 9.81]                 # gravity
friction_body_world = [0.5]                # coefficient of friction
friction_foot_world = [0.5; 0.5; 0.5; 0.5] # coefficient of friction

# inertial properties
mass_body = 1.0 
inertia_body = Array(Diagonal(ones(3)))
mass_foot = 0.1

centroidal_quadruped_ic = CentroidalQuadrupedIC(nq, nu, nw, nc,
				mass_body,
                inertia_body, 
                mass_foot,
                zeros(nq - 6),
                friction_body_world, 
                friction_foot_world, 
                gravity)

centroidal_quadruped_ic_contact_kinematics = [
    q -> q[6  .+ (1:3)],
    q -> q[9  .+ (1:3)],
    q -> q[12 .+ (1:3)],
    q -> q[15 .+ (1:3)],
    q -> q[ 0 .+ (1:3)],
]

centroidal_quadruped_ic_contact_kinematics_jacobians = [
    q -> [zeros(3, 6) I(3) zeros(3, 9)],
    q -> [zeros(3, 9) I(3) zeros(3, 6)],
    q -> [zeros(3, 12) I(3) zeros(3, 3)],
    q -> [zeros(3, 15) I(3)],
    q -> [I(3) zeros(3, 15)],
]

name(::CentroidalQuadrupedIC) = :centroidal_quadruped_ic
floating_base_dim(::CentroidalQuadrupedIC) = 6
