""" 
	planar biped 
		model inspired by Spring Flamingo 
			from "Exploiting Inherent Robustness and Natural Dynamics in the Control of Bipedal Walking Robots"
	q = (x, z, t, q1, ..., q6)
		x - lateral position 
		z - vertical position 
		t - body orientation 
		q1 - leg 1 thigh angle (absolute) 
		q2 - leg 1 calf angle  (absolute) 
		q3 - leg 1 foot angle  (absolute) 
		q4 - leg 2 thigh angle (absolute) 
		q5 - leg 2 calf angle  (absolute)
		q6 - leg 2 foot angle  (absolute) 
"""
mutable struct Biped{T}
    # dimensions
    nq::Int # generalized coordinates
    nu::Int # controls
    nw::Int # parameters
	nc::Int # contact points

    # torso
    l_torso::T
    lc_torso::T
    m_torso::T
    J_torso::T

    # leg 1
        # thigh
    l_thigh1::T
    lc_thigh1::T
    m_thigh1::T
    J_thigh1::T

        # calf
    l_calf1::T
    lc_calf1::T
    m_calf1::T
    J_calf1::T

		# foot
	lt_foot1::T # toe length
	lh_foot1::T # heel length
	m_foot1::T
	J_foot1::T

    # leg 2
        # thigh
    l_thigh2::T
    lc_thigh2::T
    m_thigh2::T
    J_thigh2::T

        # calf
    l_calf2::T
    lc_calf2::T
    m_calf2::T
    J_calf2::T

		# foot
	lt_foot2::T # toe length
	lh_foot2::T # heel length
	m_foot2::T
	J_foot2::T

	 # joint friction 
	 friction_joint::T 
	 # TODO: individual joint frictions
 
	 # environment 
	 friction_body_world::T 
	 friction_foot_world::T 
	 gravity::T
end

function kinematics_body(model::Biped, q; mode=:none) 
	x = q[1] 
	z = q[2] 
	t = q[3] 

	if mode == :ee 
		l = model.l_torso 
	elseif mode == :com 
		l = model.lc_torso 
	elseif mode == :hip 
		l = 0.0
	else 
		@error "incorrect mode specified"
	end

	return [
			x - l * sin(t);
			z + l * cos(t);
		   ]
end

function kinematics_jacobian_body(model::Biped, q; mode=:none) 
	x = q[1] 
	z = q[2] 
	t = q[3] 

	if mode == :ee 
		l = model.l_torso 
	elseif mode == :com 
		l = model.lc_torso 
	elseif mode == :hip 
		l = 0.0
	else 
		@error "incorrect mode specified"
	end

	jac = zeros(eltype(q), 2, model.nq)

	jac[1, 1] = 1.0 
	jac[1, 3] = -l * cos(t) 
	jac[2, 2] = 1.0 
	jac[2, 3] = -l * sin(t)

	return jac
end

function kinematics_thigh(model::Biped, q; leg=:none, mode=:none)
	x = q[1] 
	z = q[2] 
	t = q[3]

	if leg == :leg1 
		le = model.l_thigh1 
		lc = model.lc_thigh1
		idx = 4
	elseif leg == :leg2 
		le = model.l_thigh2
		lc = model.lc_thigh2
		idx = 7
	else
		@error "incorrect leg specified"
	end

	if mode == :ee 
		l = le 
	elseif mode == :com 
		l = lc 
	else 
		@error "incorrect mode specified"
	end


	return [
			x + l * sin(q[idx]);
			z - l * cos(q[idx]);
		   ]
end

function kinematics_jacobian_thigh(model::Biped, q; leg=:none, mode=:none)
	x = q[1] 
	z = q[2] 
	t = q[3]

	if leg == :leg1 
		le = model.l_thigh1 
		lc = model.lc_thigh1
		idx = 4
	elseif leg == :leg2 
		le = model.l_thigh2
		lc = model.lc_thigh2
		idx = 7
	else
		@error "incorrect leg specified"
	end

	if mode == :ee 
		l = le 
	elseif mode == :com 
		l = lc 
	else 
		@error "incorrect mode specified"
	end

	jac = zeros(eltype(q), 2, model.nq)
	jac[1, 1] = 1.0 
	jac[1, idx] = l * cos(q[idx])
	jac[2, 2] = 1.0 
	jac[2, idx] = l * sin(q[idx]) 

	return jac
end

function kinematics_calf(model::Biped, q; leg=:none, mode=:none)
	x = q[1] 
	z = q[2] 
	t = q[3]

	if leg == :leg1 
		le = model.l_calf1 
		lc = model.lc_calf1
		idx = 5
	elseif leg == :leg2 
		le = model.l_calf2
		lc = model.lc_calf2
		idx = 8
	else
		@error "incorrect leg specified"
	end

	if mode == :ee 
		l = le 
	elseif mode == :com 
		l = lc 
	else 
		@error "incorrect mode specified"
	end

	k_thigh = kinematics_thigh(model, q, leg=leg, mode=:ee)
	return k_thigh + [l * sin(q[idx]); -l * cos(q[idx])]
end

function kinematics_jacobian_calf(model::Biped, q; leg=:none, mode=:none)
	x = q[1] 
	z = q[2] 
	t = q[3]

	if leg == :leg1 
		le = model.l_calf1 
		lc = model.lc_calf1
		idx = 5
	elseif leg == :leg2 
		le = model.l_calf2
		lc = model.lc_calf2
		idx = 8
	else
		@error "incorrect leg specified"
	end

	if mode == :ee 
		l = le 
	elseif mode == :com 
		l = lc 
	else 
		@error "incorrect mode specified"
	end

	jac = kinematics_jacobian_thigh(model, q; leg=leg, mode=:ee) 

	jac[1, idx] += l * cos(q[idx])
	jac[2, idx] += l * sin(q[idx]) 

	return jac
end

function kinematics_foot(model::Biped, q; leg=:none, mode=:none)
	x = q[1] 
	z = q[2] 
	t = q[3]

	if leg == :leg1 
		lt = model.lt_foot1 
		lh = model.lh_foot1
		idx = 6
	elseif leg == :leg2 
		lt = model.lt_foot2
		lh = model.lh_foot2
		idx = 9
	else
		@error "incorrect leg specified"
	end

	if mode == :toe
		l = lt
	elseif mode == :heel 
		l = -lh 
	elseif mode == :com 
		l = 0.0
	else 
		@error "incorrect mode specified"
	end

	k_calf = kinematics_calf(model, q, leg=leg, mode=:ee)
	return k_calf + [l * cos(q[idx]); l * sin(q[idx])]
end

function kinematics_jacobian_foot(model::Biped, q; leg=:none, mode=:none)
	x = q[1] 
	z = q[2] 
	t = q[3]

	if leg == :leg1 
		lt = model.lt_foot1 
		lh = model.lh_foot1
		idx = 6
	elseif leg == :leg2 
		lt = model.lt_foot2
		lh = model.lh_foot2
		idx = 9
	else
		@error "incorrect leg specified"
	end

	if mode == :toe
		l = lt
	elseif mode == :heel 
		l = -lh 
	elseif mode == :com 
		l = 0.0
	else 
		@error "incorrect mode specified"
	end

	jac = kinematics_jacobian_calf(model, q, leg=leg, mode=:ee)
	jac[1, idx] += -l * sin(q[idx])
	jac[2, idx] += l * cos(q[idx])

	return jac
end

function kinematics(model::Biped, q)
	# foot
	p_toe_1 = kinematics_foot(model, q, leg=:leg1, mode=:toe) 
	p_heel_1 = kinematics_foot(model, q, leg=:leg1, mode=:heel) 
	p_toe_2 = kinematics_foot(model, q, leg=:leg2, mode=:toe) 
	p_heel_2 = kinematics_foot(model, q, leg=:leg2, mode=:heel) 

	# knee 
	p_knee_1 = kinematics_thigh(model, q, leg=:leg1, mode=:ee) 
	p_knee_2 = kinematics_thigh(model, q, leg=:leg2, mode=:ee) 

	# body
	p_body_1 = q[1:2] 
	p_body_2 = kinematics_body(model, q, mode=:ee)

	[
	 p_toe_1; p_heel_1; p_toe_2; p_heel_2;
	 p_knee_1; p_knee_2;
	 p_body_1; p_body_2;
	]
end

# Lagrangian
function lagrangian(model::Biped, q, q̇)
	L = 0.0

	# torso
	p_torso = kinematics_body(model, q, mode=:com)
	J_torso = kinematics_jacobian_body(model, q, mode=:com)
	v_torso = J_torso * q̇

	L += 0.5 * model.m_torso * transpose(v_torso) * v_torso
	L += 0.5 * model.J_torso * q̇[3]^2.0
	L -= model.m_torso * model.gravity * p_torso[2]

	# thigh 1
	p_thigh_1 = kinematics_thigh(model, q, leg=:leg1, mode=:com)
	J_thigh_1 = kinematics_jacobian_thigh(model, q, leg=:leg1, mode=:com)
	v_thigh_1 = J_thigh_1 * q̇

	L += 0.5 * model.m_thigh1 * transpose(v_thigh_1) * v_thigh_1
	L += 0.5 * model.J_thigh1 * q̇[4]^2.0
	L -= model.m_thigh1 * model.gravity * p_thigh_1[2]

	# calf 1
	p_calf_1 = kinematics_calf(model, q, leg=:leg1, mode=:com)
	J_calf_1 = kinematics_jacobian_calf(model, q, leg=:leg1, mode=:com)
	v_calf_1 = J_calf_1 * q̇

	L += 0.5 * model.m_calf1 * transpose(v_calf_1) * v_calf_1
	L += 0.5 * model.J_calf1 * q̇[5]^2.0
	L -= model.m_calf1 * model.gravity * p_calf_1[2]

	# foot 1
	p_foot_1 = kinematics_foot(model, q, leg=:leg1, mode=:com)
	J_foot_1 = kinematics_jacobian_foot(model, q, leg=:leg1, mode=:com)
	v_foot_1 = J_foot_1 * q̇

	L += 0.5 * model.m_foot1 * transpose(v_foot_1) * v_foot_1
	L += 0.5 * model.J_foot1 * q̇[6]^2.0
	L -= model.m_foot1 * model.gravity * p_foot_1[2]

	# thigh 2
	p_thigh_2 = kinematics_thigh(model, q, leg=:leg2, mode=:com)
	J_thigh_2 = kinematics_jacobian_thigh(model, q, leg=:leg2, mode=:com)
	v_thigh_2 = J_thigh_2 * q̇

	L += 0.5 * model.m_thigh2 * transpose(v_thigh_2) * v_thigh_2
	L += 0.5 * model.J_thigh2 * q̇[7]^2.0
	L -= model.m_thigh2 * model.gravity * p_thigh_2[2]

	# calf 2
	p_calf_2 = kinematics_calf(model, q, leg=:leg2, mode=:com)
	J_calf_2 = kinematics_jacobian_calf(model, q, leg=:leg2, mode=:com)
	v_calf_2 = J_calf_2 * q̇

	L += 0.5 * model.m_calf2 * transpose(v_calf_2) * v_calf_2
	L += 0.5 * model.J_calf2 * q̇[8]^2.0
	L -= model.m_calf2 * model.gravity * p_calf_2[2]

	# foot 2
	p_foot_2 = kinematics_foot(model, q, leg=:leg2, mode=:com)
	J_foot_2 = kinematics_jacobian_foot(model, q, leg=:leg2, mode=:com)
	v_foot_2 = J_foot_2 * q̇

	L += 0.5 * model.m_foot2 * transpose(v_foot_2) * v_foot_2
	L += 0.5 * model.J_foot2 * q̇[9]^2.0
	L -= model.m_foot2 * model.gravity * p_foot_2[2]

	return L
end

function signed_distance(model::Biped, q)
	k = kinematics(model, q) 

	return [k[2], k[4], k[6], k[8], k[10], k[12], k[14], k[16]]
end

function input_jacobian()
	[0.0  0.0 -1.0  1.0  0.0  0.0  0.0  0.0  0.0;
	 0.0  0.0  0.0 -1.0  1.0  0.0  0.0  0.0  0.0;
	 0.0  0.0  0.0  0.0 -1.0  1.0  0.0  0.0  0.0;
	 0.0  0.0 -1.0  0.0  0.0  0.0  1.0  0.0  0.0;
	 0.0  0.0  0.0  0.0  0.0  0.0 -1.0  1.0  0.0;
	 0.0  0.0  0.0  0.0  0.0  0.0  0.0 -1.0  1.0]
end

function contact_jacobian(model::Biped, q)
	# foot
	J_toe_1 = kinematics_jacobian_foot(model, q, leg=:leg1, mode=:toe) 
	J_heel_1 = kinematics_jacobian_foot(model, q, leg=:leg1, mode=:heel) 
	J_toe_2 = kinematics_jacobian_foot(model, q, leg=:leg2, mode=:toe) 
	J_heel_2 = kinematics_jacobian_foot(model, q, leg=:leg2, mode=:heel) 

	# knee 
	J_knee_1 = kinematics_jacobian_thigh(model, q, leg=:leg1, mode=:ee) 
	J_knee_2 = kinematics_jacobian_thigh(model, q, leg=:leg2, mode=:ee) 

	# body
	J_body_1 = kinematics_jacobian_body(model, q, mode=:hip)
	J_body_2 = kinematics_jacobian_body(model, q, mode=:ee)

	[
	 J_toe_1; 
	 J_heel_1; 
	 J_toe_2; 
	 J_heel_2;
	 J_knee_1; 
	 J_knee_2;
	 J_body_1; 
	 J_body_2;
	]
end

# dimensions
nq = 3 + 2 * 3            # generalized coordinates dimension
nu = 2 * 3                # control dimension
nc = 8                    # number of contact points
nw = 0                    # parameter dimension

# model parameters
m_torso = 12.0
m_thigh = 0.4598
m_calf = 0.306
m_foot = 0.3466

l_torso = 0.385
l_thigh = 0.42
l_calf = 0.45
l_toe = 0.1725 # distance ankle-toe
l_heel = 0.0525 # distance ankle-heel

lc_torso = 0.20
lc_thigh = 0.5 * l_thigh
lc_calf = 0.5 * l_calf 

J_torso = 0.10
J_thigh = 0.01256
J_calf = 0.00952
J_foot = 0.0015

friction_joint = 0.1 

# world parameters
friction_body_world = 0.5 
friction_foot_world = 0.5 
gravity = 9.81

biped = Biped(nq, nu, nw, nc,
			  l_torso, lc_torso, m_torso, J_torso,
			  l_thigh, lc_thigh, m_thigh, J_thigh,
			  l_calf, lc_calf, m_calf, J_calf,
			  l_toe, l_heel, m_foot, J_foot,
			  l_thigh, lc_thigh, m_thigh, J_thigh,
			  l_calf, lc_calf, m_calf, J_calf,
			  l_toe, l_heel, m_foot, J_foot,
			  friction_joint,
			  friction_body_world,
			  friction_foot_world,
			  gravity)

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
		+ transpose(input_jacobian()) * u1 
        - model.friction_joint * [0.0; 0.0; 0.0; vm2[4:9]] # joint friction
        + λ1)
end

# @variables h[1:1] q0[1:nq] q1[1:nq] u1[1:nu] w1[1:nw] λ1[1:nq] q2[1:nq]
# d = dynamics(biped, h, q0, q1, u1, w1, λ1, q2)

function residual(model, z, θ, μ)
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
    friction_body_world = θ[2nq + nu + nw .+ (1:1)] 
    friction_foot_world = θ[2nq + nu + nw + 1 .+ (1:1)] 
    h = θ[2nq + nu + nw + 2 .+ (1:1)] 

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
     friction_foot_world[1] * γ1[2] - ψ1[2];
     friction_foot_world[1] * γ1[3] - ψ1[3];
     friction_foot_world[1] * γ1[4] - ψ1[4];
     friction_body_world[1] * γ1[5] - ψ1[5];
     friction_body_world[1] * γ1[6] - ψ1[6];
     friction_body_world[1] * γ1[7] - ψ1[7];
     friction_body_world[1] * γ1[8] - ψ1[8];
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

nz = nq + nc + nc + nc + nc + nc + nc 
ny = nc + nc + nc + nc + nc + nc
nθ = 2nq + nu + nw + 2 + 1 

@variables z[1:nz] θ[1:nθ] μ[1:1]
r = residual(biped, z, θ, μ)
rz = Symbolics.jacobian(r, z)
rθ = Symbolics.jacobian(r, θ)
rz_sp = Symbolics.sparsejacobian(r, z)
rθ_sp = Symbolics.sparsejacobian(r, θ)

r_biped! = eval(Symbolics.build_function(r, z, θ, μ)[2])
rz_biped! = eval(Symbolics.build_function(rz, z, θ)[2])
rθ_biped! = eval(Symbolics.build_function(rθ, z, θ)[2])

# rz_sp! = eval(Symbolics.build_function(rz_sp.nzval, z, θ)[2])
# rθ_sp! = eval(Symbolics.build_function(rθ_sp.nzval, z, θ)[2])

# r0 = zeros(length(r))
# rz0 = zeros(size(rz))
# rθ0 = zeros(size(rθ))

# rz0_sp = zeros(size(rz_sp.nzval))
# rθ0_sp = zeros(size(rθ_sp.nzval))

# z0 = rand(nz) 
# θ0 = rand(nθ) 

# μ0 = [1.0]

# r_biped!(r0, z0, θ0, μ0)
# rz_biped!(rz0, z0, θ0)
# rθ_biped!(rθ0, z0, θ0)
# rz_sp!(rz0_sp, z0, θ0)
# rθ_sp!(rθ0_sp, z0, θ0)

# using BenchmarkTools
# @benchmark r_biped!($r0, $z0, $θ0, $μ0)
# @benchmark rz_biped!($rz0, $z0, $θ0)
# @benchmark rθ_biped!($rθ0, $z0, $θ0)
# @benchmark rz_sp!($rz0_sp, $z0, $θ0)
# @benchmark rθ_sp!($rθ0_sp, $z0, $θ0)