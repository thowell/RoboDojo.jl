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
mutable struct Biped{T} <: Model{T}
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
	friction_joint::Vector{T} 
 
	# environment 
	friction_body_world::Vector{T} 
	friction_foot_world::Vector{T} 
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
	 p_toe_1; 
	 p_heel_1; 
	 p_toe_2; 
	 p_heel_2;
	 p_knee_1; 
	 p_knee_2;
	 p_body_1; 
	 p_body_2;
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

function input_jacobian(::Biped, q)
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

# nominal configuration 
function nominal_configuration(model::Biped)
	[0.0; 1.0; 0.01 * π; -0.5 * π; 0.0 * π; -0.01 * π; 0.5 * π; 0.0 * π; -0.0 * π]
end

# friction_coefficients
function friction_coefficients(model::Biped) 
	return [model.friction_foot_world; model.friction_body_world]
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

friction_joint = [0.1; 0.1; 0.1; 0.1; 0.1; 0.1] # q1-q6

# world parameters
friction_body_world = [0.5; 0.5; 0.5; 0.5] # knee1, knee2, body1, body2
friction_foot_world = [0.5; 0.5; 0.5; 0.5] # foot1, foot2, foot3, foot4
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

biped_contact_kinematics = [q -> kinematics_foot(biped, q, leg=:leg1, mode=:toe), 
	q -> kinematics_foot(biped, q, leg=:leg1, mode=:heel),
	q -> kinematics_foot(biped, q, leg=:leg2, mode=:toe), 
	q -> kinematics_foot(biped, q, leg=:leg2, mode=:heel), 
	q -> kinematics_thigh(biped, q, leg=:leg1, mode=:ee),
	q -> kinematics_thigh(biped, q, leg=:leg2, mode=:ee),
	q -> q[1:2],
	q -> kinematics_body(biped, q, mode=:ee)]

biped_contact_kinematics_jacobians = [q -> kinematics_jacobian_foot(biped, q, leg=:leg1, mode=:toe), 
	q -> kinematics_jacobian_foot(biped, q, leg=:leg1, mode=:heel),
	q -> kinematics_jacobian_foot(biped, q, leg=:leg2, mode=:toe), 
	q -> kinematics_jacobian_foot(biped, q, leg=:leg2, mode=:heel), 
	q -> kinematics_jacobian_thigh(biped, q, leg=:leg1, mode=:ee),
	q -> kinematics_jacobian_thigh(biped, q, leg=:leg2, mode=:ee),
	q -> kinematics_jacobian_body(biped, q, mode=:hip),
	q -> kinematics_jacobian_body(biped, q, mode=:ee)]

name(::Biped) = :biped