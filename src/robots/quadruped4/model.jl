"""
    planar quadruped 
    q = (x, z, t, q1, ..., q8) 
        x - lateral position
        z - vertical position
        t - body orientation
        q1 - leg 1 shoulder angle (absolute)
        q2 - leg 1 knee angle     (absolute)
        q3 - leg 2 shoulder angle (absolute)
        q4 - leg 2 knee angle     (absolute)
        q5 - leg 3 shoulder angle (absolute)
        q6 - leg 3 knee angle     (absolute)
        q7 - leg 4 shoulder angle (absolute)
        q8 - leg 4 knee angle     (absolute)
"""
mutable struct Quadruped4{T} <: Model{T}
    # dimensions
	nq::Int # generalized coordinates 
    nu::Int # controls 
    nw::Int # parameters
    nc::Int # contact points

    # parameters

    # torso
    l_torso1::T
    l_torso2::T
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

	# leg 3
        # thigh
    l_thigh3::T
    lc_thigh3::T
    m_thigh3::T
    J_thigh3::T
        # calf
    l_calf3::T
    lc_calf3::T
    m_calf3::T
    J_calf3::T

	# leg 4
        # thigh
    l_thigh4::T
    lc_thigh4::T
    m_thigh4::T
    J_thigh4::T
        # calf
    l_calf4::T
    lc_calf4::T
    m_calf4::T
    J_calf4::T

    # joint friction 
	friction_joint::Vector{T} 

    # environment 
    friction_body_world::Vector{T} 
    friction_foot_world::Vector{T} 
    gravity::T
end

function kinematics_hip(model::Quadruped4, q; hip=:none)
	x = q[1]
	z = q[2]
    t_torso = q[3]

    if hip == :hip1 
        l = model.l_torso1 
    elseif hip == :hip2 
        l = -model.l_torso2 
    end

    return [
            x + l * cos(t_torso); 
            z + l * sin(t_torso)
           ]
end

function kinematics_jacobian_hip(model::Quadruped4, q; hip=:none)
	x = q[1]
	z = q[2]
    t_torso = q[3]

    if hip == :hip1 
        l = model.l_torso1 
    elseif hip == :hip2 
        l = -model.l_torso2 
    end

    jac = zeros(eltype(q), 2, model.nq)

    jac[1, 1] = 1.0 
    jac[1, 3] = -l * sin(t_torso)
    jac[2, 2] = 1.0 
    jac[2, 3] = l * cos(t_torso) 

    return jac
end

function kinematics_thigh(model::Quadruped4, q; leg=:none, mode = :ee)
	x = q[1]
	z = q[2]
    t_torso = q[3]

    if leg == :leg1 
        t_shoulder = q[4]
        l_torso = model.l_torso1
        le_thigh = model.l_thigh1
        lc_thigh = model.lc_thigh1
        hip = :hip1
    elseif leg == :leg2 
        t_shoulder = q[6]
        l_torso = model.l_torso1
        le_thigh = model.l_thigh2
        lc_thigh = model.lc_thigh2
        hip = :hip1
    elseif leg == :leg3 
        t_shoulder = q[8]
        l_torso = -model.l_torso2
        le_thigh = model.l_thigh3
        lc_thigh = model.lc_thigh3
        hip = :hip2
    elseif leg == :leg4 
        t_shoulder = q[10]
        l_torso = -model.l_torso2
        le_thigh = model.l_thigh4
        lc_thigh = model.lc_thigh4
        hip = :hip2
    else
        @error "incorrect leg specified"
    end

    if mode == :ee 
        l_thigh = le_thigh
    elseif mode == :com 
        l_thigh = lc_thigh
    else 
        @error "incorrect mode specified"
    end
    k_hip = kinematics_hip(model, q, hip=hip)

    return k_hip + [l_thigh * sin(t_shoulder); 
                   -l_thigh * cos(t_shoulder)]
end

function kinematics_jacobian_thigh(model::Quadruped4, q; leg=:none, mode = :ee)
	x = q[1]
	z = q[2]
    t_torso = q[3]

    if leg == :leg1 
        idx = 4
        l_torso = model.l_torso1
        le_thigh = model.l_thigh1
        lc_thigh = model.lc_thigh1
        hip = :hip1
    elseif leg == :leg2 
        idx = 6
        l_torso = model.l_torso1
        le_thigh = model.l_thigh2
        lc_thigh = model.lc_thigh2
        hip = :hip1
    elseif leg == :leg3 
        idx = 8
        l_torso = -model.l_torso2
        le_thigh = model.l_thigh3
        lc_thigh = model.lc_thigh3
        hip = :hip2
    elseif leg == :leg4 
        idx = 10
        l_torso = -model.l_torso2
        le_thigh = model.l_thigh4
        lc_thigh = model.lc_thigh4
        hip = :hip2
    else
        @error "incorrect leg specified"
    end

    t_shoulder = q[idx]

    if mode == :ee 
        l_thigh = le_thigh
    elseif mode == :com 
        l_thigh = lc_thigh
    else 
        @error "incorrect mode specified"
    end

    jac = kinematics_jacobian_hip(model, q, hip=hip)

    jac[1, idx] += l_thigh * cos(t_shoulder) 
    jac[2, idx] += l_thigh * sin(t_shoulder) 

    return jac
end

function kinematics_calf(model::Quadruped4, q; leg=:none, mode = :none)
	x = q[1]
	z = q[2]
    t_torso = q[3]

    if leg == :leg1 
        idx = 4
        l_torso = model.l_torso1
        l_thigh = model.l_thigh1
        le_calf = model.l_calf1
        lc_calf = model.lc_calf1
    elseif leg == :leg2 
        idx = 6
        l_torso = model.l_torso1
        l_thigh = model.l_thigh2
        le_calf = model.l_calf2
        lc_calf = model.lc_calf2
    elseif leg == :leg3 
        idx = 8
        l_torso = -model.l_torso2
        l_thigh = model.l_thigh3
        le_calf = model.l_calf3
        lc_calf = model.lc_calf3
    elseif leg == :leg4 
        idx = 10
        l_torso = -model.l_torso2
        l_thigh = model.l_thigh4
        le_calf = model.l_calf4
        lc_calf = model.lc_calf4
    else
        @error "incorrect leg specified"
    end

    t_shoulder = q[idx]
    t_calf = q[idx + 1]

    if mode == :ee 
        l_calf = le_calf
    elseif mode == :com 
        l_calf = lc_calf
    else 
        @error "incorrect mode specified"
    end

    k_thigh = kinematics_thigh(model, q, leg=leg, mode=:ee)

    return k_thigh + [l_calf * sin(t_calf); 
                     -l_calf * cos(t_calf)]
end

function kinematics_jacobian_calf(model::Quadruped4, q; leg=:none, mode = :none)
	x = q[1]
	z = q[2]
    t_torso = q[3]

    if leg == :leg1 
        idx = 4
        l_torso = model.l_torso1
        l_thigh = model.l_thigh1
        le_calf = model.l_calf1
        lc_calf = model.lc_calf1
    elseif leg == :leg2 
        idx = 6
        l_torso = model.l_torso1
        l_thigh = model.l_thigh2
        le_calf = model.l_calf2
        lc_calf = model.lc_calf2
    elseif leg == :leg3 
        idx = 8
        l_torso = -model.l_torso2
        l_thigh = model.l_thigh3
        le_calf = model.l_calf3
        lc_calf = model.lc_calf3
    elseif leg == :leg4 
        idx = 10
        l_torso = -model.l_torso2
        l_thigh = model.l_thigh4
        le_calf = model.l_calf4
        lc_calf = model.lc_calf4
    else
        @error "incorrect leg specified"
    end

    t_shoulder = q[idx]
    t_calf = q[idx + 1]

    if mode == :ee 
        l_calf = le_calf
    elseif mode == :com 
        l_calf = lc_calf
    else 
        @error "incorrect mode specified"
    end

    jac = kinematics_jacobian_thigh(model, q, leg=leg, mode=:ee)
    jac[1, idx + 1] += l_calf * cos(t_calf)
    jac[2, idx + 1] += l_calf * sin(t_calf)

    return jac
end

function kinematics(model::Quadruped4, q)
	p_foot_1 = kinematics_calf(model, q, leg=:leg1, mode=:ee) # foot 1
	p_foot_2 = kinematics_calf(model, q, leg=:leg2, mode=:ee) # foot 2 
	p_foot_3 = kinematics_calf(model, q, leg=:leg3, mode=:ee) # foot 3 
	p_foot_4 = kinematics_calf(model, q, leg=:leg4, mode=:ee) # foot 4
    
    # p_knee_1 = kinematics_thigh(model, q, leg=:leg1, mode=:ee) # knee 1
	# p_knee_2 = kinematics_thigh(model, q, leg=:leg2, mode=:ee) # knee 2 
	# p_knee_3 = kinematics_thigh(model, q, leg=:leg3, mode=:ee) # knee 3 
	# p_knee_4 = kinematics_thigh(model, q, leg=:leg4, mode=:ee) # knee 4

    # p_hip_1 = kinematics_hip(model, q, hip=:hip1) # hip 1 
	# p_hip_2 = kinematics_hip(model, q, hip=:hip2) # hip 2

	[
     p_foot_1; p_foot_2; p_foot_3; p_foot_4; 
    #  p_knee_1; p_knee_2; p_knee_3; p_knee_4; 
    #  p_hip_1; p_hip_2
    ]
end

# Lagrangian
function lagrangian(model::Quadruped4, q, q̇)
	L = 0.0

	# torso
	p_torso = q[1:2]
	v_torso = q̇[1:2]
	
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

	# thigh 2
	p_thigh_2 = kinematics_thigh(model, q, leg=:leg2, mode=:com)
	J_thigh_2 = kinematics_jacobian_thigh(model, q, leg=:leg2, mode=:com)
	v_thigh_2 = J_thigh_2 * q̇

	L += 0.5 * model.m_thigh2 * transpose(v_thigh_2) * v_thigh_2
	L += 0.5 * model.J_thigh2 * q̇[6]^2.0
	L -= model.m_thigh2 * model.gravity * p_thigh_2[2]

	# calf 2
	p_calf_2 = kinematics_calf(model, q, leg=:leg2, mode=:com)
	J_calf_2 = kinematics_jacobian_calf(model, q, leg=:leg2, mode=:com)
	v_calf_2 = J_calf_2 * q̇

	L += 0.5 * model.m_calf2 * transpose(v_calf_2) * v_calf_2
	L += 0.5 * model.J_calf2 * q̇[7]^2.0
	L -= model.m_calf2 * model.gravity * p_calf_2[2]

	# thigh 3
	p_thigh_3 = kinematics_thigh(model, q, leg=:leg3, mode=:com)
	J_thigh_3 = kinematics_jacobian_thigh(model, q, leg=:leg3, mode=:com)
	v_thigh_3 = J_thigh_3 * q̇

	L += 0.5 * model.m_thigh3 * transpose(v_thigh_3) * v_thigh_3
	L += 0.5 * model.J_thigh3 * q̇[8]^2.0
	L -= model.m_thigh3 * model.gravity * p_thigh_3[2]

	# calf 3
	p_calf_3 = kinematics_calf(model, q, leg=:leg3, mode=:com)
	J_calf_3 = kinematics_jacobian_calf(model, q, leg=:leg3, mode=:com)
	v_calf_3 = J_calf_3 * q̇

	L += 0.5 * model.m_calf3 * transpose(v_calf_3) * v_calf_3
	L += 0.5 * model.J_calf3 * q̇[9]^2.0
	L -= model.m_calf3 * model.gravity * p_calf_3[2]

	# thigh 4
	p_thigh_4 = kinematics_thigh(model, q, leg=:leg4, mode=:com)
	J_thigh_4 = kinematics_jacobian_thigh(model, q, leg=:leg4, mode=:com)
	v_thigh_4 = J_thigh_4 * q̇

	L += 0.5 * model.m_thigh4 * transpose(v_thigh_4) * v_thigh_4
	L += 0.5 * model.J_thigh4 * q̇[10]^2.0
	L -= model.m_thigh4 * model.gravity * p_thigh_4[2]

	# calf 4
	p_calf_4 = kinematics_calf(model, q, leg=:leg4, mode=:com)
	J_calf_4 = kinematics_jacobian_calf(model, q, leg=:leg4, mode=:com)
	v_calf_4 = J_calf_4 * q̇

	L += 0.5 * model.m_calf4 * transpose(v_calf_4) * v_calf_4
	L += 0.5 * model.J_calf4 * q̇[11]^2.0
	L -= model.m_calf4 * model.gravity * p_calf_4[2]

	return L
end

function signed_distance(model::Quadruped4, q)
	k = kinematics(model, q) 
    return [
            k[2], k[4], k[6], k[8]
           ]
        #    , k[10], k[12], k[14], k[16], k[18], k[20]]
end

function input_jacobian(::Quadruped4, q)
    [0.0  0.0 -1.0  1.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0;
     0.0  0.0  0.0 -1.0  1.0  0.0  0.0  0.0  0.0  0.0  0.0;
     0.0  0.0 -1.0  0.0  0.0  1.0  0.0  0.0  0.0  0.0  0.0;
     0.0  0.0  0.0  0.0  0.0 -1.0  1.0  0.0  0.0  0.0  0.0;
     0.0  0.0 -1.0  0.0  0.0  0.0  0.0  1.0  0.0  0.0  0.0;
     0.0  0.0  0.0  0.0  0.0  0.0  0.0 -1.0  1.0  0.0  0.0;
     0.0  0.0 -1.0  0.0  0.0  0.0  0.0  0.0  0.0  1.0  0.0;
     0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0 -1.0  1.0]
end

function contact_jacobian(model::Quadruped4, q) 
    J1 = kinematics_jacobian_calf(model, q, leg=:leg1, mode=:ee)
	J2 = kinematics_jacobian_calf(model, q, leg=:leg2, mode=:ee)
	J3 = kinematics_jacobian_calf(model, q, leg=:leg3, mode=:ee)
	J4 = kinematics_jacobian_calf(model, q, leg=:leg4, mode=:ee)

    # J5 = kinematics_jacobian_thigh(model, q, leg=:leg1, mode=:ee)
	# J6 = kinematics_jacobian_thigh(model, q, leg=:leg2, mode=:ee)
	# J7 = kinematics_jacobian_thigh(model, q, leg=:leg3, mode=:ee)
	# J8 = kinematics_jacobian_thigh(model, q, leg=:leg4, mode=:ee)

    # J9  = kinematics_jacobian_hip(model, q, hip=:hip1) 
    # J10 = kinematics_jacobian_hip(model, q, hip=:hip2) 

    [
     J1; 
     J2;
     J3;
     J4;
    #  J5; 
    #  J6; 
    #  J7; 
    #  J8; 
    #  J9; 
    #  J10
    ]
end

# nominal configuration 
function nominal_configuration(model::Quadruped4) 
    [0.0; 0.5; 0.0 * π; 0.25 * π; 0.5 * π; 0.1 * π; 0.3 * π; -0.25 * π; 0.1 * π; -0.5 * π; -0.1 * π]
end

# friction coefficients 
function friction_coefficients(model::Quadruped4) 
    [model.friction_foot_world; model.friction_body_world]
end

# dimensions
nq = 3 + 2 * 4            # generalized coordinates
nu = 2 * 4                # controls
nw = 0                    # parameters
nc = 10                   # contact points 

# parameters
gravity = 9.81            # gravity
friction_body_world = [0.5; 0.5; 0.5; 0.5; 0.5; 0.5] # coefficient of friction
friction_foot_world = [0.5; 0.5; 0.5; 0.5] # coefficient of friction
friction_joint = [0.1; 0.1; 0.1; 0.1; 0.1; 0.1; 0.1; 0.1]

# similar to Unitree A1
m_torso = 4.713 + 4 * 0.696
m_thigh = 1.013
m_calf = 0.166

J_torso = 0.01683 + 4 * 0.696 * 0.183^2.0
J_thigh = 0.00552
J_calf = 0.00299

l_torso1 = 0.183
l_torso2 = 0.183
l_thigh = 0.2
l_calf = 0.2

lc_thigh = 0.5 * l_thigh - 0.00323
lc_calf = 0.5 * l_calf - 0.006435

quadruped4 = Quadruped4(nq, nu, nw, nc,
				l_torso1, l_torso2, m_torso, J_torso,
				l_thigh, lc_thigh, m_thigh, J_thigh,
				l_calf, lc_calf, m_calf, J_calf,
				l_thigh, lc_thigh, m_thigh, J_thigh,
				l_calf, lc_calf, m_calf, J_calf,
				l_thigh, lc_thigh, m_thigh, J_thigh,
				l_calf, lc_calf, m_calf, J_calf,
				l_thigh, lc_thigh, m_thigh, J_thigh,
				l_calf, lc_calf, m_calf, J_calf,
	            friction_joint,
                friction_body_world, 
                friction_foot_world, 
                gravity)

quadruped4_contact_kinematics = [
    q -> kinematics_calf(quadruped4, q, leg=:leg1, mode=:ee),
    q -> kinematics_calf(quadruped4, q, leg=:leg2, mode=:ee),
    q -> kinematics_calf(quadruped4, q, leg=:leg3, mode=:ee),
    q -> kinematics_calf(quadruped4, q, leg=:leg4, mode=:ee),
    # q -> kinematics_thigh(quadruped, q, leg=:leg1, mode=:ee),
    # q -> kinematics_thigh(quadruped, q, leg=:leg2, mode=:ee),
    # q -> kinematics_thigh(quadruped, q, leg=:leg3, mode=:ee),
    # q -> kinematics_thigh(quadruped, q, leg=:leg4, mode=:ee),
    # q -> kinematics_hip(quadruped, q, hip=:hip1),
    # q -> kinematics_hip(quadruped, q, hip=:hip2)
]

quadruped4_contact_kinematics_jacobians = [
    q -> kinematics_jacobian_calf(quadruped4, q, leg=:leg1, mode=:ee),
    q -> kinematics_jacobian_calf(quadruped4, q, leg=:leg2, mode=:ee),
    q -> kinematics_jacobian_calf(quadruped4, q, leg=:leg3, mode=:ee),
    q -> kinematics_jacobian_calf(quadruped4, q, leg=:leg4, mode=:ee),
    # q -> kinematics_jacobian_thigh(quadruped, q, leg=:leg1, mode=:ee),
    # q -> kinematics_jacobian_thigh(quadruped, q, leg=:leg2, mode=:ee),
    # q -> kinematics_jacobian_thigh(quadruped, q, leg=:leg3, mode=:ee),
    # q -> kinematics_jacobian_thigh(quadruped, q, leg=:leg4, mode=:ee),
    # q -> kinematics_jacobian_hip(quadruped, q, hip=:hip1),
    # q -> kinematics_jacobian_hip(quadruped, q, hip=:hip2)
]

name(::Quadruped4) = :quadruped4