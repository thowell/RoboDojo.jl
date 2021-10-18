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
mutable struct Quadruped{T} 
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
	friction_joint::T 
    # TODO: individual joint frictions

    # environment 
    friction_body_world::T 
    friction_foot_world::T 
    gravity::T
end

function kinematics_hip(model::Quadruped, q; hip=:none)
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

function kinematics_jacobian_hip(model::Quadruped, q; hip=:none)
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

function kinematics_thigh(model::Quadruped, q; leg=:none, mode = :ee)
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

function kinematics_jacobian_thigh(model::Quadruped, q; leg=:none, mode = :ee)
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

function kinematics_calf(model::Quadruped, q; leg=:none, mode = :none)
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

function kinematics_jacobian_calf(model::Quadruped, q; leg=:none, mode = :none)
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

function kinematics(model::Quadruped, q)
	p_foot_1 = kinematics_calf(model, q, leg=:leg1, mode=:ee) # foot 1
	p_foot_2 = kinematics_calf(model, q, leg=:leg2, mode=:ee) # foot 2 
	p_foot_3 = kinematics_calf(model, q, leg=:leg3, mode=:ee) # foot 3 
	p_foot_4 = kinematics_calf(model, q, leg=:leg4, mode=:ee) # foot 4
    
    p_knee_1 = kinematics_thigh(model, q, leg=:leg1, mode=:ee) # knee 1
	p_knee_2 = kinematics_thigh(model, q, leg=:leg2, mode=:ee) # knee 2 
	p_knee_3 = kinematics_thigh(model, q, leg=:leg3, mode=:ee) # knee 3 
	p_knee_4 = kinematics_thigh(model, q, leg=:leg4, mode=:ee) # knee 4

    p_hip_1 = kinematics_hip(model, q, hip=:hip1) # hip 1 
	p_hip_2 = kinematics_hip(model, q, hip=:hip2) # hip 2

	[
     p_foot_1; p_foot_2; p_foot_3; p_foot_4; 
     p_knee_1; p_knee_2; p_knee_3; p_knee_4; 
     p_hip_1; p_hip_2
     ]
end

# Lagrangian
function lagrangian(model::Quadruped, q, q̇)
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

function signed_distance(model::Quadruped, q)
	k = kinematics(model, q) 
    return [k[2], k[4], k[6], k[8], k[10], k[12], k[14], k[16], k[18], k[20]]
end

function input_jacobian()
    [0.0  0.0 -1.0  1.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0;
     0.0  0.0  0.0 -1.0  1.0  0.0  0.0  0.0  0.0  0.0  0.0;
     0.0  0.0 -1.0  0.0  0.0  1.0  0.0  0.0  0.0  0.0  0.0;
     0.0  0.0  0.0  0.0  0.0 -1.0  1.0  0.0  0.0  0.0  0.0;
     0.0  0.0 -1.0  0.0  0.0  0.0  0.0  1.0  0.0  0.0  0.0;
     0.0  0.0  0.0  0.0  0.0  0.0  0.0 -1.0  1.0  0.0  0.0;
     0.0  0.0 -1.0  0.0  0.0  0.0  0.0  0.0  0.0  1.0  0.0;
     0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0 -1.0  1.0]
end

function contact_jacobian(model::Quadruped, q) 
    J1 = kinematics_jacobian_calf(model, q, leg=:leg1, mode=:ee)
	J2 = kinematics_jacobian_calf(model, q, leg=:leg2, mode=:ee)
	J3 = kinematics_jacobian_calf(model, q, leg=:leg3, mode=:ee)
	J4 = kinematics_jacobian_calf(model, q, leg=:leg4, mode=:ee)

    J5 = kinematics_jacobian_thigh(model, q, leg=:leg1, mode=:ee)
	J6 = kinematics_jacobian_thigh(model, q, leg=:leg2, mode=:ee)
	J7 = kinematics_jacobian_thigh(model, q, leg=:leg3, mode=:ee)
	J8 = kinematics_jacobian_thigh(model, q, leg=:leg4, mode=:ee)

    J9  = kinematics_jacobian_hip(model, q, hip=:hip1) 
    J10 = kinematics_jacobian_hip(model, q, hip=:hip2) 

    [
     J1; 
     J2;
     J3;
     J4;
     J5; 
     J6; 
     J7; 
     J8; 
     J9; 
     J10
    ]
end

# dimensions
nq = 3 + 2 * 4            # generalized coordinates
nu = 2 * 4                # controls
nw = 0                    # parameters
nc = 10                   # contact points 

# parameters
gravity = 9.81            # gravity
friction_body_world = 0.5 # coefficient of friction
friction_foot_world = 0.5 # coefficient of friction
friction_joint = 0.1

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

quadruped = Quadruped(nq, nu, nw, nc,
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

# variables
@variables q[1:nq] q̇[1:nq]

# Lagrangian
L = lagrangian(quadruped, q, q̇);

ddL = Symbolics.hessian(L, [q; q̇])#, simplify=true)
dLq = Symbolics.gradient(L, q)#, simplify=true)
ddLq̇q = ddL[nq .+ (1:nq), 1:nq]

M = ddL[nq .+ (1:nq), nq .+ (1:nq)]
C = ddLq̇q * q̇ - dLq

mass_matrix_sym = eval(Symbolics.build_function(M, q)[1])
dynamics_bias_sym = eval(Symbolics.build_function(C, q, q̇)[1])

# ϕ = signed_distance(quadruped, q)
# dϕ = Symbolics.jacobian(ϕ, q)
# signed_distance_jacobian = eval(Symbolics.build_function(dϕ, q)[1])

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
        - model.friction_joint * [0.0; 0.0; 0.0; vm2[4:11]] # joint friction
        + λ1)
end

# @variables h[1:1] q0[1:nq] q1[1:nq] u1[1:nu] w1[1:nw] λ1[1:nq] q2[1:nq]
# d = dynamics(quadruped, h, q0, q1, u1, w1, λ1, q2)

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

    rotation_foot1 = vector_rotation_matrix([0.0; 1.0])
    rotation_foot2 = vector_rotation_matrix([0.0; 1.0])
    rotation_foot3 = vector_rotation_matrix([0.0; 1.0])
    rotation_foot4 = vector_rotation_matrix([0.0; 1.0])

    rotation_knee1 = vector_rotation_matrix([0.0; 1.0])
    rotation_knee2 = vector_rotation_matrix([0.0; 1.0])
    rotation_knee3 = vector_rotation_matrix([0.0; 1.0])
    rotation_knee4 = vector_rotation_matrix([0.0; 1.0])

    rotation_hip1 = vector_rotation_matrix([0.0; 1.0])
    rotation_hip2 = vector_rotation_matrix([0.0; 1.0])

    # tangential velocity
    v1 = (q2 - q1) / h[1] 

    vT_foot1 = (rotation_foot1 * kinematics_jacobian_calf(model, q2, leg=:leg1, mode=:ee) * v1)[1]
    vT_foot2 = (rotation_foot2 * kinematics_jacobian_calf(model, q2, leg=:leg2, mode=:ee) * v1)[1]
    vT_foot3 = (rotation_foot3 * kinematics_jacobian_calf(model, q2, leg=:leg3, mode=:ee) * v1)[1]
    vT_foot4 = (rotation_foot4 * kinematics_jacobian_calf(model, q2, leg=:leg4, mode=:ee) * v1)[1]

    vT_knee1 = (rotation_knee1 * kinematics_jacobian_thigh(model, q2, leg=:leg1, mode=:ee) * v1)[1]
    vT_knee2 = (rotation_knee2 * kinematics_jacobian_thigh(model, q2, leg=:leg2, mode=:ee) * v1)[1]
    vT_knee3 = (rotation_knee3 * kinematics_jacobian_thigh(model, q2, leg=:leg3, mode=:ee) * v1)[1]
    vT_knee4 = (rotation_knee4 * kinematics_jacobian_thigh(model, q2, leg=:leg4, mode=:ee) * v1)[1]

    vT_hip1 = (rotation_hip1 * kinematics_jacobian_hip(model, q2, hip=:hip1) * v1)[1]
    vT_hip2 = (rotation_hip2 * kinematics_jacobian_hip(model, q2, hip=:hip2) * v1)[1]

    # contact forces
    J = contact_jacobian(model, q2)
    λ1 = transpose(J) * [
                         transpose(rotation_foot1) * [b1[1]; γ1[1]];
                         transpose(rotation_foot2) * [b1[2]; γ1[2]];
                         transpose(rotation_foot3) * [b1[3]; γ1[3]];
                         transpose(rotation_foot4) * [b1[4]; γ1[4]];
                         transpose(rotation_knee1) * [b1[5]; γ1[5]];
                         transpose(rotation_knee2) * [b1[6]; γ1[6]];
                         transpose(rotation_knee3) * [b1[7]; γ1[7]];
                         transpose(rotation_knee4) * [b1[8]; γ1[8]];
                         transpose(rotation_hip1)  * [b1[9]; γ1[9]];
                         transpose(rotation_hip2)  * [b1[10]; γ1[10]];
                         ]

	[
     dynamics(model, h, q0, q1, u1, w1, λ1, q2);
	 signed_distance(model, q2) - sγ1;
     vT_foot1 - sb1[1];
     vT_foot2 - sb1[2];
     vT_foot3 - sb1[3];
     vT_foot4 - sb1[4];
     vT_knee1 - sb1[5];
     vT_knee2 - sb1[6];
     vT_knee3 - sb1[7];
     vT_knee4 - sb1[8];
     vT_hip1  - sb1[9]; 
     vT_hip2  - sb1[10];
     friction_foot_world[1] * γ1[1] - ψ1[1];
     friction_foot_world[1] * γ1[2] - ψ1[2];
     friction_foot_world[1] * γ1[3] - ψ1[3];
     friction_foot_world[1] * γ1[4] - ψ1[4];
     friction_body_world[1] * γ1[5] - ψ1[5];
     friction_body_world[1] * γ1[6] - ψ1[6];
     friction_body_world[1] * γ1[7] - ψ1[7];
     friction_body_world[1] * γ1[8] - ψ1[8];
     friction_body_world[1] * γ1[9] - ψ1[9];
     friction_body_world[1] * γ1[10] - ψ1[10];
     γ1 .* sγ1 .- μ[1];
	 cone_product([ψ1[1]; b1[1]], [sψ1[1]; sb1[1]]) - [μ[1]; 0.0]; 
     cone_product([ψ1[2]; b1[2]], [sψ1[2]; sb1[2]]) - [μ[1]; 0.0]; 
     cone_product([ψ1[3]; b1[3]], [sψ1[3]; sb1[3]]) - [μ[1]; 0.0]; 
     cone_product([ψ1[4]; b1[4]], [sψ1[4]; sb1[4]]) - [μ[1]; 0.0]; 
     cone_product([ψ1[5]; b1[5]], [sψ1[5]; sb1[5]]) - [μ[1]; 0.0]; 
     cone_product([ψ1[6]; b1[6]], [sψ1[6]; sb1[6]]) - [μ[1]; 0.0]; 
     cone_product([ψ1[7]; b1[7]], [sψ1[7]; sb1[7]]) - [μ[1]; 0.0]; 
     cone_product([ψ1[8]; b1[8]], [sψ1[8]; sb1[8]]) - [μ[1]; 0.0]; 
     cone_product([ψ1[9]; b1[9]], [sψ1[9]; sb1[9]]) - [μ[1]; 0.0]; 
     cone_product([ψ1[10]; b1[10]], [sψ1[10]; sb1[10]]) - [μ[1]; 0.0]; 
    ]
end

nz = nq + nc + nc + nc + nc + nc + nc 
ny = nc + nc + nc + nc + nc + nc
nθ = 2nq + nu + nw + 2 + 1 

@variables z[1:nz] θ[1:nθ] μ[1:1]
r = residual(quadruped, z, θ, μ)
rz = Symbolics.jacobian(r, z)
rθ = Symbolics.jacobian(r, θ)
rz_sp = Symbolics.sparsejacobian(r, z)
rθ_sp = Symbolics.sparsejacobian(r, θ)

r_quadruped! = eval(Symbolics.build_function(r, z, θ, μ)[2])
rz_quadruped! = eval(Symbolics.build_function(rz, z, θ)[2])
rθ_quadruped! = eval(Symbolics.build_function(rθ, z, θ)[2])

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

# r_quadruped!(r0, z0, θ0, μ0)
# rz_quadruped!(rz0, z0, θ0)
# rθ_quadruped!(rθ0, z0, θ0)
# rz_sp!(rz0_sp, z0, θ0)
# rθ_sp!(rθ0_sp, z0, θ0)

# using BenchmarkTools
# @benchmark r_quadruped!($r0, $z0, $θ0, $μ0)
# @benchmark rz_quadruped!($rz0, $z0, $θ0)
# @benchmark rθ_quadruped!($rθ0, $z0, $θ0)
# @benchmark rz_sp!($rz0_sp, $z0, $θ0)
# @benchmark rθ_sp!($rθ0_sp, $z0, $θ0)