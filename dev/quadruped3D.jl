nq = 3 + 4 + 3 # position + quaternion + leg joints

function quaternion_rotation_matrix(q)
	r, i, j, k  = q

	r11 = 1.0 - 2.0 * (j^2.0 + k^2.0)
	r12 = 2.0 * (i * j - k * r)
	r13 = 2.0 * (i * k + j * r)

	r21 = 2.0 * (i * j + k * r)
	r22 = 1.0 - 2.0 * (i^2.0 + k^2.0)
	r23 = 2.0 * (j * k - i * r)

	r31 = 2.0 * (i * k - j * r)
	r32 = 2.0 * (j * k + i * r)
	r33 = 1.0 - 2.0 * (i^2.0 + j^2.0)

	[r11 r12 r13;
     r21 r22 r23;
     r31 r32 r33]
end

function rotation_matrix_x(a) 
    [1.0 0.0    0.0
     0.0 cos(a) -sin(a)
     0.0 sin(a)  cos(a)]
end

leg1_offset = [1.0; 1.0; 0.0]
l_shoulder = 1.0
l_thigh = 1.0 
l_calf = 1.0 
function kinematics_thigh_calf(a1, a2) 
    [
     l_thigh * sin(a1) + l_calf * sin(a2);
     0.0;
    -l_thigh * cos(a1) - l_calf * cos(a2);
    ]
end

function foot_kinematics(q) 
    pos = q[1:3] 
    quat = q[4:7] 
    joint = q[8:10] 

    k_thigh_calf = kinematics_thigh_calf(q[9], q[10])
    k_shoulder_thigh_calf = rotation_matrix_x(q[8]) * (k_thigh_calf + [0.0; l_shoulder; 0.0])
    k_offset_shoulder_thigh_calf = k_shoulder_thigh_calf + leg1_offset 
    k_body_offset_shoulder_thigh_calf = quaternion_rotation_matrix(quat) * k_offset_shoulder_thigh_calf 
    k = k_body_offset_shoulder_thigh_calf + pos 

    return k 
end

q0 = [0.0; 0.0; 0.0; 1.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0] 
foot_kinematics(q0)

@variables q[1:nq] q̇[1:nq]

k = foot_kinematics(q)
j = Symbolics.jacobian(k, q)

k_func = eval(Symbolics.build_function(k, q)[1])
j_func = eval(Symbolics.build_function(j, q)[1])

function lagrangian(q, q̇) 
    L = 0.0 
    mass = 1.0 
    inertia = 1.0 
    gravity = 9.81
    k_foot = k_func(q)
    v_foot = j_func(q) * q̇ 
    L += mass * transpose(v_foot) * v_foot 
    L -= mass * gravity * k_foot[3]

    return L 
end

L = lagrangian(q, q̇)

ddL = Symbolics.hessian(L, [q; q̇])#, simplify=true)
dLq = Symbolics.gradient(L, q)#, simplify=true)
ddLq̇q = ddL[nq .+ (1:nq), 1:nq]

M = ddL[nq .+ (1:nq), nq .+ (1:nq)]
C = ddLq̇q * q̇ - dLq

mass_matrix_sym = eval(Symbolics.build_function(M, q)[1])
dynamics_bias_sym = eval(Symbolics.build_function(C, q, q̇)[1])

function lagrangian_derivatives(q, v)
	D1L = -1.0 * dynamics_bias_sym(q, v)
    D2L = mass_matrix_sym(q) * v
	return D1L, D2L
end

function dynamics(h, q0, q1, u1, w1, λ1, q2)
	# evalutate at midpoint
	qm1 = 0.5 * (q0 + q1)
    vm1 = (q1 - q0) / h[1]
    qm2 = 0.5 * (q1 + q2)
    vm2 = (q2 - q1) / h[1]

	D1L1, D2L1 = lagrangian_derivatives(qm1, vm1)
	D1L2, D2L2 = lagrangian_derivatives(qm2, vm2)

	return (0.5 * h[1] * D1L1 + D2L1 + 0.5 * h[1] * D1L2 - D2L2 + λ1)
end

@variables h[1:1] q0[1:nq] q1[1:nq] u1[1:0] w1[1:0] λ1[1:nq] q2[1:nq]
d = dynamics(h, q0, q1, u1, w1, λ1, q2)

