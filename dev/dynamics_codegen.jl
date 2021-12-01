using Symbolics, LinearAlgebra

n = 2
m1 = 1.0 
m2 = 1.0 
gr = 9.81 
l1 = 1.0 
l2 = 1.0
k1(q) = [0.5 * l1 * cos(q[1]); 0.5 * l1 * sin(q[1])] 
k2(q) = k1(q) + [0.5 * l2 * cos(q[1] + q[2]); 0.5 * l1 * sin(q[1] + q[2])]

@variables q[1:n] q̇[1:n]

# kinematics 
k1s = k1(q) 
k2s = k2(q)
ks = [k1(q); k2(q)]

# kinematics Jacobians 
j1s = Symbolics.jacobian(k1s, q) 
j2s = Symbolics.jacobian(k2s, q) 
js = Symbolics.jacobian(ks, q)
j1 = eval(Symbolics.build_function(j1s, q)[1])
j2 = eval(Symbolics.build_function(j2s, q)[1])
j = eval(Symbolics.build_function(js, q)[1])

# direct approach 
@variables λ1[1:n] λ2[1:n] λ3[1:n] x[1:n]

a = js * λ1
∂a∂q = Symbolics.jacobian(a, q)

# b = transpose(\lamba    ∂a∂q) * λ2
# ∂b∂q = Symbolics.jacobian(b, q)

# mass matrix 
mass_matrix_max = Diagonal([m1; m1; m2; m2])
Ms = transpose(js) * mass_matrix_max * js
mass_matrix = eval(Symbolics.build_function(Ms, q)[1])

# mass matrix derivative
Mx = Ms * x
∂Mx∂q = Symbolics.jacobian(Mx, q)
mass_matrix_deriv = eval(Symbolics.build_function(∂Mx∂q, q, x)[1])

# dynamics bias 
Cs = ∂Mx∂q * x
dynamics_bias = eval(Symbolics.build_function(Cs, q, x)[1])

# dynamics bias derivative 
∂C∂q̇ = 2.0 * transpose(js) * mass_matrix_max * ∂a∂q 
∂C∂q = transpose(∂a∂q) * mass_matrix_max * ∂a∂q


# differentiation approach
function lagrangian(q, q̇) 
    L = 0.0 
    
    v1 = j1(q) * q̇ 
    v2 = j2(q) * q̇ 

    L += 0.5 * m1 * transpose(v1) * v1 
    # L -= m1 * gr * k1(q)[2]

    L += 0.5 * m2 * transpose(v2) * v2 
    # L -= m2 * gr * k2(q)[2] 

    return L 
end

L = lagrangian(q, q̇)

ddL = Symbolics.hessian(L, [q; q̇])
dLq = Symbolics.gradient(L, q)
ddLq̇q = ddL[n .+ (1:n), 1:n]

M = ddL[n .+ (1:n), n .+ (1:n)]
C = ddLq̇q * q̇ - dLq

mass_matrix_lag = eval(Symbolics.build_function(M, q)[1])
dynamics_bias_lag = eval(Symbolics.build_function(C, q, q̇)[1])


# test 
q0 = rand(n) 
q̇0 = rand(n)

norm(mass_matrix(q0) - mass_matrix_lag(q0), Inf)
norm(dynamics_bias(q0, q̇0) - dynamics_bias_lag(q0, q̇0), Inf)

dynamics_bias(q0, q̇0)
dynamics_bias_lag(q0, q̇0)