# if CODEGEN == :load 
#     @load path r_model rz_model rθ_model 
# else
# dimensions
nq = model.nq
nf = length(friction_coefficients(model))

# variables
@variables q[1:nq] q̇[1:nq]

# Lagrangian
L = lagrangian(model, q, q̇)

ddL = Symbolics.hessian(L, [q; q̇])
dLq = Symbolics.gradient(L, q)
ddLq̇q = ddL[nq .+ (1:nq), 1:nq]

M = ddL[nq .+ (1:nq), nq .+ (1:nq)]
C = ddLq̇q * q̇ - dLq

# mass matrix and dynamics bias
mass_matrix = eval(Symbolics.build_function(M, q)[1])
dynamics_bias = eval(Symbolics.build_function(C, q, q̇)[1])

# dimensions
nz = num_var(model)
nθ = num_data(model, nf=nf)

# residual
@variables z[1:nz] θ[1:nθ] μ[1:1]
r = residual(model, mass_matrix, dynamics_bias, contact_kinematics, contact_kinematics_jacobians, z, θ, μ)
rz = Symbolics.jacobian(r, z)
rθ = Symbolics.jacobian(r, θ)

r_model = Symbolics.build_function(r, z, θ, μ)[2]
rz_model = Symbolics.build_function(rz, z, θ)[2]
rθ_model = Symbolics.build_function(rθ, z, θ)[2]

@save path r_model rz_model rθ_model

# end

# RESIDUAL_EXPR[String(name(model)) * "_r"] = eval(r_model)
# RESIDUAL_EXPR[String(name(model)) * "_rz"] = eval(rz_model)
# RESIDUAL_EXPR[String(name(model)) * "_rθ"] = eval(rθ_model)

# using BenchmarkTools
# @benchmark r_model!($r0, $z0, $θ0, $μ0)
# @benchmark rz_model!($rz0, $z0, $θ0)
# @benchmark rθ_model!($rθ0, $z0, $θ0)