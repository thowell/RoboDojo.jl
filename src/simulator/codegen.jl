function codegen_residual(model, mass_matrix, dynamics_bias, contact_kinematics, contact_kinematics_jacobians)

    # dimensions
    nz = num_var(model)
    nθ = num_data(model, nf=length(friction_coefficients(model)))

    # residual
    @variables z[1:nz] θ[1:nθ] μ[1:1]
    r = residual(model, mass_matrix, dynamics_bias, contact_kinematics, contact_kinematics_jacobians, z, θ, μ)
    rz = Symbolics.jacobian(r, z)
    rθ = Symbolics.jacobian(r, θ)

    # methods
    r_model = Symbolics.build_function(r, z, θ, μ)[2]
    rz_model = Symbolics.build_function(rz, z, θ)[2]
    rθ_model = Symbolics.build_function(rθ, z, θ)[2]

    return r_model, rz_model, rθ_model 
end