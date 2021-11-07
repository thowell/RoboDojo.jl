function lagrangian_derivatives(mass_matrix, dynamics_bias, q, v)
    D1L = -1.0 * dynamics_bias(q, v)
    D2L = mass_matrix(q) * v
    return D1L, D2L
end

function dynamics(model::Model, mass_matrix, dynamics_bias, h, q0, q1, u1, w1, λ1, q2)
    # evalutate at midpoint
    qm1 = 0.5 * (q0 + q1)
    vm1 = (q1 - q0) / h[1]
    qm2 = 0.5 * (q1 + q2)
    vm2 = (q2 - q1) / h[1]

    D1L1, D2L1 = lagrangian_derivatives(mass_matrix, dynamics_bias, qm1, vm1)
    D1L2, D2L2 = lagrangian_derivatives(mass_matrix, dynamics_bias, qm2, vm2)

    d = 0.5 * h[1] * D1L1 + D2L1 + 0.5 * h[1] * D1L2 - D2L2 # variational integrator (midpoint)
    d .+= transpose(input_jacobian(model, qm2)) * u1             # control inputs
    d .+= λ1                                                # contact impulses
    d[(floating_base_dim(model) + 1):end] .-= model.friction_joint .* vm2[(floating_base_dim(model) + 1):end] # joint friction

    return d
end