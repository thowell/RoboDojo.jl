function codegen_dynamics(model)
    # dimensions
    nq = model.nq

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

    return mass_matrix, dynamics_bias
end