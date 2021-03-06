function dynamics(s::Simulator{T}, y::AbstractVector{T}, x::AbstractVector{T},
        u::AbstractVector{T}, w::AbstractVector{T};
        verbose=false) where T

    status = dynamics(s, x, u, w, diff_sol=false, verbose=verbose)
    # Fill the next state
    nq = s.model.nq
    y[1:nq] .= s.traj.q[3] # q3
    y[nq .+ (1:nq)] .= s.traj.v[2] # v25
    return status
end

function dynamics_jacobian_state(s::Simulator{T}, dx::AbstractMatrix{T}, x::AbstractVector{T},
        u::AbstractVector{T}, w::AbstractVector{T};
        verbose=false) where T

    status = dynamics(s, x, u, w, diff_sol=true, verbose=verbose)
    # Fill the jacobian
    # there are weird chain rules because RoboDojo uses q3 = fct(q1, q2)
    # and not (q3, v25) = fct(q2, v15)
    nq = s.model.nq
    @views dx[1:nq, 1:nq] .= s.grad.∂q3∂q2[1] # ∂q3∂q2
    @views dx[1:nq, 1:nq] .+= s.grad.∂q3∂q1[1] # ∂q3∂q2

    @views dx[nq .+ (1:nq), 1:nq] .= dx[1:nq, 1:nq] ./ s.h # ∂v25∂q2
    @views dx[nq .+ (1:nq), 1:nq][1:nq+1:nq^2] .+= - 1/s.h # ∂v25∂q2

    @views dx[1:nq, nq .+ (1:nq)] .= -s.h .* s.grad.∂q3∂q1[1] # ∂q3∂v15

    @views dx[nq .+ (1:nq), nq .+ (1:nq)] .= dx[1:nq, nq .+ (1:nq)] ./ s.h # ∂v25∂v15
    return status
end

function dynamics_jacobian_input(s::Simulator{T}, du::AbstractMatrix{T}, x::AbstractVector{T},
        u::AbstractVector{T}, w::AbstractVector{T};
        verbose=false) where T

    status = dynamics(s, x, u, w, diff_sol=true, verbose=verbose)
    # Fill the jacobian
    nq = s.model.nq
    nu = s.model.nu
    @views du[1:nq, 1:nu] .= s.grad.∂q3∂u1[1] # ∂q3∂u1
    @views du[nq .+ (1:nq), 1:nu] .= s.grad.∂q3∂u1[1] ./ s.h# ∂v25∂u1
    du .*= s.h # because we use force inputs instead of impulse inputs
    return status
end

function dynamics(s::Simulator{T}, x::AbstractVector{T},
        u::AbstractVector{T}, w::AbstractVector{T};
        diff_sol=false, verbose=false) where T

    s.ip.opts.diff_sol = diff_sol
    nq = s.model.nq
    status = false
    traj = s.traj
    # One single simulation step
    @assert length(traj.u) == 1

    # reset trajectory
    reset!(traj)
    reset!(s.grad)

    # set initial state
    @views q2 = x[1:nq]
    @views v15 = x[nq .+ (1:nq)]
    set_state!(s, q2, v15, 1)

    # set control input
    traj.u[1] .= u .* s.h # u is a force and u * h is an impulse
    # set disturbances
    traj.w[1] .= w
    # take a step
    status = step!(s, 1, verbose=verbose)
    s.ip.opts.diff_sol = false
    return status
end
