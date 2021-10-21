abstract type Policy{T} end

function policy(p::Policy, traj::Trajectory, t::T) where T
    @warn "policy not defined"
    return nothing
end

"""
    no policy
"""
struct EmptyPolicy{T} <: Policy{T}
    u::Vector{T}
end

function empty_policy(model)
    EmptyPolicy(zeros(model.nu))
end

function policy(p::EmptyPolicy, traj::Trajectory, t)
    return p.u
end

"""
    control saturation
"""
control_saturation(u, uL, uU) = min.(max.(uL, u), uU) #TODO: make allocation free
