abstract type Policy end

function policy(p::Policy, traj::ContactTraj, t::T) where T
    @warn "policy not defined"
    return nothing
end

"""
    no policy
"""
struct NoPolicy{T} <: Policy
    u::Vector{T}
end

function no_policy(model::ContactModel)
    NoPolicy(zeros(model.dim.u))
end

function policy(p::NoPolicy, x, traj::ContactTraj, t)
    return p.u
end

"""
    control saturation
"""

control_saturation(u, uL, uU) = min.(max.(uL, u), uU)
