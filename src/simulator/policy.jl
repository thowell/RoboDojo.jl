abstract type Policy{T} end

"""
    empty policy
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