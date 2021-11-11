abstract type Disturbances{T} end

"""
    empty disturbances
"""
struct EmptyDisturbances{T} <: Disturbances{T}
    w::Vector{T}
end

function empty_disturbances(model)
    EmptyDisturbances(zeros(model.nw))
end

function disturbances(d::EmptyDisturbances, x, t)
    return d.w
end
