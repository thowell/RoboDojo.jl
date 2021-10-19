abstract type Disturbances end

function disturbances(d::Disturbances, x, t)
    @warn "disturbances not defined"
    return nothing
end

"""
    empty disturbances
"""
struct EmptyDisturbances{nw,T} <: Disturbances
    w::SVector{nw,T}
end

function empty_disturbances(model)
    EmptyDisturbances(@SVector zeros(model.nw))
end

function disturbances(d::EmptyDisturbances, x, t)
    return d.w
end
