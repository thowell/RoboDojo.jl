module RobotDojo

using LinearAlgebra
using Symbolics
using SparseArrays
using BenchmarkTools
using MeshCat
using GeometryBasics
using Colors
using CoordinateTransformations
using Rotations
using Parameters
using JLD2

# Utilities
include("utils.jl")

# Solver
include("solver/lu.jl")
include("solver/cones.jl")
include("solver/indices.jl")
include("solver/interior_point.jl")

# Simulator
include(joinpath("../src/simulator/model.jl"))
include(joinpath("../src/simulator/dimensions.jl"))
include(joinpath("../src/simulator/indices.jl"))
include(joinpath("../src/simulator/disturbances.jl"))
include(joinpath("../src/simulator/trajectory.jl"))
include(joinpath("../src/simulator/policy.jl"))
include(joinpath("../src/simulator/simulator.jl"))
include(joinpath("../src/simulator/residual.jl"))
include(joinpath("../src/simulator/codegen.jl"))

# Robots
include(joinpath("../src/robots/visual_utils.jl"))

include(joinpath("../src/robots/hopper/model.jl"))
include(joinpath("../src/robots/hopper/simulator.jl"))
include(joinpath("../src/robots/hopper/visuals.jl"))

include(joinpath("../src/robots/biped/model.jl"))
include(joinpath("../src/robots/biped/visuals.jl"))

include(joinpath("../src/robots/quadruped/model.jl"))
include(joinpath("../src/robots/quadruped/visuals.jl"))

include(joinpath("../src/robots/box/model.jl"))
include(joinpath("../src/robots/box/visuals.jl"))

include(joinpath("../src/robots/integrator.jl"))
include(joinpath("../src/robots/codegen.jl"))

include(joinpath("../src/robots/load.jl"))

# load simulation environments

# RESIDUAL_EXPR = Dict{String, Any}()
# residual_expr(model::Model) = RESIDUAL_EXPR[String(name(model)) * "_r"]
# jacobian_var_expr(model::Model) = RESIDUAL_EXPR[String(name(model)) * "_rz"]
# jacobian_data_expr(model::Model) = RESIDUAL_EXPR[String(name(model)) * "_rθ"]

# model = hopper 
# path = joinpath(@__DIR__, "robots", String(name(model)), "expr/expr.jld2")
# @load path r_model rz_model rθ_model
# RESIDUAL_EXPR[String(name(model)) * "_r"] = eval(r_model)
# RESIDUAL_EXPR[String(name(model)) * "_rz"] = eval(rz_model)
# RESIDUAL_EXPR[String(name(model)) * "_rθ"] = eval(rθ_model)

# model = biped 
# path = joinpath(@__DIR__, "robots", String(name(model)), "expr/expr.jld2")
# @load path r_model rz_model rθ_model
# RESIDUAL_EXPR[String(name(model)) * "_r"] = eval(r_model)
# RESIDUAL_EXPR[String(name(model)) * "_rz"] = eval(rz_model)
# RESIDUAL_EXPR[String(name(model)) * "_rθ"] = eval(rθ_model)

# model = quadruped 
# path = joinpath(@__DIR__, "robots", String(name(model)), "expr/expr.jld2")
# @load path r_model rz_model rθ_model
# RESIDUAL_EXPR[String(name(model)) * "_r"] = eval(r_model)
# RESIDUAL_EXPR[String(name(model)) * "_rz"] = eval(rz_model)
# RESIDUAL_EXPR[String(name(model)) * "_rθ"] = eval(rθ_model)

# Policy

export
    hopper, biped, quadruped, 
    nominal_configuration,
    simulate!, Simulator, 
    visualize!, Visualizer,
    friction_coefficients,
    lagrangian, 
    num_var, num_data, 
    residual

end # module
