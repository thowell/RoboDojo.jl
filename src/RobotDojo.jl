module RobotDojo

using LinearAlgebra
using Symbolics
using SparseArrays
using StaticArrays
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
include("solver/linear_solver.jl")
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

HOPPER_CODEGEN = :load
include(joinpath(@__DIR__, "../src/robots/hopper/codegen.jl"))

BIPED_CODEGEN = :load
include(joinpath(@__DIR__, "../src/robots/biped/codegen.jl"))

QUADRUPED_CODEGEN = :load
include(joinpath(@__DIR__, "../src/robots/quadruped/codegen.jl"))

# Policy
include(joinpath("../src/policy/raibert.jl"))

export
    hopper, biped, quadruped, 
    nominal_configuration,
    simulate!, Simulator, 
    visualize!, Visualizer

end # module
