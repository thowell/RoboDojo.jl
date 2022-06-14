module RoboDojo

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
using Scratch

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
include(joinpath("../src/simulator/dynamics.jl"))

# Robots
include(joinpath("../src/robots/visual_utils.jl"))

include(joinpath("../src/robots/hopper/model.jl"))
include(joinpath("../src/robots/hopper/simulator.jl"))
include(joinpath("../src/robots/hopper/visuals.jl"))

include(joinpath("../src/robots/hopper1/model.jl"))
include(joinpath("../src/robots/hopper1/visuals.jl"))

include(joinpath("../src/robots/biped/model.jl"))
include(joinpath("../src/robots/biped/visuals.jl"))

include(joinpath("../src/robots/quadruped/model.jl"))
include(joinpath("../src/robots/quadruped/visuals.jl"))

include(joinpath("../src/robots/quadruped4/model.jl"))
include(joinpath("../src/robots/quadruped4/visuals.jl"))

include(joinpath("../src/robots/quadruped_spring/model.jl"))
include(joinpath("../src/robots/quadruped_spring/visuals.jl"))

include(joinpath("../src/robots/halfquadruped/model.jl"))
include(joinpath("../src/robots/halfquadruped/visuals.jl"))

include(joinpath("../src/robots/halfcheetah/model.jl"))
include(joinpath("../src/robots/halfcheetah/visuals.jl"))

include(joinpath("../src/robots/halfcheetah4/model.jl"))
include(joinpath("../src/robots/halfcheetah4/visuals.jl"))

include(joinpath("../src/robots/box/model.jl"))
include(joinpath("../src/robots/box/visuals.jl"))

include(joinpath("../src/robots/particle/model.jl"))
include(joinpath("../src/robots/particle/simulator.jl"))
include(joinpath("../src/robots/particle/visuals.jl"))

include(joinpath("../src/robots/centroidal_quadruped/model.jl"))
include(joinpath("../src/robots/centroidal_quadruped/simulator.jl"))
include(joinpath("../src/robots/centroidal_quadruped/visuals.jl"))

include(joinpath("../src/robots/centroidal_quadruped_ic/model.jl"))
include(joinpath("../src/robots/centroidal_quadruped_ic/simulator.jl"))
include(joinpath("../src/robots/centroidal_quadruped_ic/visuals.jl"))

include(joinpath("../src/robots/point_foot_quadruped/model.jl"))
include(joinpath("../src/robots/point_foot_quadruped/simulator.jl"))
include(joinpath("../src/robots/point_foot_quadruped/visuals.jl"))

include(joinpath("../src/robots/integrator.jl"))
include(joinpath("../src/robots/codegen.jl"))

include(joinpath("../src/robots/load.jl"))

export
    hopper, hopper1, biped, quadruped, quadruped4, halfquadruped,
    halfcheetah, halfcheetah4, box, particle, centroidal_quadruped, centroidal_quadruped_ic, point_foot_quadruped,
    nominal_configuration,
    step!,
    simulate!,
    Simulator,
    visualize!,
    Visualizer,
    friction_coefficients,
    lagrangian,
    num_var,
    num_data,
    residual,
    render,
    open

end # module
