append!(empty!(LOAD_PATH), Base.DEFAULT_LOAD_PATH)
using Pkg

################################################################################
# Generate notebooks
################################################################################
exampledir = joinpath(@__DIR__, "..", "examples")
Pkg.activate(exampledir)
Pkg.instantiate()
include(joinpath(exampledir, "generate_notebooks.jl"))

################################################################################
# Build simulation environments
################################################################################
pkgdir = joinpath(@__DIR__, "..")
Pkg.activate(pkgdir)
Pkg.instantiate()

using JLD2 
using Symbolics
using LinearAlgebra

include(joinpath("../src/utils.jl"))
include(joinpath("../src/simulator/model.jl"))
include(joinpath("../src/simulator/dimensions.jl"))
include(joinpath("../src/simulator/residual.jl"))
include(joinpath("../src/robots/integrator.jl"))

include(joinpath("../src/solver/cones.jl"))

path_codegen = joinpath(@__DIR__,"..", "src/simulator/codegen.jl")

include(joinpath(@__DIR__,"..", "src/robots/hopper/model.jl"))
include(joinpath(@__DIR__,"..", "src/robots/hopper/simulator.jl"))
model = hopper 
contact_kinematics = hopper_contact_kinematics
contact_kinematics_jacobians = hopper_contact_kinematics_jacobians 
path = joinpath(@__DIR__, "..", "src/robots", String(name(model)), "expr/expr.jld2")
include(path_codegen)

include(joinpath(@__DIR__,"..", "src/robots/biped/model.jl"))
model = biped 
contact_kinematics = biped_contact_kinematics
contact_kinematics_jacobians = biped_contact_kinematics_jacobians 
path = joinpath(@__DIR__, "..", "src/robots", String(name(model)), "expr/expr.jld2")
include(path_codegen)

include(joinpath(@__DIR__,"..", "src/robots/quadruped/model.jl"))
model = quadruped 
contact_kinematics = quadruped_contact_kinematics
contact_kinematics_jacobians = quadruped_contact_kinematics_jacobians 
path = joinpath(@__DIR__, "..", "src/robots", String(name(model)), "expr/expr.jld2")
include(path_codegen)