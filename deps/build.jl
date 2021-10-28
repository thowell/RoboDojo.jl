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
include(joinpath("../src/solver/cones.jl"))

HOPPER_CODEGEN = :generate
include(joinpath(@__DIR__, "../src/robots/hopper/model.jl"))
include(joinpath(@__DIR__, "../src/robots/hopper/simulator.jl"))
include(joinpath(@__DIR__, "../src/robots/hopper/codegen.jl"))

BIPED_CODEGEN = :generate
include(joinpath(@__DIR__, "../src/robots/biped/model.jl"))
include(joinpath(@__DIR__, "../src/robots/biped/codegen.jl"))

QUADRUPED_CODEGEN = :generate
include(joinpath(@__DIR__, "../src/robots/quadruped/model.jl"))
include(joinpath(@__DIR__, "../src/robots/quadruped/codegen.jl"))