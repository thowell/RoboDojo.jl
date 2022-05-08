using Test
using JLD2
using Symbolics
using LinearAlgebra
using Random
using BenchmarkTools
using FiniteDiff
using RoboDojo

# Utils
# include("utils.jl")
#
# # Solver
# include("solver/lu.jl")
# include("solver/random_qp.jl")
# include("solver/soc.jl")
#
# # Robots
# include("robots/hopper.jl")
# include("robots/biped.jl")
# include("robots/quadruped.jl")
#
# # Simulator
# include("simulator/utils.jl")
include("simulator/simulate.jl")
include("simulator/codegen.jl")
include("simulator/dynamics.jl")
