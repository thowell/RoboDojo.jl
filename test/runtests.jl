using Test
using ForwardDiff
using JLD2
using Symbolics
using StaticArrays
using LinearAlgebra
using Random
using BenchmarkTools
using RobotDojo

# Solver
include("solver/lu.jl")
include("solver/random_qp.jl")
include("solver/soc.jl")

# Simulator 
include("simulator/utils.jl")




