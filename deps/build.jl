append!(empty!(LOAD_PATH), Base.DEFAULT_LOAD_PATH)
using Pkg

################################################################################
# Generate notebooks
# ################################################################################
exampledir = joinpath(@__DIR__, "..", "examples")
Pkg.activate(exampledir)
Pkg.instantiate()
include(joinpath(exampledir, "generate_notebooks.jl"))

################################################################################
# Build simulation environments
################################################################################
pkgdir = joinpath(@__DIR__, "..")
Pkg.activate(pkgdir)

using JLD2 
using Symbolics
using LinearAlgebra
using Scratch 

include(joinpath("../src/utils.jl"))
include(joinpath("../src/simulator/model.jl"))
include(joinpath("../src/simulator/dimensions.jl"))
include(joinpath("../src/simulator/residual.jl"))
include(joinpath("../src/simulator/codegen.jl"))

include(joinpath("../src/robots/integrator.jl"))
include(joinpath("../src/robots/codegen.jl"))

include(joinpath("../src/solver/cones.jl"))

include(joinpath(@__DIR__,"..", "src/robots/hopper/model.jl"))  # hopper special case
include(joinpath(@__DIR__,"..", "src/robots/hopper/simulator.jl"))

path_robots = @get_scratch!("robots")

robots = [:hopper, :biped, :quadruped]

for robot in robots
    # robot model
    include(joinpath(@__DIR__,"..", "src/robots", String(robot), "model.jl"))

    # expr path
    path_expr = joinpath(path_robots, String(robot) * ".jld2")

    # kinematics
    contact_kinematics = eval(Symbol(String(robot) * "_contact_kinematics"))
    contact_kinematics_jacobians = eval(Symbol(String(robot) * "_contact_kinematics_jacobians"))

    # codegen
    mass_matrix, dynamics_bias = codegen_dynamics(eval(robot))
    r_model, rz_model, rθ_model = codegen_residual(eval(robot), mass_matrix, dynamics_bias, contact_kinematics, contact_kinematics_jacobians,
        save=true, path=path_expr)
end
