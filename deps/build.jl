append!(empty!(LOAD_PATH), Base.DEFAULT_LOAD_PATH)
using Pkg

# ################################################################################
# # Generate notebooks
# # ################################################################################
# exampledir = joinpath(@__DIR__, "..", "examples")
# Pkg.activate(exampledir)
# Pkg.instantiate()
# include(joinpath(exampledir, "generate_notebooks.jl"))

################################################################################
# Build simulation environments
################################################################################
robotdir = joinpath(@__DIR__, "..", "src/robots")
Pkg.activate(robotdir)
Pkg.instantiate()

using JLD2 
using Symbolics
using LinearAlgebra
using Scratch 

include("../src/utils.jl")
include("../src/simulator/model.jl")
include("../src/simulator/dimensions.jl")
include("../src/simulator/residual.jl")
include("../src/simulator/codegen.jl")

include("../src/robots/integrator.jl")
include("../src/robots/codegen.jl")

include("../src/solver/cones.jl")

include("../src/robots/hopper/model.jl")  # hopper special case
include("../src/robots/hopper/simulator.jl")

include("../src/robots/particle/model.jl")  # particle special case
include("../src/robots/particle/simulator.jl")

include("../src/robots/centroidal_quadruped/model.jl")  # centroidal quadruped special case
include("../src/robots/centroidal_quadruped/simulator.jl")

pkgdir = joinpath(@__DIR__, "..")
Pkg.activate(pkgdir)
path_robots = @get_scratch!("robots")

robots = [:hopper, :biped, :quadruped, :box, :particle, :centroidal_quadruped]

for robot in robots
    # robot model
    include(joinpath("../src/robots", String(robot), "model.jl"))

    # expr path
    path_expr = joinpath(path_robots, String(robot) * ".jld2")

    # kinematics
    contact_kinematics = eval(Symbol(String(robot) * "_contact_kinematics"))
    contact_kinematics_jacobians = eval(Symbol(String(robot) * "_contact_kinematics_jacobians"))

    # codegen
    r_model, rz_model, rθ_model = codegen_residual(eval(robot), codegen_dynamics(eval(robot))..., contact_kinematics, contact_kinematics_jacobians)

    @save path_expr r_model rz_model rθ_model
end

# # robot model
# robot = :box
# include(joinpath("../src/robots", String(robot), "model.jl"))

# # expr path
# path_expr = joinpath(path_robots, String(robot) * ".jld2")

# # kinematics
# contact_kinematics = eval(Symbol(String(robot) * "_contact_kinematics"))
# contact_kinematics_jacobians = eval(Symbol(String(robot) * "_contact_kinematics_jacobians"))

# # codegen
# mass_matrix, dynamics_bias = codegen_dynamics(eval(robot))
# r_model, rz_model, rθ_model = codegen_residual(eval(robot), mass_matrix, dynamics_bias, contact_kinematics, contact_kinematics_jacobians)

# @save path_expr r_model rz_model rθ_model

# # robot model
# robot = :centroidal_quadruped
# include(joinpath("../src/robots", String(robot), "model.jl"))
# include(joinpath("../src/robots", String(robot), "simulator.jl"))

# # expr path
# path_expr = joinpath(path_robots, String(robot) * ".jld2")

# # kinematics
# contact_kinematics = eval(Symbol(String(robot) * "_contact_kinematics"))
# contact_kinematics_jacobians = eval(Symbol(String(robot) * "_contact_kinematics_jacobians"))

# # codegen
# r_model, rz_model, rθ_model = codegen_residual(eval(robot), 
#     q -> mass_matrix(eval(robot), q), 
#     (q, q̇) -> dynamics_bias(eval(robot), q, q̇), 
#     contact_kinematics, contact_kinematics_jacobians)

# @save path_expr r_model rz_model rθ_model


