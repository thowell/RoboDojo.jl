@testset "Simulator: codegen" begin
    path_src = joinpath(@__DIR__, "..", "..", "src")

    CODEGEN_SAVE = false 

    path_codegen = joinpath(path_src, "simulator/codegen.jl")
    @show path_codegen
    model = hopper 
    contact_kinematics = RobotDojo.hopper_contact_kinematics
    contact_kinematics_jacobians = RobotDojo.hopper_contact_kinematics_jacobians 
    status = include(path_codegen)
    # @test status 

    model = biped 
    contact_kinematics = RobotDojo.biped_contact_kinematics
    contact_kinematics_jacobians = RobotDojo.biped_contact_kinematics_jacobians 
    # status = include(path_codegen)
    # @test status 

    model = quadruped 
    contact_kinematics = RobotDojo.quadruped_contact_kinematics
    contact_kinematics_jacobians = RobotDojo.quadruped_contact_kinematics_jacobians 
    # status = include(path_codegen) 
    # @test status
end