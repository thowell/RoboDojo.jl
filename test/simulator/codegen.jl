@testset "Simulator: codegen" begin

    robots = [hopper, biped, quadruped]

    for robot in robots 
        # kinematics
        contact_kinematics = RoboDojo.eval(Symbol(String(RoboDojo.name(robot)) * "_contact_kinematics"))
        contact_kinematics_jacobians = RoboDojo.eval(Symbol(String(RoboDojo.name(robot)) * "_contact_kinematics_jacobians"))

        # codegen
        mass_matrix, dynamics_bias = RoboDojo.codegen_dynamics(robot)
        r_model, rz_model, rθ_model = RoboDojo.codegen_residual(robot, mass_matrix, dynamics_bias, contact_kinematics, contact_kinematics_jacobians)

        # evaluate methods
        nz = num_var(robot) 
        nθ = num_data(robot, nf=length(friction_coefficients(robot)))
        r0 = zeros(nz) 
        z0 = rand(nz) 
        θ0 = rand(nθ) 
        μ0 = ones(1) 
        rz0 = rand(nz, nz) 
        rθ0 = rand(nz, nθ) 

        eval(r_model)(r0, z0, θ0, μ0)
        eval(rz_model)(rz0, z0, θ0)
        eval(rθ_model)(rθ0, z0, θ0)

        # test methods
        @test sum(r0) != 0.0
        @test sum(rz0) != 0.0
        @test sum(rθ0) != 0.0
    end
end
