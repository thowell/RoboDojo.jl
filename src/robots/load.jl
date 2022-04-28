RESIDUAL_EXPR = Dict{String, Any}()
residual_expr(model::Model) = RESIDUAL_EXPR[String(name(model)) * "_r"]
jacobian_var_expr(model::Model) = RESIDUAL_EXPR[String(name(model)) * "_rz"]
jacobian_data_expr(model::Model) = RESIDUAL_EXPR[String(name(model)) * "_rθ"]

path_robots = @get_scratch!("robots")

robots = [hopper, biped, quadruped, halfquadruped, halfcheetah, box, particle, centroidal_quadruped]

for robot in robots
    path = joinpath(path_robots, String(name(robot)) * ".jld2")

    if !isfile(path)
        # kinematics
        contact_kinematics = eval(Symbol(String(name(robot)) * "_contact_kinematics"))
        contact_kinematics_jacobians = eval(Symbol(String(name(robot)) * "_contact_kinematics_jacobians"))

        # codegen
        r_model, rz_model, rθ_model = codegen_residual(eval(robot), codegen_dynamics(eval(robot))..., contact_kinematics, contact_kinematics_jacobians)
        @save path r_model rz_model rθ_model
    else
        @load path r_model rz_model rθ_model
    end

    RESIDUAL_EXPR[String(name(robot)) * "_r"] = eval(r_model)
    RESIDUAL_EXPR[String(name(robot)) * "_rz"] = eval(rz_model)
    RESIDUAL_EXPR[String(name(robot)) * "_rθ"] = eval(rθ_model)
end
